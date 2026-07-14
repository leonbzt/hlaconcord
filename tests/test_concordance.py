from pathlib import Path

from hlaconcord import parsers
from hlaconcord.concordance import (
    ConcordanceFlag,
    ConcordanceStatus,
    ConsensusRule,
    concordance,
)
from hlaconcord.model import AlleleCall, GenotypeCall
from hlaconcord.nomenclature import ReductionBasis
from hlaconcord.normalize import normalize

FIXTURES = Path(__file__).parent / "fixtures"


def ac(reduced, *, normalized=None, resolution="FIELD_2", is_null=False, accession=None):
    return AlleleCall(
        raw=reduced,
        tool="",
        reduced=reduced,
        normalized=normalized or reduced,
        resolution=resolution,
        is_null=is_null,
        accession=accession,
    )


def gc(tool, locus, keys, sample="s1"):
    alleles = [k if isinstance(k, AlleleCall) else ac(k) for k in keys]
    return GenotypeCall(sample_id=sample, tool=tool, locus=locus, alleles=alleles)


def _one(results, locus="A"):
    return next(r for r in results if r.locus == locus)


# -- basic agreement ----------------------------------------------------------

def test_concordant_heterozygous():
    calls = [
        gc("optitype", "A", ["A*01:01", "A*02:01"]),
        gc("arcasHLA", "A", ["A*02:01", "A*01:01"]),  # order is arbitrary
    ]
    r = _one(concordance(calls))
    assert r.status is ConcordanceStatus.CONCORDANT
    assert r.consensus == ("A*01:01", "A*02:01")
    assert ConcordanceFlag.DISCORDANCE not in r.flags
    assert r.agreement("A*02:01") == "2/2"


def test_homozygous_reported_inconsistently_is_still_concordant():
    calls = [
        gc("optitype", "A", ["A*02:01", "A*02:01"]),  # homozygous, listed twice
        gc("arcasHLA", "A", ["A*02:01"]),  # homozygous, listed once
    ]
    r = _one(concordance(calls))
    assert r.status is ConcordanceStatus.CONCORDANT
    assert r.consensus == ("A*02:01",)
    assert ConcordanceFlag.SINGLETON not in r.flags


def test_discordant_with_singletons_and_majority_consensus():
    calls = [
        gc("optitype", "A", ["A*02:01", "A*01:01"]),
        gc("arcasHLA", "A", ["A*02:01", "A*03:01"]),
    ]
    r = _one(concordance(calls))
    assert r.status is ConcordanceStatus.DISCORDANT
    assert ConcordanceFlag.DISCORDANCE in r.flags
    assert ConcordanceFlag.SINGLETON in r.flags
    # only the shared allele survives majority (support*2 > n_tools)
    assert r.consensus == ("A*02:01",)


def test_single_tool_is_not_discordant():
    r = _one(concordance([gc("optitype", "A", ["A*02:01", "A*01:01"])]))
    assert r.status is ConcordanceStatus.SINGLE_TOOL
    assert ConcordanceFlag.DISCORDANCE not in r.flags
    assert r.n_tools == 1


# -- consensus rules ----------------------------------------------------------

def test_majority_vs_unanimous():
    calls = [
        gc("optitype", "A", ["A*02:01"]),
        gc("arcasHLA", "A", ["A*02:01"]),
        gc("hla-la", "A", ["A*11:01"]),  # 2/3 for A*02:01
    ]
    maj = _one(concordance(calls, rule=ConsensusRule.MAJORITY))
    assert maj.consensus == ("A*02:01",)  # 2 of 3 > half
    una = _one(concordance(calls, rule=ConsensusRule.UNANIMOUS))
    assert una.consensus == ()  # nothing called by all three


# -- flags --------------------------------------------------------------------

def test_null_allele_flagged():
    calls = [
        gc("optitype", "A", [ac("A*02:11N", is_null=True), ac("A*01:01")]),
        gc("arcasHLA", "A", [ac("A*02:11N", is_null=True), ac("A*01:01")]),
    ]
    assert ConcordanceFlag.NULL_ALLELE in _one(concordance(calls)).flags


def test_resolution_conflict_on_underresolved_call():
    # A 1-field call is prefix-compatible with a 2-field call: same allele? can't confirm.
    calls = [
        gc("optitype", "A", ["A*02"]),
        gc("arcasHLA", "A", ["A*02:01"]),
    ]
    flags = _one(concordance(calls)).flags
    assert ConcordanceFlag.RESOLUTION_CONFLICT in flags


def test_clean_depth_difference_is_not_a_conflict():
    # Same lgx key from different depths (2-field vs 3-field) is NOT a conflict or skew.
    calls = [
        gc("optitype", "A", [ac("A*02:01", normalized="A*02:01", resolution="FIELD_2")]),
        gc("arcasHLA", "A", [ac("A*02:01", normalized="A*02:01:01", resolution="FIELD_3")]),
    ]
    flags = _one(concordance(calls)).flags
    assert ConcordanceFlag.RESOLUTION_CONFLICT not in flags
    assert ConcordanceFlag.VERSION_SKEW_RESOLVED not in flags


def test_version_skew_resolved_needs_shared_accession():
    # Two tools (different DB versions) emit different names that resolve to one accession.
    calls = [
        gc("optitype", "A", [ac("A*02:01", normalized="A*02:01:01", accession="HLA00005")]),
        gc("arcasHLA", "A", [ac("A*02:01", normalized="A*02:01:99", accession="HLA00005")]),
    ]
    flags = _one(concordance(calls)).flags
    assert ConcordanceFlag.VERSION_SKEW_RESOLVED in flags


def test_distinct_alleles_collapsing_to_one_key_is_not_skew():
    # Different accessions reducing to one lgx key = lgx collapse, not version skew.
    calls = [
        gc("optitype", "A", [ac("A*02:01", normalized="A*02:01:01", accession="HLA00005")]),
        gc("arcasHLA", "A", [ac("A*02:01", normalized="A*02:01:05", accession="HLA99999")]),
    ]
    assert ConcordanceFlag.VERSION_SKEW_RESOLVED not in _one(concordance(calls)).flags


def test_invalid_call_flagged():
    bad = AlleleCall(raw="junk", tool="", parse_error="unparseable")  # reduced stays None
    calls = [
        gc("optitype", "A", [bad, ac("A*01:01")]),
        gc("arcasHLA", "A", ["A*01:01"]),
    ]
    assert ConcordanceFlag.INVALID_CALL in _one(concordance(calls)).flags


# -- integration through the real parsers -------------------------------------

def test_integration_optitype_vs_arcashla(nom):
    calls = parsers.parse("optitype", FIXTURES / "optitype" / "s1_result.tsv")
    calls += parsers.parse("arcasHLA", FIXTURES / "arcashla" / "s1.genotype.json")
    normalize(calls, nom, ReductionBasis.LGX)
    results = {r.locus: r for r in concordance(calls, basis=ReductionBasis.LGX)}
    # A/B/C typed by both tools and agree under lgx; DRB1 only by arcasHLA.
    assert results["A"].status is ConcordanceStatus.CONCORDANT
    assert results["A"].consensus == ("A*01:01", "A*02:01")
    assert results["B"].status is ConcordanceStatus.CONCORDANT
    assert results["DRB1"].status is ConcordanceStatus.SINGLE_TOOL


def test_integration_all_four_tools_harmonize(nom):
    # OptiType (2-field), arcasHLA (3-field), HLA-LA (G-group), HLA-HD (HLA- prefix, 3-field)
    # — four resolutions of the same genotype must reconcile to one lgx key per allele.
    calls = parsers.parse("optitype", FIXTURES / "optitype" / "s1_result.tsv")
    calls += parsers.parse("arcasHLA", FIXTURES / "arcashla" / "s1.genotype.json")
    calls += parsers.parse("hla-la", FIXTURES / "hla_la" / "s1_bestguess_G.txt")
    calls += parsers.parse("hla-hd", FIXTURES / "hla_hd" / "s1_final.result.txt")
    normalize(calls, nom, ReductionBasis.LGX)
    results = {r.locus: r for r in concordance(calls, basis=ReductionBasis.LGX)}

    a = results["A"]
    assert a.n_tools == 4 and a.status is ConcordanceStatus.CONCORDANT
    assert a.consensus == ("A*01:01", "A*02:01")
    assert a.agreement("A*02:01") == "4/4"
    for locus in ("B", "C"):
        assert results[locus].status is ConcordanceStatus.CONCORDANT
    # class II typed by arcasHLA + HLA-LA + HLA-HD (OptiType is class I only)
    assert results["DRB1"].n_tools == 3
    assert results["DRB1"].status is ConcordanceStatus.CONCORDANT
    assert results["DRB1"].consensus == ("DRB1*15:01",)
