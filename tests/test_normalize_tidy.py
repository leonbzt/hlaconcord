from pathlib import Path

from hlaharm import parsers
from hlaharm.model import AlleleCall
from hlaharm.nomenclature import ReductionBasis, ValidationStatus
from hlaharm.normalize import normalize, normalize_allele
from hlaharm.tidy import COLUMNS, tidy_rows

FIXTURES = Path(__file__).parent / "fixtures"


def _by_locus(calls):
    return {c.locus: c for c in calls}


# -- normalization enriches through the facade --------------------------------

def test_normalize_optitype_against_reference(nom):
    calls = parsers.parse("optitype", FIXTURES / "optitype" / "s1_result.tsv")
    normalize(calls, nom, ReductionBasis.G)
    a1 = _by_locus(calls)["A"].alleles[0]
    assert a1.raw == "A*01:01"
    assert a1.normalized == "A*01:01"
    assert a1.resolution == "FIELD_2"
    assert a1.validation is ValidationStatus.VALID_REDUCTION
    assert a1.reduced == "A*01:01:01G"  # 2-field call projected to its G-group


def test_normalize_arcashla_reductions(nom):
    calls = parsers.parse("arcasHLA", FIXTURES / "arcashla" / "s1.genotype.json")
    normalize(calls, nom, ReductionBasis.G)
    reduced = {c.locus: c.alleles[0].reduced for c in calls}
    assert reduced["A"] == "A*01:01:01G"
    assert reduced["B"] == "B*07:02:01G"
    assert reduced["DRB1"] == "DRB1*15:01:01G"
    # C*07:02:01:01 is a G-group singleton in the fixture DB -> allele kept unchanged
    assert reduced["C"] == "C*07:02:01"


def test_normalize_flags_unparseable_without_coercing(nom):
    call = AlleleCall(raw="not-an-allele", tool="optitype")
    normalize_allele(call, nom, ReductionBasis.G)
    assert call.parse_error is not None
    assert call.normalized is None and call.reduced is None
    assert call.raw == "not-an-allele"  # raw always retained


def test_normalize_preserves_null(nom):
    call = AlleleCall(raw="A*02:11N", tool="arcasHLA")
    normalize_allele(call, nom, ReductionBasis.G)
    assert call.is_null is True
    assert call.reduced == "A*02:11N"  # null survives reduction


# -- tidy table ---------------------------------------------------------------

def test_tidy_rows_shape_and_content(nom):
    calls = parsers.parse("optitype", FIXTURES / "optitype" / "s1_result.tsv")
    normalize(calls, nom, ReductionBasis.G)
    rows = tidy_rows(calls)
    assert len(rows) == 6  # 3 loci x 2 alleles
    assert all(set(r) == set(COLUMNS) for r in rows)
    a1 = next(r for r in rows if r["locus"] == "A" and r["raw"] == "A*01:01")
    assert a1["sample_id"] == "s1"
    assert a1["tool"] == "optitype"
    assert a1["reduced"] == "A*01:01:01G"
    assert a1["validation"] == "valid_reduction"
    assert a1["is_null"] == "false"
