from hlaconcord import gl
from hlaconcord.concordance import concordance
from hlaconcord.model import AlleleCall, GenotypeCall


def ac(reduced):
    return AlleleCall(raw=reduced, tool="", reduced=reduced, normalized=reduced)


def gc(tool, locus, keys, sample="s1"):
    return GenotypeCall(
        sample_id=sample, tool=tool, locus=locus, alleles=[ac(k) for k in keys]
    )


# -- primitives ---------------------------------------------------------------

def test_with_prefix_is_idempotent():
    assert gl.with_prefix("A*02:01") == "HLA-A*02:01"
    assert gl.with_prefix("HLA-A*02:01") == "HLA-A*02:01"


def test_locus_block_joins_gene_copies():
    assert gl.locus_block(["A*01:01", "A*02:01"]) == "HLA-A*01:01+HLA-A*02:01"


def test_genotype_string_orders_loci_and_skips_empty():
    out = gl.genotype_string([("DRB1", ["DRB1*15:01"]), ("A", ["A*01:01"]), ("B", [])])
    # A before DRB1 (canonical order); B omitted (no consensus)
    assert out == "HLA-A*01:01^HLA-DRB1*15:01"


def test_unknown_locus_sorts_after_known():
    out = gl.genotype_string([("XYZ", ["XYZ*01:01"]), ("A", ["A*01:01"])])
    assert out == "HLA-A*01:01^HLA-XYZ*01:01"


# -- from concordance results -------------------------------------------------

def test_consensus_gl_string_from_results():
    calls = [
        gc("optitype", "A", ["A*01:01", "A*02:01"]),
        gc("arcasHLA", "A", ["A*01:01", "A*02:01"]),
        gc("optitype", "B", ["B*07:02", "B*07:02"]),
        gc("arcasHLA", "B", ["B*07:02"]),
    ]
    results = concordance(calls)
    assert gl.consensus_gl_string(results) == "HLA-A*01:01+HLA-A*02:01^HLA-B*07:02"


def test_gl_strings_by_sample_groups():
    calls = [
        gc("optitype", "A", ["A*01:01", "A*02:01"], sample="s1"),
        gc("arcasHLA", "A", ["A*01:01", "A*02:01"], sample="s1"),
        gc("optitype", "A", ["A*03:01", "A*11:01"], sample="s2"),
        gc("arcasHLA", "A", ["A*03:01", "A*11:01"], sample="s2"),
    ]
    by_sample = gl.gl_strings_by_sample(concordance(calls))
    assert by_sample["s1"] == "HLA-A*01:01+HLA-A*02:01"
    assert by_sample["s2"] == "HLA-A*03:01+HLA-A*11:01"


def test_no_consensus_locus_yields_empty_string():
    # unanimous rule with a lone singleton -> empty consensus -> empty GL
    calls = [
        gc("optitype", "A", ["A*01:01"]),
        gc("arcasHLA", "A", ["A*02:01"]),
    ]
    from hlaconcord.concordance import ConsensusRule

    results = concordance(calls, rule=ConsensusRule.UNANIMOUS)
    assert gl.consensus_gl_string(results) == ""
