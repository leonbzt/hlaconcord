from pathlib import Path

from hlaconcord import parsers
from hlaconcord.parsers import hla_hd, hla_la

FIXTURES = Path(__file__).parent / "fixtures"


# -- HLA-LA -------------------------------------------------------------------

def test_hla_la_pairs_chromosome_rows_per_locus():
    calls = hla_la.parse(FIXTURES / "hla_la" / "s1_bestguess_G.txt")
    by_locus = {c.locus: c for c in calls}
    assert set(by_locus) == {"A", "B", "C", "DRB1"}
    assert [a.raw for a in by_locus["A"].alleles] == ["A*01:01:01G", "A*02:01:01G"]
    assert all(c.sample_id == "s1" and c.tool == "hla-la" for c in calls)


def test_hla_la_captures_quality():
    calls = hla_la.parse(FIXTURES / "hla_la" / "s1_bestguess_G.txt")
    b = next(c for c in calls if c.locus == "B")
    assert b.alleles[0].quality == 1.0
    assert b.alleles[1].quality == 0.98


def test_hla_la_sample_id_falls_back_to_dir_for_generic_filename(tmp_path):
    d = tmp_path / "S9" / "hla"
    d.mkdir(parents=True)
    f = d / "R1_bestguess_G.txt"
    f.write_text("Locus\tChromosome\tAllele\tQ1\nA\t1\tA*01:01:01G\t1\nA\t2\tA*02:01:01G\t1\n")
    assert all(c.sample_id == "S9" for c in hla_la.parse(f))


# -- HLA-HD -------------------------------------------------------------------

def test_hla_hd_parses_prefixed_alleles_and_skips_absent():
    calls = hla_hd.parse(FIXTURES / "hla_hd" / "s1_final.result.txt")
    by_locus = {c.locus: [a.raw for a in c.alleles] for c in calls}
    assert set(by_locus) == {"A", "B", "C", "DRB1"}  # DQA1 "-" and DQB1 "Not typed" dropped
    assert by_locus["A"] == ["HLA-A*01:01:01", "HLA-A*02:01:01"]
    assert all(r.startswith("HLA-") for rs in by_locus.values() for r in rs)


def test_hla_hd_homozygous_kept_as_two():
    calls = hla_hd.parse(FIXTURES / "hla_hd" / "s1_final.result.txt")
    b = next(c for c in calls if c.locus == "B")
    assert [a.raw for a in b.alleles] == ["HLA-B*07:02:01", "HLA-B*07:02:01"]


# -- registry -----------------------------------------------------------------

def test_registry_knows_all_four_tools():
    assert set(parsers.TOOLS) == {"optitype", "arcashla", "hla-la", "hla-hd"}
    assert parsers.parse("HLA-LA", FIXTURES / "hla_la" / "s1_bestguess_G.txt")
    assert parsers.parse("hla-hd", FIXTURES / "hla_hd" / "s1_final.result.txt")
