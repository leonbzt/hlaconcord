from pathlib import Path

import pytest

from hlaconcord import parsers
from hlaconcord.parsers import arcashla, optitype

FIXTURES = Path(__file__).parent / "fixtures"


# -- OptiType -----------------------------------------------------------------

def test_optitype_parses_loci_and_sample_id():
    calls = optitype.parse(FIXTURES / "optitype" / "s1_result.tsv")
    by_locus = {c.locus: c for c in calls}
    assert set(by_locus) == {"A", "B", "C"}
    assert all(c.sample_id == "s1" and c.tool == "optitype" for c in calls)
    assert [a.raw for a in by_locus["A"].alleles] == ["A*01:01", "A*02:01"]


def test_optitype_keeps_homozygous_as_two_calls():
    calls = optitype.parse(FIXTURES / "optitype" / "s1_result.tsv")
    b = next(c for c in calls if c.locus == "B")
    assert [a.raw for a in b.alleles] == ["B*07:02", "B*07:02"]


# -- arcasHLA -----------------------------------------------------------------

def test_arcashla_parses_class_i_and_ii():
    calls = arcashla.parse(FIXTURES / "arcashla" / "s1.genotype.json")
    by_locus = {c.locus: [a.raw for a in c.alleles] for c in calls}
    assert by_locus["A"] == ["A*01:01:01", "A*02:01:01"]
    assert by_locus["DRB1"] == ["DRB1*15:01:01", "DRB1*15:01:01"]
    assert all(c.sample_id == "s1" and c.tool == "arcasHLA" for c in calls)


def test_arcashla_homozygous_single_entry_is_one_call():
    calls = arcashla.parse(FIXTURES / "arcashla" / "s1.genotype.json")
    c = next(c for c in calls if c.locus == "C")
    assert [a.raw for a in c.alleles] == ["C*07:02:01"]


# -- registry -----------------------------------------------------------------

def test_registry_dispatch_is_case_insensitive():
    assert parsers.parse("OptiType", FIXTURES / "optitype" / "s1_result.tsv")
    assert parsers.parse("arcasHLA", FIXTURES / "arcashla" / "s1.genotype.json")


def test_registry_rejects_unknown_tool():
    with pytest.raises(ValueError):
        parsers.parse("nope", FIXTURES / "optitype" / "s1_result.tsv")
