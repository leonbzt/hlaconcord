from pathlib import Path

import pytest

from hlaconcord.concordance import ConcordanceStatus, ConsensusRule
from hlaconcord.nomenclature import ReductionBasis
from hlaconcord.pipeline import InputSpec, run

FIXTURES = Path(__file__).parent / "fixtures"


# -- InputSpec.parse_arg ------------------------------------------------------

def test_parse_arg_splits_tool_and_path():
    spec = InputSpec.parse_arg("optitype:/data/s1_result.tsv")
    assert spec.tool == "optitype"
    assert spec.path == Path("/data/s1_result.tsv")


def test_parse_arg_allows_colon_in_path():
    spec = InputSpec.parse_arg("hla-la:C:/win/s1_bestguess_G.txt")
    assert spec.tool == "hla-la"
    assert spec.path == Path("C:/win/s1_bestguess_G.txt")


@pytest.mark.parametrize("bad", ["optitype", ":path", "optitype:", ""])
def test_parse_arg_rejects_malformed(bad):
    with pytest.raises(ValueError, match="expected 'tool:path'"):
        InputSpec.parse_arg(bad)


# -- run ----------------------------------------------------------------------

def _specs():
    return [
        InputSpec("optitype", FIXTURES / "optitype" / "s1_result.tsv"),
        InputSpec("arcasHLA", FIXTURES / "arcashla" / "s1.genotype.json"),
        InputSpec("hla-la", FIXTURES / "hla_la" / "s1_bestguess_G.txt"),
        InputSpec("hla-hd", FIXTURES / "hla_hd" / "s1_final.result.txt"),
    ]


def test_run_harmonizes_four_tools(nom):
    result = run(_specs(), nom, db_version="3.55.0")
    assert result.db_version == "3.55.0"
    assert result.basis is ReductionBasis.LGX
    assert result.rule is ConsensusRule.MAJORITY
    assert result.samples == ["s1"]
    assert result.tools == ["arcasHLA", "hla-hd", "hla-la", "optitype"]
    by_locus = {r.locus: r for r in result.concordance}
    assert by_locus["A"].status is ConcordanceStatus.CONCORDANT
    assert by_locus["A"].consensus == ("A*01:01", "A*02:01")
    assert by_locus["DRB1"].n_tools == 3  # OptiType is class I only


def test_run_sample_override(nom):
    specs = [InputSpec("optitype", FIXTURES / "optitype" / "s1_result.tsv", sample="PATIENT_X")]
    result = run(specs, nom, db_version="3.55.0")
    assert result.samples == ["PATIENT_X"]
    assert all(gc.sample_id == "PATIENT_X" for gc in result.calls)


def test_run_source_db_version_stamped(nom):
    specs = [
        InputSpec("optitype", FIXTURES / "optitype" / "s1_result.tsv", source_db_version="3.53.0")
    ]
    result = run(specs, nom, db_version="3.55.0")
    assert all(gc.source_db_version == "3.53.0" for gc in result.calls)


def test_run_basis_and_rule_propagate(nom):
    result = run(
        _specs(), nom, db_version="3.55.0",
        basis=ReductionBasis.G, rule=ConsensusRule.UNANIMOUS,
    )
    assert result.basis is ReductionBasis.G
    assert result.rule is ConsensusRule.UNANIMOUS
    assert all(r.basis == "g" for r in result.concordance)
