import json
from pathlib import Path

import pytest

from hlaconcord import db
from hlaconcord.cli import main

FIXTURES = Path(__file__).parent / "fixtures"
FIXTURE_DB = FIXTURES / "db_3550"


def _run_inputs(data_root, out=None, extra=None):
    args = [
        "run",
        "--inputs",
        f"optitype:{FIXTURES / 'optitype' / 's1_result.tsv'}",
        f"arcasHLA:{FIXTURES / 'arcashla' / 's1.genotype.json'}",
        f"hla-la:{FIXTURES / 'hla_la' / 's1_bestguess_G.txt'}",
        f"hla-hd:{FIXTURES / 'hla_hd' / 's1_final.result.txt'}",
        "--db", "3.55.0",
        "--data-dir", str(data_root),
    ]
    if out is not None:
        args += ["-o", str(out)]
    if extra:
        args += extra
    return main(args)


# -- run ----------------------------------------------------------------------

def test_run_writes_all_outputs(data_root, tmp_path):
    out = tmp_path / "out"
    code = _run_inputs(data_root, out=out, extra=["--gl"])
    assert code == 0
    for name in ("tidy.tsv", "concordance.tsv", "concordance.json", "gl_strings.tsv"):
        assert (out / name).exists(), name
    payload = json.loads((out / "concordance.json").read_text())
    assert payload["db_version"] == "3.55.0"
    assert payload["basis"] == "lgx"
    assert payload["samples"]["s1"]["loci"]["A"]["status"] == "concordant"
    assert payload["samples"]["s1"]["gl_string"].startswith("HLA-A*01:01+HLA-A*02:01")


def test_run_stdout_only(data_root, capsys):
    code = _run_inputs(data_root)
    assert code == 0
    out = capsys.readouterr().out
    assert "sample s1" in out
    assert "OK A" in out


def test_run_discordant_exits_1(data_root, tmp_path, capsys):
    opti = tmp_path / "s2_result.tsv"
    opti.write_text(
        "\tA1\tA2\tB1\tB2\tC1\tC2\tReads\tObjective\n"
        "0\tA*01:01\tA*02:01\tB*07:02\tB*07:02\tC*07:02\tC*07:02\t9\t9\n"
    )
    arca = tmp_path / "s2.genotype.json"
    arca.write_text('{"A": ["A*01:01", "A*11:01"]}')  # disagrees on the 2nd A allele
    code = main([
        "run",
        "--inputs", f"optitype:{opti}", f"arcasHLA:{arca}",
        "--sample", "s2", "--db", "3.55.0", "--data-dir", str(data_root),
    ])
    assert code == 1
    assert "XX A" in capsys.readouterr().out


def test_run_missing_db_exits_2(data_root, capsys):
    code = _run_inputs(data_root, extra=["--db", "9.99.0"])
    # override the --db already in _run_inputs by appending a second one wins in argparse
    assert code == 2
    assert "not installed" in capsys.readouterr().err


def test_run_samplesheet(data_root, tmp_path):
    sheet = tmp_path / "samples.csv"
    sheet.write_text(
        "sample,tool,path,db_version\n"
        f"s1,optitype,{FIXTURES / 'optitype' / 's1_result.tsv'},3.55.0\n"
        f"s1,arcasHLA,{FIXTURES / 'arcashla' / 's1.genotype.json'},\n"
    )
    out = tmp_path / "out"
    code = main([
        "run", "--samplesheet", str(sheet), "--db", "3.55.0",
        "--data-dir", str(data_root), "-o", str(out), "--quiet",
    ])
    assert code == 0
    payload = json.loads((out / "concordance.json").read_text())
    assert payload["samples"]["s1"]["loci"]["A"]["status"] == "concordant"


def test_run_samplesheet_missing_column(data_root, tmp_path):
    sheet = tmp_path / "bad.csv"
    sheet.write_text("sample,tool\ns1,optitype\n")  # no path column
    code = main([
        "run", "--samplesheet", str(sheet), "--db", "3.55.0", "--data-dir", str(data_root),
    ])
    assert code == 2


def test_run_relative_samplesheet_paths_resolve_against_sheet_dir(data_root, tmp_path):
    # copy an input next to the sheet and reference it by bare filename
    opti = tmp_path / "s1_result.tsv"
    opti.write_text((FIXTURES / "optitype" / "s1_result.tsv").read_text())
    sheet = tmp_path / "samples.csv"
    sheet.write_text("sample,tool,path\ns1,optitype,s1_result.tsv\n")
    out = tmp_path / "out"
    code = main([
        "run", "--samplesheet", str(sheet), "--db", "3.55.0",
        "--data-dir", str(data_root), "-o", str(out), "--quiet",
    ])
    assert code == 0


# -- validate / normalize -----------------------------------------------------

def test_validate_exit_codes(data_root, capsys):
    code = main(["validate", "--db", "3.55.0", "--data-dir", str(data_root), "A*01:01:01:01"])
    assert code == 0
    assert "exact" in capsys.readouterr().out

    code = main(["validate", "--db", "3.55.0", "--data-dir", str(data_root), "A*99:99"])
    assert code == 1  # unknown -> nonzero for scripting


def test_validate_unparseable(data_root, capsys):
    code = main(["validate", "--db", "3.55.0", "--data-dir", str(data_root), "not-an-allele"])
    assert code == 1
    assert "unparseable" in capsys.readouterr().out


def test_normalize_preserves_null_and_shows_reduction(data_root, capsys):
    code = main([
        "normalize", "--db", "3.55.0", "--data-dir", str(data_root),
        "A*0201", "A*29:112N",
    ])
    assert code == 0
    out = capsys.readouterr().out
    assert "A*02:01" in out
    assert "A*29:112N" in out  # null suffix survives into the reduced column
    assert "yes" in out  # is-null column


# -- db -----------------------------------------------------------------------

def test_db_list_and_path(data_root, capsys):
    assert main(["db", "list", "--data-dir", str(data_root)]) == 0
    out = capsys.readouterr().out
    assert "3.55.0" in out

    assert main(["db", "path", "--data-dir", str(data_root)]) == 0
    assert str(data_root) in capsys.readouterr().out


def test_db_update_and_pin(tmp_path, monkeypatch, capsys):
    def fetch(url, dest):
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes((FIXTURE_DB / dest.name).read_bytes())

    monkeypatch.setattr(db, "_fetch_to", fetch)
    root = tmp_path / "cache"
    assert main(["db", "update", "3.55.0", "--data-dir", str(root)]) == 0
    assert "installed 3.55.0" in capsys.readouterr().out
    assert db.pinned_version(root) == "3.55.0"

    # re-pin explicitly
    assert main(["db", "pin", "3.55.0", "--data-dir", str(root)]) == 0
    assert "pinned 3.55.0" in capsys.readouterr().out


def test_db_pin_missing_exits_2(data_root, capsys):
    assert main(["db", "pin", "9.99.0", "--data-dir", str(data_root)]) == 2
    assert "cannot pin" in capsys.readouterr().err


# -- top-level ----------------------------------------------------------------

def test_no_subcommand_errors():
    with pytest.raises(SystemExit):
        main([])
