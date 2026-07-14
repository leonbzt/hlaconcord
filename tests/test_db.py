from pathlib import Path

import pytest

from hlaconcord import db
from hlaconcord.nomenclature import ReferenceData

FIXTURE_DB = Path(__file__).parent / "fixtures" / "db_3550"
_RELEASE_FILES = (
    "Allelelist.txt",
    "Allelelist_history.txt",
    "hla_nom_g.txt",
    "hla_nom_p.txt",
)


# -- version <-> branch -------------------------------------------------------

@pytest.mark.parametrize(
    "version,branch",
    [("3.55.0", "3550"), ("3.51.0", "3510"), ("3.59.0", "3590")],
)
def test_branch_of_and_version_round_trip(version, branch):
    # Modern IPD convention: major + 2-digit minor + patch (3.55.0 -> 3550).
    assert db.branch_of(version) == branch
    assert db.version_of(branch) == version


def test_branch_of_passes_through_branch_form():
    assert db.branch_of("3550") == "3550"


def test_version_of_passes_through_version_and_nondigit():
    assert db.version_of("3.55.0") == "3.55.0"
    assert db.version_of("main") == "main"


# -- discovery ----------------------------------------------------------------

def test_find_release_accepts_version_and_branch_dirs(tmp_path):
    (tmp_path / "3550").mkdir()
    for name in _RELEASE_FILES:
        (tmp_path / "3550" / name).write_text("x")
    assert db.find_release(tmp_path, "3.55.0") == tmp_path / "3550"
    assert db.find_release(tmp_path, "3550") == tmp_path / "3550"
    assert db.find_release(tmp_path, "3.99.0") is None


def test_incomplete_dir_is_not_a_release(tmp_path):
    (tmp_path / "3.55.0").mkdir()
    (tmp_path / "3.55.0" / "Allelelist.txt").write_text("x")  # missing nom files
    assert db.find_release(tmp_path, "3.55.0") is None
    assert db.list_releases(tmp_path) == []


def test_list_releases_newest_first(data_root):
    (data_root / "3.51.0").mkdir()
    for name in _RELEASE_FILES:
        (data_root / "3.51.0" / name).write_text("x")
    assert db.list_releases(data_root) == ["3.55.0", "3.51.0"]


# -- resolve_version ----------------------------------------------------------

def test_resolve_single_release_is_auto_selected(data_root):
    assert db.resolve_version(data_root, None) == "3.55.0"


def test_resolve_explicit_missing_raises(data_root):
    with pytest.raises(db.DatabaseError, match="not installed"):
        db.resolve_version(data_root, "3.99.0")


def test_resolve_none_installed_raises(tmp_path):
    with pytest.raises(db.DatabaseError, match="no IPD-IMGT/HLA release"):
        db.resolve_version(tmp_path, None)


def test_resolve_ambiguous_without_pin_raises(data_root):
    (data_root / "3.51.0").mkdir()
    for name in _RELEASE_FILES:
        (data_root / "3.51.0" / name).write_text("x")
    with pytest.raises(db.DatabaseError, match="none pinned"):
        db.resolve_version(data_root, None)


def test_pin_then_resolve_honours_pin(data_root):
    (data_root / "3.51.0").mkdir()
    for name in _RELEASE_FILES:
        (data_root / "3.51.0" / name).write_text("x")
    db.pin(data_root, "3.51.0")
    assert db.pinned_version(data_root) == "3.51.0"
    assert db.resolve_version(data_root, None) == "3.51.0"
    # explicit --db still overrides the pin
    assert db.resolve_version(data_root, "3.55.0") == "3.55.0"


def test_pin_missing_release_raises(data_root):
    with pytest.raises(db.DatabaseError, match="cannot pin"):
        db.pin(data_root, "3.99.0")


# -- load_reference -----------------------------------------------------------

def test_load_reference_returns_versioned_data(data_root):
    ref = db.load_reference(data_root, "3.55.0")
    assert isinstance(ref, ReferenceData)
    assert ref.version == "3.55.0"
    assert ref.is_allele("A*01:01:01:01")


# -- update (fetch monkeypatched to copy the fixture DB) ----------------------

def _fake_fetch(monkeypatch):
    def fetch(url, dest):
        name = dest.name
        src = FIXTURE_DB / name
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(src.read_bytes())
    monkeypatch.setattr(db, "_fetch_to", fetch)


def test_update_installs_and_pins(tmp_path, monkeypatch):
    _fake_fetch(monkeypatch)
    result = db.update(tmp_path, "3.55.0")
    assert result.version == "3.55.0"
    assert result.directory == tmp_path / "3.55.0"
    assert not result.skipped_history
    assert db.find_release(tmp_path, "3.55.0") is not None
    assert db.pinned_version(tmp_path) == "3.55.0"  # pinned by default


def test_update_no_history_skips_and_no_pin(tmp_path, monkeypatch):
    _fake_fetch(monkeypatch)
    result = db.update(tmp_path, "3.55.0", with_history=False, set_pin=False)
    assert result.skipped_history
    assert "Allelelist_history.txt" not in result.files
    assert not (tmp_path / "3.55.0" / "Allelelist_history.txt").exists()
    assert db.pinned_version(tmp_path) is None
    # still a usable release (history is optional)
    assert db.find_release(tmp_path, "3.55.0") is not None


def test_update_uses_branch_url(tmp_path, monkeypatch):
    seen = []

    def fetch(url, dest):
        seen.append(url)
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text("x")

    monkeypatch.setattr(db, "_fetch_to", fetch)
    db.update(tmp_path, "3.55.0")
    assert any("/IMGTHLA/3550/Allelelist.txt" in u for u in seen)
    assert any("/IMGTHLA/3550/wmda/hla_nom_g.txt" in u for u in seen)
