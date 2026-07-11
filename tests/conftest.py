import shutil
from pathlib import Path

import pytest

from hlaharm.nomenclature import InHouseNomenclature, ReferenceData

FIXTURES = Path(__file__).parent / "fixtures"
FIXTURE_DB = FIXTURES / "db_3550"
_RELEASE_FILES = (
    "Allelelist.txt",
    "Allelelist_history.txt",
    "hla_nom_g.txt",
    "hla_nom_p.txt",
)


@pytest.fixture(scope="session")
def reference() -> ReferenceData:
    return ReferenceData.from_dir(FIXTURE_DB, version="3.55.0")


@pytest.fixture()
def nom(reference: ReferenceData) -> InHouseNomenclature:
    return InHouseNomenclature(reference)


@pytest.fixture()
def data_root(tmp_path: Path) -> Path:
    """A cache root laid out like a real install, with the fixture DB as 3.55.0."""
    root = tmp_path / "hlaharm-data"
    release = root / "3.55.0"
    release.mkdir(parents=True)
    for name in _RELEASE_FILES:
        shutil.copy(FIXTURE_DB / name, release / name)
    return root
