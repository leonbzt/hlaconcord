"""Reference-database management (PLAN.md §9).

hlaharm validates and reduces against a pinned IPD-IMGT/HLA release. Rather than
redistribute the database (licensing is deferred, PLAN.md §15.5), releases are
fetched on demand into a local cache and every output records the version used.

Layout of the cache root::

    <root>/config.json         # {"pinned": "3.55.0"}
    <root>/3.55.0/             # a release directory
        Allelelist.txt
        Allelelist_history.txt
        hla_nom_g.txt
        hla_nom_p.txt

The root defaults to ``$HLAHARM_DATA_DIR`` then ``$XDG_DATA_HOME/hlaharm`` then
``~/.local/share/hlaharm``. A release directory named in either version form
(``3.55.0``) or IPD branch form (``3550``) is recognised, so a directory prepared
by hand or by ``git clone`` of ANHIG/IMGTHLA works without renaming.
"""

from __future__ import annotations

import json
import os
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path

from .nomenclature import ReferenceData

# IPD-IMGT/HLA reference files we consume, mapped to their path in an ANHIG/IMGTHLA
# release tree. Only the small name/group files — never the sequence files (§2.5).
RELEASE_FILES: dict[str, str] = {
    "Allelelist.txt": "Allelelist.txt",
    "Allelelist_history.txt": "Allelelist_history.txt",
    "hla_nom_g.txt": "wmda/hla_nom_g.txt",
    "hla_nom_p.txt": "wmda/hla_nom_p.txt",
}
# Files required for a directory to count as a usable release. History is optional:
# without it cross-version reconciliation degrades gracefully (§2.4) but validation
# and reduction still work.
REQUIRED_FILES = ("Allelelist.txt", "hla_nom_g.txt", "hla_nom_p.txt")

_RAW_BASE = "https://raw.githubusercontent.com/ANHIG/IMGTHLA/{branch}/{path}"
_CONFIG_NAME = "config.json"


class DatabaseError(RuntimeError):
    """A requested release is missing, malformed, or could not be fetched."""


# -- version <-> IPD branch conversion ---------------------------------------

def branch_of(version: str) -> str:
    """IPD branch/tag for a version string: ``3.55.0`` -> ``3550``.

    A string that is already in branch form (no dots) is returned unchanged.
    """
    if "." not in version:
        return version
    parts = version.split(".")
    major = int(parts[0])
    minor = int(parts[1]) if len(parts) > 1 else 0
    patch = int(parts[2]) if len(parts) > 2 else 0
    return f"{major}{minor:02d}{patch}"


def version_of(name: str) -> str:
    """Version string for a directory/branch name: ``3550`` -> ``3.55.0``.

    A name already containing dots is assumed to be a version string and returned
    unchanged; a 4+ digit branch is split major / 2-digit minor / remaining patch.
    """
    if "." in name or not name.isdigit() or len(name) < 4:
        return name
    return f"{int(name[0])}.{int(name[1:3])}.{int(name[3:])}"


# -- cache root & config ------------------------------------------------------

def default_root() -> Path:
    env = os.environ.get("HLAHARM_DATA_DIR")
    if env:
        return Path(env).expanduser()
    xdg = os.environ.get("XDG_DATA_HOME")
    base = Path(xdg).expanduser() if xdg else Path.home() / ".local" / "share"
    return base / "hlaharm"


def _config_path(root: Path) -> Path:
    return root / _CONFIG_NAME


def read_config(root: Path) -> dict:
    path = _config_path(root)
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return {}


def _write_config(root: Path, config: dict) -> None:
    root.mkdir(parents=True, exist_ok=True)
    _config_path(root).write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")


def pinned_version(root: Path) -> str | None:
    return read_config(root).get("pinned")


def pin(root: Path, version: str) -> str:
    """Record ``version`` as the default release. Requires it to be present."""
    if find_release(root, version) is None:
        raise DatabaseError(
            f"cannot pin {version!r}: not installed in {root}. "
            f"Run `hlaharm db update {version}` first."
        )
    config = read_config(root)
    config["pinned"] = version_of(version)
    _write_config(root, config)
    return config["pinned"]


# -- release discovery --------------------------------------------------------

def _is_release_dir(path: Path) -> bool:
    return path.is_dir() and all((path / f).exists() for f in REQUIRED_FILES)


def find_release(root: Path, version: str) -> Path | None:
    """Return the directory for ``version`` (accepting version or branch form)."""
    for name in _dedup(version, version_of(version), branch_of(version)):
        candidate = root / name
        if _is_release_dir(candidate):
            return candidate
    return None


def list_releases(root: Path) -> list[str]:
    """Installed release versions, newest first (by IPD branch ordering)."""
    if not root.exists():
        return []
    versions = {version_of(p.name) for p in root.iterdir() if _is_release_dir(p)}
    return sorted(versions, key=_version_sort_key, reverse=True)


def resolve_version(root: Path, requested: str | None) -> str:
    """Pick the release to use: explicit request, else pinned, else the only one.

    Raises :class:`DatabaseError` with actionable guidance when the choice is
    unavailable or ambiguous.
    """
    installed = list_releases(root)
    if requested is not None:
        if find_release(root, requested) is None:
            raise DatabaseError(
                f"release {requested!r} is not installed in {root}. "
                f"Available: {', '.join(installed) or 'none'}. "
                f"Run `hlaharm db update {requested}`."
            )
        return version_of(requested)
    pinned = pinned_version(root)
    if pinned and find_release(root, pinned) is not None:
        return version_of(pinned)
    if len(installed) == 1:
        return installed[0]
    if not installed:
        raise DatabaseError(
            f"no IPD-IMGT/HLA release installed in {root}. "
            "Run `hlaharm db update <version>` (e.g. 3.55.0) to fetch one."
        )
    raise DatabaseError(
        f"multiple releases installed ({', '.join(installed)}) and none pinned. "
        "Pass --db <version> or run `hlaharm db pin <version>`."
    )


def load_reference(root: Path, version: str) -> ReferenceData:
    release = find_release(root, version)
    if release is None:
        raise DatabaseError(f"release {version!r} is not installed in {root}")
    return ReferenceData.from_dir(release, version=version_of(version))


# -- fetching -----------------------------------------------------------------

@dataclass
class UpdateResult:
    version: str
    directory: Path
    files: list[str]
    skipped_history: bool


def _fetch_to(url: str, dest: Path) -> None:
    """Download ``url`` to ``dest`` (module-level so tests can monkeypatch it)."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(url) as response:
        tmp.write_bytes(response.read())
    tmp.replace(dest)


def update(
    root: Path,
    version: str,
    *,
    with_history: bool = True,
    set_pin: bool = True,
) -> UpdateResult:
    """Fetch a release into ``<root>/<version>/`` from ANHIG/IMGTHLA.

    ``with_history=False`` skips the large ``Allelelist_history.txt`` (~20 MB); the
    tool still validates and reduces but cannot reconcile cross-version renames.
    Pins the release as default unless ``set_pin`` is false.
    """
    branch = branch_of(version)
    version = version_of(version)
    target = root / version
    fetched: list[str] = []
    skipped_history = False
    for local_name, repo_path in RELEASE_FILES.items():
        if local_name == "Allelelist_history.txt" and not with_history:
            skipped_history = True
            continue
        url = _RAW_BASE.format(branch=branch, path=repo_path)
        try:
            _fetch_to(url, target / local_name)
        except (urllib.error.URLError, OSError) as exc:
            raise DatabaseError(f"failed to fetch {url}: {exc}") from exc
        fetched.append(local_name)
    if not _is_release_dir(target):
        raise DatabaseError(
            f"fetched files for {version} are incomplete in {target} "
            f"(need {', '.join(REQUIRED_FILES)})"
        )
    if set_pin:
        pin(root, version)
    return UpdateResult(
        version=version, directory=target, files=fetched, skipped_history=skipped_history
    )


# -- helpers ------------------------------------------------------------------

def _dedup(*items: str) -> list[str]:
    seen: dict[str, None] = {}
    for item in items:
        seen.setdefault(item, None)
    return list(seen)


def _version_sort_key(version: str) -> tuple[int, ...]:
    parts: list[int] = []
    for chunk in version.replace("-", ".").split("."):
        parts.append(int(chunk) if chunk.isdigit() else 0)
    return tuple(parts)
