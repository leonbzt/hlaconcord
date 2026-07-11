"""Loaders for IPD-IMGT/HLA reference files (one release).

Consumes the same public files py-ard consumes, from an ANHIG/IMGTHLA release:

* ``Allelelist.txt``          -> accession id + full allele names (validation truth)
* ``hla_nom_g.txt``           -> G-group definitions
* ``hla_nom_p.txt``           -> P-group definitions
* ``Allelelist_history.txt``  -> accession -> name across releases (cross-version)

Only the small name/group files are needed for name validation and reduction;
sequence files (hla.dat, fastas) are out of scope.
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path


def _iter_data_lines(path: Path) -> Iterator[str]:
    with path.open(encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            yield line


@dataclass
class ReferenceData:
    """In-memory index of one IPD-IMGT/HLA release."""

    version: str
    alleles: set[str] = field(default_factory=set)  # full names, no HLA- prefix
    prefixes: set[str] = field(default_factory=set)  # valid field-truncations of alleles
    g_group_of: dict[str, str] = field(default_factory=dict)  # allele -> G-group name
    p_group_of: dict[str, str] = field(default_factory=dict)  # allele -> P-group name
    g_groups: set[str] = field(default_factory=set)
    p_groups: set[str] = field(default_factory=set)
    history: dict[str, dict[str, str]] = field(default_factory=dict)  # accession -> {ver: name}
    _name_to_accession: dict[str, str] = field(default_factory=dict)

    # -- validation helpers ---------------------------------------------------
    def is_allele(self, name: str) -> bool:
        return name in self.alleles

    def is_valid_reduction(self, name: str) -> bool:
        return name in self.prefixes

    def accession_of(self, name: str) -> str | None:
        return self._name_to_accession.get(name)

    # -- loading --------------------------------------------------------------
    @classmethod
    def from_dir(cls, data_dir: str | Path, version: str) -> ReferenceData:
        d = Path(data_dir)
        ref = cls(version=str(version))
        ref._load_allelelist(d / "Allelelist.txt")
        ref._load_groups(d / "hla_nom_g.txt", kind="g")
        ref._load_groups(d / "hla_nom_p.txt", kind="p")
        history = d / "Allelelist_history.txt"
        if history.exists():
            ref._load_history(history)
        return ref

    def _load_allelelist(self, path: Path) -> None:
        for line in _iter_data_lines(path):
            if line.lower().startswith("alleleid"):  # column header
                continue
            parts = line.split(",")
            if len(parts) < 2:
                continue
            accession, name = parts[0].strip(), parts[1].strip()
            self.alleles.add(name)
            self._name_to_accession.setdefault(name, accession)
            self._index_prefixes(name)

    def _index_prefixes(self, name: str) -> None:
        gene, _, digits = name.partition("*")
        core = digits
        if core and core[-1] in "NLSQCA":  # drop expression suffix before splitting
            core = core[:-1]
        segs = core.split(":")
        for i in range(1, len(segs)):
            self.prefixes.add(f"{gene}*{':'.join(segs[:i])}")

    def _load_groups(self, path: Path, *, kind: str) -> None:
        if not path.exists():
            return
        mapping = self.g_group_of if kind == "g" else self.p_group_of
        groups = self.g_groups if kind == "g" else self.p_groups
        for line in _iter_data_lines(path):
            # Format: "A*;01:01:01:01/01:01:01:02N/...;01:01:01G"
            parts = line.split(";")
            if len(parts) < 3:
                continue
            locus, alleles_str, group_str = parts[0].strip(), parts[1], parts[2].strip()
            if not group_str:
                continue  # singleton: allele forms its own group, no G/P designation
            group_name = f"{locus}{group_str}"
            groups.add(group_name)
            for allele in alleles_str.split("/"):
                allele = allele.strip()
                if allele:
                    mapping[f"{locus}{allele}"] = group_name

    def _load_history(self, path: Path) -> None:
        # Format: header "HLA_ID,3550,3540,..."; rows "HLA00005,A*02:01:01:01,A*02:01:01,..."
        header: list[str] | None = None
        for line in _iter_data_lines(path):
            parts = [p.strip() for p in line.split(",")]
            if header is None:
                header = parts
                continue
            accession = parts[0]
            per_version: dict[str, str] = {}
            # Rows may carry fewer columns than the header; stop at the shorter.
            for ver, value in zip(header[1:], parts[1:], strict=False):
                if value and value.upper() != "NA":
                    per_version[ver] = value
                    self._name_to_accession.setdefault(value, accession)
            self.history[accession] = per_version
