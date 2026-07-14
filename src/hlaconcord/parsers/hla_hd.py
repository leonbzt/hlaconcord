"""HLA-HD adapter.

HLA-HD writes ``<sample>_final.result.txt``: no header, one line per locus as
``Locus<tab>allele1<tab>allele2``. Alleles carry the ``HLA-`` prefix
(``HLA-A*01:01:01``), up to 3-field; absent alleles are ``-`` or ``Not typed``.
The ``HLA-`` prefix is handled downstream by the nomenclature parser, not here.
Parsing only (PLAN.md §3, §5).
"""

from __future__ import annotations

import csv
from pathlib import Path

from ..model import AlleleCall, GenotypeCall

TOOL = "hla-hd"
_SUFFIX = "_final.result.txt"
_EMPTY = {"", "-", "not typed", "not_typed", "couldn't read", "-,-"}


def _sample_id(path: Path) -> str:
    name = path.name
    return name[: -len(_SUFFIX)] if name.endswith(_SUFFIX) else path.stem


def parse(path: str | Path) -> list[GenotypeCall]:
    path = Path(path)
    sample_id = _sample_id(path)
    calls: list[GenotypeCall] = []
    with path.open(newline="", encoding="utf-8") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if not row or not row[0].strip():
                continue
            locus = row[0].strip()
            alleles = [
                AlleleCall(raw=value, tool=TOOL)
                for cell in row[1:]
                if (value := cell.strip()) and value.lower() not in _EMPTY
            ]
            if alleles:
                calls.append(
                    GenotypeCall(sample_id=sample_id, tool=TOOL, locus=locus, alleles=alleles)
                )
    return calls
