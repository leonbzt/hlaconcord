"""OptiType adapter.

OptiType writes ``<sample>_result.tsv`` with a leading (blank-header) index column
then A1,A2,B1,B2,C1,C2,Reads,Objective. Class I only, 2-field, no ``HLA-`` prefix.
One solution row per file. Parsing only — no normalization here (PLAN.md §3).
"""

from __future__ import annotations

import csv
from pathlib import Path

from ..model import AlleleCall, GenotypeCall

TOOL = "optitype"
_LOCUS_COLUMNS = {"A": ("A1", "A2"), "B": ("B1", "B2"), "C": ("C1", "C2")}
_EMPTY = {"", "-", "no", "not typed"}


def _sample_id(path: Path) -> str:
    name = path.name
    return name[: -len("_result.tsv")] if name.endswith("_result.tsv") else path.stem


def parse(path: str | Path) -> list[GenotypeCall]:
    path = Path(path)
    sample_id = _sample_id(path)
    with path.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        return []
    row = rows[0]  # OptiType emits a single optimal solution
    calls: list[GenotypeCall] = []
    for locus, columns in _LOCUS_COLUMNS.items():
        alleles = [
            AlleleCall(raw=value, tool=TOOL)
            for col in columns
            if (value := (row.get(col) or "").strip()) and value.lower() not in _EMPTY
        ]
        if alleles:
            calls.append(GenotypeCall(sample_id=sample_id, tool=TOOL, locus=locus, alleles=alleles))
    return calls
