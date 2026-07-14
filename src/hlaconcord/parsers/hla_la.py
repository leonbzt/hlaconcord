"""HLA-LA adapter.

HLA-LA writes ``<sample>/hla/R1_bestguess_G.txt``: a header then one row **per
chromosome** — columns Locus, Chromosome, Allele, Q1 (+ coverage stats). Alleles
are in **G-group** notation (``A*02:01:01G``); many loci incl. non-classical.
Two rows per locus become one genotype. Q1 is captured as per-allele quality.
Parsing only (PLAN.md §3, §5).
"""

from __future__ import annotations

import csv
from pathlib import Path

from ..model import AlleleCall, GenotypeCall

TOOL = "hla-la"
_SUFFIXES = ("_bestguess_G.txt", "_bestguess.txt")
_GENERIC = {"R1", "R2", "READS"}  # HLA-LA names the file after the read group, not the sample
_EMPTY = {"", "?", "na", "not typed"}


def _sample_id(path: Path) -> str:
    name = path.name
    for suffix in _SUFFIXES:
        if name.endswith(suffix):
            base = name[: -len(suffix)]
            if base and base.upper() not in _GENERIC:
                return base
            # layout is <sample>/hla/R1_bestguess_G.txt — sample is the dir above "hla"
            for parent in path.parents:
                if parent.name and parent.name != "hla":
                    return parent.name
            return base or path.stem
    return path.stem


def _quality(value: str) -> float | None:
    try:
        return float(value)
    except ValueError:
        return None


def parse(path: str | Path) -> list[GenotypeCall]:
    path = Path(path)
    sample_id = _sample_id(path)
    by_locus: dict[str, list[AlleleCall]] = {}
    with path.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            locus = (row.get("Locus") or "").strip()
            allele = (row.get("Allele") or "").strip()
            if not locus or not allele or allele.lower() in _EMPTY:
                continue
            quality = _quality((row.get("Q1") or "").strip())
            by_locus.setdefault(locus, []).append(
                AlleleCall(raw=allele, tool=TOOL, quality=quality)
            )
    return [
        GenotypeCall(sample_id=sample_id, tool=TOOL, locus=locus, alleles=alleles)
        for locus, alleles in by_locus.items()
    ]
