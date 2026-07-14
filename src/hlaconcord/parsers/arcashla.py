"""arcasHLA adapter.

arcasHLA writes ``<sample>.genotype.json``: a locus -> [alleles] map, up to
3-field, no ``HLA-`` prefix. Homozygous loci may list one allele or two. Class I
and II. Parsing only (PLAN.md §3).
"""

from __future__ import annotations

import json
from pathlib import Path

from ..model import AlleleCall, GenotypeCall

TOOL = "arcasHLA"
_SUFFIX = ".genotype.json"


def _sample_id(path: Path) -> str:
    name = path.name
    return name[: -len(_SUFFIX)] if name.endswith(_SUFFIX) else path.stem


def parse(path: str | Path) -> list[GenotypeCall]:
    path = Path(path)
    sample_id = _sample_id(path)
    data = json.loads(path.read_text(encoding="utf-8"))
    calls: list[GenotypeCall] = []
    for locus, alleles in data.items():
        parsed = [AlleleCall(raw=a.strip(), tool=TOOL) for a in alleles if a and a.strip()]
        if parsed:
            calls.append(GenotypeCall(sample_id=sample_id, tool=TOOL, locus=locus, alleles=parsed))
    return calls
