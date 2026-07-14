"""Tidy long-table output (pipeline stage [5], PLAN.md §11).

One row per allele call — the auditable spine everything else builds on. Each row
carries both the raw string and the normalized/validated/reduced values, so a
reviewer can always trace a designation back to what the tool actually emitted.
"""

from __future__ import annotations

import csv
from collections.abc import Iterable
from typing import TextIO

from .model import GenotypeCall

COLUMNS = [
    "sample_id",
    "tool",
    "locus",
    "raw",
    "normalized",
    "resolution",
    "reduced",
    "validation",
    "accession",
    "is_null",
    "quality",
    "source_db_version",
    "parse_error",
]


def tidy_rows(calls: Iterable[GenotypeCall]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for genotype in calls:
        for allele in genotype.alleles:
            rows.append(
                {
                    "sample_id": genotype.sample_id,
                    "tool": genotype.tool,
                    "locus": genotype.locus,
                    "raw": allele.raw,
                    "normalized": allele.normalized or "",
                    "resolution": allele.resolution or "",
                    "reduced": allele.reduced or "",
                    "validation": allele.validation.value if allele.validation else "",
                    "accession": allele.accession or "",
                    "is_null": "true" if allele.is_null else "false",
                    "quality": "" if allele.quality is None else f"{allele.quality:g}",
                    "source_db_version": genotype.source_db_version or "",
                    "parse_error": allele.parse_error or "",
                }
            )
    return rows


def write_tidy_tsv(calls: Iterable[GenotypeCall], fh: TextIO) -> None:
    writer = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t")
    writer.writeheader()
    writer.writerows(tidy_rows(calls))
