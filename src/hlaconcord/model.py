"""Call-level data model — the objects that flow through the pipeline (PLAN.md §4).

Distinct from :class:`hlaconcord.nomenclature.Allele`, which is the parsed *name*
value object. An :class:`AlleleCall` wraps one raw string emitted by a typer plus
everything normalization derives from it, always keeping ``raw`` alongside.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .nomenclature import Allele, ValidationStatus


@dataclass
class AlleleCall:
    """One allele as a tool reported it, enriched in place by normalization."""

    raw: str
    tool: str
    quality: float | None = None  # per-allele confidence when the tool provides it (HLA-LA)
    # Filled by normalization (see hlaconcord.normalize); None until then.
    allele: Allele | None = None
    normalized: str | None = None  # canonical name, no HLA- prefix
    resolution: str | None = None
    validation: ValidationStatus | None = None
    reduced: str | None = None  # designation at the chosen comparison basis
    accession: str | None = None
    is_null: bool = False
    parse_error: str | None = None  # set (raw kept) when the string will not parse


@dataclass
class GenotypeCall:
    """A tool's call for one locus of one sample: usually two alleles."""

    sample_id: str
    tool: str
    locus: str  # gene symbol, e.g. "A", "DRB1"
    alleles: list[AlleleCall] = field(default_factory=list)
    source_db_version: str | None = None  # IPD release the tool was built on, if known
