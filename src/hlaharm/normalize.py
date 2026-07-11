"""Normalization + validation pass (pipeline stages [2]-[3], PLAN.md §3).

Takes raw :class:`GenotypeCall`s from the parsers and enriches every
:class:`AlleleCall` through the :class:`Nomenclature` facade — never a
nomenclature library directly. Non-coercive: an unparseable/unknown string is
*flagged*, and ``raw`` is always retained.
"""

from __future__ import annotations

from collections.abc import Iterable

from .model import AlleleCall, GenotypeCall
from .nomenclature import Nomenclature, ParseError, ReductionBasis


def normalize_allele(call: AlleleCall, nom: Nomenclature, basis: ReductionBasis) -> AlleleCall:
    try:
        allele = nom.parse(call.raw)
    except ParseError as exc:
        call.parse_error = str(exc)  # flagged, never coerced; raw kept
        return call
    call.allele = allele
    call.normalized = allele.name()
    call.resolution = allele.resolution.name
    call.is_null = allele.is_null
    call.validation = nom.validate(allele)
    call.reduced = nom.reduce(allele, basis)
    call.accession = nom.accession(allele)
    return call


def normalize(
    calls: Iterable[GenotypeCall], nom: Nomenclature, basis: ReductionBasis
) -> list[GenotypeCall]:
    """Enrich every allele in every call in place; returns the same calls."""
    calls = list(calls)
    for genotype in calls:
        for allele in genotype.alleles:
            normalize_allele(allele, nom, basis)
    return calls
