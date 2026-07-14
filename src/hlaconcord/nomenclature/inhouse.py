"""In-house implementation of the nomenclature facade.

Built directly on the IPD-IMGT/HLA reference files (see :mod:`reference`). This is
the shipped implementation; py-ard is used only as a test-time oracle to check it
(PLAN.md §6).
"""

from __future__ import annotations

from .base import Nomenclature, ReductionBasis, ValidationStatus
from .models import Allele, parse_allele
from .reference import ReferenceData


class InHouseNomenclature(Nomenclature):
    def __init__(self, reference: ReferenceData):
        self.ref = reference

    # -- facade ---------------------------------------------------------------
    def parse(self, text: str) -> Allele:
        return parse_allele(text)

    def validate(self, allele: Allele) -> ValidationStatus:
        name = allele.name()
        if allele.group == "G":
            return (
                ValidationStatus.G_GROUP if name in self.ref.g_groups else ValidationStatus.UNKNOWN
            )
        if allele.group == "P":
            return (
                ValidationStatus.P_GROUP if name in self.ref.p_groups else ValidationStatus.UNKNOWN
            )
        if self.ref.is_allele(name):
            return ValidationStatus.EXACT
        if self.ref.is_valid_reduction(name):
            return ValidationStatus.VALID_REDUCTION
        if self.ref.accession_of(name) is not None:
            return ValidationStatus.RENAMED
        return ValidationStatus.UNKNOWN

    def accession(self, allele: Allele) -> str | None:
        return self.ref.accession_of(allele.name())

    def reduce(self, allele: Allele, basis: ReductionBasis) -> str:
        if basis is ReductionBasis.TWO_FIELD:
            return self._two_field_name(allele)
        if basis is ReductionBasis.G:
            return self._group(allele, self.ref.g_group_of, marker="G")
        if basis is ReductionBasis.P:
            return self._group(allele, self.ref.p_group_of, marker="P")
        if basis is ReductionBasis.LGX:
            g = self._group(allele, self.ref.g_group_of, marker="G")
            return self._lgx_from_group(g)
        raise ValueError(f"unknown reduction basis {basis!r}")

    # -- reduction internals --------------------------------------------------
    def _two_field_name(self, allele: Allele) -> str:
        """Truncate to two fields, *preserving* any expression suffix.

        Dropping the suffix here silently discards null/low-expression status
        (e.g. ``A*02:11N``), the exact failure the brief forbids, so a null
        allele stays null even at 2-field resolution.
        """
        core = ":".join(allele.fields[: min(2, len(allele.fields))])
        return f"{allele.gene}*{core}{allele.expression or ''}"

    def _group(self, allele: Allele, mapping: dict[str, str], *, marker: str) -> str:
        """Reduce to a G- or P-group.

        Order of preference: the allele is already that group -> itself; an exact
        map hit on the full name; a truncation that maps unambiguously to one
        group. If IPD-IMGT/HLA assigns no group, the allele is its own ARD-level
        representative and is returned *unchanged* -- we never invent a broader
        grouping (that would over-merge distinct ARD sequences) and never drop an
        expression suffix. This matches py-ard's behaviour for ungrouped alleles.
        """
        if allele.group == marker:
            return allele.name()
        name = allele.name()
        group = mapping.get(name) or self._group_via_prefix(name, mapping)
        return group if group else name

    @staticmethod
    def _group_via_prefix(name: str, mapping: dict[str, str]) -> str | None:
        """Return the group for a truncated name iff every fuller allele agrees on one."""
        prefix = name + ":"
        found = {g for full, g in mapping.items() if full == name or full.startswith(prefix)}
        return next(iter(found)) if len(found) == 1 else None

    def _lgx_from_group(self, group_or_allele: str) -> str:
        # A named G-group ("A*02:01:01G") -> its 2-field stem; the group subsumes
        # expression, so no suffix. An ungrouped allele ("A*29:112N") -> 2-field
        # with its expression suffix preserved.
        if group_or_allele and group_or_allele[-1] in "GP":
            gene, _, digits = group_or_allele[:-1].partition("*")
            return f"{gene}*{':'.join(digits.split(':')[:2])}"
        return self._two_field_name(parse_allele(group_or_allele))
