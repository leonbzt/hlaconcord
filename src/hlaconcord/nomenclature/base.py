"""The nomenclature facade.

Every other module in hlaconcord talks to HLA nomenclature *only* through this
interface, never a third-party library directly. That makes the reducer
implementation swappable in one place: :class:`~hlaconcord.nomenclature.inhouse.InHouseNomenclature`
ships by default; a ``PyArdNomenclature`` backed by py-ard can be dropped in
behind the same interface for benchmarking or as a runtime fallback (see PLAN.md
§6). No caller needs to change.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from enum import Enum

from .models import Allele


class ReductionBasis(Enum):
    """The target basis a call is reduced to for comparison/consensus."""

    TWO_FIELD = "2field"  # plain field truncation to protein level
    G = "g"  # G-group: identical antigen-recognition-domain nucleotide sequence
    LGX = "lgx"  # 2-field representative of the G-group (G-group without the trailing g)
    P = "p"  # P-group: identical protein sequence in the ARD


# Default comparison basis for concordance. lgx (not g) because a 2-field call is
# ambiguous at G-group but collapses cleanly at lgx across every resolution the MVP
# tools emit — so a resolution difference never becomes a false discordance (§7, M1).
DEFAULT_BASIS = ReductionBasis.LGX


class ValidationStatus(Enum):
    """Outcome of validating a name against a specific IPD-IMGT/HLA release."""

    EXACT = "exact"  # a full allele name present in this release
    VALID_REDUCTION = "valid_reduction"  # a real field-prefix of >=1 allele
    G_GROUP = "g_group"  # a valid G-group designation
    P_GROUP = "p_group"  # a valid P-group designation
    RENAMED = "renamed"  # unknown here but present in allele history (renamed/deleted)
    UNKNOWN = "unknown"  # not recognised — flagged, never coerced


class Nomenclature(ABC):
    """Parse / validate / reduce / accession-resolve HLA allele designations."""

    @abstractmethod
    def parse(self, text: str) -> Allele:
        """Parse a raw allele string into an :class:`Allele`."""

    @abstractmethod
    def validate(self, allele: Allele) -> ValidationStatus:
        """Classify an allele against the reference release."""

    @abstractmethod
    def reduce(self, allele: Allele, basis: ReductionBasis) -> str:
        """Reduce an allele to the given comparison basis, returning its designation."""

    @abstractmethod
    def accession(self, allele: Allele) -> str | None:
        """Return the stable IPD-IMGT/HLA accession id for this name, if known."""
