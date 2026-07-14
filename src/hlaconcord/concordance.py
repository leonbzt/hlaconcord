"""Concordance / consensus engine (PLAN.md §7, pipeline stage [4]).

Consumes calls already normalized to a comparison basis (``AlleleCall.reduced``)
and, per sample+locus, decides whether the tools agree, builds a consensus
genotype under a configurable rule, and surfaces the flags a reviewer needs.

Design choices that keep it honest:

* **Unordered, zygosity-tolerant matching.** Allele-1/allele-2 order is arbitrary,
  and tools report homozygotes inconsistently (``A*02:01`` once vs twice). Agreement
  is decided on the *set* of distinct reduced keys, so a resolution/zygosity
  reporting quirk never fabricates discordance; the as-called multiset is retained
  for display.
* **Reconciliation before comparison.** Because keys are already version-reconciled
  and reduced, distinct raw/normalized names that collapse to one key are surfaced
  as ``VERSION_SKEW_RESOLVED`` rather than counted as disagreement.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from enum import Enum

from .model import GenotypeCall
from .nomenclature import DEFAULT_BASIS, ReductionBasis


class ConsensusRule(Enum):
    UNANIMOUS = "unanimous"  # allele must be called by every typing tool
    MAJORITY = "majority"  # allele must be called by more than half the typing tools


class ConcordanceStatus(Enum):
    CONCORDANT = "concordant"  # every typing tool reports the same allele set
    DISCORDANT = "discordant"  # typing tools disagree on the allele set
    SINGLE_TOOL = "single_tool"  # only one tool typed the locus — nothing to concord


class ConcordanceFlag(Enum):
    DISCORDANCE = "discordance"
    SINGLETON = "singleton"  # an allele only one tool called (with >=2 tools typing)
    NULL_ALLELE = "null_allele"  # a null (N) allele is present — clinically significant
    RESOLUTION_CONFLICT = "resolution_conflict"  # keys prefix-compatible but unequal
    VERSION_SKEW_RESOLVED = "version_skew_resolved"  # distinct names, one accession
    INVALID_CALL = "invalid_call"  # a tool emitted an unparseable/unknown allele here


@dataclass
class AlleleSupport:
    """Who called one reduced key at a locus, and how it was expressed."""

    key: str
    tools: list[str]
    normalized_names: set[str] = field(default_factory=set)
    resolutions: set[str] = field(default_factory=set)
    accession_names: dict[str, set[str]] = field(default_factory=dict)  # accession -> names
    is_null: bool = False

    @property
    def support(self) -> int:
        return len(self.tools)

    @property
    def reconciled_by_accession(self) -> bool:
        """True if >=2 distinct names were proven the same allele by a shared accession."""
        return any(len(names) >= 2 for names in self.accession_names.values())


@dataclass
class LocusConcordance:
    sample_id: str
    locus: str
    basis: str
    n_tools: int  # tools that produced >=1 valid call for this locus
    per_tool: dict[str, tuple[str, ...]]  # tool -> reduced keys as called (multiset)
    support: dict[str, AlleleSupport]  # reduced key -> who called it
    consensus: tuple[str, ...]  # consensus genotype under the rule
    status: ConcordanceStatus
    flags: frozenset[ConcordanceFlag]

    @property
    def concordant(self) -> bool:
        return self.status is ConcordanceStatus.CONCORDANT

    def agreement(self, key: str) -> str:
        """Human-facing support, e.g. ``3/4`` for a key called by 3 of 4 tools."""
        return f"{self.support[key].support}/{self.n_tools}"


def _has_resolution_conflict(support: dict[str, AlleleSupport]) -> bool:
    """True if one reduced key is a strict field-prefix of another.

    That means a lower-resolution call could be the same allele as a higher-res one
    but can't be confirmed equal — a resolution-driven ambiguity, distinct from clean
    discordance. Under the lgx default keys are uniform 2-field, so this stays silent;
    it fires for genuinely under-resolved input (e.g. a 1-field call).
    """
    keys = list(support)
    for i, a in enumerate(keys):
        for b in keys[i + 1 :]:
            lo, hi = sorted((a, b), key=len)
            if hi != lo and hi.startswith(lo + ":"):
                return True
    return False


def _consensus(
    support: dict[str, AlleleSupport], n_tools: int, rule: ConsensusRule
) -> tuple[str, ...]:
    if rule is ConsensusRule.UNANIMOUS:
        keep = [k for k, s in support.items() if s.support == n_tools]
    elif rule is ConsensusRule.MAJORITY:
        keep = [k for k, s in support.items() if s.support * 2 > n_tools]
    else:  # pragma: no cover - exhaustive
        raise ValueError(f"unknown consensus rule {rule!r}")
    return tuple(sorted(keep))


def _locus_concordance(
    sample_id: str,
    locus: str,
    calls: list[GenotypeCall],
    rule: ConsensusRule,
    basis: ReductionBasis,
) -> LocusConcordance:
    per_tool: dict[str, tuple[str, ...]] = {}
    per_tool_set: dict[str, frozenset[str]] = {}
    support: dict[str, AlleleSupport] = {}
    any_invalid = False

    for gc in calls:
        keys: list[str] = []
        for allele in gc.alleles:
            if allele.reduced is None:  # unparseable/unknown — flagged, not compared
                any_invalid = True
                continue
            keys.append(allele.reduced)
            info = support.setdefault(allele.reduced, AlleleSupport(key=allele.reduced, tools=[]))
            if gc.tool not in info.tools:
                info.tools.append(gc.tool)
            if allele.normalized:
                info.normalized_names.add(allele.normalized)
            if allele.resolution:
                info.resolutions.add(allele.resolution)
            if allele.accession:
                info.accession_names.setdefault(allele.accession, set()).add(
                    allele.normalized or allele.reduced
                )
            info.is_null = info.is_null or allele.is_null
        if keys:
            per_tool[gc.tool] = tuple(sorted(keys))
            per_tool_set[gc.tool] = frozenset(keys)

    n_tools = len(per_tool)
    for info in support.values():
        info.tools.sort()

    if n_tools <= 1:
        status = ConcordanceStatus.SINGLE_TOOL
    elif len(set(per_tool_set.values())) == 1:
        status = ConcordanceStatus.CONCORDANT
    else:
        status = ConcordanceStatus.DISCORDANT

    consensus = _consensus(support, n_tools, rule) if n_tools else ()

    flags: set[ConcordanceFlag] = set()
    if status is ConcordanceStatus.DISCORDANT:
        flags.add(ConcordanceFlag.DISCORDANCE)
    if n_tools >= 2 and any(s.support == 1 for s in support.values()):
        flags.add(ConcordanceFlag.SINGLETON)
    if any(s.is_null for s in support.values()):
        flags.add(ConcordanceFlag.NULL_ALLELE)
    if _has_resolution_conflict(support):
        flags.add(ConcordanceFlag.RESOLUTION_CONFLICT)
    # Genuine version skew: >=2 differently-named calls proven one allele by a shared
    # accession (§2.4). A depth difference (A*02:01 vs A*02:01:01), or one lone name that
    # happens to resolve while others don't, does NOT qualify.
    if any(s.reconciled_by_accession for s in support.values()):
        flags.add(ConcordanceFlag.VERSION_SKEW_RESOLVED)
    if any_invalid:
        flags.add(ConcordanceFlag.INVALID_CALL)

    return LocusConcordance(
        sample_id=sample_id,
        locus=locus,
        basis=basis.value,
        n_tools=n_tools,
        per_tool=per_tool,
        support=support,
        consensus=consensus,
        status=status,
        flags=frozenset(flags),
    )


def concordance(
    calls: list[GenotypeCall],
    rule: ConsensusRule = ConsensusRule.MAJORITY,
    basis: ReductionBasis = DEFAULT_BASIS,
) -> list[LocusConcordance]:
    """Group normalized calls by sample+locus and assess concordance per group.

    ``basis`` is metadata for the report only — reduction already happened in
    normalization; the engine compares ``AlleleCall.reduced`` verbatim.
    """
    grouped: dict[tuple[str, str], list[GenotypeCall]] = defaultdict(list)
    for gc in calls:
        grouped[(gc.sample_id, gc.locus)].append(gc)
    return [
        _locus_concordance(sample_id, locus, group, rule, basis)
        for (sample_id, locus), group in sorted(grouped.items())
    ]
