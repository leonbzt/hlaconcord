"""GL String (Genotype List) export (PLAN.md §2.3, §11).

GL Strings are the standardized interchange format for HLA genotypes. We *emit*
them from the harmonized consensus; we do not accept them as input (the MVP typers
don't produce them). Delimiters, in increasing scope:

    /   allele ambiguity      (this position is one of several alleles)
    ~   in-phase haplotype
    +   gene copy / genotype  (the two copies at a locus)
    |   genotype ambiguity
    ^   locus separator

We use ``/`` (within a position), ``+`` (between the two gene copies), and ``^``
(between loci). Names carry the ``HLA-`` prefix, per the GL String convention.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence

from .concordance import LocusConcordance

# Canonical reporting order for HLA loci; anything unlisted sorts after, alphabetically.
_LOCUS_ORDER = [
    "A", "B", "C", "E", "F", "G",
    "DRA", "DRB1", "DRB3", "DRB4", "DRB5",
    "DQA1", "DQB1", "DPA1", "DPB1",
    "MICA", "MICB",
]


def _locus_sort_key(locus: str) -> tuple[int, str]:
    try:
        return (_LOCUS_ORDER.index(locus), "")
    except ValueError:
        return (len(_LOCUS_ORDER), locus)


def with_prefix(name: str) -> str:
    """Ensure the ``HLA-`` prefix GL Strings expect (idempotent)."""
    return name if name.startswith("HLA-") else f"HLA-{name}"


def locus_block(alleles: Sequence[str]) -> str:
    """One locus' gene copies joined with ``+``, e.g. ``HLA-A*01:01+HLA-A*02:01``.

    Each entry is emitted once in the order given (the consensus is already a sorted
    set of distinct keys, so a homozygous or single-copy call yields one term — its
    zygosity is genuinely unknown and is not fabricated as ``X+X``).
    """
    return "+".join(with_prefix(a) for a in alleles)


def genotype_string(loci: Iterable[tuple[str, Sequence[str]]]) -> str:
    """Build a GL String from ``(locus, alleles)`` pairs, canonical locus order.

    Loci with no alleles (no consensus) are omitted — a GL String asserts genotypes,
    so an empty locus is left unstated rather than encoded as a false call.
    """
    blocks = [
        locus_block(alleles)
        for locus, alleles in sorted(loci, key=lambda item: _locus_sort_key(item[0]))
        if alleles
    ]
    return "^".join(blocks)


def consensus_gl_string(results: Iterable[LocusConcordance]) -> str:
    """GL String for one sample's consensus across its loci.

    Pass the :class:`LocusConcordance` results for a single sample; the per-locus
    consensus genotype becomes each locus block.
    """
    return genotype_string((r.locus, r.consensus) for r in results)


def gl_strings_by_sample(results: Iterable[LocusConcordance]) -> dict[str, str]:
    """Map each sample id to its consensus GL String."""
    by_sample: dict[str, list[LocusConcordance]] = {}
    for result in results:
        by_sample.setdefault(result.sample_id, []).append(result)
    return {sample: consensus_gl_string(items) for sample, items in by_sample.items()}
