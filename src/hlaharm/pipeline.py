"""End-to-end orchestration: parse -> normalize -> concord (PLAN.md §3).

The importable entry point behind the CLI. Given a set of ``(tool, path)`` inputs
and a nomenclature backend, it runs the whole pipeline and returns both the
enriched calls (for the tidy table) and the per-locus concordance (for the report),
so callers integrating hlaharm into a larger workflow never have to touch the CLI.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from . import parsers
from .concordance import ConsensusRule, LocusConcordance, concordance
from .model import GenotypeCall
from .nomenclature import DEFAULT_BASIS, Nomenclature, ReductionBasis
from .normalize import normalize


@dataclass
class InputSpec:
    """One typer output to ingest."""

    tool: str
    path: Path
    sample: str | None = None  # override the sample id parsed from the filename
    source_db_version: str | None = None  # IPD release the tool itself was built on

    @classmethod
    def parse_arg(cls, spec: str) -> InputSpec:
        """Parse a ``tool:path`` CLI argument (path may itself contain colons)."""
        tool, sep, path = spec.partition(":")
        if not sep or not tool or not path:
            raise ValueError(
                f"invalid input {spec!r}: expected 'tool:path' "
                f"(known tools: {', '.join(parsers.TOOLS)})"
            )
        return cls(tool=tool, path=Path(path))


@dataclass
class PipelineResult:
    calls: list[GenotypeCall]
    concordance: list[LocusConcordance]
    db_version: str
    basis: ReductionBasis
    rule: ConsensusRule
    samples: list[str] = field(default_factory=list)
    tools: list[str] = field(default_factory=list)


def run(
    inputs: list[InputSpec],
    nom: Nomenclature,
    *,
    db_version: str,
    basis: ReductionBasis = DEFAULT_BASIS,
    rule: ConsensusRule = ConsensusRule.MAJORITY,
) -> PipelineResult:
    """Parse every input, normalize through ``nom``, and compute concordance."""
    calls: list[GenotypeCall] = []
    for spec in inputs:
        parsed = parsers.parse(spec.tool, spec.path)
        for gc in parsed:
            if spec.sample is not None:
                gc.sample_id = spec.sample
            if spec.source_db_version is not None:
                gc.source_db_version = spec.source_db_version
        calls.extend(parsed)

    normalize(calls, nom, basis)
    results = concordance(calls, rule=rule, basis=basis)

    return PipelineResult(
        calls=calls,
        concordance=results,
        db_version=db_version,
        basis=basis,
        rule=rule,
        samples=sorted({gc.sample_id for gc in calls}),
        tools=sorted({gc.tool for gc in calls}),
    )
