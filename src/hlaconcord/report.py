"""Concordance report emitters (pipeline stage [5], PLAN.md §11).

A tidy one-row-per-locus table for downstream analysis, plus a compact
human-readable summary for eyeballing where tools agree and disagree.
"""

from __future__ import annotations

import csv
import json
from collections.abc import Iterable, Mapping
from typing import TextIO

from .concordance import ConcordanceStatus, LocusConcordance

CONCORDANCE_COLUMNS = [
    "sample_id",
    "locus",
    "basis",
    "n_tools",
    "status",
    "consensus",
    "agreement",
    "flags",
    "per_tool",
]


def _consensus_str(result: LocusConcordance) -> str:
    return " + ".join(result.consensus)


def _agreement_str(result: LocusConcordance) -> str:
    return "; ".join(f"{key} {result.agreement(key)}" for key in result.consensus)


def _per_tool_str(result: LocusConcordance) -> str:
    return " | ".join(
        f"{tool}: {'/'.join(keys)}" for tool, keys in sorted(result.per_tool.items())
    )


def _flags_str(result: LocusConcordance) -> str:
    return ",".join(sorted(flag.value for flag in result.flags))


def concordance_rows(results: Iterable[LocusConcordance]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for result in results:
        rows.append(
            {
                "sample_id": result.sample_id,
                "locus": result.locus,
                "basis": result.basis,
                "n_tools": str(result.n_tools),
                "status": result.status.value,
                "consensus": _consensus_str(result),
                "agreement": _agreement_str(result),
                "flags": _flags_str(result),
                "per_tool": _per_tool_str(result),
            }
        )
    return rows


def write_concordance_tsv(results: Iterable[LocusConcordance], fh: TextIO) -> None:
    writer = csv.DictWriter(fh, fieldnames=CONCORDANCE_COLUMNS, delimiter="\t")
    writer.writeheader()
    writer.writerows(concordance_rows(results))


# -- JSON (programmatic consumers, PLAN.md §11) -------------------------------

def _locus_json(result: LocusConcordance) -> dict:
    return {
        "status": result.status.value,
        "n_tools": result.n_tools,
        "consensus": list(result.consensus),
        "agreement": {key: result.agreement(key) for key in result.consensus},
        "flags": sorted(flag.value for flag in result.flags),
        "per_tool": {tool: list(keys) for tool, keys in sorted(result.per_tool.items())},
    }


def concordance_json(
    results: Iterable[LocusConcordance],
    *,
    meta: Mapping[str, object] | None = None,
    gl_strings: Mapping[str, str] | None = None,
) -> dict:
    """Structured, sample-grouped concordance for programmatic consumers.

    ``meta`` (db version, basis, tool set, ...) is merged at the top level;
    ``gl_strings`` maps sample id -> consensus GL String when supplied.
    """
    samples: dict[str, dict] = {}
    for result in results:
        sample = samples.setdefault(result.sample_id, {"loci": {}})
        sample["loci"][result.locus] = _locus_json(result)
    if gl_strings:
        for sample_id, gl in gl_strings.items():
            if sample_id in samples:
                samples[sample_id]["gl_string"] = gl
    return {**dict(meta or {}), "samples": samples}


def write_concordance_json(
    results: Iterable[LocusConcordance],
    fh: TextIO,
    *,
    meta: Mapping[str, object] | None = None,
    gl_strings: Mapping[str, str] | None = None,
) -> None:
    results = list(results)
    json.dump(concordance_json(results, meta=meta, gl_strings=gl_strings), fh, indent=2)
    fh.write("\n")


_MARK = {
    ConcordanceStatus.CONCORDANT: "OK ",
    ConcordanceStatus.DISCORDANT: "XX ",
    ConcordanceStatus.SINGLE_TOOL: "-- ",
}


def format_concordance(results: Iterable[LocusConcordance]) -> str:
    """Compact per-locus summary, grouped by sample."""
    lines: list[str] = []
    current: str | None = None
    for result in results:
        if result.sample_id != current:
            current = result.sample_id
            lines.append(f"# sample {current} (basis={result.basis})")
        mark = _MARK[result.status]
        consensus = _consensus_str(result) or "(no consensus)"
        line = f"  {mark}{result.locus:6s} {consensus}"
        flags = _flags_str(result)
        if flags:
            line += f"   [{flags}]"
        lines.append(line)
    return "\n".join(lines)
