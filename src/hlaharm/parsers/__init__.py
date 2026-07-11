"""Per-tool parser adapters. Each exposes ``parse(path) -> list[GenotypeCall]``.

Adding a typer means adding an adapter here — core logic is untouched (PLAN.md §5).
"""

from __future__ import annotations

from pathlib import Path

from ..model import GenotypeCall
from . import arcashla, hla_hd, hla_la, optitype

# Keyed by lower-cased tool name so callers can pass "OptiType", "arcasHLA", etc.
PARSERS = {
    optitype.TOOL: optitype.parse,
    arcashla.TOOL.lower(): arcashla.parse,
    hla_la.TOOL: hla_la.parse,
    hla_hd.TOOL: hla_hd.parse,
}

TOOLS = sorted(PARSERS)


def parse(tool: str, path: str | Path) -> list[GenotypeCall]:
    """Dispatch to the named tool's adapter."""
    try:
        adapter = PARSERS[tool.lower()]
    except KeyError:
        raise ValueError(f"unknown tool {tool!r}; known: {', '.join(TOOLS)}") from None
    return adapter(path)


__all__ = ["PARSERS", "TOOLS", "parse", "optitype", "arcashla", "hla_la", "hla_hd"]
