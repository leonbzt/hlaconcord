"""Core HLA allele model and parser.

The parser is deliberately lenient about *reading* an allele string but never
coerces an unrecognised or malformed name into a valid-looking one: malformed
input raises :class:`ParseError`, and semantic validity (does this name exist in a
given IPD-IMGT/HLA release?) is decided later by the validator (see
``reference.py`` / ``inhouse.py``). This keeps the "never silently wrong" invariant
that matters for a clinically-adjacent tool.

Allele name anatomy (per IPD-IMGT/HLA)::

    HLA-A*02:01:01:01L
    └┬┘ └┬┘ └────┬────┘└┬┘
    pfx gene   fields   expression suffix

Fields carry decreasing significance: field 1 = allele group, field 2 = protein
(the clinical high-resolution level), field 3 = synonymous coding, field 4 =
non-coding. Expression suffixes ``N L S Q C A`` are significant and preserved; a
trailing ``G`` / ``P`` denotes a G-group / P-group designation, not an allele.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, replace
from enum import Enum

# Expression suffixes: Null, Low, Secreted, Questionable, Cytoplasm, Aberrant.
EXPRESSION_SUFFIXES = frozenset("NLSQCA")
# Group designation markers.
GROUP_SUFFIXES = frozenset("GP")

_ALLELE_RE = re.compile(
    r"^(?P<prefix>HLA-)?"
    r"(?P<gene>[A-Za-z0-9]+)"
    r"\*"
    r"(?P<digits>\d+(?::\d+)*)"
    r"(?P<suffix>[A-Za-z])?$"
)


class ParseError(ValueError):
    """Raised when a string is not a well-formed HLA allele designation."""


class Resolution(Enum):
    FIELD_1 = "1-field"
    FIELD_2 = "2-field"
    FIELD_3 = "3-field"
    FIELD_4 = "4-field"
    G_GROUP = "G-group"
    P_GROUP = "P-group"


_FIELD_RESOLUTION = {
    1: Resolution.FIELD_1,
    2: Resolution.FIELD_2,
    3: Resolution.FIELD_3,
    4: Resolution.FIELD_4,
}


@dataclass(frozen=True)
class Allele:
    """A parsed HLA allele designation.

    ``raw`` preserves the exact source string so a normalized value never has to
    stand in for the original in a record.
    """

    gene: str
    fields: tuple[str, ...]
    expression: str | None = None
    group: str | None = None
    raw: str | None = None
    had_prefix: bool = False
    legacy: bool = False

    @property
    def locus(self) -> str:
        return self.gene

    @property
    def is_null(self) -> bool:
        return self.expression == "N"

    @property
    def resolution(self) -> Resolution:
        if self.group == "G":
            return Resolution.G_GROUP
        if self.group == "P":
            return Resolution.P_GROUP
        return _FIELD_RESOLUTION[len(self.fields)]

    def name(self, *, prefix: bool = False) -> str:
        """Canonical designation, e.g. ``A*02:01:01:01L`` (or ``HLA-A*…`` if requested)."""
        core = f"{self.gene}*{':'.join(self.fields)}"
        if self.group:
            core += self.group
        elif self.expression:
            core += self.expression
        return f"HLA-{core}" if prefix else core

    def truncated(self, n_fields: int) -> Allele:
        """Return this allele reduced to its first ``n_fields`` fields.

        Drops any expression suffix and group marker: the result is a
        lower-resolution *reduction*, not a specific allele. Null-ness is not
        preserved here by design — callers that must not lose it read
        :attr:`is_null` on the original before reducing.
        """
        if not 1 <= n_fields <= len(self.fields):
            raise ValueError(
                f"cannot reduce {self.name()} ({len(self.fields)} fields) to {n_fields}"
            )
        return replace(self, fields=self.fields[:n_fields], expression=None, group=None)

    def __str__(self) -> str:  # pragma: no cover - trivial
        return self.name()


def _split_legacy(digits: str) -> tuple[str, ...]:
    """Split a colon-less legacy run (e.g. ``0201``) into 2-digit fields.

    Best-effort only: odd-length runs are genuinely ambiguous without the
    database and are refused rather than guessed.
    """
    if len(digits) % 2 != 0:
        raise ParseError(
            f"ambiguous legacy allele digits {digits!r}: odd length needs the "
            "IPD-IMGT/HLA database to resolve field boundaries"
        )
    return tuple(digits[i : i + 2] for i in range(0, len(digits), 2))


def parse_allele(text: str) -> Allele:
    """Parse an HLA allele string into an :class:`Allele`.

    Raises :class:`ParseError` on anything not well-formed.
    """
    s = text.strip()
    if not s:
        raise ParseError("empty allele string")
    m = _ALLELE_RE.match(s)
    if not m:
        raise ParseError(f"not a well-formed HLA allele: {text!r}")

    gene = m["gene"].upper()
    if gene == "HLA":
        raise ParseError(f"missing gene/locus in {text!r}")

    digits = m["digits"]
    if ":" in digits:
        fields = tuple(digits.split(":"))
        legacy = False
    elif len(digits) <= 3:
        # Modern low-resolution single field, e.g. A*02.
        fields = (digits,)
        legacy = False
    else:
        fields = _split_legacy(digits)
        legacy = True

    if not 1 <= len(fields) <= 4:
        raise ParseError(f"{text!r} has {len(fields)} fields (expected 1-4)")

    expression = group = None
    suffix = m["suffix"]
    if suffix:
        su = suffix.upper()
        if su in GROUP_SUFFIXES:
            group = su
        elif su in EXPRESSION_SUFFIXES:
            expression = su
        else:
            raise ParseError(f"unknown allele suffix {suffix!r} in {text!r}")

    return Allele(
        gene=gene,
        fields=fields,
        expression=expression,
        group=group,
        raw=text,
        had_prefix=bool(m["prefix"]),
        legacy=legacy,
    )
