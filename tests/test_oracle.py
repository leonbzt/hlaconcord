"""M0 oracle gate: in-house reducer vs py-ard on real IPD-IMGT/HLA data.

py-ard (NMDP, LGPL-3.0) is a *test-only* oracle -- never a runtime dependency
(PLAN.md §6). The M0 spike (PLAN.md §12) established what "agreement" means here,
because py-ard is not a pure ground truth:

* **G basis (our default).** In-house reproduces py-ard exactly *except* three
  G-groups where py-ard relabels the group to a different member. Membership is
  identical; only the group's name differs, and in-house emits the name printed
  in IPD-IMGT/HLA's own ``hla_nom_g.txt`` -- the more IMGT-faithful choice. These
  are allow-listed below; anything else is a regression.

* **lgx basis.** py-ard additionally collapses ~0.5% of alleles that IPD assigns
  to *no* published G-group, using exon-2/3 sequence identity we deliberately do
  not ingest (PLAN.md §6, sequence files out of scope). Those are rare alleles
  the MVP typers do not emit. We accept and *count* them, but assert no other
  lgx disagreement appears (i.e. none where in-house did assign a G-group).

Skipped unless both py-ard and a release are available::

    pip install -e ".[oracle]"
    export HLACONCORD_IMGT_DIR=/path/to/IMGTHLA/release   # Allelelist.txt, hla_nom_g.txt, ...
    export HLACONCORD_IMGT_VERSION=3550                    # py-ard version string
    pytest tests/test_oracle.py -s
"""

from __future__ import annotations

import os
import random

import pytest

from hlaconcord.nomenclature import (
    InHouseNomenclature,
    ReductionBasis,
    ReferenceData,
    parse_allele,
)

pyard = pytest.importorskip("pyard", reason="oracle extra not installed")

IMGT_DIR = os.environ.get("HLACONCORD_IMGT_DIR")
IMGT_VERSION = os.environ.get("HLACONCORD_IMGT_VERSION")

pytestmark = pytest.mark.skipif(
    not (IMGT_DIR and IMGT_VERSION),
    reason="set HLACONCORD_IMGT_DIR and HLACONCORD_IMGT_VERSION to run the oracle gate",
)

SAMPLE_SIZE = int(os.environ.get("HLACONCORD_ORACLE_SAMPLE", "5000"))
MVP_LOCI = ("A*", "B*", "C*", "DRB1*", "DQB1*")

# G-group name pairs (in-house == IMGT wmda file, py-ard == its relabel) that share
# identical membership. Verified against release 3.55.0.
KNOWN_G_LABEL_DELTAS = frozenset({
    ("A*02:17:01G", "A*02:17:02G"),
    ("C*02:10:01G", "C*02:02:37G"),
    ("DRB3*03:22:01G", "DRB3*02:171:01G"),
})


def _two_field_of_group(group_name: str) -> str:
    gene, _, digits = group_name[:-1].partition("*")
    return f"{gene}*{':'.join(digits.split(':')[:2])}"


# Same relabels expressed at lgx (2-field) resolution, keeping only pairs that differ.
KNOWN_LGX_DELTAS = frozenset(
    pair
    for ours, theirs in KNOWN_G_LABEL_DELTAS
    if (pair := (_two_field_of_group(ours), _two_field_of_group(theirs)))[0] != pair[1]
)


@pytest.fixture(scope="module")
def oracle():
    return pyard.init(IMGT_VERSION)


@pytest.fixture(scope="module")
def inhouse():
    ref = ReferenceData.from_dir(IMGT_DIR, version=IMGT_VERSION)
    return InHouseNomenclature(ref)


@pytest.fixture(scope="module")
def sample(inhouse):
    # sort first: ref.alleles is a set and CPython randomizes str hashing per process,
    # so sorting is what makes the seeded draw reproducible across runs (a CI gate must be).
    names = sorted(n for n in inhouse.ref.alleles if n.startswith(MVP_LOCI))
    random.seed(0)
    return random.sample(names, min(SAMPLE_SIZE, len(names)))


def _reduce_pairs(inhouse, oracle, sample, basis, pyard_mode):
    """Yield (name, ours, theirs) for every allele py-ard could reduce."""
    for name in sample:
        allele = inhouse.parse(name)
        ours = inhouse.reduce(allele, basis)
        try:
            theirs = oracle.redux(name, pyard_mode)
        except Exception:  # py-ard raises on a few alleles it cannot reduce; skip them
            continue
        yield name, ours, theirs


def test_g_basis_matches_pyard_modulo_known_label_deltas(inhouse, oracle, sample):
    unexpected = []
    label_deltas = 0
    for name, ours, theirs in _reduce_pairs(inhouse, oracle, sample, ReductionBasis.G, "G"):
        if ours == theirs:
            continue
        if (ours, theirs) in KNOWN_G_LABEL_DELTAS:
            label_deltas += 1
        else:
            unexpected.append((name, ours, theirs))
    print(f"\n[G] known label deltas: {label_deltas}; unexpected: {len(unexpected)}")
    for name, ours, theirs in unexpected[:25]:
        print(f"  {name}: in-house={ours!r} pyard={theirs!r}")
    assert not unexpected, f"{len(unexpected)} unexpected G-basis disagreements"


def test_lgx_divergence_is_bounded_and_understood(inhouse, oracle, sample):
    unexpected = []
    label_deltas = 0
    ars_extra = 0  # py-ard reduced an allele IPD assigns to no G-group (sequence-based)
    for name, ours, theirs in _reduce_pairs(inhouse, oracle, sample, ReductionBasis.LGX, "lgx"):
        if ours == theirs:
            continue
        if (ours, theirs) in KNOWN_LGX_DELTAS:
            label_deltas += 1
        elif not inhouse.reduce(inhouse.parse(name), ReductionBasis.G).endswith("G"):
            # in-house found no published G-group; py-ard's stronger collapse is the
            # accepted, out-of-scope ARS capability gap.
            ars_extra += 1
        else:
            unexpected.append((name, ours, theirs))
    total = sum(1 for _ in sample)
    print(
        f"\n[lgx] label deltas: {label_deltas}; ars-extra (accepted): {ars_extra} "
        f"({ars_extra / max(total, 1):.2%}); unexpected: {len(unexpected)}"
    )
    for name, ours, theirs in unexpected[:25]:
        print(f"  {name}: in-house={ours!r} pyard={theirs!r}")
    assert not unexpected, f"{len(unexpected)} unexpected lgx disagreements"


def _truncate(name: str, n_fields: int) -> str | None:
    """Clean field-truncation (drops any suffix), as OptiType/arcasHLA emit."""
    allele = parse_allele(name)
    if len(allele.fields) < n_fields:
        return None
    return f"{allele.gene}*{':'.join(allele.fields[:n_fields])}"


@pytest.mark.parametrize("n_fields", [2, 3])
def test_lgx_on_truncated_inputs_matches_pyard(inhouse, oracle, sample, n_fields):
    """lgx is our default basis; MVP tools emit 2-/3-field calls, not full names.

    A 2-field call is ambiguous at G-group but must collapse cleanly at lgx (§7),
    so this guards the OptiType/HLA-HD path the full-name sampling above misses.
    """
    inputs = {t for name in sample if (t := _truncate(name, n_fields))}
    unexpected = []
    label_deltas = ars_extra = 0
    for text in inputs:
        allele = inhouse.parse(text)
        ours = inhouse.reduce(allele, ReductionBasis.LGX)
        try:
            theirs = oracle.redux(text, "lgx")
        except Exception:
            continue
        if ours == theirs:
            continue
        if (ours, theirs) in KNOWN_LGX_DELTAS:
            label_deltas += 1
        elif not inhouse.reduce(allele, ReductionBasis.G).endswith("G"):
            ars_extra += 1
        else:
            unexpected.append((text, ours, theirs))
    print(
        f"\n[lgx/{n_fields}-field] inputs: {len(inputs)}; label deltas: {label_deltas}; "
        f"ars-extra (accepted): {ars_extra}; unexpected: {len(unexpected)}"
    )
    for text, ours, theirs in unexpected[:25]:
        print(f"  {text}: in-house={ours!r} pyard={theirs!r}")
    assert not unexpected, f"{len(unexpected)} unexpected {n_fields}-field lgx disagreements"
