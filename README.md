# hlaconcord

[![CI](https://github.com/leonbzt/hlaconcord/actions/workflows/ci.yml/badge.svg)](https://github.com/leonbzt/hlaconcord/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/hlaconcord.svg)](https://pypi.org/project/hlaconcord/)
[![Python](https://img.shields.io/pypi/pyversions/hlaconcord.svg)](https://pypi.org/project/hlaconcord/)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

**Harmonize and validate HLA typing output across multiple typers — and see, at a
glance, where they truly agree and where they don't.**

Running several HLA typers and cross-checking their calls is best practice, especially
for clinically-adjacent work — it's how you drive down the error rate. But every typer
speaks a slightly different dialect: different file formats, different resolutions
(2‑field vs. 3/4‑field vs. G‑group), and — the subtle one — they're built against
**different IPD‑IMGT/HLA database versions**, where allele *names* change over time.
So in practice, cross-checking means renaming alleles by hand and building concordance
tables in a spreadsheet.

`hlaconcord` does that for you, reproducibly:

1. **Parses** each typer's native output.
2. **Normalizes** every call to canonical IPD‑IMGT/HLA nomenclature (preserving the
   things that matter clinically — like the null‑allele `N` suffix — and never coercing
   an unrecognized name into a valid-looking one).
3. **Validates** each name against a pinned release of the official allele database.
4. **Reconciles** database-version skew: names that changed between releases are matched
   by their stable accession id, so a version difference never masquerades as a real
   disagreement.
5. **Reports** a per-locus concordance/consensus view with the flags a reviewer needs.

The result: heterogeneous typer outputs in, one validated, comparable answer out — with
genuine disagreements surfaced and formatting/version noise filtered away.

## Supported typers

| Typer | Input | Output parsed | Resolution |
|---|---|---|---|
| [OptiType](https://github.com/FRED-2/OptiType) | DNA/RNA FASTQ | `*_result.tsv` | 2‑field, class I (A/B/C) |
| [arcasHLA](https://github.com/RabadanLab/arcasHLA) | RNA‑seq BAM | `*.genotype.json` | up to 3‑field, class I + II |
| [HLA‑LA](https://github.com/DiltheyLab/HLA-LA) | WGS BAM/CRAM | `*_bestguess_G.txt` | G‑group, many loci |
| [HLA‑HD](https://www.genome.med.kyoto-u.ac.jp/HLA-HD/) | FASTQ | `*_final.result.txt` | up to 3‑field, many loci |

Adding a typer is an isolated parser — the harmonization core is untouched.

## Install

```bash
pip install hlaconcord
```

The IPD‑IMGT/HLA reference database is **fetched on demand, not bundled** (only the small
name/group files, never the sequences). Install a release once before your first run:

```bash
hlacc db update 3.55.0     # cached under your local data dir; recorded in every output
```

No runtime dependencies — the nomenclature core is in-house and standard-library only.
Requires Python 3.10+.

## A clear application

You ran four typers on one sample and want a single, trustworthy genotype — plus an
honest account of where the tools disagree. Point `hlaconcord` at the four raw files:

```bash
hlacc run \
  --inputs optitype:examples/s1/s1_result.tsv \
           arcasHLA:examples/s1/s1.genotype.json \
           hla-la:examples/s1/s1_bestguess_G.txt \
           hla-hd:examples/s1/s1_final.result.txt \
  --db 3.55.0 --gl -o out/
```

The four tools reported locus A in four *different* dialects:

```
optitype  A*02:01          (2-field)
arcasHLA  A*02:01:01        (3-field)
hla-la    A*02:01:01G       (G-group)
hla-hd    HLA-A*02:01:01    (HLA- prefixed)
```

`hlaconcord` reduces all four to one comparison key and reports agreement:

```
# sample s1 (basis=lgx)
  OK A      A*01:01 + A*02:01
  OK B      B*07:02
  OK C      C*07:02
  OK DRB1   DRB1*15:01

# GL s1
HLA-A*01:01+HLA-A*02:01^HLA-B*07:02^HLA-C*07:02^HLA-DRB1*15:01
```

Four formats, four resolutions → one concordant answer per locus, plus a standard
[GL String](https://glstring.org/) for downstream systems. (This example is in
[`examples/`](examples/) and is runnable end-to-end.)

### What it catches that a spreadsheet misses

The value shows up when the raw names *look* like they disagree. Consider a donor where
the tools' raw calls differ at two loci — here's an illustrative report:

```
# sample donor_07 (basis=lgx)
  OK   A      A*01:01 + A*02:01
  OK   B      B*08:01 + B*44:02     [version_skew_resolved]
  XX   DQB1   DQB1*03:01            [discordance]
```

- **B — reconciled, not a disagreement.** One typer was built on an older IPD‑IMGT/HLA
  release and emitted an allele name that has since been renamed. `hlaconcord` matched it
  to the newer name by its stable accession id and flagged it `version_skew_resolved` —
  the manual workflow would have logged a false discordance here.
- **DQB1 — a real discordance.** The tools genuinely disagree, and it's flagged `XX` so a
  human looks at exactly the locus that needs attention, not all of them.

`run` exits `0` when every locus concords, `1` if any locus is discordant, and `2` on a
configuration/database error — so it fits straight into a pipeline or `make` target.

## Outputs

With `-o DIR`, `run` writes:

- **`tidy.tsv`** — one row per allele call, the raw string kept beside the
  normalized / validated / reduced values and the accession. The auditable spine: you can
  always trace a designation back to what the tool actually emitted.
- **`concordance.tsv`** — one row per sample×locus: status, consensus genotype, per‑tool
  agreement (e.g. `A*02:01 4/4`), and flags: discordance, singleton (one‑tool‑only call),
  null allele, resolution conflict, version‑skew‑resolved, invalid call.
- **`concordance.json`** — the same, sample‑grouped and metadata‑stamped (db version,
  basis, tools), for programmatic consumers; includes each sample's consensus GL String.
- **`gl_strings.tsv`** — consensus GL String per sample (with `--gl`).

Without `-o`, `run` prints the human‑readable summary to stdout.

## Command-line reference

```bash
# harmonize a sample (or a batch)
hlacc run --inputs TOOL:PATH ... [--sample S] [--compare lgx|g|p|2field]
          [--consensus majority|unanimous] [--db VER] [--gl] [-o DIR]
hlacc run --samplesheet samples.csv -o out/    # columns: sample,tool,path[,db_version]

# quick nomenclature utilities
hlacc validate  HLA-A*02:01 A*99:99            # classify names against a release
hlacc normalize A*0201 A*29:112N               # canonical form + accession + reduction

# reference database management
hlacc db list | update <ver> | pin <ver> | path
```

**Comparison basis** (`--compare`, default **`lgx`**) is the 2‑field representative of an
allele's G‑group. It collapses every resolution the typers emit to one key, so a
resolution difference never becomes a false discordance; `g` is available as a stricter
antigen‑recognition‑domain lens. See [`PLAN.md`](PLAN.md) §7.

## Python API

The CLI is a thin layer over an importable pipeline:

```python
from hlaconcord import db
from hlaconcord.nomenclature import InHouseNomenclature
from hlaconcord.pipeline import InputSpec, run

nom = InHouseNomenclature(db.load_reference(db.default_root(), "3.55.0"))
result = run(
    [InputSpec("optitype", "s1_result.tsv"),
     InputSpec("arcasHLA", "s1.genotype.json")],
    nom, db_version="3.55.0",
)
for locus in result.concordance:
    print(locus.locus, locus.status.value, locus.consensus)
```

## How it stays correct

Getting HLA nomenclature subtly wrong in a clinically-adjacent tool is the failure mode
that matters, so correctness is designed in, not assumed:

- The shipped nomenclature reducer/validator is **in-house and built directly on the
  official IPD‑IMGT/HLA files**. [`py-ard`](https://github.com/nmdp-bioinformatics/py-ard)
  (NMDP) is used **only as a test-time oracle** to cross-check it — never a runtime
  dependency. Continuous integration reduces thousands of real alleles through both and
  asserts agreement (≥99.4% at `lgx`; the handful of documented divergences are
  allow-listed and understood).
- Normalization is **non‑lossy and non‑coercive**: the raw call is always retained, an
  unrecognized name is *flagged* rather than rewritten, and expression suffixes — including
  the clinically critical null `N` — are never silently dropped.

See [`PLAN.md`](PLAN.md) for the full design, the domain model, and the empirical results.

## Project status

Early but real: **v0.1.0 on PyPI**, four typers supported, the full
parse → normalize → validate → reconcile → concord → report pipeline implemented and
covered by 113 tests plus the py‑ard oracle gate in CI. Not yet clinically validated —
treat output as a decision aid, not a diagnostic. Feedback and typer requests welcome via
[issues](https://github.com/leonbzt/hlaconcord/issues).

## Development

```bash
pip install -e ".[dev]"      # add ",oracle" to also run the py-ard cross-check
ruff check src tests
pytest                       # 113 tests; set HLACONCORD_IMGT_DIR + install .[oracle] for the oracle gate
```

Releasing is documented in [`RELEASING.md`](RELEASING.md).

## License and data

`hlaconcord` is MIT-licensed. It reads and fetches reference data from the
[IPD‑IMGT/HLA database](https://www.ebi.ac.uk/ipd/imgt/hla/), which carries its own terms
and citation requirements — that data is neither bundled with nor covered by this license.
Please cite IPD‑IMGT/HLA if you use `hlaconcord` in published work.
