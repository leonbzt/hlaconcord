# hlaharm

Harmonize and validate HLA typing output across multiple typers.

Ingest the raw output of N HLA typers (OptiType, arcasHLA, HLA-LA, HLA-HD, …),
normalize every call to canonical IPD-IMGT/HLA nomenclature at a chosen resolution,
validate each allele name against the official allele database, and emit one tidy
table plus a concordance/consensus view showing where the tools agree and disagree.

The field's own best practice is to run several typers and cross-check them; today
that means renaming alleles by hand and building concordance tables in a spreadsheet.
hlaharm makes that workflow reproducible — and does the one thing the manual version
gets wrong, reconciling **cross-database-version** name skew via stable accession ids
so a version difference never masquerades as a real disagreement.

See [`PLAN.md`](PLAN.md) for the full design and rationale.

## Status

**M0–M5 complete** — nomenclature core, all four MVP parsers, concordance/consensus
engine, CLI, and packaging/CI hardening. The engine
(parse → normalize → validate → reconcile → concord → emit) is implemented and tested
(113 tests + an optional py-ard oracle gate).

## Install

```bash
pip install -e ".[dev]"          # add ",oracle" to also run the py-ard cross-check
```

The IPD-IMGT/HLA reference database is **fetched, not bundled** (only the small
name/group files, never the sequences). Install a release before your first run:

```bash
hlaharm db update 3.55.0         # into a local cache (XDG data dir by default)
```

## Quickstart

```bash
hlaharm run \
  --inputs optitype:examples/s1/s1_result.tsv \
           arcasHLA:examples/s1/s1.genotype.json \
           hla-la:examples/s1/s1_bestguess_G.txt \
           hla-hd:examples/s1/s1_final.result.txt \
  --db 3.55.0 --gl -o out/
```

```
# sample s1 (basis=lgx)
  OK A      A*01:01 + A*02:01
  OK B      B*07:02
  OK C      C*07:02
  OK DRB1   DRB1*15:01
```

Four different formats/resolutions of the same genotype — 2-field, 3-field, G-group,
and `HLA-`-prefixed — reduce to one comparison key per allele and concord across loci.
See [`examples/`](examples/) for the full walkthrough.

## CLI

```bash
hlaharm run --inputs TOOL:PATH ... [--sample S] [--compare lgx|g|p|2field]
            [--consensus majority|unanimous] [--db VER] [--gl] [-o DIR]
hlaharm run --samplesheet samples.csv -o out/   # batch: sample,tool,path[,db_version]

hlaharm validate  HLA-A*02:01 A*99:99            # classify names against a release
hlaharm normalize A*0201 A*29:112N              # canonical form + accession + reduction

hlaharm db list | update <ver> | pin <ver> | path
```

`run` exits `0` when every locus concords, `1` if any locus is discordant, and `2` on
a configuration or database error — so it drops into a pipeline or `make` target.

Comparison basis (`--compare`) defaults to **`lgx`** — the 2-field representative of an
allele's G-group — because it collapses every resolution the MVP tools emit to one key
without turning a resolution difference into a false discordance (`g` is available as a
stricter ARD lens; see [`PLAN.md`](PLAN.md) §7).

## Outputs

With `-o DIR`, `run` writes:

- `tidy.tsv` — one row per allele call, raw string kept alongside the
  normalized / validated / reduced values and the accession (the auditable spine).
- `concordance.tsv` — one row per sample×locus: status, consensus, per-tool agreement,
  and flags (discordance, singleton, null allele, resolution conflict,
  version-skew-resolved, invalid call).
- `concordance.json` — the same, sample-grouped and metadata-stamped, for
  programmatic consumers; includes each sample's consensus GL string.
- `gl_strings.tsv` — consensus [GL String](https://glstring.org/) per sample (with `--gl`).

## Library

The CLI is a thin layer over an importable pipeline:

```python
from hlaharm import db
from hlaharm.nomenclature import InHouseNomenclature
from hlaharm.pipeline import InputSpec, run

root = db.default_root()
nom = InHouseNomenclature(db.load_reference(root, "3.55.0"))
result = run(
    [InputSpec("optitype", "s1_result.tsv"), InputSpec("arcasHLA", "s1.genotype.json")],
    nom, db_version="3.55.0",
)
for locus in result.concordance:
    print(locus.locus, locus.status.value, locus.consensus)
```

Package layout:

- `hlaharm.nomenclature` — allele model + parser, the `Nomenclature` facade, the
  in-house reducer/validator (`InHouseNomenclature`), and the reference loaders.
- `hlaharm.parsers` — OptiType, arcasHLA, HLA-LA, HLA-HD adapters behind a registry.
- `hlaharm.normalize` / `hlaharm.tidy` — enrich calls through the facade; tidy table.
- `hlaharm.concordance` / `hlaharm.report` / `hlaharm.gl` — per sample×locus
  concordance and consensus, TSV/JSON/human reports, and GL-String export.
- `hlaharm.db` — reference-release cache management.
- `hlaharm.pipeline` / `hlaharm.cli` — orchestration and the command-line surface.

## Nomenclature design note

The in-house reducer is the shipped implementation. `py-ard` (NMDP, LGPL-3.0) is used
only as a **test-time oracle** to cross-check it — never a runtime dependency, so no
LGPL code ships in the distributed artifact. The M0 spike verified agreement on real
IPD-IMGT/HLA 3.55.0 data (≥99.4% at `lgx`; residual divergences are documented and
allow-listed). See [`PLAN.md`](PLAN.md) §6, §12.

## Development

```bash
pip install -e ".[dev]"
ruff check src tests
pytest                       # 113 tests; add ".[oracle]" + HLAHARM_IMGT_DIR for the oracle gate
```
