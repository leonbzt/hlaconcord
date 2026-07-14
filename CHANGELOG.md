# Changelog

All notable changes to hlaconcord are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/); this project uses semantic
versioning once it reaches 1.0.

## [Unreleased]

## [0.1.0] — 2026-07-15

First tagged release. Package `hlaconcord`, CLI command `hlacc`.

### Added
- **CLI (`hlacc`)** — `run`, `validate`, `normalize`, and `db` subcommands
  (M4). `run` ingests `--inputs tool:path …` or a `--samplesheet`, harmonizes
  across tools, and writes a tidy table, concordance report, JSON, and (with
  `--gl`) consensus GL strings. Exit code is `1` on any discordant locus, `2` on
  a configuration/database error.
- **Reference-database management** (`hlaconcord.db`) — on-demand fetch of an
  IPD-IMGT/HLA release into a local cache (`db update`), release discovery and
  pinning (`db list`, `db pin`, `db path`), and version↔IPD-branch conversion.
- **GL String export** (`hlaconcord.gl`) — standardized Genotype List strings from
  the harmonized consensus.
- **Pipeline entry point** (`hlaconcord.pipeline.run`) — importable parse →
  normalize → concord orchestration behind the CLI.
- **JSON report** (`hlaconcord.report.concordance_json`) — sample-grouped,
  metadata-stamped output for programmatic consumers.
- Packaging metadata, `ruff` lint configuration, an `examples/` validation set,
  and a GitHub Actions CI workflow (fast suite + optional py-ard oracle gate).

### Notes
- The IPD-IMGT/HLA database is **fetched, not redistributed** — `run` requires a
  release installed via `hlacc db update <version>` first.

## Milestones delivered before this changelog was started
- **M0** — nomenclature core (in-house reducer/validator + accession history),
  validated against the py-ard oracle on real release 3.55.0.
- **M1** — OptiType + arcasHLA parsers, call model, normalization pass, tidy table.
- **M2** — concordance/consensus engine + report.
- **M3** — HLA-LA + HLA-HD parsers; all four MVP typers harmonize on real data.

See `PLAN.md` §12 for the full milestone record and empirical results.
