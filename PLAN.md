# hlaconcord — HLA typing output harmonizer + validator

**Status:** planning draft · **Date:** 2026-07-02

A Python tool that ingests the raw outputs of multiple HLA typers, normalizes every
call to canonical IPD-IMGT/HLA nomenclature at a chosen resolution, validates each
allele name against the official allele database, and emits one tidy table plus a
concordance/consensus view showing where the tools agree and disagree.

Inspector + validator + converter, combined, for one problem: making the field's own
recommended "run several typers and cross-check" workflow reproducible instead of hand-built.

---

## 0. Decisions assumed for this draft

The following were not confirmed and use recommended defaults. Flagged here so they are
easy to revisit:

| Decision | Assumed default | Alternatives to reconsider |
|---|---|---|
| MVP typers | OptiType, arcasHLA, HLA-LA, HLA-HD | All 8+ up front; or minimal 2-tool proof |
| Reference DB | Bundle one pinned IPD-IMGT/HLA release + `update` command | Fetch-at-runtime; or version-selectable cache |
| Interface | CLI-first, importable library underneath | Library-first; or both equal from day one |
| Nomenclature core | **In-house slim reducer/validator; `py-ard` as a test-only oracle** — see §6 | Adopt `py-ard` at runtime (behind the facade) if the M0 spike finds bad edge cases |

---

## 1. Why this is the well-founded gap

- Best practice in the field is to run multiple typers and compare, especially for
  clinical use, to drive down the error rate.
- But tools take different inputs (FASTQ / BAM / CRAM), emit different formats and
  resolutions (2-field vs 3/4-field vs G-group), and — the real problem — are built
  against **different IPD-IMGT/HLA database versions** that update at different rates.
- So today, cross-validation means manually renaming alleles and hand-building
  concordance tables. No standard harmonizer does this.
- Not tied to a single vendor's moving target → lower churn than a vendor-format chaser.

The novelty is the **multi-tool parsing + cross-version reconciliation + concordance**
layer. It is explicitly *not* to reinvent HLA nomenclature reduction — that logic is
subtle, clinically load-bearing, and already implemented in maintained libraries we
should evaluate leaning on (§6).

---

## 2. The domain model (this is the crux — get it exactly right)

Getting HLA nomenclature subtly wrong in a clinically-adjacent tool is the bad failure
mode. The whole design is organized around never doing that silently.

### 2.1 Allele name structure

```
HLA-A*02:01:01:01L
└┬┘ └┬┘ └─────┬─────┘└┬┘
gene locus   fields   expression suffix
```

- **Fields** are colon-separated and carry decreasing biological significance:
  - Field 1 — **allele group** (often maps to a serological antigen), e.g. `A*02`
  - Field 2 — **specific protein**: differs by ≥1 non-synonymous (amino-acid) change,
    e.g. `A*02:01`. **This 2-field level is the clinical high-resolution standard.**
  - Field 3 — synonymous (silent) coding differences, e.g. `A*02:01:01`
  - Field 4 — non-coding differences (introns/UTRs), e.g. `A*02:01:01:01`
- **Expression suffixes** (must be preserved, never dropped):
  `N` null (not expressed), `L` low surface expression, `S` secreted/soluble,
  `Q` questionable expression, `C` cytoplasm-only, `A` aberrant. **`N` (null) is
  clinically critical** — a silently dropped `N` is a serious defect.
- **Legacy nomenclature** (pre-2010) used no colons: `A*0201`. Some tools/DBs still
  emit or embed these; the normalizer must convert them correctly.

### 2.2 Resolution levels and ARD-equivalence groups

- 1-field (allele group) → 2-field (protein, clinical hi-res) → 3-field → 4-field.
- **G groups** (`A*02:01:01G`): alleles with identical sequence across the exons
  encoding the antigen-recognition domain (exons 2+3 for class I, exon 2 for class II).
  Many NGS typers can only resolve to ARD level, so they report G-groups.
- **P groups** (`A*02:01P`): alleles with identical *protein* sequence in the ARD.
- G/P group memberships come from IPD-IMGT/HLA files `hla_nom_g.txt` / `hla_nom_p.txt`.
  Concordance across tools frequently has to happen at G- or P-group level, not just by
  truncating fields.

### 2.3 Ambiguity encodings the tool should understand

- **MAC / NMDP codes** (Multiple Allele Codes), e.g. `A*02:AB` — a compressed code for a
  set of alleles sharing the first field. Expansion table is maintained by NMDP.
- **GL Strings** (Genotype List strings, standardized): the interchange format we should
  be able to *emit*. Delimiters by increasing scope:
  `/` allele ambiguity · `~` in-phase haplotype · `+` gene-copy/genotype ·
  `|` genotype ambiguity · `^` locus separator.

### 2.4 The cross-version problem, and the key insight

Allele **names can change between IPD-IMGT/HLA releases** (renames, splits, deletions).
The stable identifier across releases is the **accession number** (e.g. `HLA00005`), not
the name. `Allelelist_history.txt` maps accession → name per release.

**Design consequence:** to compare calls from tools built on different DB versions,
resolve each name to its accession (where possible), then project to the target version's
name *before* declaring agreement or disagreement. This is what prevents a version skew
from masquerading as a real discordance — the exact failure the manual workflow suffers.

### 2.5 Reference files we consume (from the ANHIG/IMGTHLA release)

- `Allelelist.txt` — canonical allele IDs + names for a release (validation source of truth).
- `Allelelist_history.txt` — accession → name across releases (cross-version reconciliation).
- `hla_nom_g.txt`, `hla_nom_p.txt` — G/P group definitions.
- Release/version string like `3.55.0` (roughly quarterly).
- Sequence files (`hla.dat`, fastas) are **not needed** for name validation → out of MVP scope.

---

## 3. Pipeline / architecture

```
inputs (per tool)
    │
    ▼
[1] Ingest / parse        per-tool adapters → raw GenotypeCall objects
    │
    ▼
[2] Normalize             legacy→canonical, prefix handling, suffix preservation,
    │                     resolution tagging, accession resolution
    ▼
[3] Validate              against chosen DB version: real allele? valid reduction?
    │                     valid G/P group? deleted/renamed? → status per allele
    ▼
[4] Reconcile + concord   project to common DB version + comparison resolution;
    │                     set-match genotypes per locus; consensus + flags
    ▼
[5] Emit                  tidy long table · concordance/consensus report · JSON · GL strings
```

Each stage is independently testable. Parsing is deliberately separated from
normalization so a new typer only adds an adapter, not core logic.

The nomenclature logic in stages [2]–[3] (parse, reduce, validate, accession-resolve)
lives behind a single **`hlaconcord.nomenclature` facade** (§6). Every other module calls that
interface and never a third-party library directly, so the reducer implementation can be
swapped (in-house ↔ py-ard) in one file without touching the parsers or the concordance
engine.

---

## 4. Internal data model

- **`Allele`**: `raw`, `gene`, `fields` (tuple), `suffix`, `source_tool`,
  `source_db_version`, `resolution`, `accession?`, `normalized_name`,
  `validation_status`, `g_group?`, `p_group?`. Always keeps `raw` alongside normalized.
- **`GenotypeCall`**: `locus`, `alleles` (usually 2; handle homozygous & "not typed"),
  `tool`, `quality?`/`confidence?` when the tool provides it.
- **`Sample`**: `id`, mapping `tool → [GenotypeCall]`.
- **`ConcordanceResult`**: per locus — tool×allele matrix, consensus genotype,
  agreement stats (e.g. 3/4), and flags (discordance, singleton call, null allele,
  resolution conflict, version-skew-resolved).
- **`Nomenclature` (facade)**: the single interface for parse / normalize / reduce /
  validate / accession-resolve. `InHouseNomenclature` is the default, shipped
  implementation; `PyArdNomenclature` is a swappable alternative behind the same
  interface (§6). No other module imports a nomenclature library directly.

Invariant: normalization is **non-lossy and non-coercive** — it never rewrites an
unrecognized name into a valid-looking one; it flags it. Raw is always retained.

---

## 5. Per-tool parsers (MVP set)

Enumerating the real format differences the adapters must absorb:

| Tool | Input type | Output artifact | Resolution emitted | Prefix / quirks | Loci |
|---|---|---|---|---|---|
| **OptiType** | DNA/RNA FASTQ | `*_result.tsv` (cols A1,A2,B1,B2,C1,C2, Reads, Objective) | 2-field | `A*02:01`, no `HLA-` prefix | Class I only (A/B/C) |
| **arcasHLA** | RNA-seq BAM | `genotype.json` (locus→[alleles]) + `genotypes.tsv` | up to 3-field | ties to the IMGT version it was built on | Class I + II |
| **HLA-LA** | WGS BAM/CRAM | `hla/R1_bestguess_G.txt` (Locus, Chromosome, Allele, Q…) | **G-group** (`A*02:01:01G`) | per-chromosome rows; quality column | Many, incl. non-classical |
| **HLA-HD** | FASTQ | `*_final.result.txt` (per locus, 2 alleles/line) | up to 3-field | `HLA-` prefix; `-`/`Not typed` for absent | Many loci |

Adapter contract: `parse(path) -> list[GenotypeCall]` with `tool` and (where knowable)
`source_db_version` populated. Optional format auto-detection, but the CLI always allows
explicit `tool:path` to avoid guessing. Golden-file fixtures per tool (§8).

---

## 6. Nomenclature normalizer + validator (highest-risk component)

**Decision: build a slim in-house core; use `py-ard` (NMDP, LGPL-3.0) as a test-only
oracle, not a shipped runtime dependency.**

Rationale: the four MVP typers emit plain allele names and G-groups only — no MAC codes,
no serologic names, no ambiguous GL-string *input*. The large, hard parts of py-ard (MAC
expansion, GL ambiguity algebra, serologic conversion) are therefore out of scope, and the
slice we actually need — parse a name, reduce to G-group / 2-field, validate — is a small
module built directly on the same IPD-IMGT/HLA files py-ard itself consumes
(`Allelelist.txt`, `hla_nom_g.txt`, `hla_nom_p.txt`). Building it in-house removes an LGPL
dependency from the distributed (possibly for-profit) artifact, removes py-ard's
network/ephemeral-cache and API-coupling risks, and costs us nothing on cross-version
handling — py-ard does **not** do cross-version reconciliation, so that layer is ours either
way. LGPL note: py-ard used only in tests does not propagate into the shipped artifact.

The one real downside of in-house is **correctness** in a clinically-adjacent domain. It is
neutralized by using py-ard purely as a **test oracle**: CI runs a large allele set through
both the in-house reducer and `pyard` and asserts agreement at the reductions we rely on.

Everything sits behind one facade (§3, §4):

```
hlaconcord.nomenclature.Nomenclature        # the interface every other module calls
├── InHouseNomenclature   (default, shipped)
└── PyArdNomenclature     (optional; behind the same interface, for fallback/benchmark)
```

The in-house implementation must:

- Parse any well-formed allele string (with/without `HLA-` prefix, legacy no-colon form).
- Preserve expression suffixes; treat null (`N`) as significant and never silently drop.
- Tag the resolution of each call.
- Reduce to a target basis on request: field truncation (2-field) **and** G/P-group
  projection, handling alleles that have no assigned G-group correctly.
- Resolve name → accession using `Allelelist_history.txt` for cross-version work.
- Validate a name as one of: exact allele in DB version · valid lower-resolution reduction
  (a real prefix of ≥1 full allele) · valid G/P group · deleted/renamed (flag with the
  mapping) · **unknown** (flag, never coerce).

If the M0 spike (§12) surfaces edge cases the in-house reducer gets wrong, adopting
`PyArdNomenclature` at runtime is a one-file swap behind the facade — reversible and cheap.

---

## 7. Concordance / consensus engine

- Choose a **comparison basis** (`--compare`, default **`lgx`**; also `g`, `p`, `2field`),
  because tools report at different depths and we must not turn a resolution difference into
  a false discordance. **Default revised from `g` to `lgx` on M0/M1 real-data evidence:** a
  2-field call (OptiType, often HLA-HD) is *ambiguous* at G-group — `A*02:01:*` spans several
  ARD sequences, so `A*02:01` maps to no single G-group (py-ard itself returns an ambiguity
  list there). `lgx` collapses every emitted resolution of the same allele to one 2-field key
  (`A*02:01`, `A*02:01:01`, `A*02:01:01G` → `A*02:01`); the in-house reducer matches py-ard
  ≥99.4% across 2-/3-field/G inputs. `g` stays available as a **stricter ARD lens** for when
  all calls resolve to ≥3 fields. Internal comparison key = the chosen basis; human-facing
  report = 2-field. (`DEFAULT_BASIS` in `hlaconcord.nomenclature`.)
- Reduce every call to that comparison basis, after version reconciliation (§2.4).
- Per locus, compare genotypes as **unordered multisets** (allele-1/allele-2 ordering is
  arbitrary; handle homozygous correctly).
- Consensus rule is configurable: unanimous / majority / weighted; report agreement level
  (e.g. `A*02:01` called by 3/4 tools).
- Flags surfaced: real discordance, singleton (one-tool-only) calls, null alleles,
  resolution conflicts, and **version-skew reconciled** (would have looked discordant on
  raw names but are the same allele by accession).

---

## 8. Testing strategy (non-negotiable given clinical adjacency)

- **Golden fixtures**: real/public example outputs from each MVP tool → parser tests.
- **Nomenclature unit tests** across tricky cases: legacy no-colon, every suffix
  (`N/L/S/Q/C/A`), G and P groups, homozygous, "not typed", MAC codes, GL strings.
- **Property tests**: normalization is idempotent; reduce-then-validate holds; round-trip.
- **Cross-version test**: an allele renamed between two releases resolves via accession to
  the same identity → *not* flagged discordant.
- **Known-answer concordance**: constructed multi-tool inputs with predetermined
  agreement/disagreement.
- **Explicit failure-mode tests**: never drop a null suffix; never coerce an invalid name;
  reject/flag unknown genes and malformed strings.

---

## 9. Reference database management

- Bundle one pinned IPD-IMGT/HLA release for reproducible, citable output.
- `hlacc db list | update | pin <version>` to fetch newer releases into a local cache.
- Always record the DB version used in every output row and report header.
- Only fetch the small name/group files (`Allelelist*`, `hla_nom_g/p`), not sequences.

---

## 10. CLI (surface sketch)

```
hlacc run \
  --inputs optitype:S1.tsv arcashla:S1.json hla-la:S1_bestguess_G.txt hla-hd:S1_final.result.txt \
  --sample S1 --resolution 2 --db 3.55.0 -o out/

hlacc run --samplesheet samples.csv -o out/     # batch: sample,tool,path[,db_version]
hlacc validate  HLA-A*02:01  A*99:99            # quick name check
hlacc normalize A*0201                          # show canonical form + accession
hlacc db list|update|pin
```

Backed by an importable package (`hlaconcord.parse`, `hlaconcord.nomenclature`,
`hlaconcord.validate`, `hlaconcord.concordance`) for pipeline integration.

## 11. Outputs

- **Tidy long table** (CSV/TSV/Parquet): `sample, locus, allele_slot, tool, raw_call,
  normalized_call, resolution, validation_status, accession, g_group, p_group, db_version`.
- **Concordance report**: per sample×locus consensus genotype, per-tool agreement,
  discordance flags — CSV + human-readable summary (optional HTML).
- **JSON** for programmatic consumers; **GL string** export option.

---

## 12. Milestones

- **M0 — domain core + spike. ✓ DONE.** Built the `Nomenclature` facade and the in-house
  reducer (parse, reduce to G-group / 2-field / lgx / P, validate) on the IPD files, plus the
  accession-history layer. Spike ran against real release 3.55.0 (5,000 A/B/C/DRB1/DQB1
  alleles) with `pyard.init('3550')` as oracle. **Outcome — build in-house, confirmed:**
    - **G basis (our default): in-house matches py-ard exactly** except three G-groups where
      py-ard *relabels* the group to a different member (same membership; in-house emits the
      name in IPD's `hla_nom_g.txt`, the more faithful choice): `A*02:17:01G`, `C*02:10:01G`,
      `DRB3*03:22:01G`. Allow-listed in `tests/test_oracle.py`.
    - **lgx basis: 0 unexpected**; py-ard additionally collapses ~0.58% of alleles that IPD
      puts in *no* published G-group, via exon-2/3 sequence identity we deliberately don't
      ingest (§6). Rare alleles the MVP typers don't emit; accepted and counted.
    - The spike caught one real bug — reductions dropped the expression suffix, silently
      voiding null (`N`) status — now fixed; every basis preserves it (regression test added).
    - **Key reframing:** py-ard is *not* a pure ground truth (it relabels 3 groups vs the IPD
      file), so the CI oracle asserts "no divergence beyond the allow-listed, understood set,"
      not byte-equality. py-ard stays a test-only oracle; the shipped path is in-house.
    - **Cross-version churn (3.51.0 → 3.55.0):** of 30,453 MVP-locus alleles common to both
      releases, 72 were renamed (0.24%). The G-group basis alone collapses 36% of that churn;
      accession-history reconciliation then recovers **100% (72/72)** — every old name resolves
      to the correct accession in the newer DB, so version skew never masquerades as
      discordance. Confirms §2.4 + §7.
  Test-heavy throughout. M0 complete.
- **M1 — first two parsers + tidy table. ✓ DONE.** OptiType (`_result.tsv`) + arcasHLA
  (`genotype.json`) adapters behind a registry (`parse(path) -> list[GenotypeCall]`);
  call-level model (`AlleleCall`/`GenotypeCall`); normalization pass enriching every call
  through the facade (parse→validate→reduce→accession, non-coercive); tidy long-table emitter.
  Verified end-to-end on real 3.55.0 data (OptiType 2-field + arcasHLA 3-field harmonize to a
  common lgx key). Revised the default comparison basis `g`→`lgx` on this evidence (§7).
- **M2 — concordance/consensus engine + report. ✓ DONE.** Groups normalized calls by
  sample+locus; unordered, zygosity-tolerant set matching (homozygous reporting quirks never
  fabricate discordance); configurable consensus (unanimous/majority) with per-allele
  agreement (`A*02:01 3/4`); flags: discordance, singleton, null allele, resolution conflict
  (prefix-compatible-but-unequal — silent under lgx by design), version-skew-resolved
  (accession-proven, not mere depth differences), invalid call. Tidy concordance table +
  human summary emitters. Known-answer tests + real-data integration.
- **M3 — remaining MVP parsers. ✓ DONE.** HLA-LA (`_bestguess_G.txt`, per-chromosome rows
  paired per locus, G-group alleles, Q1 quality captured; sample-id falls back to the sample
  dir for the generic `R1_` filename) + HLA-HD (`_final.result.txt`, `HLA-` prefix, `-`/`Not
  typed` dropped). Verified all four tools harmonize on real 3.55.0 data: 2-field, 3-field,
  G-group, and prefixed calls reduce to one lgx key per allele and concord across loci.
- **M4 — CLI + reference-DB management. ✓ DONE.** `hlacc` CLI (`hlaconcord.cli:main`)
  with `run` / `validate` / `normalize` / `db`. `run` ingests `--inputs tool:path …` or a
  `--samplesheet` (sample,tool,path[,db_version]), harmonizes, and writes the tidy table,
  concordance TSV, sample-grouped JSON, and (with `--gl`) consensus GL strings; exit code
  `1` on any discordance, `2` on a config/DB error. Importable `hlaconcord.pipeline.run`
  behind it. Reference-DB layer (`hlaconcord.db`): on-demand fetch of a release's small
  name/group files from ANHIG/IMGTHLA into a local cache (`db update`), discovery + pinning
  (`db list|pin|path`), version↔IPD-branch conversion. GL-String export (`hlaconcord.gl`) from
  the consensus. Verified end-to-end on real 3.55.0 data; ~49 new tests.
- **M5 — hardening. ✓ DONE.** Packaging metadata (classifiers, URLs, entry point) + `ruff`
  lint config (clean); MIT `LICENSE` (with an IPD-IMGT/HLA data-terms notice) and
  `CHANGELOG.md`; expanded `README.md` (install / quickstart / CLI / outputs / library);
  an `examples/` validation set (four-typer sample + expected output); GitHub Actions CI —
  a fast matrix (py3.10–3.12: ruff + the 113-test suite) plus a py-ard **oracle gate** job
  that bootstraps its reference data with `hlacc db update` and runs the M0 cross-check.

---

## 13. Risks & failure modes

| Risk | Mitigation |
|---|---|
| Silent wrong normalization (the big one) | Validate everything; never coerce; keep raw alongside normalized; in-house reducer checked against the `py-ard` oracle in CI (allow-listing only the understood divergences — M0 §12) |
| G-group *label* differs between sources | py-ard relabels 3 IPD groups (same membership); we emit the IPD `hla_nom_g.txt` name and pin the deltas in the oracle test, so a py-ard release change surfaces as a test diff, not a silent shift |
| DB version skew → false discordance | Accession-history reconciliation before comparison (§2.4); G-group comparison basis (§7) collapses most field-3/4 churn |
| MAC / GL-string subtlety | Out of MVP scope (the MVP tools don't emit them); if added later, validate against the `py-ard` oracle |
| Dropped null (`N`) allele | Treat suffix as significant; explicit tests |
| Non-classical loci / coverage differences between tools | Model loci explicitly; concordance handles missing loci as "not called", not disagreement |

---

## 14. Maintenance (per the brief)

- **Reference DB currency** — regular, mechanical (`db update`); track ANHIG/IMGTHLA releases.
- **New parser when a typer appears** — occasional; isolated to an adapter.
- **Upstream nomenclature lib** — if we adopt `py-ard`, track its releases.

Lower churn than a vendor-format chaser. Clinical relevance → citations and sticky adoption.

---

## 15. Decisions & open questions

**Resolved:**
1. MVP typers: OptiType, arcasHLA, HLA-LA, HLA-HD. ✓
2. Nomenclature core: in-house slim reducer/validator; `py-ard` as a **test-only oracle**,
   not a runtime dependency; behind a swappable facade (§6). **Confirmed by the M0 spike (§12):**
   in-house matches py-ard on the G basis modulo 3 allow-listed group *relabels* (where in-house
   is the more IMGT-faithful), and py-ard's only extra reach is sequence-based lgx on rare
   ungrouped alleles that MVP typers don't emit. ✓
3. Comparison basis: default **`lgx`**, report at 2-field; configurable
   (`--compare lgx|g|p|2field`) (§7). **Revised from `g` on M1 real-data evidence** — a
   2-field call (OptiType/HLA-HD) is ambiguous at G-group and would false-flag as discordant;
   `lgx` collapses all emitted resolutions cleanly and matches py-ard ≥99.4%. `g` retained as
   a stricter lens. **Confirmed.** ✓

**Still open:**
4. Whether MAC expansion and GL-string *input* (not just output) are in scope for v1
   (currently out — the MVP tools don't emit them).

**Resolved since:**
5. Licensing / distribution of the IPD-IMGT/HLA data. **Sidestepped for now by not
   redistributing it:** hlaconcord fetches a release's small name/group files on demand
   (`hlacc db update`) into a per-user cache and records the version in every output,
   rather than bundling the database. `LICENSE` carries an IPD-IMGT/HLA data-terms/citation
   notice. If we later choose to bundle a pinned release for turnkey/offline use, the
   attribution/citation terms must be revisited then. (§9 "bundle a release" is thus
   deferred, not done — the fetch path supersedes it for the MVP.)
