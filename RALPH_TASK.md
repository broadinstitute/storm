---
task: Build STORM v0.2 (Integrated-callset-first + overlays + backends)
---

## Goal (one sentence)
Build STORM as a Rust crate with a Python (Jupyter-friendly) front-end that ingests an integrated LR SV VCF as default genotype authority, supports ordered overlays (TRGT now; HLA/KIR later), materializes an Arrow/Parquet cache, and runs association via pluggable backends including an internal StormGLM (no mixed models).

## Non-goals (this milestone)
Do NOT implement mixed-model association (no SAIGE/REGENIE execution).
Do NOT implement Hail integration beyond export-friendly formats.
Do NOT optimize for full biobank-scale throughput yet (correctness > speed).

## Core invariants (must be explicit in code)
Separate “presence/zygosity” (primary) from “allele value” (overlays can refine).
All model/encoding/backend selection is plan-driven and deterministic.
Every resolved genotype has provenance: presence_source, allele_source, QC gating.
No silent assumptions for missing genotypes (missing must be handled explicitly).

---

# A. Rust CLI skeleton + project structure
[ ] Add/confirm CLI entrypoint(s) for:
  [ ] `storm cache build`
  [ ] `storm cache verify`
  [ ] `storm explain`
  [ ] `storm assoc run`
  [ ] `storm assoc export` (export-only stubs)
[ ] Establish Rust module layout (suggested):
  [ ] `vcf/` (parsing)
  [ ] `catalog/` (TRExplorer ingest + interval index)
  [ ] `test_units/` (construction/grouping)
  [ ] `resolver/` (primary+overlay resolution)
  [ ] `encoding/` (S/M/D/G/T/Bins)
  [ ] `plan/` (YAML schema + rule engine)
  [ ] `assoc/` (backend trait + StormGLM)
  [ ] `io/` (Arrow/Parquet)
  [ ] `provenance/`
[ ] Add `.ralph/progress.md` and `.ralph/guardrails.md` (if not present).

---

# B. VCF ingestion (primary: integrated callset)
[ ] Implement integrated SV VCF parser:
  [ ] Parse CHROM, POS, ID, REF, ALT, FILTER.
  [ ] Parse INFO fields (at least SVTYPE, SVLEN, END if present).
  [ ] Parse per-sample GT (and optionally GQ/DP if present).
  [ ] Emit normalized `SvRecord` table + per-sample genotype table.
[ ] Define stable `sv_id` semantics:
  [ ] Prefer VCF ID when present; else derive stable hash from (chrom,pos,ref,alt).
[ ] Ensure missing GT handling is explicit (no “absent = 0/0” assumptions).

---

# C. TRExplorer catalog ingestion (BED + EH JSON)
[ ] Implement catalog ingest:
  [ ] Read BED intervals (chr,start,end,motif).
  [ ] Read EH JSON annotations (ReferenceRegion, LocusId, Motif, CanonicalMotif, NumRepeatsInReference, etc.).
  [ ] Join BED+JSON into internal `CatalogLocus` objects keyed by JSON `LocusId`.
  [ ] Build interval index for overlap queries.
[ ] Expose catalog helpers:
  [ ] motif length (bp)
  [ ] reference repeat count (NumRepeatsInReference)
  [ ] reference tract length (repeat units and/or bp)

---

# D. TestUnit construction (SV + RepeatProxySV + Repeat)
[ ] Define internal `VariantKind`:
  [ ] `SV`
  [ ] `RepeatProxySV` (derived from SV alleles overlapping catalog loci/clusters)
  [ ] `Repeat` (TRGT-derived true allele lengths)
  [ ] placeholders for `HLA`, `KIR`
[ ] Implement SV→RepeatProxySV mapping:
  [ ] For each SV record, overlap with catalog loci.
  [ ] If overlaps: tag with `catalog_locus_id` (and `cluster_id` if supported).
[ ] Group SV records into RepeatProxy TestUnits:
  [ ] `test_unit_id = catalog_locus_id` (or cluster-based id later)
  [ ] maintain list of contributing `sv_id`s.

---

# E. TRGT overlay ingestion (VCF)
[ ] Implement TRGT VCF parser:
  [ ] Parse site fields: TRID, END, MOTIFS/STRUC where available.
  [ ] Parse per-sample fields: GT and AL (e.g., "39,47").
  [ ] Parse minimal QC fields (configurable; do not overfit).
[ ] Map TRGT calls to catalog loci:
  [ ] Use interval overlap tolerance + motif match where possible.
  [ ] Record mapping confidence / ambiguity.
[ ] Represent TRGT allele lengths in repeat units (preferred) and/or bp.

---

# F. Resolver (primary + ordered overlays)
[ ] Implement resolver core:
  [ ] Input: primary SV genotypes + 0..N overlays (ordered).
  [ ] Output: resolved genotype record per (test_unit_id, sample_id) with provenance.
[ ] Presence layer default:
  [ ] from integrated callset for SV and RepeatProxySV.
[ ] Allele layer for repeat-like units:
  [ ] Prefer TRGT AL if QC-pass.
  [ ] Else fallback to RepeatProxySV reconstruction (if valid).
  [ ] Else fallback to carrier-only (presence) with allele fields null.
[ ] Record provenance fields:
  [ ] `presence_source`
  [ ] `allele_source`
  [ ] QC flags / override reasons.

---

# G. RepeatProxySV diploid reconstruction (multi-allelic reversal)
[ ] Implement grouping of allele rows per catalog locus:
  [ ] INS and DEL rows included.
  [ ] Use SVLEN sign rule: INS +SVLEN bp; DEL −abs(SVLEN) bp.
[ ] Implement diploid inference per sample:
  [ ] Start from baseline (reference tract length from catalog).
  [ ] Apply deltas based on GT per allele row:
    [ ] 0/0: none
    [ ] 0/1: one haplotype
    [ ] 1/1: both haplotypes
  [ ] Handle compound het: two different 0/1 alleles → one on each haplotype.
  [ ] If >2 distinct non-ref alleles present → mark sample ambiguous and drop for that TestUnit.
[ ] Canonicality check:
  [ ] If delta bp divisible by motif length → convert to repeat units.
  [ ] Else mark non-canonical; configurable behavior (exclude from repeat-proxy vs keep bp).

---

# H. Feature computation (phenotype-agnostic)
[ ] Compute per-TestUnit features from resolved genotypes:
  [ ] call rate
  [ ] carrier count
  [ ] allele spectrum summaries (min/median/max)
  [ ] tail metrics for max allele
  [ ] multimodality heuristic (simple v0.2)
[ ] Write features to `features.parquet`.

---

# I. Arrow/Parquet cache
[ ] Implement `storm cache build`:
  [ ] Write `test_units.parquet`
  [ ] Write `genotypes.parquet` (long form: test_unit_id, sample_id, presence/zyg, L1,L2, qc, sources)
  [ ] Write `catalog.parquet` (or compact form)
  [ ] Write `features.parquet`
  [ ] Write `provenance.json` (input hashes, versions, git commit)
[ ] Implement `storm cache verify`:
  [ ] Validate schemas and expected row counts.

---

# J. Plan (YAML) + rule engine (“executable pre-registration”)
[ ] Implement plan schema + parser:
  [ ] ordered rules
  [ ] match predicates on VariantKind + features + catalog tags
  [ ] action chooses encoding + model + backend preference
[ ] Implement rare-event ladder:
  [ ] <20 carriers → BinomiRare
  [ ] 20–200 → Firth
  [ ] >200 → logistic
[ ] Record chosen `rule_id`, encoding, model, backend for every test result.

---

# K. Encodings (predictor builders)
[ ] Implement encoders:
  [ ] S = L1 + L2
  [ ] M = max(L1, L2)
  [ ] D = |L1 − L2|
  [ ] G = tail indicator (q=0.995 on M; cohort or stratum-specific)
  [ ] T = RINT(S) (VNTR option)
  [ ] Bins(S) / categorical bins (narrow-window option)
[ ] Ensure encoders produce:
  [ ] predictor vector X
  [ ] sample mask (valid samples)
  [ ] encoding metadata (thresholds/bins used)

---

# L. Association backends
[ ] Define backend trait:
  [ ] `run(test_unit_id, X, phenotype, covariates, mask, model_spec) -> ResultRow`
[ ] Implement internal backend `StormGLM` (no mixed models):
  [ ] Linear regression (quantitative)
  [ ] Logistic regression (binary)
  [ ] Factor regression (multi-df LRT) for categorical predictors
  [ ] BinomiRare for rare binary predictors
  [ ] Firth logistic (feature-flagged OK if needed; but implement if possible)
  [ ] Covariates matrix (incl PCs)
[ ] Population stratification support (internal):
  [ ] accept PCs as covariates
  [ ] optional ancestry-stratified runs + fixed-effect meta-analysis
[ ] External backend export stubs:
  [ ] Export SAIGE inputs (phenotype/covariates/predictors)
  [ ] Export REGENIE inputs (phenotype/covariates/predictors)
  [ ] No execution required; export-only.

---

# M. Results writing (backend-agnostic)
[ ] Write results to Parquet with stable schema:
  [ ] ids: test_unit_id, variant_kind, catalog_locus_id, contributing_sv_ids
  [ ] model meta: encoding, model, backend, rule_id, plan_hash
  [ ] stats: beta/OR, SE, CI, p, df, convergence flags
  [ ] counts: N, carriers, call_rate
  [ ] provenance: presence_source, allele_source, qc flags

---

# N. Explainability
[ ] Implement `storm explain <test_unit_id> [--sample <id>]`:
  [ ] contributing SV records (IDs, SVTYPE/SVLEN)
  [ ] overlay TRGT record(s) used (if any)
  [ ] final resolved genotype (presence, L1,L2)
  [ ] chosen encoding/model/backend/rule_id
  [ ] provenance + QC flags

---

# O. Python front-end (Jupyter)
[ ] Provide Python package:
  [ ] build cache (shell out or bindings)
  [ ] load cache Parquet into pandas/polars
  [ ] run association (StormGLM) and load results
  [ ] `explain()` wrapper
[ ] Add a small notebook showing:
  [ ] load VCFs + catalog
  [ ] build cache
  [ ] run one association
  [ ] explain one repeat locus

---

# Acceptance criteria
[ ] Unit tests exist for:
  [ ] integrated VCF parsing
  [ ] TRGT VCF parsing (AL parsing)
  [ ] catalog BED+JSON ingest + locus ID mapping
  [ ] RepeatProxy grouping + diploid reconstruction (incl compound het and ambiguity)
  [ ] overlay precedence (TRGT overrides allele; presence from integrated)
  [ ] StormGLM models on tiny synthetic matrices
[ ] End-to-end run on minimal inputs produces:
  [ ] cache Parquet outputs with stable schema
  [ ] results Parquet output with stable schema
  [ ] deterministic `storm explain` output
[ ] When complete, write `DONE` in `.ralph/progress.md`.