# STORM v0.2: Integrated-callset–first SV/TR association framework

## One-liner
Build STORM as a Rust crate with a Python (Jupyter-friendly) front-end that:
1) uses an integrated long-read SV callset as the default authority for variant presence/zygosity,
2) supports ordered overlay datasets (TRGT; later HLA/KIR) that refine or override allele values in specific regions/loci,
3) materializes an internal Arrow/Parquet cache,
4) runs association tests via pluggable backends, including a minimal internal GLM backend.

STORM = Structural & Tandem-Repeat Optimized Regression Models.

---

## Non-goals (for this milestone)
- No mixed models (SAIGE/REGENIE execution not required yet).
- No Hail integration beyond export-friendly formats.
- No full genome-wide scan optimizations; correctness > speed.

---

## Core concepts (must be explicit in code)
- **Primary genotype authority**: integrated SV callset (presence/zygosity).
- **Overlay datasets**: refine allele values (e.g., TRGT lengths; HLA/KIR alleles).
- **TestUnit**: the atomic unit of association (SV, repeat locus/cluster, HLA gene, etc.).
- **Resolver**: deterministic engine that combines primary + overlays into a final analysis genotype with provenance.
- **Plan-driven**: encoding, model, and backend are chosen by a pre-specified plan, never ad hoc.

---

## Deliverables

### A. Rust core (`storm` crate)

#### A1. VCF ingestion — integrated callset (primary)
- Parse site fields: CHROM, POS, ID, REF, ALT, FILTER.
- Parse INFO: SVTYPE, SVLEN, END (if present).
- Parse per-sample FORMAT: GT (plus optional GQ/DP).
- Emit normalized `SvRecord` + per-sample genotype table.
- Treat missing GTs as invalid (no silent assumptions).

#### A2. VCF ingestion — TRGT overlay
- Parse site fields: TRID, END, MOTIFS/STRUC.
- Parse per-sample fields: GT, AL (allele lengths, e.g. "39,47").
- Gate overlay usage on simple QC (configurable).
- Map TRGT calls to TRExplorer loci by interval overlap + motif.

#### A3. Catalog ingestion — TRExplorer
- Read BED (interval backbone) and EH JSON (annotations).
- Build interval index for overlap queries.
- Use JSON `LocusId` as stable `catalog_locus_id`.
- Expose motif length and reference repeat count.

---

### B. TestUnit construction

- VariantKind:
  - `SV` (standard biallelic / CN SVs)
  - `RepeatProxySV` (SV alleles overlapping a TRExplorer locus/cluster)
  - `Repeat` (true repeat genotypes from TRGT)
  - (placeholders for `HLA`, `KIR`)

- Group SV records overlapping the same `catalog_locus_id` (or cluster) into one RepeatProxy TestUnit.
- Keep explicit links to contributing SV IDs.

---

### C. Resolver (primary + overlays)

- Deterministic, ordered overlay application.
- Separate layers:
  - **Presence layer** (default: integrated callset).
  - **Allele layer** (TRGT preferred; SV-proxy fallback).
- Record provenance:
  - presence_source
  - allele_source
  - overlay rule applied
  - QC gating decisions.

#### C1. RepeatProxySV reconstruction
- For each TestUnit with multiple INS/DEL alleles:
  - Baseline = reference repeat length (from catalog).
  - Apply deltas per allele row using GT:
    - INS → +SVLEN bp
    - DEL → −abs(SVLEN) bp
  - Convert bp → repeat units iff divisible by motif length.
  - Reconstruct diploid (L1, L2):
    - 0/1 on different alleles ⇒ compound heterozygote.
  - If ambiguous (>2 non-ref alleles), flag sample as invalid for this TestUnit.

---

### D. Arrow / Parquet cache

- `storm cache build` writes:
  - `test_units.parquet`
  - `genotypes.parquet`
    - test_unit_id, sample_id
    - presence, zygosity
    - L1, L2 (or null)
    - qc flags
    - presence_source, allele_source
  - `catalog.parquet`
  - `features.parquet` (phenotype-agnostic features)
  - `provenance.json` (input hashes, git commit, timestamps)

- `storm cache verify` validates schema + row counts.

---

### E. Feature computation (phenotype-agnostic)
For each TestUnit:
- call rate
- carrier count
- allele spectrum summaries
- tail metrics (for derived indicators)
- multimodality flags (simple heuristic)

---

### F. Plan + rule engine (“executable pre-registration”)

- YAML schema specifying ordered rules:
  - match on VariantKind + features + catalog tags
  - choose:
    - encoding (S, M, tail indicator, categorical bins, etc.)
    - model (linear, logistic, binomirare, firth, factor)
    - backend preference
- Implement rare-event ladder:
  - <20 carriers → BinomiRare
  - 20–200 → Firth
  - >200 → logistic
- First-match-wins, deterministic.
- Must emit rule_id used.

---

### G. Association backends

#### G1. Backend interface
Define a trait like:
- `run(test_unit, predictor, phenotype, covariates, mask) -> ResultRow`

#### G2. Internal backend: StormGLM
Support:
- Linear regression (quantitative traits).
- Logistic regression (binary traits).
- Categorical (factor) regression with multi-df LRT.
- BinomiRare exact test for rare binary predictors.
- Firth logistic (feature-flagged but strongly recommended).
- Covariates (incl. PCs).
- Optional ancestry-stratified runs + fixed-effect meta-analysis.

No mixed models.

#### G3. External backend stubs (export-only)
- Export inputs for SAIGE / REGENIE:
  - phenotype file
  - covariates
  - predictor matrices (continuous, binary, dummy-coded)
- No execution required yet.

---

### H. Results schema (backend-agnostic)
Write results to Parquet with:
- identifiers: test_unit_id, variant_kind, catalog_locus_id
- model metadata: encoding, model, backend, rule_id, plan_hash
- statistics: beta/OR, SE, CI, p, df
- counts: N, carriers, call_rate
- provenance: presence_source, allele_source, flags

---

### I. Explainability
- `storm explain <test_unit_id> [--sample id]`
- Must print:
  - contributing SV records
  - overlays applied (TRGT, etc.)
  - final resolved genotype
  - encoding, model, backend, rule_id
  - provenance + QC notes

---

### J. Python front-end
- Python package to:
  - build/read cache
  - run association via StormGLM
  - load results into pandas/polars
- Include a small Jupyter demo notebook.

---

## Acceptance Criteria
- Unit tests for:
  - VCF parsing (integrated + TRGT).
  - Catalog ingestion (BED + JSON).
  - RepeatProxy grouping + diploid reconstruction.
  - Overlay precedence (TRGT overrides allele, not presence).
  - Each StormGLM model on tiny synthetic matrices.
- End-to-end run on fixtures:
  - build cache
  - run association
  - deterministic `explain` output.
- Stable Parquet schemas + provenance recorded.

---

## Plan of attack (PR-sized steps)
1) CLI skeleton + integrated VCF parser.
2) TRExplorer catalog ingestion + overlap mapping.
3) TestUnit construction + RepeatProxy reconstruction.
4) TRGT overlay ingestion + resolver precedence.
5) Arrow/Parquet cache writer + verifier.
6) Plan/rule engine.
7) StormGLM backend + results writer.
8) Explain command + Python demo.

---

## Completion condition
When:
- cache builds from fixtures,
- StormGLM runs at least one association,
- `storm explain` shows SV vs TRGT resolution clearly,
write **DONE** in `.ralph/progress.md` and stop.