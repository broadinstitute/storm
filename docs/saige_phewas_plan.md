# SAIGE PheWAS Integration Plan (Storm SV + TR)

This document turns the SAIGE integration into incremental, testable phases.

## Goal

Run a single genome-wide PheWAS framework with SAIGE using two clean test strata,
preserving TR-specific repeat logic where it adds genuine modeling value:

- standard SV/non-repeat dosage tests
- repeat-aware quantitative tests using a dosage-weighted site-level TR proxy

## Scope and assumptions

- Input MatrixTable is `mt1` from `storm.annotate_svs(...)`.
- `mt1` has row fields `allele_id`, `svlen`, and optional `tr` struct.
- SAIGE is installed and loadable in R.
- One combined results table downstream, with stratum labels and TR metadata
  annotations preserved for locus-specific follow-up.

## Analysis strata

Two strata:

1. `standard_sv`
   - Variants without usable repeat-length encoding, OR TR variants where
     `repeat_units_estimate` variance across samples is below a dispersion
     threshold (i.e., effectively biallelic in length).
   - Predictor: genotype dosage (`0/1/2`).

2. `tr_quantitative`
   - Variants with defined `repeat_units_estimate`.
   - Predictor: `dosage * repeat_units_estimate` (optionally standardized).
   - Interpretation: dosage-weighted site-level repeat proxy, not a true
     diploid per-sample repeat count.
   - Catalog locus membership (`tr_locus_id`) and rule applicability
     (`rule_applicable`) carried as metadata columns, not as separate test types.

## Encoding caveat (accepted tradeoff)

Current `tr_quantitative` encoding uses row-level `repeat_units_estimate`
projected to entries via genotype dosage:

- `0/0` -> `0`
- `0/1` -> `1 * repeat_units_estimate`
- `1/1` -> `2 * repeat_units_estimate`

This does not capture within-genotype heterogeneity in true repeat composition
(e.g., interruptions, per-sample allele sequence variation) that may be blurred
in a large merged multi-caller joint SV representation. We accept this tradeoff
for the primary scan and reserve allele-resolved TR quantitation for follow-up.

## Phase plan

### Phase 1 (implemented now): feature inventory + strata assignment

Deliverables:

- Derive a per-feature table from `mt1.rows()`.
- Assign `feature_class` (`standard_sv` or `tr_quantitative`) from
  availability of finite `repeat_units_estimate`.
- Include key metadata: `tr_locus_id`, `rule_applicable`, `motif`, `svtype`.
- Print counts by class for chr22; extrapolate to genome-wide expectation.

Purpose:

- Validate feature population and strata balance before building SAIGE marker files.
- Validate that two-strata assignment is stable and reproducible.

### Phase 2: export SAIGE-ready marker matrices

Deliverables:

- Stratum-specific marker exports with consistent `sample_id` order.
- Quantitative TR encoding for `tr_quantitative`:
  `predictor = dosage * repeat_units_estimate` (standardization optional).
- Standard dosage marker export for `standard_sv`.

Implementation (Storm):

- **In-memory path**: `storm.build_dense_saige_marker_matrix(mt, stratum=...)` then
  `storm.align_matrix_cols_to_manifest(mt, manifest_ht)`; or
  `storm.export_saige_stratum_vcfs(mt, manifest_ht, out_dir=..., prefix=...)`.
- **Sparse staging replay**: `storm.build_feature_vcf_row_lookup(mt)` plus
  `storm.dense_marker_matrix_from_long(long_ht, row_lookup)` (e.g. long tables
  from `storm.build_long_predictor_tables` or from exported `.ht` / TSV).
- **SAIGE Step 2**: `storm.export_saige_dosage_vcf(mt, "…vcf.bgz")` writes FORMAT
  `DS`; use `--vcfField=DS` with `--vcfFile` / `--vcfFileIndex` as in SAIGE docs.
  `tr_quantitative` values can exceed 2.0 (dosage-weighted repeat proxy); ensure
  your SAIGE build tolerates that for continuous encodings.

Checks:

- No missing/invalid marker IDs.
- No sample ID/order drift across exports.
- Distribution plots of `tr_quantitative` predictors to confirm signal range.

### Phase 3: phenotype/covariate + null model wiring

Deliverables:

- Stable phenotype/covariate table(s): age, sex, sequencing depth, ancestry PCs.

Storm helpers (synthetic cohort for pipeline tests):

- `storm.build_synthetic_saige_pheno_covar_table(manifest_ht, n_pcs=5)` and
  `storm.export_saige_phenotype_covariate_tsv` for a tab-delimited Step 1 file;
  `storm.saige_synthetic_covar_col_list()` for `--covarColList`.
- R driver: `scripts/saige/run_saige.R` + `scripts/saige/example_chr22.config` — runs
  SAIGE `step1_fitNULLGLMM.R` then `step2_SPAtests.R` for both DS VCF strata
  (`--vcfField=DS`). Set `SAIGE_EXTDATA` or install the **SAIGE** R package; fill in
  `plink_prefix` and/or sparse GRM paths in the config before Step 1.
- SAIGE step 1 command templates and outputs (`.rda`, variance ratios) — one null
  model per ancestry stratum per phenotype.
- GRM construction strategy documented (sparse vs. full, and rationale).
- Sample inclusion parity checks against marker exports.
- Binary vs. quantitative phenotype handling documented separately (different
  SAIGE modes).

### Phase 4: step 2 runs + unified result assembly

Deliverables:

- SAIGE step 2 per stratum, per ancestry, per phenotype.
- METAL meta-analysis across ancestry strata.
- Merged results table with columns:
  - `feature_id`, `feature_class`, `test_type`, `phenotype`, `beta`, `se`, `pval`
  - TR metadata: `tr_locus_id`, `rule_applicable`, `motif`
  - FDR columns: stratum-specific FDR + global FDR
- Secondary dosage re-test at genome-wide significant `tr_quantitative` hits,
  to confirm the continuous encoding was doing real work (anticipates reviewer
  requests).

## Operational guidance

- Start with one phenotype on chr22 to de-risk the full pipeline end-to-end.
- Scale to full genome only after chr22 smoke run validates all four phases.
- Use class-specific QC thresholds: `tr_quantitative` markers will be sparser
  and noisier than `standard_sv` markers.
- Keep all derived files versioned with sidecar `schema_version`, `ruleset_version`,
  and SAIGE command/config hash.
- Locus-specific follow-up (catalog rule-bin models, motif purity, inheritance
  mode, pathogenicity threshold testing) is deferred to post-discovery and
  triggered only by genome-wide significant hits in the primary scan.

## Next implementation target

Phase 3: phenotype/covariate tables, SAIGE Step 1 null models, and GRM strategy
(documented command templates and sample inclusion parity with marker VCFs).