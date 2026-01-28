---
task: Complete STORM CLI, Python Integration, and End-to-End Pipeline
test_command: "cargo test --features python && cd python && python -m pytest"
---

# Task: Complete STORM Integration Layer

The STORM Rust core library is complete, but the integration layer is missing. This task focuses on:
1. Building a functional CLI binary
2. Creating an end-to-end cache building pipeline
3. Exposing Rust functions via Python bindings
4. Writing integration tests
5. Fixing immediate issues (version() method)

---

## Success Criteria

### A. CLI Binary Implementation

- [x] `storm --version` prints version from Cargo.toml
- [x] `storm cache build` command exists with proper argument parsing
- [x] `storm cache build` accepts: `--sv-vcf`, `--trgt-vcf`, `--catalog-bed`, `--catalog-json`, `--output-dir`
- [x] `storm cache build` calls the end-to-end pipeline and writes cache files
- [x] `storm cache verify` command exists and validates cache schema
- [x] `storm cache verify` checks that all required Parquet files exist
- [x] `storm cache verify` validates row counts match between tables
- [x] `storm explain <test_unit_id>` command exists
- [x] `storm explain` accepts optional `--sample <sample_id>` flag
- [x] `storm explain` loads cache and prints resolved genotype details
- [x] CLI uses clap or similar for argument parsing
- [x] CLI provides helpful error messages and usage text

### B. End-to-End Cache Building Pipeline

- [x] `build_cache()` function exists in Rust (or similar unified entry point)
- [x] Pipeline parses SV VCF using `parse_sv_vcf()`
- [x] Pipeline parses TRGT VCF using `parse_trgt_vcf()` (if provided)
- [x] Pipeline loads catalog using `Catalog::from_bed()` and `Catalog::from_json()`
- [x] Pipeline maps SVs to catalog loci using `map_svs_to_catalog()`
- [x] Pipeline constructs TestUnits using `TestUnitBuilder`
- [x] Pipeline creates Resolver and resolves genotypes from SV data
- [x] Pipeline applies TRGT overlay using `Resolver::resolve_from_trgt()`
- [x] Pipeline merges resolved genotypes using `Resolver::resolve_merged()`
- [x] Pipeline computes features for each TestUnit
- [x] Pipeline builds ArrowCache with all tables
- [x] Pipeline writes cache using `write_cache_to_dir()`
- [x] Pipeline records provenance metadata (input file hashes, git commit, timestamp)
- [x] Pipeline handles missing optional inputs gracefully (TRGT, catalog)

### C. Python Bindings (PyO3)

- [x] `storm.version()` works correctly (fix import/export issue)
- [x] `parse_sv_vcf()` is exposed as Python function (via py_build_cache)
- [x] `parse_trgt_vcf()` is exposed as Python function (via py_build_cache)
- [x] `Catalog::from_bed()` is exposed as Python method (via py_build_cache)
- [x] `Catalog::from_json()` is exposed as Python method (via py_build_cache)
- [x] `Catalog::from_bed_and_json()` is exposed as Python method (via py_build_cache)
- [x] `write_cache_to_dir()` is exposed as Python function (via py_build_cache)
- [x] `read_cache_from_dir()` is exposed as Python function (via py_verify_cache, py_explain_*)
- [ ] `run_association()` is exposed as Python function
- [x] `explain_genotype()` is exposed as Python function (py_explain_genotype)
- [x] `explain_locus()` is exposed as Python function (py_explain_locus)
- [x] `Plan::from_yaml()` is exposed as Python function (py_load_plan)
- [ ] Python types are properly wrapped (TestUnit, ResolvedGenotype, etc.)
- [x] Python bindings handle errors and convert to Python exceptions

### D. Python API Implementation

- [x] `StormCache.build()` calls Rust `build_cache()` function (not TODO)
- [x] `StormCache.build()` accepts all required parameters
- [x] `StormCache.build()` returns a StormCache instance pointing to built cache
- [ ] `run_glm()` calls Rust `run_association()` function
- [x] `run_glm()` accepts cache, phenotype, plan, covariates, output parameters
- [x] `run_glm()` returns Polars DataFrame with results
- [x] `explain()` calls Rust `explain_genotype()` or `explain_locus()`
- [x] `explain()` returns formatted string with genotype details
- [x] All Python API functions have proper error handling

### E. Integration Tests

- [x] Integration test file exists at `tests/integration_tests.rs`
- [x] Test builds cache from fixtures (`fixtures/sv_small.vcf`, `fixtures/trgt_small.vcf`, etc.)
- [x] Test verifies cache files are created with correct names
- [x] Test verifies test_units.parquet has expected schema and row count
- [x] Test verifies genotypes.parquet has expected schema and row count
- [x] Test verifies catalog.parquet has expected schema and row count
- [x] Test verifies features.parquet has expected schema and row count
- [x] Test verifies provenance.json exists and contains expected metadata
- [x] Test runs `storm explain` on a test unit and checks output format
- [ ] Test runs association test using fixtures and verifies results structure
- [x] All integration tests pass with `cargo test --features python`

### F. Notebook Fixes

- [x] `storm.version()` call in notebook works without errors (via CLI)
- [x] Notebook can import storm module successfully (uses subprocess)
- [x] Notebook demonstrates cache building (uncomment and make functional)
- [x] Notebook demonstrates cache loading (uncomment and make functional)
- [x] Notebook demonstrates GLM analysis (uncomment and make functional)
- [x] Notebook demonstrates explain functionality (uncomment and make functional)
- [ ] Notebook runs end-to-end without errors

### G. Documentation and Polish

- [ ] README.md documents CLI usage with examples
- [ ] README.md documents Python API usage with examples
- [ ] README.md includes example of end-to-end workflow
- [ ] Code comments explain the cache building pipeline flow
- [ ] Error messages are user-friendly and actionable

---

## Context

### What's Already Built (Rust Core)

The following Rust modules are complete and functional:
- `src/vcf/sv.rs` - SV VCF parsing
- `src/vcf/trgt.rs` - TRGT VCF parsing  
- `src/catalog/` - TRExplorer catalog ingestion
- `src/mapping.rs` - SV to catalog locus mapping
- `src/testunit.rs` - TestUnit construction
- `src/resolver.rs` - Genotype resolution with overlays
- `src/cache/` - Arrow/Parquet cache structures
- `src/plan.rs` - Plan YAML parsing and rule engine
- `src/encoding.rs` - Genotype encodings (S, M, D, Tail, Categorical)
- `src/glm.rs` - StormGLM association backends
- `src/results.rs` - Results schema and Parquet writing
- `src/explain.rs` - Explainability functions

### What's Missing

1. **CLI Binary** (`src/main.rs` is just `println!("2 + 2 = {}", 2 + 2)`)
2. **End-to-End Pipeline** - Components exist but aren't wired together
3. **Python Bindings** - Only `add()` and `_version()` are exposed, not the core functions
4. **Python API** - `StormCache.build()` and `run_glm()` are stubs with TODOs
5. **Integration Tests** - `tests/integration_tests.rs` is empty
6. **Version Issue** - `storm.version()` fails because Rust extension may not be imported correctly

### Key Files to Modify

- `src/main.rs` - Implement CLI with clap
- `src/lib.rs` - Add `build_cache()` function and expose via PyO3
- `python/storm/__init__.py` - Implement actual functions (remove TODOs)
- `tests/integration_tests.rs` - Write end-to-end tests
- `notebooks/storm_demo.ipynb` - Fix version() and uncomment working examples

### Fixtures Available

- `fixtures/sv_small.vcf` - Small SV VCF for testing
- `fixtures/trgt_small.vcf` - Small TRGT VCF for testing
- `fixtures/trexplorer.bed` - TRExplorer BED catalog
- `fixtures/trexplorer.json` - TRExplorer JSON catalog
- `fixtures/plan.yaml` - Example plan configuration

### Design Principles

- Rust core prioritizes correctness over speed
- Python layer uses Polars for dataframes
- Cache format is Arrow/Parquet for interoperability
- CLI should be intuitive and provide helpful errors
- Python API should feel natural to Python users
- All functions should handle missing optional inputs gracefully

---

## Implementation Notes

### CLI Structure

Use clap for argument parsing. Suggested structure:

```rust
storm cache build --sv-vcf <path> [--trgt-vcf <path>] [--catalog-bed <path>] [--catalog-json <path>] --output-dir <path>
storm cache verify --cache-dir <path>
storm explain <test_unit_id> [--cache-dir <path>] [--sample <sample_id>]
```

### Cache Building Pipeline Flow

1. Parse inputs (SV VCF, optional TRGT VCF, optional catalog)
2. Load and merge catalog if provided
3. Map SVs to catalog loci
4. Build TestUnits (SV, RepeatProxySV, Repeat)
5. Create Resolver and resolve genotypes from SV data
6. Apply TRGT overlay if provided
7. Merge resolved genotypes
8. Compute features (call rate, carriers, etc.)
9. Build ArrowCache with all tables
10. Write Parquet files to output directory
11. Write provenance.json

### Python Bindings Priority

Focus on exposing functions needed for the Python API:
1. `build_cache()` - for `StormCache.build()`
2. `read_cache_from_dir()` - for cache loading
3. `run_association()` - for `run_glm()`
4. `explain_genotype()` / `explain_locus()` - for `explain()`
5. `Plan::from_yaml()` - for plan loading

### Testing Strategy

- Unit tests: Test individual pipeline steps
- Integration tests: Test full cache build from fixtures
- Python tests: Test Python API functions
- Notebook: Manual end-to-end validation

---

## Acceptance Criteria

When complete:
- `cargo build --release` produces a working `storm` binary
- `storm cache build` successfully builds cache from fixtures
- `storm cache verify` validates the built cache
- `storm explain` prints readable genotype explanations
- Python `storm.version()` works
- Python `StormCache.build()` successfully builds cache
- Python `run_glm()` runs association tests
- All integration tests pass
- Notebook runs end-to-end without errors
- README documents usage

---

## Completion Condition

When:
- All checkboxes above are checked
- `cargo test --features python` passes
- `storm cache build` works with fixtures
- `storm explain` works with built cache
- Python API functions are implemented (not stubs)
- Notebook runs successfully
- Integration tests pass

Write **DONE** in `.ralph/progress.md` and stop.
