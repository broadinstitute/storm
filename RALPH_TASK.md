---
task: Support Real Data — Integrated SV BCF and Per-Sample TRGT VCF.gz
test_command: "cargo test --features python && cd python && python -m pytest -W error"
---

# Task: Support Real Data — Integrated SV BCF and Per-Sample TRGT VCF.gz

The codebase currently works with small text fixtures (plain VCF, one SV file, one TRGT file with few samples). Real development data is:

- **Integrated SV BCF** (~13k samples): `scratch/chr6.31803187_32050925.bcf`
- **Per-sample TRGT VCF.gz** (~10k files): `scratch/trgt/*.vcf.gz` (one sample per file)

This task is to identify and plan all changes needed so we can build caches and run association workflows on this real data. Use a couple of real TRGT files (and the BCF header) to inform the plan; do not assume the full ~10k TRGT set is loaded at once in early implementation.

---

## Success Criteria

### A. BCF input for SV calls

- [x] SV pipeline can read BCF (binary) as well as plain VCF
- [x] BCF path is accepted by CLI and Python `build_cache` (no format rename required)
- [x] Parsed SV records (SVTYPE, SVLEN, END, GT) are equivalent whether source is VCF or BCF
- [x] Existing fixture-based tests still pass (plain VCF unchanged)
- [x] Optional: document or test with `scratch/chr6.31803187_32050925.bcf` (manual or CI if available)

### B. Compressed VCF input (gzip / BGZF)

- [ ] TRGT parser can read `.vcf.gz` (gzip or BGZF) in addition to plain `.vcf`
- [ ] SV parser can read `.vcf.gz` if we ever pass compressed VCF for SV
- [ ] Detection by extension or magic; no requirement to rename files
- [ ] Fixture tests still pass with plain VCF; add test(s) with a small compressed VCF if feasible

### C. Multiple TRGT files (per-sample VCFs)

- [ ] Cache build accepts **multiple** TRGT paths (not just one)
- [ ] CLI: accept multiple `--trgt-vcf` and/or `--trgt-dir` (with glob) and/or `--trgt-list` (file of paths)
- [ ] Lib API: `build_cache` (and Python) accept a list/slice of TRGT paths
- [ ] Merge strategy: one sample per file; merge by locus (TRID/chrom/pos/end) so each `TrgtRecord` has genotypes from all samples across the provided files
- [ ] Sample IDs: taken from VCF header (real files use numeric ID, e.g. `1524337`); no change to parser contract
- [ ] Works with a small subset (e.g. 2–3 TRGT files) for tests and incremental validation

### D. Sample alignment and catalog

- [ ] Document that BCF sample IDs are numeric (e.g. `1000234`) and TRGT sample ID is the single sample in each file (same numeric ID in real data)
- [ ] Joint SV + TRGT: same sample set or explicit subset; document ordering/subset requirements
- [ ] Catalog: real region uses chr6:31803187–32050925; catalog format (BED/JSON) unchanged; TRID in real TRGT is coordinate-style (e.g. `6-31803583-31803598-T`) — code treats TRID as opaque string; ensure catalog matching or test-unit ID logic does not assume fixture-style IDs (e.g. `TR001`)

### E. Performance and scale (planning only)

- [ ] Plan documents expected behavior: BCF ~13k samples, TRGT ~10k files; no requirement to load all 10k TRGT in memory at once in v1 if we can specify a subset
- [ ] Consider: parallel parsing of TRGT files, streaming/blocked BCF reading, or chunked merge — captured in plan; implementation can be incremental

### F. Documentation and testing

- [ ] Update `DEV_DATA.md` (or equivalent) with: paths, format summary, and how to run cache build with real data (BCF + subset of TRGT)
- [ ] RALPH_TASK.md contains the full development plan (this file)
- [ ] Unit/integration tests still pass; add tests for: multiple TRGT paths, merge-by-locus, and (if feasible) one BCF and one .vcf.gz in tests

---

## Context

### Current behavior

- **SV**: `parse_sv_vcf(path)` in `src/vcf/sv.rs` uses `File::open` + `BufReader` + `reader.lines()`. **Plain text only**; BCF is binary and will fail if opened as text.
- **TRGT**: `parse_trgt_vcf(path)` in `src/vcf/trgt.rs` uses the same pattern. **Plain text only**; `.vcf.gz` is gzip-compressed and will fail.
- **Cache build**: `build_cache(sv_vcf_path, trgt_vcf_path, ...)` in `src/lib.rs` takes **one** optional `trgt_vcf_path`. It parses one SV file and optionally one TRGT file; TRGT records are expected to have multiple samples in `genotypes` (fixture: 2 samples in one file). Real data: **one sample per TRGT file**, ~10k files.
- **CLI**: `storm cache build --sv-vcf <path> [--trgt-vcf <path>]` — single TRGT path only.
- **Catalog**: BED + JSON (e.g. `fixtures/trexplorer.bed`, `fixtures/trexplorer.json`). Fixture uses chr1/chr2/chr3 and IDs like `TR001`. Real region is chr6; TRGT INFO TRID is like `6-31803583-31803598-T`. Parser and resolver use TRID as string; no code change for ID format.

### Real data (from FIXTURE_COMPARISON.md and inspection)

- **SV BCF** (`scratch/chr6.31803187_32050925.bcf`): ~12,680 samples; format `GT:FT:SQ:GQ:PS:NE:DP:AD:KS`; INFO includes SVTYPE, SVLEN, END. Parser only needs GT + SVTYPE/SVLEN/END — compatible once we can read BCF.
- **TRGT** (`scratch/trgt/*.vcf.gz`): One sample per file; header sample name is numeric (e.g. `1524337`); format `GT:AL:ALLR:SD:MC:MS:AP:AM`; INFO TRID, END, MOTIFS, STRUC. Parser uses GT and AL — compatible; TRID format differs from fixture but is still a string.

### Key files to touch

| Area | Files | Change |
|------|--------|--------|
| SV input | `src/vcf/sv.rs` | Add BCF reading path (e.g. noodles-bcf) or reader abstraction that dispatches on format; keep existing text VCF path. |
| TRGT input | `src/vcf/trgt.rs` | Add decompression for `.vcf.gz` (e.g. flate2 or noodles bgzf); keep existing text path. |
| Multi-TRGT | `src/lib.rs` | `build_cache`: accept `Option<&[String]>` or `Vec<&str>` for TRGT paths; loop over paths, parse each, merge TrgtRecords by locus (trid/chrom/pos/end), then run existing TrueRepeat/test-unit logic. |
| CLI | `src/main.rs` | Add multiple `--trgt-vcf` and/or `--trgt-dir` / `--trgt-list`; pass list into `build_cache`. |
| Python | `python/storm/__init__.py`, `src/lib.rs` (py_build_cache) | Extend `build_cache(sv_vcf, trgt_vcf=None, ...)` to accept `trgt_vcf: Optional[Union[str, List[str]]]`; pass list to Rust. |
| Tests | `src/vcf/sv.rs`, `src/vcf/trgt.rs`, `tests/integration_tests.rs`, `python/tests/test_storm.py` | Add tests for BCF (or stub), gzip VCF, multiple TRGT merge; keep fixture tests. |
| Docs | `DEV_DATA.md`, `README.md` | Document real data paths, BCF + multi-TRGT usage, sample alignment. |

---

## Development Plan

### Phase 1: Input format support (BCF + gzip)

1. **BCF for SV**
   - Add noodles BCF support: enable `bcf` feature in `Cargo.toml` for noodles (or use separate `noodles-bcf` if versioning differs). Implement a BCF reader path that yields the same logical data as `parse_sv_vcf`: sample names from header; for each record, chrom, pos, id, ref, alt, INFO SVTYPE/SVLEN/END, and per-sample GT. Keep `parse_sv_vcf` for plain VCF; add `read_sv_vcf_or_bcf(path)` (or internal helper) that chooses by extension or magic (e.g. first bytes) and calls either the existing parser (on a `Read` impl) or the BCF path. Ensure `build_cache` uses this so a single `sv_vcf` path can be BCF or VCF.
   - Tests: existing `parse_sv_vcf` tests unchanged; add test that BCF reader (or a small BCF fixture) produces expected SV records if feasible; otherwise document manual test with `scratch/...bcf`.

2. **Gzip/BGZF for VCF**
   - In `sv.rs` and `trgt.rs`, support reading from a decompressing reader when the path ends in `.gz` or the stream is BGZF/gzip. Options: (a) `flate2::read::GzDecoder` for plain gzip, (b) noodles `bgzf` for BGZF. Prefer one abstraction (e.g. `open_vcf_read(path) -> impl BufRead`) used by both SV and TRGT parsers so that `parse_sv_vcf` and `parse_trgt_vcf` accept `.vcf.gz` without changing their function signatures (they take a path). Refactor current `File::open` + `BufReader` to use this reader so plain `.vcf` still works.
   - Tests: add a small `.vcf.gz` fixture (or generate one from existing fixture) and test that parsing yields the same records as plain VCF.

### Phase 2: Multiple TRGT paths and merge

3. **Lib API: multiple TRGT paths**
   - Change `build_cache` to accept `trgt_vcf_paths: Option<&[impl AsRef<str>]>` (or `Option<Vec<&str>>`) instead of `trgt_vcf_path: Option<&str>`. Call sites: `src/main.rs`, `src/lib.rs` (py_build_cache), Python `StormCache.build(..., trgt_vcf=...)`.
   - Single path: treat as one-element list for backward compatibility.

4. **Merge logic**
   - For each path in `trgt_vcf_paths`: call `parse_trgt_vcf(path)` (which now supports .gz), get `(sample_names, records)`. Real data: one sample per file, so `sample_names.len() == 1`. For each `TrgtRecord` in `records`, use a canonical locus key (e.g. `(rec.chrom.clone(), rec.pos, rec.end, rec.trid.clone())`) and either insert a new entry into a `HashMap<LocusKey, TrgtRecord>` or merge the single sample’s genotype into the existing record’s `genotypes` map. After processing all files, collect `TrgtRecord`s from the map (values). Use this merged list as the current “trgt_records” in the rest of `build_cache` (test units, resolution, features, cache write). Ensure `all_samples` includes every sample ID seen across TRGT files (and SV) and is sorted as today.

5. **CLI**
   - `storm cache build`: add `--trgt-vcf <path>` (multiple) and/or `--trgt-dir <dir>` (expand to `dir/*.vcf.gz` or `dir/*.vcf`) and/or `--trgt-list <file>` (one path per line). Pass the collected list of TRGT paths to `build_cache`.

6. **Python**
   - `StormCache.build(..., trgt_vcf=None, ...)`: allow `trgt_vcf` to be `str | List[str] | None`. If `str`, wrap in a one-element list; if `List[str]`, use as-is. Pass to Rust `py_build_cache(sv_vcf, trgt_vcf_paths, ...)`. Extend Rust `py_build_cache` to accept a list (e.g. `Vec<String>` from Python).

### Phase 3: Catalog and sample alignment

7. **Catalog and TRID**
   - No change to catalog format. For real chr6 region, user may provide a BED/JSON with entries overlapping chr6:31803187–32050925; IDs can be coordinate-based or symbolic. Code already treats TRID as opaque in TRGT; test units from TRGT use `rec.trid` as unit id. Document that real TRGT TRIDs look like `6-31803583-31803598-T`.

8. **Sample alignment**
   - Document in `DEV_DATA.md`: BCF sample order/IDs (numeric); TRGT one sample per file, same numeric ID in header; for joint SV+TRGT, use the intersection or a defined subset and consistent ordering. No code change required if we merge TRGT by locus and use `all_samples` as today.

### Phase 4: Tests and documentation

9. **Tests**
   - Unit: multiple TRGT paths → merge produces one record per locus with N samples’ genotypes when N files each have one sample; fixture with 2 samples in one file still works (one path, two samples in that file).
   - Integration: `build_cache` with two TRGT fixture files (e.g. split fixture into two single-sample VCFs) and verify cache contents (test units, genotype count). Optional: BCF fixture or skip BCF in CI and document manual test; optional: small .vcf.gz in tests.
   - Python: `build_cache` with `trgt_vcf=[path1, path2]` and with `trgt_vcf=path1`; both succeed and produce consistent caches.

10. **Documentation**
    - Update `DEV_DATA.md`: paths to BCF and TRGT dir; example commands for building cache with BCF and a subset of TRGT (e.g. first 100 files); note on sample overlap and ordering.
    - README: mention BCF and .vcf.gz support and multiple TRGT files where the “Usage” or “Input” section is described.

---

## Acceptance Criteria

When complete:

- Cache can be built from `scratch/chr6.31803187_32050925.bcf` (SV) and a subset of `scratch/trgt/*.vcf.gz` (e.g. 2–10 files) without converting BCF to VCF or decompressing TRGT on disk.
- CLI and Python accept multiple TRGT paths; merge-by-locus is correct (one TrgtRecord per locus, genotypes from all provided files).
- All existing fixture-based tests pass; new tests cover multi-TRGT merge and (if feasible) compressed VCF and BCF.
- `DEV_DATA.md` (or equivalent) describes real data and how to run the pipeline; RALPH_TASK.md contains this plan.

---

## Completion Condition

When:

- All checkboxes in Success Criteria are checked.
- `cargo test --features python` passes.
- `cd python && python -m pytest -W error` passes.
- Cache build works with real BCF and a small set of real TRGT .vcf.gz files (validated manually or in CI).
- Documentation is updated.

Write **DONE** in `.ralph/progress.md` and stop.
