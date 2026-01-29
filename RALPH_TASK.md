---
task: Real-Data Section — Load All TRGT Files and Add Combination Stats
test_command: "cargo test --features python && cd python && python -m pytest -W error"
---

# Task: Real-Data Section — Load All TRGT Files and Add Combination Stats

The notebook's real-data section currently runs successfully but uses only a subset of TRGT files (e.g. 200) for speed, and it does not clearly show how the integrated callset (BCF) and TRGT calls are combined. This task updates the real-data section to:

1. **Load all TRGT files** for the real-data workflow (no subset cap).
2. **Add combination statistics** so we can see how BCF and TRGT are combined (samples, test units by source, etc.).

---

## Success Criteria

### A. Load All TRGT Files in Real-Data Section

- [x] Real-data cache build uses **all** TRGT files in `scratch/trgt/` (e.g. `trgt_files` with no `N_TRGT_FILES` limit, or `N_TRGT_FILES = len(trgt_files)`).
- [x] Remove or repurpose the "subset for demo speed" logic so the default real-data run is full cohort.
- [x] Optionally: add a short note or variable (e.g. `USE_TRGT_SUBSET = False` or `N_TRGT_FILES = None` meaning "all") so users can still cap TRGT for quick tests if desired.
- [x] Document expected runtime/memory for full TRGT load (~10k files).

### B. Combination Statistics (How BCF + TRGT Are Combined)

- [x] **Notebook:** Add a "Combination stats" subsection (or expand the existing cache-stats section) that clearly reports:
  - **Input counts:** Number of TRGT files loaded; note that BCF is one file with many samples.
  - **Samples:** Total number of unique samples in the cache (union of BCF and TRGT samples).
  - **Test units by source:** Counts by `unit_type` (e.g. Sv vs TrueRepeat), so we see how many loci come from the integrated callset vs TRGT.
  - **Genotypes:** Total genotype records (or a short summary), so we see scale of the combined callset.
- [x] **Backend (optional but recommended):** Extend `build_cache` (and provenance) to expose combination counts so the notebook can display them:
  - **Option 1 (minimal):** Return or write to provenance: `num_samples_sv` (samples in BCF), `num_samples_trgt` (samples in TRGT files), `num_samples_both` (samples present in both). Then the notebook can show "Samples: BCF only N, TRGT only M, both K, total unique T".
  - **Option 2 (current only):** If we do not extend the backend, the notebook still reports: total samples, test units by type, TRGT file count, total genotypes, and a short markdown explanation that "samples = union of BCF and TRGT".
- [x] Combination stats are printed or displayed in the notebook right after cache build (or in the same section).

### C. Documentation and Polish

- [x] Markdown in the notebook explains what "combination" means (union of samples; test units from SV vs TRGT).
- [x] Expected runtime for full TRGT load is noted (e.g. "~5–10 min for ~10k TRGT files" or similar from prior runs).
- [x] All cells still run successfully with real data present; no regressions.

---

## Context

### Current Behavior

- **Notebook:** Real-data section sets `N_TRGT_FILES = 200` and builds cache with `trgt_subset = trgt_files[:N_TRGT_FILES]`, so only 200 TRGT files are used.
- **Cache stats shown:** Build stats (num_test_units, num_samples, num_genotypes) and test units by `unit_type` (e.g. Sv only if no catalog). There is no explicit "samples from BCF vs TRGT" or "how we combine" summary.
- **Backend:** `build_cache` returns `CacheBuildStats { num_test_units, num_samples, num_genotypes, num_catalog_entries }`. It does not return per-source sample counts (BCF-only, TRGT-only, overlap).

### Desired Behavior

- Real-data section uses **all** TRGT files by default (with optional subset for quick tests).
- A clear "Combination stats" block shows: TRGT files loaded, total samples, test units by type (Sv, TrueRepeat), total genotypes, and ideally sample overlap (BCF-only, TRGT-only, both) if the backend is extended.

### Key Files

| File | Change |
|------|--------|
| `notebooks/storm_demo.ipynb` | Use all TRGT files; add combination-stats subsection and optional subset switch. |
| `src/lib.rs` | Optionally extend `CacheBuildStats` and build logic to track/return `num_samples_sv`, `num_samples_trgt`, `num_samples_both`. |
| `src/cache/` (provenance) | Optionally add combination fields to provenance JSON. |
| Python `storm/__init__.py`, `py_build_cache` | Optionally expose new stats from Rust. |

---

## Implementation Plan

### Phase 1: Notebook — Use All TRGT Files

1. In the real-data cache-build cell, stop using a fixed subset (e.g. remove `N_TRGT_FILES = 200` or set it to use all files).
   - Use `trgt_subset = [str(f) for f in trgt_files]` (all files), or introduce `N_TRGT_FILES = None` (meaning all) and `trgt_subset = trgt_files[:N_TRGT_FILES] if N_TRGT_FILES else trgt_files`, then `trgt_subset = [str(f) for f in trgt_subset]`.
2. Update the print message to say "Building cache with N TRGT files" where N = len(trgt_files).
3. Add a short markdown or comment: expected runtime for full TRGT load (~5–10 min) and optional note on how to cap TRGT for a quick test (e.g. set `N_TRGT_FILES = 100`).

### Phase 2: Combination Statistics

4. **Notebook:** Add a "Combination stats" subsection after cache build (or merge into existing stats):
   - TRGT files loaded: `len(trgt_files)` (or len of list actually passed to build).
   - Total samples: from build stats or `real_cache.genotypes["sample_id"].n_unique()`.
   - Test units by type: `real_cache.test_units.group_by("unit_type").agg(pl.len().alias("count"))` (or equivalent).
   - Total genotypes: from build stats or `len(real_cache.genotypes)`.
5. **Backend (optional):** In `build_cache`, after merging `all_samples`, compute:
   - `num_samples_sv` = sv_samples.len()
   - `num_samples_trgt` = number of unique samples that came from TRGT (from parse_and_merge_trgt_vcfs or by tracking which IDs were added from TRGT).
   - `num_samples_both` = number of sample IDs that appear in both sv_samples and the TRGT sample set.
   Add these to `CacheBuildStats` and to provenance JSON. Expose in Python `py_build_cache` return value / `_build_stats` so the notebook can print "Samples: BCF only X, TRGT only Y, both Z, total T".

### Phase 3: Polish

6. Add markdown explaining that the cache combines the integrated callset (BCF) and TRGT by taking the union of samples and building test units from both SV records and TRGT loci.
7. Verify all cells run successfully with real data.

---

## Acceptance Criteria

When complete:

- Real-data section builds the cache using **all** TRGT files in `scratch/trgt/` by default.
- Combination statistics are visible in the notebook (TRGT file count, total samples, test units by type, total genotypes; optionally BCF-only / TRGT-only / both sample counts).
- Notebook still runs to completion with real data; no regressions.

---

## Completion Condition

When:

- All checkboxes above are checked.
- `cargo test --features python` passes.
- `cd python && python -m pytest -W error` passes.
- Real-data section uses all TRGT files and displays combination stats.

Write **DONE** in `.ralph/progress.md` and stop.
