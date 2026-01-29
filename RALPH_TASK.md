---
task: Real-Data Association Testing Workflow in Jupyter Notebook
test_command: "cargo test --features python && cd python && python -m pytest -W error"
---

# Task: Real-Data Association Testing Workflow in Jupyter Notebook

The infrastructure to read real data (BCF, gzip VCF, multiple TRGT files) is complete. Now we need to update the Jupyter notebook to demonstrate a full association testing workflow using the real development data.

**Real data:**
- SV BCF: `scratch/chr6.31803187_32050925.bcf` (~12,680 samples)
- TRGT: `scratch/trgt/*.vcf.gz` (~10,003 files, one sample per file)
- Sample overlap: ~10,000 samples present in both

**Goal:** The notebook should show a complete, reproducible association testing workflow that runs on real data and produces meaningful statistical results (not the "Insufficient samples" errors from the 2-sample fixtures).

---

## Success Criteria

### A. Build Cache from Real Data

- [x] Notebook has a section that builds a cache from the real BCF
- [x] Uses a manageable subset of TRGT files (e.g., 100–500 files) for demo speed
- [x] Documents how to scale up to the full ~10k TRGT files
- [x] Cache build completes in reasonable time (document expected duration)
- [x] Shows cache statistics (samples, test units, genotypes)

### B. Simulated Phenotype for Association Testing

- [x] Generate a simulated continuous phenotype for the samples in the cache
- [x] Generate a simulated binary phenotype (case/control) 
- [x] Phenotypes are reproducible (seeded random generation)
- [x] Optionally: show how to load phenotypes from an external file (CSV/Parquet)

### C. Run Association Tests with Real Data

- [x] Run `run_glm()` on the real-data cache with the simulated phenotype
- [x] Results show actual beta, SE, and p-values (not null/error)
- [x] Demonstrate linear regression (continuous phenotype)
- [x] Demonstrate logistic regression (binary phenotype)
- [x] Show how to use a plan YAML for model/encoding selection
- [x] Show covariate adjustment (at minimum: simulated age/sex or PCs)

### D. Interpret and Visualize Results

- [x] Display association results as a Polars DataFrame
- [x] Filter/sort results by p-value to show top hits
- [x] Create a simple visualization (e.g., Manhattan-style plot or p-value histogram) if matplotlib/seaborn available
- [x] Explain what the results mean (encoding, model, beta interpretation)

### E. Performance and Practical Guidance

- [x] Document expected runtime for different TRGT subset sizes (100, 500, all)
- [x] Note memory considerations for large sample counts
- [x] Provide guidance on running the full pipeline outside the notebook (CLI)

### F. Cleanup and Polish

- [x] Keep the existing fixture-based demo sections (or clearly separate them)
- [x] Notebook runs end-to-end without errors when real data is present
- [x] Graceful handling if `scratch/` data is not available (skip or use fixtures)
- [x] Clear section headers and explanatory markdown
- [x] Cleanup temporary cache directories at the end

---

## Context

### Current Notebook State

The notebook (`notebooks/storm_demo.ipynb`) currently:
- Uses `fixtures/sv_small.vcf` and `fixtures/trgt_small.vcf` (2 samples)
- Shows "Insufficient samples for regression" in association results
- Demonstrates cache build, load, explain, and verify
- Uses CLI examples that also use fixtures

### What Needs to Change

1. **Add real-data sections** after the fixture demo (or replace it):
   - Build cache from BCF + subset of TRGT `.vcf.gz`
   - Generate phenotypes for the samples
   - Run association testing
   - Show results

2. **Phenotype generation:**
   - Get sample IDs from the cache
   - Create simulated phenotypes aligned to those samples
   - Pass to `run_glm()`

3. **Results interpretation:**
   - The current notebook shows results but with nulls
   - With real data, we should have actual statistics

### Key Files

| File | Change |
|------|--------|
| `notebooks/storm_demo.ipynb` | Add real-data association workflow sections |
| `fixtures/plan.yaml` | May need adjustment for real data (motif patterns, thresholds) |

### API Reference

```python
# Build cache with multiple TRGT files
cache = storm.StormCache.build(
    sv_vcf="scratch/chr6.31803187_32050925.bcf",
    trgt_vcf=["scratch/trgt/file1.vcf.gz", "scratch/trgt/file2.vcf.gz", ...],
    output_dir="real_cache",
)

# Or use glob to get TRGT file list
from pathlib import Path
trgt_files = sorted(Path("scratch/trgt").glob("*.vcf.gz"))[:100]  # first 100
cache = storm.StormCache.build(
    sv_vcf="scratch/chr6.31803187_32050925.bcf",
    trgt_vcf=[str(f) for f in trgt_files],
    output_dir="real_cache",
)

# Get sample IDs from cache
sample_ids = cache.genotypes["sample_id"].unique().to_list()

# Create phenotype (Polars Series or dict)
import polars as pl
import numpy as np
np.random.seed(42)
phenotype = pl.Series("phenotype", np.random.randn(len(sample_ids)))

# Run association
results = storm.run_glm(
    cache=cache,
    phenotype=phenotype,
    model="linear",       # or "logistic" for binary
    encoding="S",         # S, M, D, or binary
)
```

---

## Implementation Plan

### Phase 1: Real-Data Cache Build Section

1. Add markdown explaining the real data and where it comes from
2. Check if `scratch/` exists; if not, note that fixtures will be used
3. Build cache from BCF + subset of TRGT (e.g., 100 files for quick demo)
4. Display cache statistics and sample count

### Phase 2: Phenotype Generation Section

5. Extract unique sample IDs from the cache
6. Generate a simulated continuous phenotype (seeded numpy random)
7. Generate a simulated binary phenotype (0/1 with ~50% prevalence)
8. Optionally show loading from CSV

### Phase 3: Association Testing Section

9. Run `run_glm()` with continuous phenotype and linear model
10. Run `run_glm()` with binary phenotype and logistic model
11. Show results DataFrame with actual statistics
12. Demonstrate plan-based model selection
13. Demonstrate covariate adjustment

### Phase 4: Results Interpretation

14. Sort results by p-value, show top 10
15. Explain columns (beta, se, p_value, encoding, model)
16. Optional: simple matplotlib histogram of p-values or -log10(p) plot
17. Discuss what a "hit" would mean in this context

### Phase 5: Polish

18. Ensure notebook runs cleanly with and without real data
19. Add cleanup cell at end
20. Review and improve markdown explanations
21. Test end-to-end execution

---

## Acceptance Criteria

When complete:

- Notebook demonstrates full association testing on real data
- Results show actual beta, SE, and p-values (not null/Insufficient samples)
- Both linear and logistic regression are demonstrated
- Sample sizes are meaningful (100+ samples minimum)
- Notebook runs without errors when `scratch/` data is present
- Notebook gracefully handles missing `scratch/` data

---

## Completion Condition

When:

- All checkboxes above are checked
- `cargo test --features python` passes
- `cd python && python -m pytest -W error` passes
- Notebook runs end-to-end and displays association results with real statistics
- Notebook is saved with output cells showing the workflow

Write **DONE** in `.ralph/progress.md` and stop.
