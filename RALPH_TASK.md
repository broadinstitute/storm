---
task: Fix Bugs, Complete Association Testing, and Polish STORM
test_command: "cargo test --features python && cd python && python -m pytest -W error"
---

# Task: Fix Bugs, Complete Association Testing, and Polish STORM

The STORM integration layer is complete, but there are several bugs and incomplete implementations that need attention. This task focuses on:
1. Fixing the Polars Series conversion bug in `run_glm()`
2. Implementing actual association testing in `py_run_association()`
3. Adding missing test dependencies
4. Improving error handling and test coverage
5. General polish and robustness improvements

---

## Success Criteria

### A. Fix Python API Bugs

- [ ] Fix Polars Series conversion in `run_glm()` (line 206 checks `to_dict` but should check `to_list`)
- [ ] `run_glm()` handles Polars Series correctly without warnings
- [ ] `run_glm()` handles dict-like phenotypes correctly
- [ ] `run_glm()` handles list/array phenotypes correctly
- [ ] All Python tests pass without warnings (`pytest -W error`)
- [ ] Error messages are clear when phenotype format is incorrect

### B. Implement Real Association Testing

- [ ] `py_run_association()` calls actual `run_association()` from `glm.rs` module
- [ ] `py_run_association()` loads test units from cache
- [ ] `py_run_association()` loads genotypes from cache
- [ ] `py_run_association()` extracts encoded predictor values for each test unit
- [ ] `py_run_association()` applies plan rules to select encoding/model
- [ ] `py_run_association()` runs association test for each test unit
- [ ] `py_run_association()` returns real statistics (beta, SE, p-value, etc.)
- [ ] `py_run_association()` handles covariates if provided
- [ ] `py_run_association()` respects plan configuration (rare-variant ladder, etc.)
- [ ] Results include proper metadata (encoding, model, rule_id, etc.)

### C. Test Infrastructure

- [ ] Add `pytest` to `dev-requirements.txt` or `python/pyproject.toml`
- [ ] Add `polars` to test dependencies (if not already present)
- [ ] All Python tests pass with `pytest -W error` (no warnings)
- [ ] Add test for error cases (invalid phenotype format, missing cache, etc.)
- [ ] Add test for association testing with real data
- [ ] Add test for plan-based model selection
- [ ] Add test for covariates handling
- [ ] CI runs Python tests automatically (check `.github/workflows/ci.yml`)

### D. Error Handling and Robustness

- [ ] `run_glm()` provides clear error message when cache is missing required tables
- [ ] `run_glm()` provides clear error message when phenotype length doesn't match samples
- [ ] `run_glm()` handles missing sample IDs gracefully
- [ ] `py_run_association()` handles empty cache gracefully
- [ ] `py_run_association()` handles invalid plan YAML gracefully
- [ ] All error messages are user-friendly and actionable
- [ ] Error handling is consistent across Python API functions

### E. Code Quality and Documentation

- [ ] Remove placeholder comments in `py_run_association()` (line 862-863)
- [ ] Add docstrings explaining phenotype format requirements
- [ ] Add docstrings explaining plan-based model selection
- [ ] Add examples in docstrings for common use cases
- [ ] Update README with association testing examples
- [ ] Document expected phenotype formats (Polars Series, dict, list)
- [ ] Code follows consistent style (run formatters if needed)

### F. Verification and Testing

- [ ] Run full test suite: `cargo test --features python && cd python && python -m pytest -W error`
- [ ] Verify Python extension builds correctly with `maturin develop`
- [ ] Verify notebook runs without warnings
- [ ] Verify CLI commands work end-to-end
- [ ] Verify association results are reasonable (not all zeros/ones)
- [ ] Verify plan-based model selection works correctly

---

## Context

### Current Status

**What Works:**
- CLI binary (`storm cache build`, `storm cache verify`, `storm explain`)
- End-to-end cache building pipeline
- Python bindings for most functions
- Python API for cache building and loading
- Integration tests (Rust)
- Python tests (16 tests, but 1 warning)

**What's Broken/Incomplete:**

1. **Polars Series Bug** (`python/storm/__init__.py:206`)
   - Checks `hasattr(phenotype, 'to_dict')` but Polars Series has `to_list()`, not `to_dict()`
   - Falls through to `dict(phenotype)` which fails
   - Causes warning in tests: "cannot convert dictionary update sequence element #0 to a sequence"

2. **Placeholder Association Testing** (`src/lib.rs:862-863`)
   - `py_run_association()` returns dummy results (all zeros/ones)
   - Comment says "Full implementation would call run_association from glm module"
   - Needs to actually call `glm::run_association()` with proper data extraction

3. **Missing Test Dependencies**
   - `pytest` not in `dev-requirements.txt`
   - Tests work but dependency should be documented

### Key Files to Modify

- `python/storm/__init__.py` - Fix Polars Series conversion (line 206)
- `src/lib.rs` - Implement real association testing (lines 837-885)
- `dev-requirements.txt` or `python/pyproject.toml` - Add pytest
- `python/tests/test_storm.py` - Add error case tests
- `README.md` - Add association testing examples

### Implementation Notes

#### Fixing Polars Series Conversion

Change line 206 from:
```python
if hasattr(phenotype, 'to_dict'):
```

To:
```python
if isinstance(phenotype, pl.Series) or hasattr(phenotype, 'to_list'):
```

Or more robustly, check for Polars Series first:
```python
if HAS_POLARS and isinstance(phenotype, pl.Series):
    # Handle Polars Series
elif hasattr(phenotype, 'to_list'):
    # Handle other list-like objects
else:
    # Handle dict-like objects
```

#### Implementing Real Association Testing

The `py_run_association()` function needs to:

1. Load cache and extract test units
2. For each test unit:
   - Load genotypes from cache
   - Extract encoded predictor values (using plan rules)
   - Match phenotype values to samples
   - Call `glm::run_association()` with proper parameters
   - Collect results
3. Return results as JSON

Key challenge: Extracting encoded values from cache and matching to phenotype samples.

The `glm::run_association()` signature is:
```rust
pub fn run_association(
    unit_id: &str,
    predictor: &[f64],
    outcome: &[f64],
    model: &Model,
    encoding: &Encoding,
    covariates: Option<&Covariates>,
) -> Result<AssociationResult, GlmError>
```

So we need to:
- Extract predictor values from cache (features table or compute from genotypes)
- Extract outcome values matching sample order
- Select model/encoding from plan
- Build covariates if provided
- Call `run_association()` for each unit

### Testing Strategy

1. **Unit tests** for phenotype conversion logic
2. **Integration tests** for association testing with fixtures
3. **Error case tests** for invalid inputs
4. **Regression tests** to ensure no warnings

---

## Acceptance Criteria

When complete:
- All Python tests pass with `pytest -W error` (zero warnings)
- `run_glm()` handles Polars Series correctly
- `py_run_association()` returns real association results (not placeholders)
- Association results include proper statistics and metadata
- Error messages are clear and helpful
- README includes association testing examples
- All tests documented in dev-requirements.txt

---

## Completion Condition

When:
- All checkboxes above are checked
- `cargo test --features python` passes
- `cd python && python -m pytest -W error` passes (zero warnings)
- `run_glm()` works correctly with Polars Series
- `py_run_association()` returns real association results
- All error cases handled gracefully
- Documentation updated

Write **DONE** in `.ralph/progress.md` and stop.
