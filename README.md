# STORM: Structural & Tandem-Repeat Optimized Regression Models

A high-performance Rust library with Python bindings for association testing of structural variants and tandem repeats using long-read sequencing data.

## Features

- **VCF Parsing**: Parse integrated SV VCFs and TRGT VCFs
- **Catalog Integration**: Load TRExplorer BED/JSON catalogs for tandem repeat annotations
- **Genotype Resolution**: Map SVs to repeat loci and reconstruct diploid allele lengths
- **Multiple Encodings**: S (sum), M (max), D (diff), Tail, Binary, Categorical
- **Association Testing**: Linear regression, logistic regression, BinomiRare, Firth regression
- **Caching**: Arrow/Parquet-based cache for efficient data storage and access
- **Explainability**: Detailed genotype explanations for debugging and QC

## Installation

### From Source (Rust CLI)

```bash
# Clone the repository
git clone https://github.com/broadinstitute/storm.git
cd storm

# Build the release binary
cargo build --release

# The binary will be available at target/release/storm
```

### Python Package

```bash
# Install maturin for building Python extensions
pip install maturin

# Build and install the Python package
maturin develop --features python
```

## CLI Usage

### Building a Cache

Build a cache from input VCF files:

```bash
storm cache build \
    --sv-vcf input/sv.vcf.gz \
    --trgt-vcf input/trgt.vcf.gz \
    --catalog-bed catalog/trexplorer.bed \
    --catalog-json catalog/trexplorer.json \
    --output-dir storm_cache
```

Options:
- `--sv-vcf`: Path to integrated SV VCF (required)
- `--trgt-vcf`: Path to TRGT VCF (optional)
- `--catalog-bed`: Path to TRExplorer BED file (optional)
- `--catalog-json`: Path to TRExplorer JSON file (optional)
- `--output-dir`: Output directory for cache (default: `storm_cache`)

### Verifying a Cache

Validate that a cache has the correct schema and data:

```bash
storm cache verify --cache-dir storm_cache
```

### Explaining Genotypes

Get detailed information about genotypes at a locus:

```bash
# Explain all genotypes at a locus
storm explain sv_sv1 --cache-dir storm_cache

# Explain a specific sample
storm explain sv_sv1 --cache-dir storm_cache --sample SAMPLE1
```

## Python API

### Building a Cache

```python
import storm

# Build cache using the Rust backend
cache = storm.StormCache.build(
    sv_vcf="input/sv.vcf.gz",
    trgt_vcf="input/trgt.vcf.gz",
    catalog_bed="catalog/trexplorer.bed",
    catalog_json="catalog/trexplorer.json",
    output_dir="storm_cache"
)
```

### Loading and Exploring a Cache

```python
import storm
import polars as pl

# Load existing cache
cache = storm.load_cache("storm_cache")

# Access tables as Polars DataFrames
print(cache.test_units)
print(cache.genotypes)
print(cache.catalog)
print(cache.features)
```

### Running Association Tests

```python
import storm

# Run GLM association testing
results = storm.run_glm(
    cache=cache,
    phenotype=phenotype_series,
    plan="plan.yaml",
    covariates=covariates_df,
    output="results.parquet"
)

print(results.sort("p_value").head(10))
```

### Explaining Genotypes

```python
import storm

# Explain all genotypes at a locus
explanation = storm.explain(cache, "repeat_TR001")
print(explanation)

# Explain for a specific sample
explanation = storm.explain(cache, "repeat_TR001", sample_id="SAMPLE1")
print(explanation)
```

## Plan Configuration

STORM uses YAML plan files to configure analysis:

```yaml
name: my_analysis
default_encoding: s
default_model: linear

rules:
  - name: rare_binary
    condition:
      max_carriers: 10
    encoding: binary
    model: binomirare
    
  - name: cag_repeats  
    condition:
      motif_pattern: CAG
    encoding: m
    model: logistic

covariates:
  - age
  - sex
  
num_pcs: 10
```

## Cache Format

The cache consists of Parquet files:

| File | Description |
|------|-------------|
| `test_units.parquet` | Test unit definitions (SV, RepeatProxy, TrueRepeat) |
| `genotypes.parquet` | Resolved genotypes per sample |
| `catalog.parquet` | TRExplorer catalog entries |
| `features.parquet` | Computed features (S, M, D, binary) |
| `provenance.parquet` | Build metadata |
| `provenance.json` | Human-readable provenance |

## Encodings

STORM supports multiple genotype encodings:

| Encoding | Description | Formula |
|----------|-------------|---------|
| S (Sum) | Sum of allele lengths | L1 + L2 |
| M (Max) | Maximum allele length | max(L1, L2) |
| D (Diff) | Absolute difference | \|L1 - L2\| |
| Tail | Threshold indicator | 1 if max(L1,L2) > threshold |
| Binary | Carrier status | 1 if has non-ref allele |
| Categorical | Binned categories | Bin by length ranges |

## Development

### Running Tests

```bash
# Run all tests
cargo test

# Run with Python feature
cargo test --features python

# Run integration tests only
cargo test --test integration_tests
```

### Building Documentation

```bash
cd docs
make html
```

## License

MIT License - see LICENSE file for details.

## Citation

If you use STORM in your research, please cite:

```
[Citation information]
```
