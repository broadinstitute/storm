# STORM: Structural & Tandem-Repeat Optimized Regression Models

A high-performance framework for association testing of structural variants and tandem repeats using long-read sequencing data.

## Features

- **SV/TRGT Integration**: Parse and integrate structural variants with TRGT repeat genotypes
- **Genotype Resolution**: Resolve diploid allele lengths from multiple data sources
- **Multiple Encodings**: S (sum), M (max), D (diff), Tail, Categorical, Binary
- **Association Testing**: Linear, logistic, Firth, and BinomiRare models
- **Arrow/Parquet Cache**: Efficient storage for large-scale analysis
- **Plan-based Configuration**: YAML-driven analysis rules

## Installation

### Rust Binary

```bash
# Clone and build
git clone https://github.com/broadinstitute/storm.git
cd storm
cargo build --release

# The binary will be at target/release/storm
```

### Python Package

```bash
# Create virtual environment
python -m venv venv
source venv/bin/activate

# Install with maturin
pip install maturin
maturin develop
```

## CLI Usage

### Build a Cache

```bash
storm cache build \
    --sv-vcf integrated_sv.vcf \
    --trgt-vcf trgt.vcf \
    --catalog-bed trexplorer.bed \
    --catalog-json trexplorer.json \
    --output-dir my_cache
```

### Verify a Cache

```bash
storm cache verify --cache-dir my_cache
```

### Explain Genotypes

```bash
# Explain all samples at a locus
storm explain TR_123 --cache-dir my_cache

# Explain a specific sample
storm explain TR_123 --cache-dir my_cache --sample SAMPLE1
```

## Python API

### Build Cache

```python
import storm

# Build cache from input files
cache = storm.StormCache.build(
    sv_vcf="integrated_sv.vcf",
    trgt_vcf="trgt.vcf",
    catalog_bed="trexplorer.bed",
    catalog_json="trexplorer.json",
    output_dir="my_cache",
)

print(f"Test units: {cache._build_stats['num_test_units']}")
print(f"Samples: {cache._build_stats['num_samples']}")
```

### Load and Explore Cache

```python
# Load existing cache
cache = storm.load_cache("my_cache")

# Access tables as Polars DataFrames
print(cache.test_units)
print(cache.genotypes)
print(cache.features)
print(cache.catalog)
```

### Explain Genotypes

```python
# Explain all samples at a locus
explanation = storm.explain(cache, "TR_123")
print(explanation)

# Explain a specific sample
explanation = storm.explain(cache, "TR_123", sample_id="SAMPLE1")
print(explanation)
```

### Run Association Testing

```python
import polars as pl

# Create phenotype series (sample order matches cache)
phenotype = pl.Series("phenotype", [0.5, 1.2, 0.8, ...])

# Run GLM analysis
results = storm.run_glm(
    cache=cache,
    phenotype=phenotype,
)

# Filter significant results
significant = results.filter(pl.col("p_value") < 5e-8)
print(significant)
```

### Verify Cache

```python
result = storm.verify_cache("my_cache")
print(f"Valid: {result['is_valid']}")
print(f"Test units: {result['num_test_units']}")
print(f"Genotypes: {result['num_genotypes']}")
```

### Load Plan Configuration

```python
# Load plan YAML
plan = storm.load_plan("analysis_plan.yaml")
print(plan)
```

## End-to-End Workflow

```python
import storm
import polars as pl

# 1. Build cache from input files
cache = storm.StormCache.build(
    sv_vcf="data/sv.vcf",
    trgt_vcf="data/trgt.vcf",
    catalog_bed="data/catalog.bed",
    catalog_json="data/catalog.json",
    output_dir="analysis_cache",
)

# 2. Verify cache integrity
result = storm.verify_cache("analysis_cache")
assert result["is_valid"]

# 3. Load phenotype data
phenotype = pl.read_csv("phenotypes.csv")["trait"]

# 4. Run association testing
results = storm.run_glm(cache, phenotype)

# 5. Save results
results.write_parquet("association_results.parquet")

# 6. Investigate top hits
top_hits = results.sort("p_value").head(10)
for unit_id in top_hits["unit_id"]:
    print(storm.explain(cache, unit_id))
```

## Cache Structure

The cache directory contains Parquet files:

| File | Description |
|------|-------------|
| `test_units.parquet` | Test unit definitions (id, chrom, start, end, type, motif) |
| `genotypes.parquet` | Resolved genotypes (unit_id, sample_id, allele1, allele2) |
| `catalog.parquet` | Catalog entries (id, chrom, start, end, motif, metadata) |
| `features.parquet` | Computed features (unit_id, sample_id, S, M, D, binary) |
| `provenance.parquet` | Build provenance metadata |
| `provenance.json` | Human-readable provenance |

## Requirements

- Rust 1.70+ (for building)
- Python 3.8+ (for Python bindings)
- Polars (for DataFrame operations)

## Development data

Full-scale SV BCF and TRGT VCFs for development live under `scratch/` (gitignored). See **[DEV_DATA.md](DEV_DATA.md)** for paths (`scratch/chr6.31803187_32050925.bcf`, `scratch/trgt/*.vcf.gz`) and usage.

## Running Tests

```bash
# Rust tests
cargo test

# All tests including integration
cargo test --features python
```

## License

MIT License - see LICENSE file for details.
