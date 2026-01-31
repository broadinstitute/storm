# Development data (scratch)

Local, full-scale datasets for driving development. These live under `scratch/` (gitignored).

## Datasets

### Integrated SV callset (BCF)

- **Path:** `scratch/chr6.31803187_32050925.bcf`
- **Index:** `scratch/chr6.31803187_32050925.bcf.csi`
- **Region:** chr6:31,803,187–32,050,925
- **Samples:** ~12,680
- **Content:** Multi-sample SV callset (DEL, INS, DUP, INV, BND); format `GT:FT:SQ:GQ:PS:NE:DP:AD:KS`

### TRExplorer repeat catalog (optional, for RepeatProxy units)

- **BED:** `scratch/TRExplorer.repeat_catalog_v2.hg38.1_to_1000bp_motifs.bed.gz` (plain or .gz)
- **JSON (annotations):** `scratch/TRExplorer.repeat_catalog_v2.hg38.1_to_1000bp_motifs.EH.with_annotations.json.gz` (plain or .gz)
- **Note:** Use the **main** BED and the **EH** JSON for STORM’s catalog. The `.TRGT.bed.gz` file is TRGT-specific and is not used by STORM’s cache build. STORM supports gzipped catalog files (`.bed.gz`, `.json.gz`).

### TRGT calls (per-sample VCF.gz)

- **Path:** `scratch/trgt/*.vcf.gz`
- **Samples:** ~10,000 (one file per sample)
- **Region:** Same chr6 region (`chr6.31803187_32050925`) as the BCF
- **Naming:** `{sample_id}_trgt.TR_Explorer_1.01.chr6.31803187_32050925.vcf.gz`
- **Content:** TRGT TR Explorer 1.01; format `GT:AL:ALLR:SD:MC:MS:AP:AM`

## Sample overlap

- BCF has ~12,680 samples; TRGT directory has ~10,023 `.vcf.gz` files.
- Sample IDs match: BCF column names are numeric IDs (e.g. `1000234`); TRGT sample ID is the filename prefix before `_trgt`.
- **When both BCF and TRGT are provided, the cache uses the intersection of samples** (only samples present in both). This ensures genotypes align for joint SV + TRGT association testing.

## Usage

**SV only (BCF):**

```bash
storm cache build \
  --sv-vcf scratch/chr6.31803187_32050925.bcf \
  --output-dir my_cache
```

**TRGT only (many VCFs via directory):**

```bash
storm cache build \
  --sv-vcf fixtures/sv_small.vcf \
  --trgt-dir scratch/trgt \
  --output-dir my_cache
```

**Multiple specific TRGT files:**

```bash
storm cache build \
  --sv-vcf scratch/chr6.31803187_32050925.bcf \
  --trgt-vcf scratch/trgt/1000234_trgt.TR_Explorer_1.01.chr6.31803187_32050925.vcf.gz \
  --trgt-vcf scratch/trgt/1000291_trgt.TR_Explorer_1.01.chr6.31803187_32050925.vcf.gz \
  --output-dir my_cache
```

**TRGT files from a list file:**

```bash
# Create a list file with one path per line
ls scratch/trgt/*.vcf.gz | head -100 > trgt_files.txt
storm cache build \
  --sv-vcf scratch/chr6.31803187_32050925.bcf \
  --trgt-list trgt_files.txt \
  --output-dir my_cache
```

**Joint SV + TRGT:**  
The cache includes only samples that appear in both the BCF and the TRGT files (intersection), so genotypes align for association testing.

### Python API

```python
import storm

# Single TRGT file (as list)
storm.build_cache(
    sv_vcf="scratch/chr6.31803187_32050925.bcf",
    trgt_vcf=["scratch/trgt/1000234_trgt.TR_Explorer_1.01.chr6.31803187_32050925.vcf.gz"],
    output_dir="my_cache"
)

# Multiple TRGT files
trgt_files = ["scratch/trgt/1000234_...vcf.gz", "scratch/trgt/1000291_...vcf.gz"]
storm.build_cache(
    sv_vcf="scratch/chr6.31803187_32050925.bcf",
    trgt_vcf=trgt_files,
    output_dir="my_cache"
)
```

## Notes

- `scratch/` is in `.gitignore`; nothing under it is committed.
- For quick tests, use a subset of TRGT files (e.g. first 100) to keep runtime and memory manageable.
