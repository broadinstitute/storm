# Development data (scratch)

Local, full-scale datasets for driving development. These live under `scratch/` (gitignored).

## Datasets

### Integrated SV callset (BCF)

- **Path:** `scratch/chr6.31803187_32050925.bcf`
- **Index:** `scratch/chr6.31803187_32050925.bcf.csi`
- **Region:** chr6:31,803,187–32,050,925
- **Samples:** ~12,680
- **Content:** Multi-sample SV callset (DEL, INS, DUP, INV, BND); format `GT:FT:SQ:GQ:PS:NE:DP:AD:KS`

### TRGT calls (per-sample VCF.gz)

- **Path:** `scratch/trgt/*.vcf.gz`
- **Samples:** ~10,000 (one file per sample)
- **Region:** Same chr6 region (`chr6.31803187_32050925`) as the BCF
- **Naming:** `{sample_id}_trgt.TR_Explorer_1.01.chr6.31803187_32050925.vcf.gz`
- **Content:** TRGT TR Explorer 1.01; format `GT:AL:ALLR:SD:MC:MS:AP:AM`

## Sample overlap

- BCF has ~12,680 samples; TRGT directory has ~10,023 `.vcf.gz` files.
- Sample IDs match: BCF column names are numeric IDs (e.g. `1000234`); TRGT sample ID is the filename prefix before `_trgt`.
- Use the subset of samples present in both for joint SV + TRGT workflows.

## Usage

**SV only (BCF):**

```bash
storm build-cache \
  --plan fixtures/plan.yaml \
  --sv-vcf scratch/chr6.31803187_32050925.bcf \
  --catalog fixtures/trexplorer.json \
  --out-dir my_cache
```

**TRGT only (many VCFs):**  
Current API expects paths to VCF/BCF; you can pass a list of paths to `scratch/trgt/*.vcf.gz` (or a subset) where the pipeline accepts multiple TRGT inputs.

**Joint SV + TRGT:**  
Use the same sample ordering or sample subset in both the BCF and the TRGT file list so genotypes align.

## Notes

- `scratch/` is in `.gitignore`; nothing under it is committed.
- For quick tests, use a subset of TRGT files (e.g. first 100) to keep runtime and memory manageable.
