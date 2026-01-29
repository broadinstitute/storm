# Test Fixture vs Real Data Comparison

## Summary

The test fixtures **adequately resemble** the real data files for basic parsing, but there are some important differences in scale and format complexity that should be noted.

## TRGT VCF Files

### Test Fixture (`fixtures/trgt_small.vcf`)
- **Samples**: 2 (SAMPLE1, SAMPLE2)
- **Records**: 3
- **Format fields**: `GT:AL:SD` (3 fields)
- **Fields used by parser**: GT, AL ✅
- **Additional fields**: SD (not used by parser)

### Real Data (development: `scratch/trgt/*.vcf.gz`)
- **Samples**: ~10,000 (one file per sample; sample ID in filename, e.g. `1000234_trgt.TR_Explorer_1.01.chr6.31803187_32050925.vcf.gz`)
- **Records**: Many (region chr6:31803187–32050925, same as integrated BCF)
- **Format fields**: `GT:AL:ALLR:SD:MC:MS:AP:AM` (8 fields)
- **Fields used by parser**: GT, AL ✅
- **Additional fields**: ALLR, SD, MC, MS, AP, AM (not used by parser)

### Compatibility Assessment: ✅ **COMPATIBLE**
- The parser only uses GT and AL fields, which are present in both
- Additional fields in real data are ignored (gracefully handled)
- Single sample per file vs multiple samples: parser handles both correctly

## SV/BCF Files

### Test Fixture (`fixtures/sv_small.vcf`)
- **Samples**: 2 (SAMPLE1, SAMPLE2)
- **Records**: 3
- **Format fields**: `GT` (1 field)
- **INFO fields**: Basic (SVTYPE, SVLEN, END)
- **Fields used by parser**: GT, SVTYPE, SVLEN, END ✅

### Real Data (development: `scratch/chr6.31803187_32050925.bcf`)
- **Samples**: ~12,680 samples (overlap with TRGT ~10k for joint workflows)
- **Records**: Many (region chr6:31803187–32050925)
- **Format fields**: `GT:FT:SQ:GQ:PS:NE:DP:AD:KS` (9 fields)
- **INFO fields**: Extensive (50+ fields including AC, AF, DP, GQ, SCORE, SUPP_*, etc.)
- **Fields used by parser**: GT, SVTYPE, SVLEN, END ✅

### Compatibility Assessment: ✅ **COMPATIBLE**
- The parser only uses GT format field and SVTYPE/SVLEN/END INFO fields, all present in both
- Additional format fields in real data are ignored (gracefully handled)
- Extensive INFO fields in real data are parsed but only required fields are used
- Scale difference (2 vs 12,680 samples): parser handles both correctly

## Key Differences

### 1. **Scale**
- **Fixtures**: 2-3 samples, 3 records
- **Real data**: ~10k TRGT samples, ~12.7k BCF samples (see `DEV_DATA.md` for paths)
- **Impact**: None - parsers handle any number of samples/records

### 2. **Format Field Complexity**
- **Fixtures**: Minimal format fields (1-3 fields)
- **Real data**: More format fields (8-9 fields)
- **Impact**: None - parsers use `.position()` to find required fields and ignore extras

### 3. **INFO Field Complexity**
- **Fixtures**: Basic INFO fields
- **Real data**: Extensive INFO fields (50+)
- **Impact**: None - parsers only extract required fields (TRID, END, MOTIFS, STRUC for TRGT; SVTYPE, SVLEN, END for SV)

### 4. **File Format**
- **Fixtures**: Plain VCF
- **Real data**: Compressed VCF (.vcf.gz) and BCF (.bcf)
- **Impact**: None - parsers work with uncompressed streams (bcftools/gunzip handles decompression)

## Recommendations

### ✅ **Current Fixtures Are Adequate For:**
1. Testing basic parsing logic
2. Testing format field extraction
3. Testing genotype parsing
4. Unit testing individual components

### ⚠️ **Consider Adding Tests For:**
1. **Large sample counts**: Test with 100+ samples to ensure performance
2. **Missing fields**: Test behavior when optional fields are missing
3. **Compressed files**: Test direct parsing of .vcf.gz files (if supported)
4. **Edge cases**: Missing genotypes, phased vs unphased, etc.

### 🔍 **Potential Issues to Watch For:**
1. **Memory usage**: BCF has ~12,680 samples; TRGT is one file per sample (~10k files) - ensure efficient parsing
2. **Performance**: Large files may need streaming/chunked processing
3. **Missing data**: Real data may have more missing genotypes (./.) - ensure robust handling

## Conclusion

**The test fixtures adequately resemble the real data files** for the purposes of testing the parsing logic. The parsers are designed to be flexible and handle:
- Variable numbers of samples ✅
- Additional format fields ✅
- Extensive INFO fields ✅
- Missing optional fields ✅

The main difference is scale, which is expected and doesn't affect compatibility. The parsers should work correctly with the real data files.

**Development data:** Full-scale SV BCF and TRGT VCFs are in `scratch/`; see **DEV_DATA.md** for paths and usage.
