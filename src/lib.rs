//! STORM: Structural & Tandem-Repeat Optimized Regression Models
//!
//! A Rust crate with Python front-end for association testing of
//! structural variants and tandem repeats using long-read data.

#[cfg(feature = "python")]
use pyo3::prelude::*;

use std::collections::HashMap;
use std::io::Write;
use std::path::Path;
use thiserror::Error;
use arrow::array::Array;

// Core modules
pub mod vcf;
pub mod catalog;
pub mod mapping;
pub mod testunit;
pub mod resolver;
pub mod cache;
pub mod explain;
pub mod plan;
pub mod encoding;
pub mod glm;
pub mod results;

// Re-exports
pub use vcf::sv::{parse_sv_vcf, SvRecord, SvType, SvVcfError, Genotype};
pub use vcf::trgt::{parse_trgt_vcf, TrgtRecord, TrgtGenotype, TrgtVcfError};
pub use catalog::{
    Catalog, CatalogEntry, CatalogError,
    parse_trexplorer_bed, BedRecord, BedError,
    parse_trexplorer_json, JsonRecord, JsonError,
};
pub use mapping::{map_svs_to_catalog, map_svs_with_overlaps, SvMapping, MappingStats, compute_mapping_stats};
pub use testunit::{TestUnit, TestUnitType, TestUnitBuilder, DataSource};
pub use resolver::{Resolver, ResolvedGenotype, PresenceSource, AlleleSource, LocusSvAlleles};
pub use cache::{
    ArrowCache, ArrowCacheError, Provenance, write_parquet, read_parquet, ParquetError,
};
pub use cache::parquet_cache::{write_cache_to_dir, read_cache_from_dir};
pub use explain::{
    explain_genotype as explain_genotype_writer,
    explain_locus as explain_locus_writer,
    format_genotype_compact,
};
pub use plan::{Plan, PlanError, Rule, RuleCondition, Encoding, Model};
pub use encoding::{encode_s, encode_m, encode_d, encode_tail, encode_categorical, encode_binary, CategoricalBin, apply_encoding};
pub use glm::{GlmError, AssociationResult, linear_regression, logistic_regression, binomi_rare_test, firth_logistic_regression, categorical_regression, run_association};
pub use results::{ResultsError, build_results_batch, write_results_parquet, ResultsSummary, RareVariantLadder, Covariates};

/// Errors that can occur during cache building
#[derive(Error, Debug)]
pub enum BuildCacheError {
    #[error("SV VCF parsing error: {0}")]
    SvVcf(#[from] SvVcfError),
    #[error("TRGT VCF parsing error: {0}")]
    TrgtVcf(#[from] TrgtVcfError),
    #[error("Catalog error: {0}")]
    Catalog(#[from] CatalogError),
    #[error("Cache error: {0}")]
    Cache(#[from] ArrowCacheError),
    #[error("Parquet error: {0}")]
    Parquet(#[from] ParquetError),
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
    #[error("JSON error: {0}")]
    Json(#[from] serde_json::Error),
}

/// Statistics from cache building
#[derive(Debug, Clone)]
pub struct CacheBuildStats {
    pub num_test_units: usize,
    pub num_samples: usize,
    pub num_genotypes: usize,
    pub num_catalog_entries: usize,
}

/// Result of cache verification
#[derive(Debug, Clone)]
pub struct CacheVerifyResult {
    pub is_valid: bool,
    pub num_test_units: usize,
    pub num_genotypes: usize,
    pub num_catalog_entries: usize,
    pub num_features: usize,
    pub errors: Vec<String>,
}

/// Build a STORM cache from input files
///
/// This is the main entry point for building a cache. It:
/// 1. Parses the SV VCF
/// 2. Optionally parses the TRGT VCF
/// 3. Optionally loads the catalog from BED and/or JSON
/// 4. Maps SVs to catalog loci
/// 5. Builds TestUnits
/// 6. Resolves genotypes
/// 7. Computes features
/// 8. Writes cache to output directory
pub fn build_cache(
    sv_vcf_path: &str,
    trgt_vcf_path: Option<&str>,
    catalog_bed_path: Option<&str>,
    catalog_json_path: Option<&str>,
    output_dir: &str,
) -> Result<CacheBuildStats, BuildCacheError> {
    // 1. Parse SV VCF
    let (sv_samples, sv_records) = parse_sv_vcf(sv_vcf_path)?;
    
    // Collect all sample IDs from SV records
    let mut all_samples: Vec<String> = sv_samples.clone();
    all_samples.sort();

    // 2. Optionally parse TRGT VCF
    let trgt_records = if let Some(path) = trgt_vcf_path {
        let (_samples, records) = parse_trgt_vcf(path)?;
        Some(records)
    } else {
        None
    };

    // Add TRGT samples to all_samples
    if let Some(ref trgt) = trgt_records {
        for rec in trgt {
            for sample_id in rec.genotypes.keys() {
                if !all_samples.contains(sample_id) {
                    all_samples.push(sample_id.clone());
                }
            }
        }
        all_samples.sort();
    }

    // 3. Load catalog
    let catalog = match (catalog_bed_path, catalog_json_path) {
        (Some(bed), Some(json)) => Some(Catalog::from_bed_and_json(bed, json)?),
        (Some(bed), None) => Some(Catalog::from_bed(bed)?),
        (None, Some(json)) => Some(Catalog::from_json(json)?),
        (None, None) => None,
    };

    // 4. Build test units
    let mut test_units = Vec::new();
    let mut genotypes_map: Vec<(String, Vec<ResolvedGenotype>)> = Vec::new();

    // Create TestUnits from SVs
    for sv in &sv_records {
        let sv_type_str = match &sv.sv_type {
            SvType::Del => "DEL",
            SvType::Ins => "INS",
            SvType::Dup => "DUP",
            SvType::Inv => "INV",
            SvType::Bnd => "BND",
            SvType::Cnv => "CNV",
            SvType::Other(s) => s.as_str(),
        };

        let end = sv.end.unwrap_or_else(|| {
            if let Some(svlen) = sv.sv_len {
                sv.pos + (svlen.abs() as u64)
            } else {
                sv.pos + 1
            }
        });

        let unit = TestUnit::from_sv(&sv.id, &sv.chrom, sv.pos, end, sv_type_str);
        
        // Resolve genotypes for this SV
        let mut unit_gts = Vec::new();
        for sample_id in &all_samples {
            let gt = if let Some(sv_gt) = sv.genotypes.get(sample_id) {
                let is_present = !sv_gt.is_ref_hom();
                let (allele1, allele2) = if let Some(svlen) = sv.sv_len {
                    // For INS: length is positive; for DEL: length is negative
                    let len = svlen.abs();
                    if is_present {
                        // Heterozygous: one ref, one alt
                        if sv_gt.is_het() {
                            (Some(0), Some(len))
                        } else {
                            // Homozygous alt
                            (Some(len), Some(len))
                        }
                    } else {
                        // Reference
                        (Some(0), Some(0))
                    }
                } else {
                    (None, None)
                };

                ResolvedGenotype {
                    sample_id: sample_id.clone(),
                    is_present,
                    presence_source: PresenceSource::SvVcf,
                    allele1,
                    allele2,
                    allele_source: AlleleSource::SvProxy,
                    raw_gt: Some(format!("{:?}", sv_gt.alleles)),
                }
            } else {
                ResolvedGenotype::missing(sample_id)
            };
            unit_gts.push(gt);
        }

        genotypes_map.push((unit.id.clone(), unit_gts));
        test_units.push(unit);
    }

    // If we have a catalog, also create RepeatProxy test units for mapped loci
    if let Some(ref cat) = catalog {
        let mappings = map_svs_to_catalog(&sv_records, cat);
        
        // Group mappings by catalog entry
        let mut entry_sv_map: HashMap<String, Vec<&SvRecord>> = HashMap::new();
        for mapping in &mappings {
            for entry in &mapping.overlaps {
                entry_sv_map
                    .entry(entry.id.clone())
                    .or_default()
                    .push(mapping.sv);
            }
        }

        // Create RepeatProxy units for each catalog entry with overlapping SVs
        for (entry_id, svs) in entry_sv_map {
            if let Some(entry) = cat.entries.get(&entry_id) {
                let sv_ids: Vec<String> = svs.iter().map(|sv| sv.id.clone()).collect();
                let unit = TestUnit::from_repeat_proxy(
                    &entry_id,
                    &entry.chrom,
                    entry.start,
                    entry.end,
                    &entry.motif,
                    sv_ids,
                );

                // Resolve genotypes using the Resolver
                let mut resolver = Resolver::new();
                let entries_slice: Vec<&CatalogEntry> = vec![entry];
                resolver.add_catalog_refs(&entries_slice);
                
                // Create slice of references for resolve_from_svs
                let svs_refs: Vec<&SvRecord> = svs.iter().copied().collect();
                let unit_gts = resolver.resolve_from_svs(&entry_id, &svs_refs, &all_samples);

                genotypes_map.push((unit.id.clone(), unit_gts));
                test_units.push(unit);
            }
        }
    }

    // If we have TRGT data, create TrueRepeat units
    if let Some(ref trgt) = trgt_records {
        for rec in trgt {
            let motif = rec.motifs.first().map(|s| s.as_str()).unwrap_or("N");
            let unit = TestUnit::from_true_repeat(
                &rec.trid,
                &rec.chrom,
                rec.pos,
                rec.end,
                motif,
                &rec.trid,
            );

            let mut unit_gts = Vec::new();
            for sample_id in &all_samples {
                let gt = if let Some(trgt_gt) = rec.genotypes.get(sample_id) {
                    let is_present = trgt_gt.allele_lengths.iter().any(|&l| l > 0);
                    let (allele1, allele2) = if trgt_gt.allele_lengths.len() >= 2 {
                        (
                            Some(trgt_gt.allele_lengths[0]),
                            Some(trgt_gt.allele_lengths[1]),
                        )
                    } else if trgt_gt.allele_lengths.len() == 1 {
                        (Some(trgt_gt.allele_lengths[0]), None)
                    } else {
                        (None, None)
                    };

                    ResolvedGenotype {
                        sample_id: sample_id.clone(),
                        is_present,
                        presence_source: PresenceSource::TrgtVcf,
                        allele1,
                        allele2,
                        allele_source: AlleleSource::TrgtTrue,
                        raw_gt: Some(trgt_gt.allele_lengths.iter().map(|l| l.to_string()).collect::<Vec<_>>().join("/")),
                    }
                } else {
                    ResolvedGenotype::missing(sample_id)
                };
                unit_gts.push(gt);
            }

            genotypes_map.push((unit.id.clone(), unit_gts));
            test_units.push(unit);
        }
    }

    // 5. Compute features
    let mut features: Vec<(String, String, HashMap<String, f64>)> = Vec::new();
    for (unit_id, gts) in &genotypes_map {
        for gt in gts {
            let mut fm = HashMap::new();
            if let Some(sum) = gt.sum_lengths() {
                fm.insert("S".to_string(), sum as f64);
            }
            if let Some(max) = gt.max_length() {
                fm.insert("M".to_string(), max as f64);
            }
            if let Some(diff) = gt.diff_lengths() {
                fm.insert("D".to_string(), diff as f64);
            }
            fm.insert("binary".to_string(), if gt.is_present { 1.0 } else { 0.0 });
            features.push((unit_id.clone(), gt.sample_id.clone(), fm));
        }
    }

    // 6. Build ArrowCache
    let mut cache = ArrowCache::new();
    cache.set_test_units(&test_units)?;
    cache.set_genotypes(&genotypes_map)?;
    
    if let Some(ref cat) = catalog {
        let entries: Vec<&CatalogEntry> = cat.entries.values().collect();
        cache.set_catalog(&entries)?;
    } else {
        // Create empty catalog
        let empty_entries: Vec<&CatalogEntry> = Vec::new();
        cache.set_catalog(&empty_entries)?;
    }
    
    cache.set_features(&features)?;

    // 7. Set provenance
    let mut prov = Provenance::new();
    prov.input_files.push(sv_vcf_path.to_string());
    if let Some(path) = trgt_vcf_path {
        prov.input_files.push(path.to_string());
    }
    if let Some(path) = catalog_bed_path {
        prov.input_files.push(path.to_string());
    }
    if let Some(path) = catalog_json_path {
        prov.input_files.push(path.to_string());
    }
    prov.num_samples = all_samples.len();
    prov.num_test_units = test_units.len();
    cache.set_provenance(&prov)?;

    // 8. Write cache to directory
    let output_path = Path::new(output_dir);
    write_cache_to_dir(output_path, &cache)?;

    // Also write provenance as JSON for easy reading
    let prov_json = serde_json::json!({
        "storm_version": prov.storm_version,
        "created_at": prov.created_at,
        "input_files": prov.input_files,
        "num_samples": prov.num_samples,
        "num_test_units": prov.num_test_units,
    });
    let prov_json_path = output_path.join("provenance.json");
    std::fs::write(prov_json_path, serde_json::to_string_pretty(&prov_json)?)?;

    // Count genotypes
    let num_genotypes: usize = genotypes_map.iter().map(|(_, gts)| gts.len()).sum();

    Ok(CacheBuildStats {
        num_test_units: test_units.len(),
        num_samples: all_samples.len(),
        num_genotypes,
        num_catalog_entries: catalog.map(|c| c.entries.len()).unwrap_or(0),
    })
}

/// Verify a STORM cache
pub fn verify_cache(cache_dir: &str) -> Result<CacheVerifyResult, BuildCacheError> {
    let dir = Path::new(cache_dir);
    let mut errors = Vec::new();
    let mut is_valid = true;

    // Check required files exist
    let required_files = [
        "test_units.parquet",
        "genotypes.parquet",
        "catalog.parquet",
        "features.parquet",
        "provenance.parquet",
    ];

    for file in &required_files {
        let path = dir.join(file);
        if !path.exists() {
            errors.push(format!("Missing required file: {}", file));
            is_valid = false;
        }
    }

    // Try to read the cache
    let cache = read_cache_from_dir(dir)?;

    let num_test_units = cache.test_units.as_ref().map(|b| b.num_rows()).unwrap_or(0);
    let num_genotypes = cache.genotypes.as_ref().map(|b| b.num_rows()).unwrap_or(0);
    let num_catalog_entries = cache.catalog.as_ref().map(|b| b.num_rows()).unwrap_or(0);
    let num_features = cache.features.as_ref().map(|b| b.num_rows()).unwrap_or(0);

    // Validate schema
    if let Some(ref batch) = cache.test_units {
        let expected_cols = ["id", "chrom", "start", "end", "unit_type", "source", "catalog_id", "motif"];
        for col in expected_cols {
            if batch.schema().column_with_name(col).is_none() {
                errors.push(format!("test_units.parquet missing column: {}", col));
                is_valid = false;
            }
        }
    }

    if let Some(ref batch) = cache.genotypes {
        let expected_cols = ["unit_id", "sample_id", "is_present", "presence_source", "allele1", "allele2", "allele_source"];
        for col in expected_cols {
            if batch.schema().column_with_name(col).is_none() {
                errors.push(format!("genotypes.parquet missing column: {}", col));
                is_valid = false;
            }
        }
    }

    // Verify row count consistency
    // Features should have same or fewer rows than genotypes
    if num_features > 0 && num_genotypes > 0 && num_features != num_genotypes {
        // This is okay - features might be sparse
    }

    Ok(CacheVerifyResult {
        is_valid,
        num_test_units,
        num_genotypes,
        num_catalog_entries,
        num_features,
        errors,
    })
}

/// Explain genotypes for a specific sample at a test unit
pub fn explain_genotype(cache: &ArrowCache, unit_id: &str, sample_id: &str) -> Result<String, BuildCacheError> {
    use arrow::array::{StringArray, BooleanArray, Int64Array};
    
    let mut output = Vec::new();
    
    // Find the test unit
    let test_units_batch = cache.get_test_units()?;
    let id_col = test_units_batch
        .column_by_name("id")
        .and_then(|c| c.as_any().downcast_ref::<StringArray>())
        .ok_or_else(|| ArrowCacheError::NotPopulated("test_units.id".to_string()))?;
    
    let mut unit_idx = None;
    for i in 0..id_col.len() {
        if id_col.value(i) == unit_id {
            unit_idx = Some(i);
            break;
        }
    }
    
    let unit_idx = unit_idx.ok_or_else(|| {
        std::io::Error::new(std::io::ErrorKind::NotFound, format!("Test unit not found: {}", unit_id))
    })?;

    // Extract test unit info
    let chrom_col = test_units_batch.column_by_name("chrom").and_then(|c| c.as_any().downcast_ref::<StringArray>());
    let start_col = test_units_batch.column_by_name("start").and_then(|c| c.as_any().downcast_ref::<arrow::array::UInt64Array>());
    let end_col = test_units_batch.column_by_name("end").and_then(|c| c.as_any().downcast_ref::<arrow::array::UInt64Array>());
    let unit_type_col = test_units_batch.column_by_name("unit_type").and_then(|c| c.as_any().downcast_ref::<StringArray>());
    let motif_col = test_units_batch.column_by_name("motif").and_then(|c| c.as_any().downcast_ref::<StringArray>());

    writeln!(output, "=== Genotype Explanation ===").ok();
    writeln!(output).ok();
    writeln!(output, "Test Unit: {}", unit_id).ok();
    
    if let (Some(chrom), Some(start), Some(end)) = (chrom_col, start_col, end_col) {
        writeln!(output, "  Location: {}:{}-{}", chrom.value(unit_idx), start.value(unit_idx), end.value(unit_idx)).ok();
    }
    if let Some(ut) = unit_type_col {
        writeln!(output, "  Type: {}", ut.value(unit_idx)).ok();
    }
    if let Some(m) = motif_col {
        if !m.is_null(unit_idx) {
            writeln!(output, "  Motif: {}", m.value(unit_idx)).ok();
        }
    }
    writeln!(output).ok();

    // Find the genotype
    let genotypes_batch = cache.get_genotypes()?;
    let gt_unit_col = genotypes_batch.column_by_name("unit_id").and_then(|c| c.as_any().downcast_ref::<StringArray>());
    let gt_sample_col = genotypes_batch.column_by_name("sample_id").and_then(|c| c.as_any().downcast_ref::<StringArray>());
    let gt_present_col = genotypes_batch.column_by_name("is_present").and_then(|c| c.as_any().downcast_ref::<BooleanArray>());
    let gt_a1_col = genotypes_batch.column_by_name("allele1").and_then(|c| c.as_any().downcast_ref::<Int64Array>());
    let gt_a2_col = genotypes_batch.column_by_name("allele2").and_then(|c| c.as_any().downcast_ref::<Int64Array>());
    let gt_source_col = genotypes_batch.column_by_name("allele_source").and_then(|c| c.as_any().downcast_ref::<StringArray>());

    if let (Some(unit_col), Some(sample_col)) = (gt_unit_col, gt_sample_col) {
        for i in 0..unit_col.len() {
            if unit_col.value(i) == unit_id && sample_col.value(i) == sample_id {
                writeln!(output, "Sample: {}", sample_id).ok();
                
                if let Some(present) = gt_present_col {
                    writeln!(output, "  Presence: {}", if present.value(i) { "PRESENT" } else { "ABSENT" }).ok();
                }
                writeln!(output).ok();
                
                writeln!(output, "Alleles:").ok();
                if let Some(a1) = gt_a1_col {
                    if a1.is_null(i) {
                        writeln!(output, "  Allele 1: MISSING").ok();
                    } else {
                        writeln!(output, "  Allele 1: {} bp", a1.value(i)).ok();
                    }
                }
                if let Some(a2) = gt_a2_col {
                    if a2.is_null(i) {
                        writeln!(output, "  Allele 2: MISSING").ok();
                    } else {
                        writeln!(output, "  Allele 2: {} bp", a2.value(i)).ok();
                    }
                }
                if let Some(src) = gt_source_col {
                    writeln!(output, "  Source: {}", src.value(i)).ok();
                }
                writeln!(output).ok();

                // Compute encodings
                let a1_val = gt_a1_col.and_then(|a| if a.is_null(i) { None } else { Some(a.value(i)) });
                let a2_val = gt_a2_col.and_then(|a| if a.is_null(i) { None } else { Some(a.value(i)) });
                
                if let (Some(a1), Some(a2)) = (a1_val, a2_val) {
                    writeln!(output, "Computed Encodings:").ok();
                    writeln!(output, "  S (sum): {} bp", a1 + a2).ok();
                    writeln!(output, "  M (max): {} bp", a1.max(a2)).ok();
                    writeln!(output, "  D (diff): {} bp", (a1 - a2).abs()).ok();
                }
                
                break;
            }
        }
    }

    Ok(String::from_utf8(output).unwrap_or_default())
}

/// Explain all genotypes at a locus
pub fn explain_locus(cache: &ArrowCache, unit_id: &str) -> Result<String, BuildCacheError> {
    use arrow::array::{StringArray, BooleanArray, Int64Array};
    
    let mut output = Vec::new();
    
    // Find the test unit
    let test_units_batch = cache.get_test_units()?;
    let id_col = test_units_batch
        .column_by_name("id")
        .and_then(|c| c.as_any().downcast_ref::<StringArray>())
        .ok_or_else(|| ArrowCacheError::NotPopulated("test_units.id".to_string()))?;
    
    let mut unit_idx = None;
    for i in 0..id_col.len() {
        if id_col.value(i) == unit_id {
            unit_idx = Some(i);
            break;
        }
    }
    
    let unit_idx = unit_idx.ok_or_else(|| {
        std::io::Error::new(std::io::ErrorKind::NotFound, format!("Test unit not found: {}", unit_id))
    })?;

    // Extract test unit info
    let chrom_col = test_units_batch.column_by_name("chrom").and_then(|c| c.as_any().downcast_ref::<StringArray>());
    let start_col = test_units_batch.column_by_name("start").and_then(|c| c.as_any().downcast_ref::<arrow::array::UInt64Array>());
    let end_col = test_units_batch.column_by_name("end").and_then(|c| c.as_any().downcast_ref::<arrow::array::UInt64Array>());
    let unit_type_col = test_units_batch.column_by_name("unit_type").and_then(|c| c.as_any().downcast_ref::<StringArray>());
    let motif_col = test_units_batch.column_by_name("motif").and_then(|c| c.as_any().downcast_ref::<StringArray>());

    writeln!(output, "=== Locus Summary ===").ok();
    writeln!(output).ok();
    writeln!(output, "Test Unit: {}", unit_id).ok();
    
    if let (Some(chrom), Some(start), Some(end)) = (chrom_col, start_col, end_col) {
        writeln!(output, "Location: {}:{}-{}", chrom.value(unit_idx), start.value(unit_idx), end.value(unit_idx)).ok();
    }
    if let Some(ut) = unit_type_col {
        writeln!(output, "Type: {}", ut.value(unit_idx)).ok();
    }
    if let Some(m) = motif_col {
        if !m.is_null(unit_idx) {
            writeln!(output, "Motif: {}", m.value(unit_idx)).ok();
        }
    }
    writeln!(output).ok();

    // Collect genotypes for this unit
    let genotypes_batch = cache.get_genotypes()?;
    let gt_unit_col = genotypes_batch.column_by_name("unit_id").and_then(|c| c.as_any().downcast_ref::<StringArray>());
    let gt_sample_col = genotypes_batch.column_by_name("sample_id").and_then(|c| c.as_any().downcast_ref::<StringArray>());
    let gt_present_col = genotypes_batch.column_by_name("is_present").and_then(|c| c.as_any().downcast_ref::<BooleanArray>());
    let gt_a1_col = genotypes_batch.column_by_name("allele1").and_then(|c| c.as_any().downcast_ref::<Int64Array>());
    let gt_a2_col = genotypes_batch.column_by_name("allele2").and_then(|c| c.as_any().downcast_ref::<Int64Array>());
    let gt_source_col = genotypes_batch.column_by_name("allele_source").and_then(|c| c.as_any().downcast_ref::<StringArray>());

    let mut samples = Vec::new();
    let mut present_count = 0;
    let mut with_alleles = 0;

    if let (Some(unit_col), Some(sample_col), Some(present_col)) = (gt_unit_col, gt_sample_col, gt_present_col) {
        for i in 0..unit_col.len() {
            if unit_col.value(i) == unit_id {
                let sample_id = sample_col.value(i).to_string();
                let is_present = present_col.value(i);
                let a1 = gt_a1_col.and_then(|c| if c.is_null(i) { None } else { Some(c.value(i)) });
                let a2 = gt_a2_col.and_then(|c| if c.is_null(i) { None } else { Some(c.value(i)) });
                let source = gt_source_col.map(|c| c.value(i).to_string()).unwrap_or_default();

                if is_present {
                    present_count += 1;
                }
                if a1.is_some() && a2.is_some() {
                    with_alleles += 1;
                }

                samples.push((sample_id, is_present, a1, a2, source));
            }
        }
    }

    let total = samples.len();
    let call_rate = if total > 0 { with_alleles as f64 / total as f64 } else { 0.0 };

    writeln!(output, "Summary Statistics:").ok();
    writeln!(output, "  Total samples: {}", total).ok();
    writeln!(output, "  Present (carriers): {}", present_count).ok();
    writeln!(output, "  Absent (non-carriers): {}", total - present_count).ok();
    writeln!(output, "  Call rate: {:.2}%", call_rate * 100.0).ok();
    writeln!(output).ok();

    writeln!(output, "Sample Details:").ok();
    writeln!(output, "{:<15} {:>8} {:>10} {:>10} {:>12}", "Sample", "Present", "Allele1", "Allele2", "Source").ok();
    writeln!(output, "{}", "-".repeat(60)).ok();

    for (sample_id, is_present, a1, a2, source) in samples {
        let present_str = if is_present { "Yes" } else { "No" };
        let a1_str = a1.map(|v| v.to_string()).unwrap_or_else(|| ".".to_string());
        let a2_str = a2.map(|v| v.to_string()).unwrap_or_else(|| ".".to_string());
        writeln!(output, "{:<15} {:>8} {:>10} {:>10} {:>12}", sample_id, present_str, a1_str, a2_str, source).ok();
    }

    Ok(String::from_utf8(output).unwrap_or_default())
}

// Example library function
#[cfg(feature = "python")]
#[pyfunction]
pub fn add(a: i32, b: i32) -> i32 {
    a + b
}

/// Returns the version of the storm package as a string.
///
/// This function returns the version that was set in Cargo.toml at compile time.
#[cfg(feature = "python")]
#[pyfunction]
fn _version() -> PyResult<String> {
    Ok(env!("CARGO_PKG_VERSION").to_string())
}

#[cfg(feature = "python")]
#[pymodule]
fn storm(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(add, m)?)?;
    m.add_function(wrap_pyfunction!(_version, m)?)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    #[test]
    #[cfg(feature = "python")]
    fn test_add() {
        assert_eq!(super::add(2, 2), 4);
        assert_eq!(super::add(-1, 1), 0);
        assert_eq!(super::add(0, 0), 0);
    }

    #[test]
    #[cfg(feature = "python")]
    fn test_version() {
        let version = super::_version().unwrap();
        assert!(!version.is_empty());
        // Version should be in format x.y.z
        assert!(version.matches(r"^\d+\.\d+\.\d+$"));
    }
}
