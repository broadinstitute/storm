//! Integration tests for STORM
//!
//! These tests verify the end-to-end functionality of the STORM pipeline.

use arrow::array::Array;
use storm::{
    build_cache, verify_cache, explain_genotype, explain_locus,
    parse_sv_vcf, parse_trgt_vcf,
    Catalog,
    cache::read_cache_from_dir,
};

fn test_dir(name: &str) -> std::path::PathBuf {
    let dir = std::env::temp_dir().join(format!("storm_integ_{}", name));
    let _ = std::fs::remove_dir_all(&dir);
    dir
}

// ============================================================================
// VCF Parsing Tests
// ============================================================================

#[test]
fn test_sv_vcf_parsing() {
    let (samples, records) = parse_sv_vcf("fixtures/sv_small.vcf")
        .expect("Failed to parse SV VCF");
    
    assert_eq!(samples.len(), 2, "Should have 2 samples");
    assert_eq!(records.len(), 3, "Should have 3 SV records");
    
    // Check first record
    let sv = &records[0];
    assert_eq!(sv.id, "sv1");
    assert_eq!(sv.chrom, "chr1");
    assert_eq!(sv.pos, 1000);
}

#[test]
fn test_trgt_vcf_parsing() {
    let (samples, records) = parse_trgt_vcf("fixtures/trgt_small.vcf")
        .expect("Failed to parse TRGT VCF");
    
    assert_eq!(samples.len(), 2, "Should have 2 samples");
    assert_eq!(records.len(), 3, "Should have 3 TRGT records");
}

// ============================================================================
// Catalog Tests
// ============================================================================

#[test]
fn test_catalog_from_bed() {
    let catalog = Catalog::from_bed("fixtures/trexplorer.bed")
        .expect("Failed to load BED catalog");
    
    assert!(catalog.entries.len() > 0, "Should have catalog entries");
}

#[test]
fn test_catalog_from_json() {
    let catalog = Catalog::from_json("fixtures/trexplorer.json")
        .expect("Failed to load JSON catalog");
    
    assert!(catalog.entries.len() > 0, "Should have catalog entries");
}

#[test]
fn test_catalog_from_bed_and_json() {
    let catalog = Catalog::from_bed_and_json(
        "fixtures/trexplorer.bed",
        "fixtures/trexplorer.json",
    ).expect("Failed to load merged catalog");
    
    assert!(catalog.entries.len() > 0, "Should have catalog entries");
}

// ============================================================================
// Cache Building Tests
// ============================================================================

#[test]
fn test_build_cache_sv_only() {
    let output_dir = test_dir("build_sv");
    
    let stats = build_cache(
        "fixtures/sv_small.vcf",
        None,
        None,
        None,
        output_dir.to_str().unwrap(),
    ).expect("Failed to build cache");
    
    assert_eq!(stats.num_test_units, 3, "Should have 3 test units");
    assert_eq!(stats.num_samples, 2, "Should have 2 samples");
    assert_eq!(stats.num_genotypes, 6, "Should have 6 genotypes");
    
    // Verify files exist
    assert!(output_dir.join("test_units.parquet").exists(), "test_units.parquet should exist");
    assert!(output_dir.join("genotypes.parquet").exists(), "genotypes.parquet should exist");
    assert!(output_dir.join("features.parquet").exists(), "features.parquet should exist");
    assert!(output_dir.join("provenance.json").exists(), "provenance.json should exist");
    
    let _ = std::fs::remove_dir_all(&output_dir);
}

#[test]
fn test_build_cache_from_fixtures() {
    let output_dir = test_dir("build_full");
    
    let stats = build_cache(
        "fixtures/sv_small.vcf",
        Some(&["fixtures/trgt_small.vcf"][..]),
        Some("fixtures/trexplorer.bed"),
        Some("fixtures/trexplorer.json"),
        output_dir.to_str().unwrap(),
    ).expect("Failed to build full cache");
    
    // 3 SV + some repeat proxy + 3 TRGT = at least 6
    assert!(stats.num_test_units >= 6, "Should have at least 6 test units");
    assert_eq!(stats.num_samples, 2, "Should have 2 samples");
    assert!(stats.num_catalog_entries > 0, "Should have catalog entries");
    
    let _ = std::fs::remove_dir_all(&output_dir);
}

// ============================================================================
// Multi-TRGT Cache Tests
// ============================================================================

#[test]
fn test_build_cache_multi_trgt() {
    let output_dir = test_dir("build_multi_trgt");

    // Use SV and TRGT fixtures that share sample IDs so intersection is non-empty.
    // Samples = intersection(BCF, TRGT); trgt_small.vcf has SAMPLE1, SAMPLE2 (same as sv_small.vcf).
    let stats = build_cache(
        "fixtures/sv_small.vcf",
        Some(&["fixtures/trgt_small.vcf"][..]),
        None,
        None,
        output_dir.to_str().unwrap(),
    ).expect("Failed to build cache with TRGT");

    // 3 SV units + 3 TRGT units = 6 test units; samples = intersection = 2 (SAMPLE1, SAMPLE2)
    assert!(stats.num_test_units >= 6, "Should have at least 6 test units, got {}", stats.num_test_units);
    assert_eq!(stats.num_samples, 2, "Should have 2 samples (intersection of BCF and TRGT), got {}", stats.num_samples);
    
    // Verify provenance includes both TRGT files
    let prov_path = output_dir.join("provenance.json");
    let prov_content = std::fs::read_to_string(&prov_path)
        .expect("Failed to read provenance.json");
    let prov: serde_json::Value = serde_json::from_str(&prov_content)
        .expect("Failed to parse provenance.json");
    
    let input_files = prov.get("input_files").and_then(|v| v.as_array());
    assert!(input_files.is_some(), "Should have input_files in provenance");
    let files: Vec<&str> = input_files.unwrap().iter()
        .filter_map(|v| v.as_str())
        .collect();
    
    // Should have sv_small.vcf and TRGT file
    assert!(files.iter().any(|f| f.contains("sv_small.vcf")), "Should include sv_small.vcf");
    assert!(files.iter().any(|f| f.contains("trgt_small.vcf")), "Should include trgt_small.vcf");
    
    let _ = std::fs::remove_dir_all(&output_dir);
}

#[test]
fn test_build_cache_multi_trgt_merged_genotypes() {
    use storm::vcf::trgt::parse_and_merge_trgt_vcfs;
    
    // Test the merge function directly to verify genotypes are combined
    let paths = vec!["fixtures/trgt_sample1.vcf", "fixtures/trgt_sample2.vcf"];
    let (samples, records) = parse_and_merge_trgt_vcfs(&paths)
        .expect("Failed to merge TRGT files");
    
    // Should have 2 samples total
    assert_eq!(samples.len(), 2, "Should have 2 samples after merge");
    
    // Should have 3 records (same loci in both files)
    assert_eq!(records.len(), 3, "Should have 3 merged records");
    
    // Each record should have genotypes from both samples
    for rec in &records {
        assert_eq!(rec.genotypes.len(), 2, 
            "Record {} should have 2 genotypes", rec.trid);
    }
}

// ============================================================================
// Cache Verification Tests
// ============================================================================

#[test]
fn test_verify_cache() {
    let output_dir = test_dir("verify");
    
    build_cache(
        "fixtures/sv_small.vcf",
        None,
        Some("fixtures/trexplorer.bed"),
        None,
        output_dir.to_str().unwrap(),
    ).expect("Failed to build cache");
    
    let result = verify_cache(output_dir.to_str().unwrap())
        .expect("Failed to verify cache");
    
    assert!(result.is_valid, "Cache should be valid");
    assert_eq!(result.errors.len(), 0, "Should have no errors");
    assert!(result.num_test_units > 0, "Should have test units");
    
    let _ = std::fs::remove_dir_all(&output_dir);
}

// ============================================================================
// Schema Tests (using verify_cache which reads cache internally)
// ============================================================================

#[test]
fn test_cache_has_all_tables() {
    let output_dir = test_dir("tables");
    
    build_cache(
        "fixtures/sv_small.vcf",
        None,
        Some("fixtures/trexplorer.bed"),
        Some("fixtures/trexplorer.json"),
        output_dir.to_str().unwrap(),
    ).expect("Failed to build cache");
    
    let result = verify_cache(output_dir.to_str().unwrap())
        .expect("Failed to verify cache");
    
    assert!(result.is_valid, "Cache should be valid");
    assert!(result.num_test_units > 0, "Should have test units");
    assert!(result.num_genotypes > 0, "Should have genotypes");
    assert!(result.num_catalog_entries > 0, "Should have catalog entries");
    assert!(result.num_features > 0, "Should have features");
    
    let _ = std::fs::remove_dir_all(&output_dir);
}

// ============================================================================
// Provenance Tests
// ============================================================================

#[test]
fn test_provenance_metadata() {
    let output_dir = test_dir("provenance");
    
    build_cache(
        "fixtures/sv_small.vcf",
        Some(&["fixtures/trgt_small.vcf"][..]),
        None,
        None,
        output_dir.to_str().unwrap(),
    ).expect("Failed to build cache");
    
    // Read provenance.json
    let prov_path = output_dir.join("provenance.json");
    let prov_content = std::fs::read_to_string(&prov_path)
        .expect("Failed to read provenance.json");
    
    let prov: serde_json::Value = serde_json::from_str(&prov_content)
        .expect("Failed to parse provenance.json");
    
    assert!(prov.get("storm_version").is_some(), "Should have storm_version");
    assert!(prov.get("created_at").is_some(), "Should have created_at");
    assert!(prov.get("input_files").is_some(), "Should have input_files");
    assert!(prov.get("num_samples").is_some(), "Should have num_samples");
    
    let _ = std::fs::remove_dir_all(&output_dir);
}

// ============================================================================
// Read Cache and Explain Tests
// ============================================================================

#[test]
fn test_read_and_explain_locus() {
    let output_dir = test_dir("read_explain_locus");
    
    build_cache(
        "fixtures/sv_small.vcf",
        None,
        None,
        None,
        output_dir.to_str().unwrap(),
    ).expect("Failed to build cache");
    
    // Read the cache back
    let cache = read_cache_from_dir(&output_dir)
        .expect("Failed to read cache from dir");
    
    // Check that we got data
    assert!(cache.test_units.is_some(), "Should have test_units");
    assert!(cache.genotypes.is_some(), "Should have genotypes");
    
    // Test explain_locus
    let explanation = explain_locus(&cache, "sv_sv1")
        .expect("Failed to explain locus");
    
    assert!(explanation.contains("Locus Summary"), "Should have header");
    assert!(explanation.contains("sv_sv1"), "Should contain unit ID");
    
    let _ = std::fs::remove_dir_all(&output_dir);
}

#[test]
fn test_read_and_explain_genotype() {
    let output_dir = test_dir("read_explain_gt");
    
    build_cache(
        "fixtures/sv_small.vcf",
        None,
        None,
        None,
        output_dir.to_str().unwrap(),
    ).expect("Failed to build cache");
    
    // Read the cache back
    let cache = read_cache_from_dir(&output_dir)
        .expect("Failed to read cache from dir");
    
    // Test explain_genotype
    let explanation = explain_genotype(&cache, "sv_sv1", "SAMPLE1")
        .expect("Failed to explain genotype");
    
    assert!(explanation.contains("Genotype Explanation"), "Should have header");
    assert!(explanation.contains("sv_sv1"), "Should contain unit ID");
    assert!(explanation.contains("SAMPLE1"), "Should contain sample");
    
    let _ = std::fs::remove_dir_all(&output_dir);
}

// ============================================================================
// Association Testing Tests
// ============================================================================

#[test]
fn test_run_association_linear() {
    use storm::glm::run_association;
    use storm::plan::{Model, Encoding};
    
    // Simple test with synthetic data
    let x = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0];
    let y = vec![2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0];
    
    let result = run_association(
        "test_unit_assoc",
        &x,
        &y,
        &Model::Linear,
        &Encoding::S,
        None,
    ).expect("Failed to run association");
    
    assert_eq!(result.unit_id, "test_unit_assoc");
    assert!((result.beta - 2.0).abs() < 0.01, "Beta should be ~2.0");
    assert!(result.p_value < 0.001, "Should be highly significant");
    assert_eq!(result.n_samples, 10);
    assert_eq!(result.model, "Linear");
    assert_eq!(result.encoding, "S");
}

#[test]
fn test_run_association_logistic() {
    use storm::glm::run_association;
    use storm::plan::{Model, Encoding};
    
    // Binary outcome data
    let x = vec![0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0];
    let y = vec![0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0];
    
    let result = run_association(
        "test_logistic",
        &x,
        &y,
        &Model::Logistic,
        &Encoding::Binary,
        None,
    ).expect("Failed to run logistic association");
    
    assert_eq!(result.unit_id, "test_logistic");
    assert!(result.beta > 0.0, "Beta should be positive (x=1 -> higher y)");
    assert!(result.se > 0.0, "SE should be positive");
    assert!(result.p_value >= 0.0 && result.p_value <= 1.0, "P-value should be valid");
    assert_eq!(result.model, "Logistic");
}

#[test]
fn test_association_result_structure() {
    use storm::glm::run_association;
    use storm::plan::{Model, Encoding};
    
    // Use larger sample size and non-perfect correlation for realistic test
    let x = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0];
    let y = vec![2.1, 3.9, 6.1, 8.0, 10.2, 11.9, 14.1, 16.0, 17.8, 20.1];
    
    let result = run_association(
        "test_structure",
        &x,
        &y,
        &Model::Linear,
        &Encoding::M,
        None,
    ).expect("Failed to run association");
    
    // Verify all fields of AssociationResult are populated
    assert!(!result.unit_id.is_empty(), "unit_id should not be empty");
    assert!(result.beta.is_finite(), "beta should be finite");
    assert!(result.se.is_finite(), "se should be finite");
    assert!(result.statistic.is_finite(), "statistic should be finite");
    assert!(result.p_value >= 0.0 && result.p_value <= 1.0, "p_value should be in [0,1]");
    assert_eq!(result.n_samples, 10, "n_samples should match input length");
    assert!(result.call_rate >= 0.0 && result.call_rate <= 1.0, "call_rate should be in [0,1]");
    assert!(!result.model.is_empty(), "model should not be empty");
    assert!(!result.encoding.is_empty(), "encoding should not be empty");
}

// ============================================================================
// Canonical Repeat Unit Tests
// ============================================================================

#[test]
fn test_canonical_repeat_units_only() {
    use arrow::array::StringArray;
    
    let output_dir = test_dir("canonical_repeat");
    
    // Build cache with SV + TRGT + catalog
    build_cache(
        "fixtures/sv_small.vcf",
        Some(&["fixtures/trgt_small.vcf"][..]),
        Some("fixtures/trexplorer.bed"),
        None,
        output_dir.to_str().unwrap(),
    ).expect("Failed to build cache");
    
    // Read cache and verify test unit types
    let cache = read_cache_from_dir(&output_dir).expect("Failed to read cache");
    let test_units = cache.test_units.as_ref().expect("Should have test_units");
    
    let unit_type_col = test_units.column_by_name("unit_type")
        .and_then(|c| c.as_any().downcast_ref::<StringArray>())
        .expect("Should have unit_type column");
    
    let id_col = test_units.column_by_name("id")
        .and_then(|c| c.as_any().downcast_ref::<StringArray>())
        .expect("Should have id column");
    
    // Count unit types
    let mut sv_count = 0;
    let mut repeat_count = 0;
    let mut true_repeat_count = 0;
    let mut repeat_proxy_count = 0;
    
    for i in 0..unit_type_col.len() {
        let unit_type = unit_type_col.value(i);
        let id = id_col.value(i);
        
        match unit_type {
            "Sv" => sv_count += 1,
            "Repeat" => repeat_count += 1,
            "TrueRepeat" => {
                true_repeat_count += 1;
                // Fail if we see separate TrueRepeat units (should be merged into Repeat)
                panic!("Found unexpected TrueRepeat unit: {}", id);
            }
            "RepeatProxy" => {
                repeat_proxy_count += 1;
                // Fail if we see separate RepeatProxy units (should be merged into Repeat)
                panic!("Found unexpected RepeatProxy unit: {}", id);
            }
            _ => {}
        }
    }
    
    // Should have only canonical Repeat units, not TrueRepeat or RepeatProxy
    assert!(repeat_count > 0, "Should have canonical Repeat units");
    assert_eq!(true_repeat_count, 0, "Should NOT have separate TrueRepeat units");
    assert_eq!(repeat_proxy_count, 0, "Should NOT have separate RepeatProxy units");
    assert!(sv_count > 0, "Should have SV units for non-overlapping SVs");
    
    let _ = std::fs::remove_dir_all(&output_dir);
}

#[test]
fn test_sv_overlapping_catalog_no_separate_unit() {
    use arrow::array::StringArray;
    use std::io::Write;
    
    let output_dir = test_dir("sv_overlap_catalog");
    
    // Create a VCF with an SV that overlaps catalog locus TR001 (chr1:10000-10050)
    let sv_vcf_path = output_dir.parent().unwrap().join("sv_overlap_test.vcf");
    std::fs::create_dir_all(output_dir.parent().unwrap()).ok();
    let mut f = std::fs::File::create(&sv_vcf_path).expect("Failed to create VCF");
    writeln!(f, "##fileformat=VCFv4.2").ok();
    writeln!(f, "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type\">").ok();
    writeln!(f, "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Length\">").ok();
    writeln!(f, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">").ok();
    writeln!(f, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">").ok();
    writeln!(f, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\tSAMPLE2").ok();
    // This SV overlaps TR001 (chr1:10000-10050)
    writeln!(f, "chr1\t9900\tsv_inside\tN\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN=500;END=10400\tGT\t0/1\t1/1").ok();
    // This SV is outside all catalog regions
    writeln!(f, "chr1\t1000\tsv_outside\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-500;END=1500\tGT\t0/1\t0/0").ok();
    drop(f);
    
    // Build cache
    build_cache(
        sv_vcf_path.to_str().unwrap(),
        None,
        Some("fixtures/trexplorer.bed"),
        None,
        output_dir.to_str().unwrap(),
    ).expect("Failed to build cache");
    
    // Read cache and verify
    let cache = read_cache_from_dir(&output_dir).expect("Failed to read cache");
    let test_units = cache.test_units.as_ref().expect("Should have test_units");
    
    let id_col = test_units.column_by_name("id")
        .and_then(|c| c.as_any().downcast_ref::<StringArray>())
        .expect("Should have id column");
    
    let unit_type_col = test_units.column_by_name("unit_type")
        .and_then(|c| c.as_any().downcast_ref::<StringArray>())
        .expect("Should have unit_type column");
    
    // Collect unit IDs and types
    let mut found_sv_inside = false;
    let mut found_sv_outside = false;
    let mut found_repeat_tr001 = false;
    
    for i in 0..id_col.len() {
        let id = id_col.value(i);
        let unit_type = unit_type_col.value(i);
        
        if id == "sv_sv_inside" {
            found_sv_inside = true;
        }
        if id == "sv_sv_outside" {
            found_sv_outside = true;
        }
        if id == "repeat_TR001" && unit_type == "Repeat" {
            found_repeat_tr001 = true;
        }
    }
    
    // sv_inside overlaps TR001, should NOT have separate unit
    assert!(!found_sv_inside, "SV overlapping catalog should NOT have separate sv_* unit");
    // sv_outside does not overlap catalog, should have separate unit
    assert!(found_sv_outside, "SV outside catalog SHOULD have sv_* unit");
    // TR001 should exist as canonical repeat unit (with SV contributing to it)
    assert!(found_repeat_tr001, "Should have canonical repeat_TR001 unit");
    
    let _ = std::fs::remove_dir_all(&output_dir);
    let _ = std::fs::remove_file(&sv_vcf_path);
}

#[test]
fn test_trgt_merged_into_catalog_unit() {
    use arrow::array::{StringArray, BooleanArray, Int64Array};
    
    let output_dir = test_dir("trgt_merged");
    
    // Build cache with TRGT + catalog (no overlapping SVs)
    build_cache(
        "fixtures/sv_small.vcf",  // SVs don't overlap catalog
        Some(&["fixtures/trgt_small.vcf"][..]),
        Some("fixtures/trexplorer.bed"),
        None,
        output_dir.to_str().unwrap(),
    ).expect("Failed to build cache");
    
    // Read cache
    let cache = read_cache_from_dir(&output_dir).expect("Failed to read cache");
    let test_units = cache.test_units.as_ref().expect("Should have test_units");
    let genotypes = cache.genotypes.as_ref().expect("Should have genotypes");
    
    // Find repeat_TR001 unit
    let tu_id_col = test_units.column_by_name("id")
        .and_then(|c| c.as_any().downcast_ref::<StringArray>())
        .expect("Should have id column");
    let tu_source_col = test_units.column_by_name("source")
        .and_then(|c| c.as_any().downcast_ref::<StringArray>())
        .expect("Should have source column");
    
    let mut tr001_idx = None;
    for i in 0..tu_id_col.len() {
        if tu_id_col.value(i) == "repeat_TR001" {
            tr001_idx = Some(i);
            break;
        }
    }
    
    let tr001_idx = tr001_idx.expect("Should have repeat_TR001 unit");
    let tr001_source = tu_source_col.value(tr001_idx);
    
    // Source should be TrgtVcf (TRGT data merged in)
    assert_eq!(tr001_source, "TrgtVcf", "repeat_TR001 should have TrgtVcf source");
    
    // Check genotypes for TR001 have TRGT allele data
    let gt_unit_col = genotypes.column_by_name("unit_id")
        .and_then(|c| c.as_any().downcast_ref::<StringArray>())
        .expect("Should have unit_id column");
    let gt_allele_source_col = genotypes.column_by_name("allele_source")
        .and_then(|c| c.as_any().downcast_ref::<StringArray>())
        .expect("Should have allele_source column");
    let gt_allele1_col = genotypes.column_by_name("allele1")
        .and_then(|c| c.as_any().downcast_ref::<Int64Array>())
        .expect("Should have allele1 column");
    
    let mut found_trgt_genotype = false;
    for i in 0..gt_unit_col.len() {
        if gt_unit_col.value(i) == "repeat_TR001" {
            let allele_source = gt_allele_source_col.value(i);
            // Should be TrgtTrue (TRGT-derived alleles)
            assert_eq!(allele_source, "TrgtTrue", "Genotype should have TrgtTrue allele source");
            // Should have allele lengths from TRGT
            assert!(!gt_allele1_col.is_null(i), "Should have allele1 value from TRGT");
            found_trgt_genotype = true;
        }
    }
    
    assert!(found_trgt_genotype, "Should have genotypes for repeat_TR001 from TRGT");
    
    let _ = std::fs::remove_dir_all(&output_dir);
}

#[test]
fn test_comparison_mode_shadow_units() {
    use arrow::array::StringArray;
    use storm::build_cache_with_options;
    
    let output_dir = test_dir("comparison_mode");
    
    // Build cache WITH comparison mode enabled
    build_cache_with_options(
        "fixtures/sv_small.vcf",
        Some(&["fixtures/trgt_small.vcf"][..]),
        Some("fixtures/trexplorer.bed"),
        None,
        output_dir.to_str().unwrap(),
        true,  // comparison_mode = true
    ).expect("Failed to build cache with comparison mode");
    
    // Read cache
    let cache = read_cache_from_dir(&output_dir).expect("Failed to read cache");
    let test_units = cache.test_units.as_ref().expect("Should have test_units");
    
    let id_col = test_units.column_by_name("id")
        .and_then(|c| c.as_any().downcast_ref::<StringArray>())
        .expect("Should have id column");
    
    // Count shadow units
    let mut shadow_count = 0;
    let mut canonical_repeat_count = 0;
    
    for i in 0..id_col.len() {
        let id = id_col.value(i);
        if id.starts_with("shadow_") {
            shadow_count += 1;
        } else if id.starts_with("repeat_") {
            canonical_repeat_count += 1;
        }
    }
    
    // Should have both canonical repeat units AND shadow units
    assert!(canonical_repeat_count > 0, "Should have canonical repeat units");
    assert!(shadow_count > 0, "Comparison mode should emit shadow units");
    
    let _ = std::fs::remove_dir_all(&output_dir);
}

// ============================================================================
// Canonical Repeat Unit Tests
// ============================================================================

#[test]
fn test_canonical_repeat_units_sv_overlap() {
    use arrow::array::StringArray;
    
    let output_dir = test_dir("canonical_repeat_sv");
    
    // sv_with_overlap.vcf has:
    //   sv_overlap1: chr1:10010 (overlaps TR001 at chr1:10000-10050)
    //   sv_outside1: chr1:1000-1500 (outside all catalog regions)
    //   sv_overlap2: chr2:20010-20030 (overlaps TR003 at chr2:20000-20030)
    //   sv_outside2: chr3:99999 (outside all catalog regions)
    let stats = build_cache(
        "fixtures/sv_with_overlap.vcf",
        None,
        Some("fixtures/trexplorer.bed"),
        None,
        output_dir.to_str().unwrap(),
    ).expect("Failed to build cache");
    
    // Should have:
    // - 2 SV units (sv_outside1, sv_outside2) for SVs outside catalog
    // - 2 Repeat units (TR001, TR003) for catalog loci with SV overlap
    // SVs that overlap catalog (sv_overlap1, sv_overlap2) should NOT have separate units
    assert_eq!(stats.num_test_units, 4, 
        "Should have 4 test units: 2 SV (outside catalog) + 2 Repeat (catalog with overlap)");
    
    // Read cache and verify unit types
    let cache = read_cache_from_dir(&output_dir).expect("Failed to read cache");
    let test_units = cache.test_units.as_ref().expect("Should have test_units");
    
    let id_col = test_units.column_by_name("id")
        .and_then(|c| c.as_any().downcast_ref::<StringArray>())
        .expect("Should have id column");
    let type_col = test_units.column_by_name("unit_type")
        .and_then(|c| c.as_any().downcast_ref::<StringArray>())
        .expect("Should have unit_type column");
    
    let mut sv_units = Vec::new();
    let mut repeat_units = Vec::new();
    
    for i in 0..id_col.len() {
        let id = id_col.value(i);
        let unit_type = type_col.value(i);
        
        if unit_type == "Sv" {
            sv_units.push(id.to_string());
        } else if unit_type == "Repeat" {
            repeat_units.push(id.to_string());
        }
    }
    
    // SVs overlapping catalog should NOT appear as separate SV units
    assert!(!sv_units.iter().any(|id| id.contains("overlap")), 
        "SVs overlapping catalog should NOT have separate sv_<id> units");
    
    // SVs outside catalog SHOULD have sv_<id> units
    assert!(sv_units.iter().any(|id| id.contains("outside1")),
        "SV outside catalog (sv_outside1) should have sv_<id> unit");
    assert!(sv_units.iter().any(|id| id.contains("outside2")),
        "SV outside catalog (sv_outside2) should have sv_<id> unit");
    
    // Catalog loci with SV overlap should have canonical repeat units
    assert!(repeat_units.iter().any(|id| id.contains("TR001")),
        "Catalog TR001 (with SV overlap) should have repeat unit");
    assert!(repeat_units.iter().any(|id| id.contains("TR003")),
        "Catalog TR003 (with SV overlap) should have repeat unit");
    
    let _ = std::fs::remove_dir_all(&output_dir);
}

#[test]
fn test_canonical_repeat_units_trgt_merged() {
    use arrow::array::StringArray;
    
    let output_dir = test_dir("canonical_repeat_trgt");
    
    // Build cache with both SV and TRGT to test merging
    let stats = build_cache(
        "fixtures/sv_with_overlap.vcf",
        Some(&["fixtures/trgt_small.vcf"][..]),
        Some("fixtures/trexplorer.bed"),
        None,
        output_dir.to_str().unwrap(),
    ).expect("Failed to build cache");
    
    // Read cache and verify unit types
    let cache = read_cache_from_dir(&output_dir).expect("Failed to read cache");
    let test_units = cache.test_units.as_ref().expect("Should have test_units");
    
    let id_col = test_units.column_by_name("id")
        .and_then(|c| c.as_any().downcast_ref::<StringArray>())
        .expect("Should have id column");
    let type_col = test_units.column_by_name("unit_type")
        .and_then(|c| c.as_any().downcast_ref::<StringArray>())
        .expect("Should have unit_type column");
    
    let mut unit_types: std::collections::HashMap<String, String> = std::collections::HashMap::new();
    for i in 0..id_col.len() {
        unit_types.insert(id_col.value(i).to_string(), type_col.value(i).to_string());
    }
    
    // Verify canonical repeat units (not TrueRepeat or RepeatProxy)
    // When catalog exists, repeat units should be of type "Repeat"
    for (id, unit_type) in &unit_types {
        if id.starts_with("repeat_") {
            assert_eq!(unit_type, "Repeat", 
                "Unit {} should be canonical 'Repeat' type, not '{}'", id, unit_type);
        }
    }
    
    // There should be NO separate TrueRepeat units (TRGT is merged into catalog units)
    let true_repeat_count = unit_types.values()
        .filter(|t| *t == "TrueRepeat")
        .count();
    assert_eq!(true_repeat_count, 0, 
        "There should be no separate TrueRepeat units when catalog exists");
    
    // There should be NO RepeatProxy units (replaced by canonical Repeat)
    let proxy_count = unit_types.values()
        .filter(|t| *t == "RepeatProxy")
        .count();
    assert_eq!(proxy_count, 0,
        "There should be no RepeatProxy units (replaced by canonical Repeat)");
    
    let _ = std::fs::remove_dir_all(&output_dir);
}

#[test]
fn test_trgt_without_catalog_creates_repeat_units() {
    use arrow::array::StringArray;
    
    let output_dir = test_dir("trgt_no_catalog");
    
    // Build cache with TRGT but NO catalog
    let stats = build_cache(
        "fixtures/sv_small.vcf",
        Some(&["fixtures/trgt_small.vcf"][..]),
        None,  // No catalog
        None,
        output_dir.to_str().unwrap(),
    ).expect("Failed to build cache");
    
    // Should have SV units + repeat units from TRGT
    assert!(stats.num_test_units >= 6, "Should have at least 6 test units");
    
    // Read cache and verify
    let cache = read_cache_from_dir(&output_dir).expect("Failed to read cache");
    let test_units = cache.test_units.as_ref().expect("Should have test_units");
    
    let id_col = test_units.column_by_name("id")
        .and_then(|c| c.as_any().downcast_ref::<StringArray>())
        .expect("Should have id column");
    let type_col = test_units.column_by_name("unit_type")
        .and_then(|c| c.as_any().downcast_ref::<StringArray>())
        .expect("Should have unit_type column");
    
    let mut has_repeat_units = false;
    for i in 0..id_col.len() {
        let id = id_col.value(i);
        let unit_type = type_col.value(i);
        
        if id.starts_with("repeat_") {
            has_repeat_units = true;
            // Without catalog, TRGT creates canonical Repeat units
            assert_eq!(unit_type, "Repeat", 
                "TRGT units without catalog should be 'Repeat' type, got '{}'", unit_type);
        }
    }
    
    assert!(has_repeat_units, "Should have repeat units from TRGT data");
    
    let _ = std::fs::remove_dir_all(&output_dir);
}
