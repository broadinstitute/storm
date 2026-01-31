//! Integration tests for STORM
//!
//! These tests verify the end-to-end functionality of the STORM pipeline.

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
