//! Integration tests for STORM
//!
//! These tests verify the end-to-end functionality of the STORM pipeline.

use std::sync::atomic::{AtomicU64, Ordering};

use storm::{
    build_cache, verify_cache, explain_genotype, explain_locus,
    parse_sv_vcf, parse_trgt_vcf,
    Catalog,
    cache::read_cache_from_dir,
};

// Counter for unique test directories
static TEST_COUNTER: AtomicU64 = AtomicU64::new(0);

fn unique_test_dir(prefix: &str) -> std::path::PathBuf {
    let counter = TEST_COUNTER.fetch_add(1, Ordering::SeqCst);
    let pid = std::process::id();
    std::env::temp_dir().join(format!("storm_test_{}_{}_{}_{}", prefix, pid, counter, 
        std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos()))
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
    let output_dir = unique_test_dir("sv_only");
    
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
    
    // Cleanup
    let _ = std::fs::remove_dir_all(&output_dir);
}

#[test]
fn test_build_cache_from_fixtures() {
    let output_dir = unique_test_dir("full");

    let stats = build_cache(
        "fixtures/sv_small.vcf",
        Some("fixtures/trgt_small.vcf"),
        Some("fixtures/trexplorer.bed"),
        Some("fixtures/trexplorer.json"),
        output_dir.to_str().unwrap(),
    ).expect("Failed to build cache");

    // Should have SV units, proxy units, and TRGT units
    assert!(stats.num_test_units >= 6, "Should have at least 6 test units");
    assert_eq!(stats.num_samples, 2, "Should have 2 samples");
    assert!(stats.num_genotypes > 0, "Should have genotypes");

    // Verify all files were created
    assert!(output_dir.join("test_units.parquet").exists());
    assert!(output_dir.join("genotypes.parquet").exists());
    assert!(output_dir.join("catalog.parquet").exists());
    assert!(output_dir.join("features.parquet").exists());
    assert!(output_dir.join("provenance.parquet").exists());
    assert!(output_dir.join("provenance.json").exists());

    // Cleanup
    let _ = std::fs::remove_dir_all(&output_dir);
}

// ============================================================================
// Cache Verification Tests
// ============================================================================

#[test]
fn test_verify_cache() {
    let output_dir = unique_test_dir("verify");

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

    // Cleanup
    let _ = std::fs::remove_dir_all(&output_dir);
}

// ============================================================================
// Schema Tests
// ============================================================================

#[test]
fn test_test_units_schema() {
    let output_dir = unique_test_dir("schema_tu");

    build_cache(
        "fixtures/sv_small.vcf",
        None,
        None,
        None,
        output_dir.to_str().unwrap(),
    ).expect("Failed to build cache");

    let cache = read_cache_from_dir(&output_dir).expect("Failed to read cache");
    let batch = cache.test_units.expect("Should have test_units");

    // Check schema
    let schema = batch.schema();
    assert!(schema.column_with_name("id").is_some(), "Should have id column");
    assert!(schema.column_with_name("chrom").is_some(), "Should have chrom column");
    assert!(schema.column_with_name("start").is_some(), "Should have start column");
    assert!(schema.column_with_name("end").is_some(), "Should have end column");
    assert!(schema.column_with_name("unit_type").is_some(), "Should have unit_type column");

    // Check row count
    assert_eq!(batch.num_rows(), 3, "Should have 3 test units");

    // Cleanup
    let _ = std::fs::remove_dir_all(&output_dir);
}

#[test]
fn test_genotypes_schema() {
    let output_dir = unique_test_dir("schema_gt");

    build_cache(
        "fixtures/sv_small.vcf",
        None,
        None,
        None,
        output_dir.to_str().unwrap(),
    ).expect("Failed to build cache");

    let cache = read_cache_from_dir(&output_dir).expect("Failed to read cache");
    let batch = cache.genotypes.expect("Should have genotypes");

    // Check schema
    let schema = batch.schema();
    assert!(schema.column_with_name("unit_id").is_some(), "Should have unit_id column");
    assert!(schema.column_with_name("sample_id").is_some(), "Should have sample_id column");
    assert!(schema.column_with_name("is_present").is_some(), "Should have is_present column");
    assert!(schema.column_with_name("allele1").is_some(), "Should have allele1 column");
    assert!(schema.column_with_name("allele2").is_some(), "Should have allele2 column");

    // Check row count (3 units * 2 samples = 6)
    assert_eq!(batch.num_rows(), 6, "Should have 6 genotype rows");

    // Cleanup
    let _ = std::fs::remove_dir_all(&output_dir);
}

#[test]
fn test_catalog_schema() {
    let output_dir = unique_test_dir("schema_cat");

    build_cache(
        "fixtures/sv_small.vcf",
        None,
        Some("fixtures/trexplorer.bed"),
        Some("fixtures/trexplorer.json"),
        output_dir.to_str().unwrap(),
    ).expect("Failed to build cache");

    let cache = read_cache_from_dir(&output_dir).expect("Failed to read cache");
    let batch = cache.catalog.expect("Should have catalog");

    // Check schema
    let schema = batch.schema();
    assert!(schema.column_with_name("id").is_some(), "Should have id column");
    assert!(schema.column_with_name("chrom").is_some(), "Should have chrom column");
    assert!(schema.column_with_name("motif").is_some(), "Should have motif column");

    // Check row count
    assert!(batch.num_rows() > 0, "Should have catalog entries");

    // Cleanup
    let _ = std::fs::remove_dir_all(&output_dir);
}

#[test]
fn test_features_schema() {
    let output_dir = unique_test_dir("schema_feat");

    build_cache(
        "fixtures/sv_small.vcf",
        None,
        None,
        None,
        output_dir.to_str().unwrap(),
    ).expect("Failed to build cache");

    let cache = read_cache_from_dir(&output_dir).expect("Failed to read cache");
    let batch = cache.features.expect("Should have features");

    // Check schema
    let schema = batch.schema();
    assert!(schema.column_with_name("unit_id").is_some(), "Should have unit_id column");
    assert!(schema.column_with_name("sample_id").is_some(), "Should have sample_id column");

    // Check row count
    assert_eq!(batch.num_rows(), 6, "Should have 6 feature rows");

    // Cleanup
    let _ = std::fs::remove_dir_all(&output_dir);
}

// ============================================================================
// Explain Tests
// ============================================================================

#[test]
fn test_explain_locus_integration() {
    let output_dir = unique_test_dir("explain_locus");

    build_cache(
        "fixtures/sv_small.vcf",
        None,
        None,
        None,
        output_dir.to_str().unwrap(),
    ).expect("Failed to build cache");

    let cache = read_cache_from_dir(&output_dir).expect("Failed to read cache");

    let explanation = explain_locus(&cache, "sv_sv1")
        .expect("Failed to explain locus");

    assert!(explanation.contains("Locus Summary"), "Should have locus summary header");
    assert!(explanation.contains("sv_sv1"), "Should contain unit ID");
    assert!(explanation.contains("SAMPLE1"), "Should contain sample name");

    // Cleanup
    let _ = std::fs::remove_dir_all(&output_dir);
}

#[test]
fn test_explain_genotype_integration() {
    let output_dir = unique_test_dir("explain_gt");

    build_cache(
        "fixtures/sv_small.vcf",
        None,
        None,
        None,
        output_dir.to_str().unwrap(),
    ).expect("Failed to build cache");

    let cache = read_cache_from_dir(&output_dir).expect("Failed to read cache");

    let explanation = explain_genotype(&cache, "sv_sv1", "SAMPLE1")
        .expect("Failed to explain genotype");

    assert!(explanation.contains("Genotype Explanation"), "Should have genotype header");
    assert!(explanation.contains("sv_sv1"), "Should contain unit ID");
    assert!(explanation.contains("SAMPLE1"), "Should contain sample name");
    assert!(explanation.contains("Allele"), "Should contain allele info");

    // Cleanup
    let _ = std::fs::remove_dir_all(&output_dir);
}

// ============================================================================
// Provenance Tests
// ============================================================================

#[test]
fn test_provenance_metadata() {
    let output_dir = unique_test_dir("provenance");

    build_cache(
        "fixtures/sv_small.vcf",
        None,
        None,
        None,
        output_dir.to_str().unwrap(),
    ).expect("Failed to build cache");

    let prov_path = output_dir.join("provenance.json");
    let prov_str = std::fs::read_to_string(prov_path).expect("Failed to read provenance.json");
    let prov: serde_json::Value = serde_json::from_str(&prov_str).expect("Failed to parse JSON");

    assert!(prov.get("storm_version").is_some(), "Should have storm_version");
    assert!(prov.get("created_at").is_some(), "Should have created_at");
    assert!(prov.get("input_files").is_some(), "Should have input_files");
    assert!(prov.get("num_samples").is_some(), "Should have num_samples");
    assert!(prov.get("num_test_units").is_some(), "Should have num_test_units");

    // Cleanup
    let _ = std::fs::remove_dir_all(&output_dir);
}
