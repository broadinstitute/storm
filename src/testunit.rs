//! TestUnit module - core abstraction for association testing
//!
//! A TestUnit represents something we can test for association:
//! - An SV (structural variant) directly
//! - A repeat locus proxied by overlapping SVs (INS/DEL sizes)
//! - A true repeat with direct TRGT measurements

use std::collections::HashMap;

/// The type of test unit
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TestUnitType {
    /// A standalone structural variant
    Sv,
    /// A repeat locus with SV-based proxy measurements
    RepeatProxy,
    /// A repeat locus with true TRGT measurements
    TrueRepeat,
}

/// Source of the genotype/measurement data
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DataSource {
    /// From an integrated SV VCF
    SvVcf,
    /// From a TRGT VCF
    TrgtVcf,
    /// Merged from both sources
    Merged,
}

/// A test unit for association analysis
#[derive(Debug, Clone)]
pub struct TestUnit {
    /// Unique identifier for this test unit
    pub id: String,
    /// Chromosome
    pub chrom: String,
    /// Start position
    pub start: u64,
    /// End position
    pub end: u64,
    /// Type of test unit
    pub unit_type: TestUnitType,
    /// Data source
    pub source: DataSource,
    /// Associated catalog entry ID (if any)
    pub catalog_id: Option<String>,
    /// Associated SV IDs (if any)
    pub sv_ids: Vec<String>,
    /// Repeat motif (if applicable)
    pub motif: Option<String>,
    /// Additional metadata
    pub metadata: HashMap<String, String>,
}

impl TestUnit {
    /// Create a new TestUnit for an SV
    pub fn from_sv(
        sv_id: &str,
        chrom: &str,
        start: u64,
        end: u64,
        sv_type: &str,
    ) -> Self {
        let mut metadata = HashMap::new();
        metadata.insert("sv_type".to_string(), sv_type.to_string());

        TestUnit {
            id: format!("sv_{}", sv_id),
            chrom: chrom.to_string(),
            start,
            end,
            unit_type: TestUnitType::Sv,
            source: DataSource::SvVcf,
            catalog_id: None,
            sv_ids: vec![sv_id.to_string()],
            motif: None,
            metadata,
        }
    }

    /// Create a TestUnit for a repeat-proxy locus (SV-based)
    pub fn from_repeat_proxy(
        catalog_id: &str,
        chrom: &str,
        start: u64,
        end: u64,
        motif: &str,
        sv_ids: Vec<String>,
    ) -> Self {
        let mut metadata = HashMap::new();
        metadata.insert("proxy_sv_count".to_string(), sv_ids.len().to_string());

        TestUnit {
            id: format!("proxy_{}", catalog_id),
            chrom: chrom.to_string(),
            start,
            end,
            unit_type: TestUnitType::RepeatProxy,
            source: DataSource::SvVcf,
            catalog_id: Some(catalog_id.to_string()),
            sv_ids,
            motif: Some(motif.to_string()),
            metadata,
        }
    }

    /// Create a TestUnit for a true repeat locus (TRGT-based)
    pub fn from_true_repeat(
        catalog_id: &str,
        chrom: &str,
        start: u64,
        end: u64,
        motif: &str,
        trid: &str,
    ) -> Self {
        let mut metadata = HashMap::new();
        metadata.insert("trid".to_string(), trid.to_string());

        TestUnit {
            id: format!("repeat_{}", catalog_id),
            chrom: chrom.to_string(),
            start,
            end,
            unit_type: TestUnitType::TrueRepeat,
            source: DataSource::TrgtVcf,
            catalog_id: Some(catalog_id.to_string()),
            sv_ids: Vec::new(),
            motif: Some(motif.to_string()),
            metadata,
        }
    }

    /// Check if this is an SV test unit
    pub fn is_sv(&self) -> bool {
        self.unit_type == TestUnitType::Sv
    }

    /// Check if this is a repeat-proxy test unit
    pub fn is_repeat_proxy(&self) -> bool {
        self.unit_type == TestUnitType::RepeatProxy
    }

    /// Check if this is a true repeat test unit
    pub fn is_true_repeat(&self) -> bool {
        self.unit_type == TestUnitType::TrueRepeat
    }
}

/// Builder for creating test units from parsed data
pub struct TestUnitBuilder {
    units: Vec<TestUnit>,
}

impl TestUnitBuilder {
    /// Create a new builder
    pub fn new() -> Self {
        TestUnitBuilder { units: Vec::new() }
    }

    /// Add test units from SV records
    pub fn add_from_svs(&mut self, svs: &[crate::vcf::sv::SvRecord]) {
        for sv in svs {
            let end = sv.end.unwrap_or_else(|| {
                if let Some(svlen) = sv.sv_len {
                    sv.pos + svlen.unsigned_abs()
                } else {
                    sv.pos + 1
                }
            });

            let sv_type = format!("{:?}", sv.sv_type);
            self.units.push(TestUnit::from_sv(
                &sv.id,
                &sv.chrom,
                sv.pos,
                end,
                &sv_type,
            ));
        }
    }

    /// Add test units from repeat-proxy mappings
    pub fn add_from_repeat_proxies(&mut self, mappings: &[crate::mapping::SvMapping]) {
        use std::collections::HashSet;

        // Group SVs by overlapping catalog entry
        let mut locus_svs: HashMap<String, (Vec<String>, &crate::catalog::CatalogEntry)> =
            HashMap::new();

        for mapping in mappings {
            for entry in &mapping.overlaps {
                locus_svs
                    .entry(entry.id.clone())
                    .or_insert_with(|| (Vec::new(), entry))
                    .0
                    .push(mapping.sv.id.clone());
            }
        }

        // Create test units for each locus with overlapping SVs
        for (catalog_id, (sv_ids, entry)) in locus_svs {
            // Deduplicate SV IDs
            let unique_ids: Vec<String> = sv_ids
                .into_iter()
                .collect::<HashSet<_>>()
                .into_iter()
                .collect();

            self.units.push(TestUnit::from_repeat_proxy(
                &catalog_id,
                &entry.chrom,
                entry.start,
                entry.end,
                &entry.motif,
                unique_ids,
            ));
        }
    }

    /// Add test units from TRGT records
    pub fn add_from_trgt(&mut self, records: &[crate::vcf::trgt::TrgtRecord]) {
        for rec in records {
            // Use TRID as both catalog_id and trid for now
            self.units.push(TestUnit::from_true_repeat(
                &rec.trid,
                &rec.chrom,
                rec.pos,
                rec.end,
                rec.motifs.first().map(|s| s.as_str()).unwrap_or(""),
                &rec.trid,
            ));
        }
    }

    /// Build and return the test units
    pub fn build(self) -> Vec<TestUnit> {
        self.units
    }

    /// Get current count
    pub fn len(&self) -> usize {
        self.units.len()
    }

    /// Check if empty
    pub fn is_empty(&self) -> bool {
        self.units.is_empty()
    }
}

impl Default for TestUnitBuilder {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_testunit_from_sv() {
        let unit = TestUnit::from_sv("sv1", "chr1", 1000, 1500, "DEL");
        
        assert_eq!(unit.id, "sv_sv1");
        assert_eq!(unit.chrom, "chr1");
        assert_eq!(unit.start, 1000);
        assert_eq!(unit.end, 1500);
        assert!(unit.is_sv());
        assert!(!unit.is_repeat_proxy());
        assert!(!unit.is_true_repeat());
        assert_eq!(unit.source, DataSource::SvVcf);
    }

    #[test]
    fn test_testunit_from_repeat_proxy() {
        let unit = TestUnit::from_repeat_proxy(
            "TR001",
            "chr1",
            10000,
            10050,
            "CAG",
            vec!["sv1".to_string(), "sv2".to_string()],
        );
        
        assert_eq!(unit.id, "proxy_TR001");
        assert!(unit.is_repeat_proxy());
        assert_eq!(unit.catalog_id, Some("TR001".to_string()));
        assert_eq!(unit.sv_ids.len(), 2);
        assert_eq!(unit.motif, Some("CAG".to_string()));
    }

    #[test]
    fn test_testunit_from_true_repeat() {
        let unit = TestUnit::from_true_repeat("TR001", "chr1", 10000, 10050, "CAG", "TR001");
        
        assert_eq!(unit.id, "repeat_TR001");
        assert!(unit.is_true_repeat());
        assert_eq!(unit.source, DataSource::TrgtVcf);
    }

    #[test]
    fn test_builder_from_svs() {
        use crate::vcf::sv::parse_sv_vcf;

        let (_, svs) = parse_sv_vcf("fixtures/sv_small.vcf").expect("Failed to parse SV VCF");
        
        let mut builder = TestUnitBuilder::new();
        builder.add_from_svs(&svs);
        
        let units = builder.build();
        assert_eq!(units.len(), 3);
        assert!(units.iter().all(|u| u.is_sv()));
    }

    #[test]
    fn test_builder_from_trgt() {
        use crate::vcf::trgt::parse_trgt_vcf;

        let (_, records) = parse_trgt_vcf("fixtures/trgt_small.vcf").expect("Failed to parse TRGT VCF");
        
        let mut builder = TestUnitBuilder::new();
        builder.add_from_trgt(&records);
        
        let units = builder.build();
        assert_eq!(units.len(), 3);
        assert!(units.iter().all(|u| u.is_true_repeat()));
    }
}
