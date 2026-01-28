//! Arrow in-memory cache

use std::collections::HashMap;
use arrow::record_batch::RecordBatch;
use thiserror::Error;

use crate::testunit::TestUnit;
use crate::resolver::ResolvedGenotype;
use crate::catalog::CatalogEntry;
use super::{
    build_test_units_batch, build_genotypes_batch, build_catalog_batch,
    build_features_batch, Provenance,
};

/// Errors that can occur with Arrow cache
#[derive(Error, Debug)]
pub enum ArrowCacheError {
    #[error("Arrow error: {0}")]
    Arrow(#[from] arrow::error::ArrowError),
    #[error("Cache not populated: {0}")]
    NotPopulated(String),
}

/// In-memory Arrow cache
#[derive(Debug, Default)]
pub struct ArrowCache {
    /// Test units table
    pub test_units: Option<RecordBatch>,
    /// Genotypes table
    pub genotypes: Option<RecordBatch>,
    /// Catalog table
    pub catalog: Option<RecordBatch>,
    /// Features table
    pub features: Option<RecordBatch>,
    /// Provenance metadata
    pub provenance: Option<RecordBatch>,
}

impl ArrowCache {
    /// Create a new empty cache
    pub fn new() -> Self {
        Self::default()
    }

    /// Set test units from slice
    pub fn set_test_units(&mut self, units: &[TestUnit]) -> Result<(), ArrowCacheError> {
        self.test_units = Some(build_test_units_batch(units)?);
        Ok(())
    }

    /// Set genotypes
    pub fn set_genotypes(
        &mut self,
        genotypes: &[(String, Vec<ResolvedGenotype>)],
    ) -> Result<(), ArrowCacheError> {
        self.genotypes = Some(build_genotypes_batch(genotypes)?);
        Ok(())
    }

    /// Set catalog entries
    pub fn set_catalog(&mut self, entries: &[&CatalogEntry]) -> Result<(), ArrowCacheError> {
        self.catalog = Some(build_catalog_batch(entries)?);
        Ok(())
    }

    /// Set features
    pub fn set_features(
        &mut self,
        features: &[(String, String, HashMap<String, f64>)],
    ) -> Result<(), ArrowCacheError> {
        self.features = Some(build_features_batch(features)?);
        Ok(())
    }

    /// Set provenance
    pub fn set_provenance(&mut self, prov: &Provenance) -> Result<(), ArrowCacheError> {
        self.provenance = Some(prov.to_batch()?);
        Ok(())
    }

    /// Get test units batch
    pub fn get_test_units(&self) -> Result<&RecordBatch, ArrowCacheError> {
        self.test_units
            .as_ref()
            .ok_or_else(|| ArrowCacheError::NotPopulated("test_units".to_string()))
    }

    /// Get genotypes batch
    pub fn get_genotypes(&self) -> Result<&RecordBatch, ArrowCacheError> {
        self.genotypes
            .as_ref()
            .ok_or_else(|| ArrowCacheError::NotPopulated("genotypes".to_string()))
    }

    /// Get catalog batch
    pub fn get_catalog(&self) -> Result<&RecordBatch, ArrowCacheError> {
        self.catalog
            .as_ref()
            .ok_or_else(|| ArrowCacheError::NotPopulated("catalog".to_string()))
    }

    /// Get features batch
    pub fn get_features(&self) -> Result<&RecordBatch, ArrowCacheError> {
        self.features
            .as_ref()
            .ok_or_else(|| ArrowCacheError::NotPopulated("features".to_string()))
    }

    /// Get provenance batch
    pub fn get_provenance(&self) -> Result<&RecordBatch, ArrowCacheError> {
        self.provenance
            .as_ref()
            .ok_or_else(|| ArrowCacheError::NotPopulated("provenance".to_string()))
    }

    /// Check if all tables are populated
    pub fn is_complete(&self) -> bool {
        self.test_units.is_some()
            && self.genotypes.is_some()
            && self.catalog.is_some()
            && self.features.is_some()
            && self.provenance.is_some()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_arrow_cache_basic() {
        let mut cache = ArrowCache::new();
        
        // Initially not populated
        assert!(cache.get_test_units().is_err());
        assert!(!cache.is_complete());

        // Add test units
        let units = vec![TestUnit::from_sv("sv1", "chr1", 1000, 1500, "DEL")];
        cache.set_test_units(&units).expect("Failed to set test units");

        assert!(cache.get_test_units().is_ok());
        assert!(!cache.is_complete()); // Still missing other tables
    }

    #[test]
    fn test_arrow_cache_complete() {
        let mut cache = ArrowCache::new();

        // Set all tables
        let units = vec![TestUnit::from_sv("sv1", "chr1", 1000, 1500, "DEL")];
        cache.set_test_units(&units).unwrap();

        let gts = vec![(
            "sv1".to_string(),
            vec![ResolvedGenotype {
                sample_id: "S1".to_string(),
                is_present: true,
                presence_source: crate::resolver::PresenceSource::SvVcf,
                allele1: Some(30),
                allele2: Some(30),
                allele_source: crate::resolver::AlleleSource::Reference,
                raw_gt: None,
            }],
        )];
        cache.set_genotypes(&gts).unwrap();

        let catalog = crate::catalog::Catalog::from_bed("fixtures/trexplorer.bed").unwrap();
        let entries: Vec<&CatalogEntry> = catalog.entries.values().collect();
        cache.set_catalog(&entries).unwrap();

        let features: Vec<(String, String, HashMap<String, f64>)> = vec![];
        cache.set_features(&features).unwrap();

        let prov = Provenance::new();
        cache.set_provenance(&prov).unwrap();

        assert!(cache.is_complete());
    }
}
