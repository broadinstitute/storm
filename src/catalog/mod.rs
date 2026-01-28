//! Tandem repeat catalog module
//!
//! Handles ingestion of TRExplorer BED and JSON files into a unified catalog.

mod bed;
mod json;

pub use bed::{parse_trexplorer_bed, BedRecord, BedError};
pub use json::{parse_trexplorer_json, JsonRecord, JsonError};

use std::collections::HashMap;
use serde::Serialize;
use thiserror::Error;

/// Errors that can occur with catalog operations
#[derive(Error, Debug)]
pub enum CatalogError {
    #[error("BED error: {0}")]
    Bed(#[from] BedError),
    #[error("JSON error: {0}")]
    Json(#[from] JsonError),
    #[error("Missing record for ID: {0}")]
    MissingRecord(String),
}

/// A unified catalog entry combining BED and JSON data
#[derive(Debug, Clone, Serialize)]
pub struct CatalogEntry {
    /// Repeat ID
    pub id: String,
    /// Chromosome
    pub chrom: String,
    /// Start position (0-based)
    pub start: u64,
    /// End position
    pub end: u64,
    /// Repeat motif
    pub motif: String,
    /// Reference repeat count
    pub repeat_count: u32,
    /// Strand
    pub strand: char,
    /// Repeat type (simple, complex, etc.)
    pub repeat_type: String,
    /// Pathogenicity status
    pub pathogenicity: String,
    /// Associated gene (from JSON)
    pub gene: Option<String>,
    /// Associated disease (from JSON)
    pub disease: Option<String>,
    /// Normal max repeat count (from JSON)
    pub normal_max: Option<u32>,
    /// Pathogenic minimum repeat count (from JSON)
    pub pathogenic_min: Option<u32>,
}

/// Unified repeat catalog
#[derive(Debug, Clone, Default, Serialize)]
pub struct Catalog {
    /// Entries indexed by ID
    pub entries: HashMap<String, CatalogEntry>,
    /// Entries indexed by genomic interval (chrom -> sorted entries)
    pub by_chrom: HashMap<String, Vec<String>>,
}

impl Catalog {
    /// Create a new empty catalog
    pub fn new() -> Self {
        Self::default()
    }

    /// Load catalog from BED file only
    pub fn from_bed(path: &str) -> Result<Self, CatalogError> {
        let bed_records = parse_trexplorer_bed(path)?;
        let mut catalog = Catalog::new();
        
        for rec in bed_records {
            let entry = CatalogEntry {
                id: rec.id.clone(),
                chrom: rec.chrom.clone(),
                start: rec.start,
                end: rec.end,
                motif: rec.motif.clone(),
                repeat_count: rec.repeat_count,
                strand: rec.strand,
                repeat_type: rec.repeat_type.clone(),
                pathogenicity: rec.pathogenicity.clone(),
                gene: None,
                disease: None,
                normal_max: None,
                pathogenic_min: None,
            };
            
            catalog.by_chrom
                .entry(rec.chrom.clone())
                .or_default()
                .push(rec.id.clone());
            catalog.entries.insert(rec.id, entry);
        }

        // Sort entries by position within each chromosome
        for entries in catalog.by_chrom.values_mut() {
            entries.sort_by(|a, b| {
                let ea = catalog.entries.get(a).unwrap();
                let eb = catalog.entries.get(b).unwrap();
                ea.start.cmp(&eb.start)
            });
        }

        Ok(catalog)
    }

    /// Load catalog from JSON file only
    pub fn from_json(path: &str) -> Result<Self, CatalogError> {
        let json_records = parse_trexplorer_json(path)?;
        let mut catalog = Catalog::new();
        
        for rec in json_records {
            let entry = CatalogEntry {
                id: rec.id.clone(),
                chrom: rec.chrom.clone(),
                start: rec.start,
                end: rec.end,
                motif: rec.motif.clone(),
                repeat_count: rec.repeat_count,
                strand: rec.strand.chars().next().unwrap_or('+'),
                repeat_type: rec.repeat_type.clone(),
                pathogenicity: rec.pathogenicity.clone(),
                gene: rec.gene.clone(),
                disease: rec.disease.clone(),
                normal_max: rec.normal_max,
                pathogenic_min: rec.pathogenic_min,
            };
            
            catalog.by_chrom
                .entry(rec.chrom.clone())
                .or_default()
                .push(rec.id.clone());
            catalog.entries.insert(rec.id, entry);
        }

        // Sort entries by position within each chromosome
        for entries in catalog.by_chrom.values_mut() {
            entries.sort_by(|a, b| {
                let ea = catalog.entries.get(a).unwrap();
                let eb = catalog.entries.get(b).unwrap();
                ea.start.cmp(&eb.start)
            });
        }

        Ok(catalog)
    }

    /// Join BED and JSON files into a unified catalog
    /// JSON provides additional annotations (gene, disease, thresholds)
    pub fn from_bed_and_json(bed_path: &str, json_path: &str) -> Result<Self, CatalogError> {
        let mut catalog = Self::from_bed(bed_path)?;
        let json_records = parse_trexplorer_json(json_path)?;
        
        // Create a lookup map for JSON records
        let json_map: HashMap<String, JsonRecord> = json_records
            .into_iter()
            .map(|r| (r.id.clone(), r))
            .collect();

        // Merge JSON data into existing entries
        for entry in catalog.entries.values_mut() {
            if let Some(json_rec) = json_map.get(&entry.id) {
                entry.gene = json_rec.gene.clone();
                entry.disease = json_rec.disease.clone();
                entry.normal_max = json_rec.normal_max;
                entry.pathogenic_min = json_rec.pathogenic_min;
            }
        }

        Ok(catalog)
    }

    /// Get number of entries
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Check if catalog is empty
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Get entry by ID
    pub fn get(&self, id: &str) -> Option<&CatalogEntry> {
        self.entries.get(id)
    }

    /// Find overlapping entries for a genomic region
    pub fn find_overlapping(&self, chrom: &str, start: u64, end: u64) -> Vec<&CatalogEntry> {
        let mut results = Vec::new();
        
        if let Some(ids) = self.by_chrom.get(chrom) {
            for id in ids {
                if let Some(entry) = self.entries.get(id) {
                    // Check for overlap: not (end1 <= start2 or end2 <= start1)
                    if !(end <= entry.start || entry.end <= start) {
                        results.push(entry);
                    }
                }
            }
        }
        
        results
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_catalog_from_bed() {
        let catalog = Catalog::from_bed("fixtures/trexplorer.bed").expect("Failed to load BED");
        
        assert_eq!(catalog.len(), 4);
        
        let entry = catalog.get("TR001").unwrap();
        assert_eq!(entry.chrom, "chr1");
        assert_eq!(entry.start, 10000);
        assert_eq!(entry.end, 10050);
        assert_eq!(entry.motif, "CAG");
        assert_eq!(entry.repeat_count, 16);
    }

    #[test]
    fn test_catalog_from_json() {
        let catalog = Catalog::from_json("fixtures/trexplorer.json").expect("Failed to load JSON");
        
        assert_eq!(catalog.len(), 4);
        
        let entry = catalog.get("TR001").unwrap();
        assert_eq!(entry.gene, Some("HTT".to_string()));
        assert_eq!(entry.normal_max, Some(26));
        assert_eq!(entry.pathogenic_min, Some(36));
    }

    #[test]
    fn test_catalog_join() {
        let catalog = Catalog::from_bed_and_json(
            "fixtures/trexplorer.bed",
            "fixtures/trexplorer.json"
        ).expect("Failed to join");
        
        assert_eq!(catalog.len(), 4);
        
        // Check that BED data is present
        let entry = catalog.get("TR001").unwrap();
        assert_eq!(entry.motif, "CAG");
        
        // Check that JSON data is merged
        assert_eq!(entry.gene, Some("HTT".to_string()));
        assert_eq!(entry.normal_max, Some(26));
    }

    #[test]
    fn test_find_overlapping() {
        let catalog = Catalog::from_bed("fixtures/trexplorer.bed").expect("Failed to load BED");
        
        // Query that overlaps TR001
        let results = catalog.find_overlapping("chr1", 10025, 10100);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].id, "TR001");
        
        // Query that overlaps nothing
        let results = catalog.find_overlapping("chr1", 0, 1000);
        assert_eq!(results.len(), 0);
        
        // Query that overlaps both TR001 and TR002
        let results = catalog.find_overlapping("chr1", 10025, 50050);
        assert_eq!(results.len(), 2);
    }
}
