//! TRExplorer JSON file parser

use std::fs::File;
use std::io::BufReader;
use std::path::Path;
use serde::Deserialize;
use thiserror::Error;

/// Errors that can occur during JSON parsing
#[derive(Error, Debug)]
pub enum JsonError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
    #[error("JSON parse error: {0}")]
    Json(#[from] serde_json::Error),
}

/// A parsed JSON record from TRExplorer
#[derive(Debug, Clone, Deserialize)]
pub struct JsonRecord {
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
    /// Repeat count
    pub repeat_count: u32,
    /// Strand
    pub strand: String,
    /// Repeat type
    pub repeat_type: String,
    /// Pathogenicity status
    pub pathogenicity: String,
    /// Associated gene
    pub gene: Option<String>,
    /// Associated disease
    pub disease: Option<String>,
    /// Normal maximum repeat count
    pub normal_max: Option<u32>,
    /// Pathogenic minimum repeat count
    pub pathogenic_min: Option<u32>,
}

/// Parse a TRExplorer JSON file
pub fn parse_trexplorer_json<P: AsRef<Path>>(path: P) -> Result<Vec<JsonRecord>, JsonError> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let records: Vec<JsonRecord> = serde_json::from_reader(reader)?;
    Ok(records)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_trexplorer_json() {
        let records = parse_trexplorer_json("fixtures/trexplorer.json").expect("Failed to parse JSON");
        
        assert_eq!(records.len(), 4);

        let rec = &records[0];
        assert_eq!(rec.id, "TR001");
        assert_eq!(rec.chrom, "chr1");
        assert_eq!(rec.start, 10000);
        assert_eq!(rec.end, 10050);
        assert_eq!(rec.gene, Some("HTT".to_string()));
        assert_eq!(rec.disease, None);
        assert_eq!(rec.normal_max, Some(26));
        assert_eq!(rec.pathogenic_min, Some(36));

        let rec = &records[2];
        assert_eq!(rec.id, "TR003");
        assert_eq!(rec.gene, Some("FMR1".to_string()));
        assert_eq!(rec.disease, Some("Fragile X".to_string()));
    }
}
