//! TRExplorer BED file parser

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use thiserror::Error;

/// Errors that can occur during BED parsing
#[derive(Error, Debug)]
pub enum BedError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
    #[error("Parse error at line {line}: {message}")]
    ParseError { line: usize, message: String },
}

/// A parsed BED record from TRExplorer
#[derive(Debug, Clone)]
pub struct BedRecord {
    /// Chromosome
    pub chrom: String,
    /// Start position (0-based)
    pub start: u64,
    /// End position
    pub end: u64,
    /// Repeat ID
    pub id: String,
    /// Repeat motif
    pub motif: String,
    /// Repeat count
    pub repeat_count: u32,
    /// Strand (+/-)
    pub strand: char,
    /// Repeat type (simple, complex, etc.)
    pub repeat_type: String,
    /// Pathogenicity status
    pub pathogenicity: String,
}

/// Parse a TRExplorer BED file
pub fn parse_trexplorer_bed<P: AsRef<Path>>(path: P) -> Result<Vec<BedRecord>, BedError> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut records = Vec::new();

    for (line_num, line_result) in reader.lines().enumerate() {
        let line = line_result?;
        
        // Skip empty lines and comments
        if line.is_empty() || line.starts_with('#') || line.starts_with("track") || line.starts_with("browser") {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        
        if fields.len() < 9 {
            return Err(BedError::ParseError {
                line: line_num + 1,
                message: format!("Expected at least 9 fields, got {}", fields.len()),
            });
        }

        let chrom = fields[0].to_string();
        let start: u64 = fields[1].parse().map_err(|_| BedError::ParseError {
            line: line_num + 1,
            message: format!("Invalid start position: {}", fields[1]),
        })?;
        let end: u64 = fields[2].parse().map_err(|_| BedError::ParseError {
            line: line_num + 1,
            message: format!("Invalid end position: {}", fields[2]),
        })?;
        let id = fields[3].to_string();
        let motif = fields[4].to_string();
        let repeat_count: u32 = fields[5].parse().map_err(|_| BedError::ParseError {
            line: line_num + 1,
            message: format!("Invalid repeat count: {}", fields[5]),
        })?;
        let strand = fields[6].chars().next().unwrap_or('+');
        let repeat_type = fields[7].to_string();
        let pathogenicity = fields[8].to_string();

        records.push(BedRecord {
            chrom,
            start,
            end,
            id,
            motif,
            repeat_count,
            strand,
            repeat_type,
            pathogenicity,
        });
    }

    Ok(records)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_trexplorer_bed() {
        let records = parse_trexplorer_bed("fixtures/trexplorer.bed").expect("Failed to parse BED");
        
        assert_eq!(records.len(), 4);

        let rec = &records[0];
        assert_eq!(rec.chrom, "chr1");
        assert_eq!(rec.start, 10000);
        assert_eq!(rec.end, 10050);
        assert_eq!(rec.id, "TR001");
        assert_eq!(rec.motif, "CAG");
        assert_eq!(rec.repeat_count, 16);
        assert_eq!(rec.strand, '+');
        assert_eq!(rec.repeat_type, "simple");
        assert_eq!(rec.pathogenicity, "polymorphic");

        let rec = &records[2];
        assert_eq!(rec.id, "TR003");
        assert_eq!(rec.strand, '-');
    }
}
