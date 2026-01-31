//! TRExplorer BED file parser
//! Supports plain .bed and gzip/BGZF .bed.gz.

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use thiserror::Error;
use noodles::bgzf;

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

fn open_bed_reader<P: AsRef<Path>>(path: P) -> Result<Box<dyn BufRead>, BedError> {
    let path = path.as_ref();
    let file = File::open(path)?;
    let reader: Box<dyn BufRead> = if path.extension().map(|e| e == "gz").unwrap_or(false) {
        Box::new(BufReader::new(bgzf::Reader::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };
    Ok(reader)
}

/// Parse a TRExplorer BED file (plain .bed or .bed.gz)
pub fn parse_trexplorer_bed<P: AsRef<Path>>(path: P) -> Result<Vec<BedRecord>, BedError> {
    let reader = open_bed_reader(path)?;
    let mut records = Vec::new();

    for (line_num, line_result) in reader.lines().enumerate() {
        let line = line_result?;
        
        // Skip empty lines and comments
        if line.is_empty() || line.starts_with('#') || line.starts_with("track") || line.starts_with("browser") {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        
        // Require at least chrom, start, end, name (BED4). Optional columns 5–9 match extended TRExplorer format.
        if fields.len() < 4 {
            return Err(BedError::ParseError {
                line: line_num + 1,
                message: format!("Expected at least 4 fields (chrom, start, end, name), got {}", fields.len()),
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
        // Extended TRExplorer format has 9+ fields: motif, repeat_count, strand, repeat_type, pathogenicity
        let (motif, repeat_count, strand, repeat_type, pathogenicity) = if fields.len() >= 9 {
            let rc: u32 = fields[5].parse().map_err(|_| BedError::ParseError {
                line: line_num + 1,
                message: format!("Invalid repeat count: {}", fields[5]),
            })?;
            (
                fields[4].to_string(),
                rc,
                fields[6].chars().next().unwrap_or('+'),
                fields[7].to_string(),
                fields[8].to_string(),
            )
        } else {
            // BED4/BED5: no motif column (column 4 is name in BED4, score in BED5)
            (
                String::new(),
                fields.get(5).and_then(|s| s.parse().ok()).unwrap_or(0u32),
                fields.get(6).and_then(|s| s.chars().next()).unwrap_or('+'),
                fields.get(7).map(|s| s.to_string()).unwrap_or_else(|| "unknown".to_string()),
                fields.get(8).map(|s| s.to_string()).unwrap_or_else(|| "unknown".to_string()),
            )
        };

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

    #[test]
    fn test_parse_bed5_minimal() {
        // Real TRExplorer catalog may use BED5 (chrom, start, end, name, score)
        use std::io::Write;
        let dir = std::env::temp_dir().join("storm_bed5_test");
        std::fs::create_dir_all(&dir).ok();
        let path = dir.join("five_field.bed");
        let mut f = std::fs::File::create(&path).unwrap();
        writeln!(f, "chr6\t31803187\t31803200\trepeat_1\t0").unwrap();
        writeln!(f, "chr6\t32000000\t32000015\trepeat_2\t0").unwrap();
        f.sync_all().unwrap();
        drop(f);

        let records = parse_trexplorer_bed(&path).expect("Failed to parse 5-field BED");
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].chrom, "chr6");
        assert_eq!(records[0].start, 31803187);
        assert_eq!(records[0].end, 31803200);
        assert_eq!(records[0].id, "repeat_1");
        assert_eq!(records[0].motif, "");
        assert_eq!(records[0].repeat_count, 0);
        assert_eq!(records[0].strand, '+');
        assert_eq!(records[0].repeat_type, "unknown");
        assert_eq!(records[0].pathogenicity, "unknown");
    }
}
