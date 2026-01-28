//! TRGT VCF parser
//!
//! Parses TRGT VCF files extracting TRID, AL (allele lengths), and GT fields.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use thiserror::Error;

/// Errors that can occur during TRGT VCF parsing
#[derive(Error, Debug)]
pub enum TrgtVcfError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
    #[error("Missing required field: {0}")]
    MissingField(String),
    #[error("Parse error: {0}")]
    ParseError(String),
}

/// A parsed TRGT record
#[derive(Debug, Clone)]
pub struct TrgtRecord {
    /// Chromosome
    pub chrom: String,
    /// 1-based position
    pub pos: u64,
    /// End position
    pub end: u64,
    /// Tandem repeat ID (from INFO TRID)
    pub trid: String,
    /// Repeat motifs
    pub motifs: Vec<String>,
    /// Repeat structure
    pub struc: Option<String>,
    /// Reference allele
    pub ref_allele: String,
    /// Alternate alleles
    pub alt_alleles: Vec<String>,
    /// Sample allele data (sample_id -> TrgtGenotype)
    pub genotypes: HashMap<String, TrgtGenotype>,
}

/// TRGT genotype with allele lengths
#[derive(Debug, Clone, PartialEq)]
pub struct TrgtGenotype {
    /// Genotype allele indices
    pub gt: Vec<Option<u8>>,
    /// Allele lengths in bp
    pub allele_lengths: Vec<i64>,
    /// Whether genotype is phased
    pub phased: bool,
}

impl TrgtGenotype {
    /// Parse from GT and AL fields
    pub fn from_fields(gt_str: &str, al_str: &str) -> Self {
        let phased = gt_str.contains('|');
        let sep = if phased { '|' } else { '/' };
        
        let gt: Vec<Option<u8>> = gt_str
            .split(sep)
            .map(|a| {
                if a == "." {
                    None
                } else {
                    a.parse::<u8>().ok()
                }
            })
            .collect();

        let allele_lengths: Vec<i64> = al_str
            .split(',')
            .filter_map(|s| s.parse().ok())
            .collect();

        TrgtGenotype {
            gt,
            allele_lengths,
            phased,
        }
    }

    /// Get the diploid allele lengths based on GT indices
    pub fn diploid_lengths(&self) -> Option<(i64, i64)> {
        if self.allele_lengths.len() >= 2 && self.gt.len() >= 2 {
            Some((self.allele_lengths[0], self.allele_lengths[1]))
        } else {
            None
        }
    }
}

/// Parse a TRGT VCF file and return records
pub fn parse_trgt_vcf<P: AsRef<Path>>(path: P) -> Result<(Vec<String>, Vec<TrgtRecord>), TrgtVcfError> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    let mut sample_names: Vec<String> = Vec::new();
    let mut records: Vec<TrgtRecord> = Vec::new();

    for line_result in reader.lines() {
        let line = line_result?;

        // Skip comment lines but parse header for sample names
        if line.starts_with("##") {
            continue;
        }

        if line.starts_with("#CHROM") {
            // Parse header line for sample names
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() > 9 {
                sample_names = fields[9..].iter().map(|s| s.to_string()).collect();
            }
            continue;
        }

        // Parse data line
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 8 {
            continue;
        }

        let chrom = fields[0].to_string();
        let pos: u64 = fields[1]
            .parse()
            .map_err(|_| TrgtVcfError::ParseError(format!("Invalid POS: {}", fields[1])))?;
        let ref_allele = fields[3].to_string();
        let alt_alleles: Vec<String> = fields[4].split(',').map(|s| s.to_string()).collect();
        let info = fields[7];

        // Parse INFO fields
        let info_map: HashMap<&str, &str> = info
            .split(';')
            .filter_map(|kv| {
                let mut parts = kv.splitn(2, '=');
                match (parts.next(), parts.next()) {
                    (Some(k), Some(v)) => Some((k, v)),
                    _ => None,
                }
            })
            .collect();

        let trid = info_map
            .get("TRID")
            .map(|s| s.to_string())
            .ok_or_else(|| TrgtVcfError::MissingField("TRID".to_string()))?;

        let end: u64 = info_map
            .get("END")
            .and_then(|s| s.parse().ok())
            .unwrap_or(pos);

        let motifs: Vec<String> = info_map
            .get("MOTIFS")
            .map(|s| s.split(',').map(|m| m.to_string()).collect())
            .unwrap_or_default();

        let struc = info_map.get("STRUC").map(|s| s.to_string());

        // Parse genotypes
        let mut genotypes = HashMap::new();
        if fields.len() > 9 {
            let format_fields: Vec<&str> = fields[8].split(':').collect();
            let gt_index = format_fields.iter().position(|&f| f == "GT");
            let al_index = format_fields.iter().position(|&f| f == "AL");

            for (i, sample_field) in fields[9..].iter().enumerate() {
                let sample_values: Vec<&str> = sample_field.split(':').collect();
                
                let gt_str = gt_index
                    .and_then(|idx| sample_values.get(idx))
                    .copied()
                    .unwrap_or("./.");
                
                let al_str = al_index
                    .and_then(|idx| sample_values.get(idx))
                    .copied()
                    .unwrap_or("0,0");

                let sample_name = sample_names
                    .get(i)
                    .cloned()
                    .unwrap_or_else(|| format!("SAMPLE{}", i));
                
                genotypes.insert(sample_name, TrgtGenotype::from_fields(gt_str, al_str));
            }
        }

        records.push(TrgtRecord {
            chrom,
            pos,
            end,
            trid,
            motifs,
            struc,
            ref_allele,
            alt_alleles,
            genotypes,
        });
    }

    Ok((sample_names, records))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_trgt_genotype_parsing() {
        let gt = TrgtGenotype::from_fields("0/1", "9,12");
        assert_eq!(gt.gt, vec![Some(0), Some(1)]);
        assert_eq!(gt.allele_lengths, vec![9, 12]);
        assert!(!gt.phased);
        assert_eq!(gt.diploid_lengths(), Some((9, 12)));

        let gt = TrgtGenotype::from_fields("1|1", "15,15");
        assert!(gt.phased);
        assert_eq!(gt.diploid_lengths(), Some((15, 15)));
    }

    #[test]
    fn test_parse_trgt_vcf() {
        let (samples, records) = parse_trgt_vcf("fixtures/trgt_small.vcf").expect("Failed to parse TRGT VCF");
        
        assert_eq!(samples, vec!["SAMPLE1", "SAMPLE2"]);
        assert_eq!(records.len(), 3);

        // Check first record
        let rec = &records[0];
        assert_eq!(rec.chrom, "chr1");
        assert_eq!(rec.pos, 10000);
        assert_eq!(rec.trid, "TR001");
        assert_eq!(rec.end, 10009);
        assert_eq!(rec.motifs, vec!["CAG"]);
        assert_eq!(rec.struc, Some("(CAG)n".to_string()));
        
        let gt1 = rec.genotypes.get("SAMPLE1").unwrap();
        assert_eq!(gt1.gt, vec![Some(0), Some(1)]);
        assert_eq!(gt1.allele_lengths, vec![9, 12]);

        let gt2 = rec.genotypes.get("SAMPLE2").unwrap();
        assert_eq!(gt2.gt, vec![Some(1), Some(2)]);
        assert_eq!(gt2.allele_lengths, vec![12, 15]);

        // Check second record
        let rec = &records[1];
        assert_eq!(rec.trid, "TR002");
        assert_eq!(rec.motifs, vec!["AT"]);

        // Check third record
        let rec = &records[2];
        assert_eq!(rec.trid, "TR003");
    }
}
