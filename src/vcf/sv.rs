//! Structural Variant VCF parser
//!
//! Parses integrated SV VCF files extracting SVTYPE, SVLEN, and GT fields.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use serde::Serialize;
use thiserror::Error;

/// Errors that can occur during SV VCF parsing
#[derive(Error, Debug)]
pub enum SvVcfError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
    #[error("Missing required INFO field: {0}")]
    MissingInfoField(String),
    #[error("Invalid SVTYPE: {0}")]
    InvalidSvType(String),
    #[error("Parse error: {0}")]
    ParseError(String),
}

/// Structural variant types
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub enum SvType {
    Del,
    Ins,
    Dup,
    Inv,
    Bnd,
    Cnv,
    Other(String),
}

impl SvType {
    pub fn from_str(s: &str) -> Self {
        match s.to_uppercase().as_str() {
            "DEL" => SvType::Del,
            "INS" => SvType::Ins,
            "DUP" => SvType::Dup,
            "INV" => SvType::Inv,
            "BND" => SvType::Bnd,
            "CNV" => SvType::Cnv,
            other => SvType::Other(other.to_string()),
        }
    }
}

/// A parsed structural variant record
#[derive(Debug, Clone, Serialize)]
pub struct SvRecord {
    /// Chromosome
    pub chrom: String,
    /// 1-based position
    pub pos: u64,
    /// Variant ID
    pub id: String,
    /// Reference allele
    pub ref_allele: String,
    /// Alternate allele
    pub alt_allele: String,
    /// Structural variant type
    pub sv_type: SvType,
    /// Length of the SV (negative for deletions)
    pub sv_len: Option<i64>,
    /// End position
    pub end: Option<u64>,
    /// Sample genotypes (sample_id -> genotype string)
    pub genotypes: HashMap<String, Genotype>,
}

/// Parsed genotype
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct Genotype {
    /// Allele indices (0 = ref, 1+ = alt)
    pub alleles: Vec<Option<u8>>,
    /// Whether the genotype is phased
    pub phased: bool,
}

impl Genotype {
    /// Parse a genotype string like "0/1" or "1|0"
    pub fn from_str(s: &str) -> Self {
        let phased = s.contains('|');
        let sep = if phased { '|' } else { '/' };
        let alleles = s
            .split(sep)
            .map(|a| {
                if a == "." {
                    None
                } else {
                    a.parse::<u8>().ok()
                }
            })
            .collect();
        Genotype { alleles, phased }
    }

    /// Check if this is a reference homozygote (0/0)
    pub fn is_ref_hom(&self) -> bool {
        self.alleles.iter().all(|a| *a == Some(0))
    }

    /// Check if heterozygous
    pub fn is_het(&self) -> bool {
        let non_missing: Vec<_> = self.alleles.iter().filter_map(|a| *a).collect();
        non_missing.len() >= 2 && non_missing[0] != non_missing[1]
    }

    /// Check if alt homozygote
    pub fn is_alt_hom(&self) -> bool {
        let non_missing: Vec<_> = self.alleles.iter().filter_map(|a| *a).collect();
        non_missing.len() >= 2 && non_missing.iter().all(|&a| a > 0) && non_missing[0] == non_missing[1]
    }

    /// Check if carrier (has at least one alt allele)
    pub fn is_carrier(&self) -> bool {
        self.alleles.iter().any(|a| a.map(|v| v > 0).unwrap_or(false))
    }
}

/// Parse an SV VCF file and return records
pub fn parse_sv_vcf<P: AsRef<Path>>(path: P) -> Result<(Vec<String>, Vec<SvRecord>), SvVcfError> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    let mut sample_names: Vec<String> = Vec::new();
    let mut records: Vec<SvRecord> = Vec::new();

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
            .map_err(|_| SvVcfError::ParseError(format!("Invalid POS: {}", fields[1])))?;
        let id = fields[2].to_string();
        let ref_allele = fields[3].to_string();
        let alt_allele = fields[4].to_string();
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

        let sv_type = info_map
            .get("SVTYPE")
            .map(|s| SvType::from_str(s))
            .ok_or_else(|| SvVcfError::MissingInfoField("SVTYPE".to_string()))?;

        let sv_len: Option<i64> = info_map.get("SVLEN").and_then(|s| s.parse().ok());

        let end: Option<u64> = info_map.get("END").and_then(|s| s.parse().ok());

        // Parse genotypes
        let mut genotypes = HashMap::new();
        if fields.len() > 9 {
            let format_fields: Vec<&str> = fields[8].split(':').collect();
            let gt_index = format_fields.iter().position(|&f| f == "GT");

            if let Some(gt_idx) = gt_index {
                for (i, sample_field) in fields[9..].iter().enumerate() {
                    let sample_values: Vec<&str> = sample_field.split(':').collect();
                    if let Some(gt_str) = sample_values.get(gt_idx) {
                        let sample_name = sample_names
                            .get(i)
                            .cloned()
                            .unwrap_or_else(|| format!("SAMPLE{}", i));
                        genotypes.insert(sample_name, Genotype::from_str(gt_str));
                    }
                }
            }
        }

        records.push(SvRecord {
            chrom,
            pos,
            id,
            ref_allele,
            alt_allele,
            sv_type,
            sv_len,
            end,
            genotypes,
        });
    }

    Ok((sample_names, records))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_genotype_parsing() {
        let gt = Genotype::from_str("0/1");
        assert_eq!(gt.alleles, vec![Some(0), Some(1)]);
        assert!(!gt.phased);
        assert!(gt.is_het());

        let gt = Genotype::from_str("1|1");
        assert_eq!(gt.alleles, vec![Some(1), Some(1)]);
        assert!(gt.phased);
        assert!(gt.is_alt_hom());

        let gt = Genotype::from_str("0/0");
        assert!(gt.is_ref_hom());

        let gt = Genotype::from_str("./.");
        assert_eq!(gt.alleles, vec![None, None]);
    }

    #[test]
    fn test_sv_type_parsing() {
        assert_eq!(SvType::from_str("DEL"), SvType::Del);
        assert_eq!(SvType::from_str("del"), SvType::Del);
        assert_eq!(SvType::from_str("INS"), SvType::Ins);
        assert_eq!(SvType::from_str("DUP"), SvType::Dup);
        assert_eq!(SvType::from_str("INV"), SvType::Inv);
        assert_eq!(SvType::from_str("BND"), SvType::Bnd);
        assert_eq!(SvType::from_str("CNV"), SvType::Cnv);
        assert!(matches!(SvType::from_str("UNKNOWN"), SvType::Other(_)));
    }

    #[test]
    fn test_parse_sv_vcf() {
        let (samples, records) = parse_sv_vcf("fixtures/sv_small.vcf").expect("Failed to parse VCF");
        
        assert_eq!(samples, vec!["SAMPLE1", "SAMPLE2"]);
        assert_eq!(records.len(), 3);

        // Check first record (DEL)
        let rec = &records[0];
        assert_eq!(rec.chrom, "chr1");
        assert_eq!(rec.pos, 1000);
        assert_eq!(rec.id, "sv1");
        assert_eq!(rec.sv_type, SvType::Del);
        assert_eq!(rec.sv_len, Some(-500));
        assert_eq!(rec.end, Some(1500));
        
        let gt1 = rec.genotypes.get("SAMPLE1").unwrap();
        assert!(gt1.is_het());
        
        let gt2 = rec.genotypes.get("SAMPLE2").unwrap();
        assert!(gt2.is_alt_hom());

        // Check second record (INS)
        let rec = &records[1];
        assert_eq!(rec.sv_type, SvType::Ins);
        assert_eq!(rec.sv_len, Some(200));

        // Check third record (DUP)
        let rec = &records[2];
        assert_eq!(rec.sv_type, SvType::Dup);
        assert_eq!(rec.chrom, "chr2");
    }
}
