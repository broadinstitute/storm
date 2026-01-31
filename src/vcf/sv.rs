//! Structural Variant VCF parser
//!
//! Parses integrated SV VCF/BCF files extracting SVTYPE, SVLEN, and GT fields.
//! Supports both plain VCF and binary BCF formats.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;
use serde::Serialize;
use thiserror::Error;
use noodles::bcf;
use noodles::vcf::{self as noodles_vcf, variant::RecordBuf};
use noodles::vcf::variant::record_buf::info::field::Value as InfoValue;
use noodles::vcf::variant::record_buf::samples::sample::Value as SampleValue;
use noodles::vcf::variant::record::samples::series::value::Genotype as GenotypeIter;
use noodles::vcf::variant::record::samples::series::value::genotype::Phasing;
use noodles::bgzf;

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

/// Detect if a file is BCF format based on extension or magic bytes
fn is_bcf_file<P: AsRef<Path>>(path: P) -> Result<bool, SvVcfError> {
    let path = path.as_ref();
    
    // Check extension first
    if let Some(ext) = path.extension() {
        if ext == "bcf" {
            return Ok(true);
        }
    }
    
    // Check magic bytes: BCF starts with "BCF"
    let mut file = File::open(path)?;
    let mut magic = [0u8; 3];
    if file.read_exact(&mut magic).is_ok() && &magic == b"BCF" {
        return Ok(true);
    }
    
    Ok(false)
}

/// Detect if a file is gzip-compressed based on extension or magic bytes
fn is_gzip_file<P: AsRef<Path>>(path: P) -> Result<bool, SvVcfError> {
    let path = path.as_ref();
    
    // Check extension first
    if let Some(ext) = path.extension() {
        if ext == "gz" {
            return Ok(true);
        }
    }
    
    // Check magic bytes: gzip starts with 0x1f 0x8b
    let mut file = File::open(path)?;
    let mut magic = [0u8; 2];
    if file.read_exact(&mut magic).is_ok() && magic == [0x1f, 0x8b] {
        return Ok(true);
    }
    
    Ok(false)
}

/// Open a VCF file, handling both plain text and gzip/BGZF-compressed formats.
/// Uses BGZF reader for .gz so multi-block bgzip (typical for indexed .vcf.gz) is read fully.
fn open_vcf_reader<P: AsRef<Path>>(path: P) -> Result<Box<dyn BufRead>, SvVcfError> {
    let path = path.as_ref();

    if is_gzip_file(path)? {
        let file = File::open(path)?;
        let bgzf_reader = bgzf::Reader::new(file);
        Ok(Box::new(BufReader::new(bgzf_reader)))
    } else {
        let file = File::open(path)?;
        Ok(Box::new(BufReader::new(file)))
    }
}

/// Parse an SV VCF or BCF file and return records
/// Automatically detects format based on extension or magic bytes
pub fn parse_sv_vcf<P: AsRef<Path>>(path: P) -> Result<(Vec<String>, Vec<SvRecord>), SvVcfError> {
    if is_bcf_file(&path)? {
        parse_sv_bcf(path)
    } else {
        parse_sv_vcf_text(path)
    }
}

/// Parse an SV BCF file using noodles
fn parse_sv_bcf<P: AsRef<Path>>(path: P) -> Result<(Vec<String>, Vec<SvRecord>), SvVcfError> {
    let mut reader = File::open(&path)
        .map(bcf::io::Reader::new)
        .map_err(|e| SvVcfError::Io(e))?;
    
    let header = reader.read_header()
        .map_err(|e| SvVcfError::ParseError(format!("Failed to read BCF header: {}", e)))?;
    
    // Extract sample names from header
    let sample_names: Vec<String> = header.sample_names()
        .iter()
        .map(|s| s.to_string())
        .collect();
    
    let mut records: Vec<SvRecord> = Vec::new();
    
    for result in reader.record_bufs(&header) {
        let record = result
            .map_err(|e| SvVcfError::ParseError(format!("Failed to read BCF record: {}", e)))?;
        
        if let Some(sv_record) = parse_noodles_record(&record, &sample_names, &header)? {
            records.push(sv_record);
        }
    }
    
    Ok((sample_names, records))
}

/// Parse a noodles RecordBuf into an SvRecord
fn parse_noodles_record(
    record: &RecordBuf,
    sample_names: &[String],
    _header: &noodles_vcf::Header,
) -> Result<Option<SvRecord>, SvVcfError> {
    let chrom = record.reference_sequence_name().to_string();
    let pos = record.variant_start()
        .map(|p| p.get() as u64)
        .unwrap_or(0);
    
    // IDs - as_ref() to get the slice of ids
    let id = record.ids()
        .as_ref()
        .first()
        .map(|id| id.to_string())
        .unwrap_or_else(|| ".".to_string());
    
    let ref_allele = record.reference_bases().to_string();
    
    // Alt alleles - as_ref() to get the slice
    let alt_allele = record.alternate_bases()
        .as_ref()
        .first()
        .map(|a| a.to_string())
        .unwrap_or_else(|| ".".to_string());
    
    // Parse INFO fields - use get() with string key
    let info = record.info();
    
    // Get SVTYPE
    let sv_type = match info.get("SVTYPE") {
        Some(Some(InfoValue::String(s))) => SvType::from_str(s),
        _ => return Err(SvVcfError::MissingInfoField("SVTYPE".to_string())),
    };
    
    // Get SVLEN
    let sv_len: Option<i64> = match info.get("SVLEN") {
        Some(Some(InfoValue::Integer(n))) => Some(*n as i64),
        Some(Some(InfoValue::Array(arr))) => {
            use noodles::vcf::variant::record_buf::info::field::value::Array;
            match arr {
                Array::Integer(vals) => vals.first().copied().flatten().map(|n| n as i64),
                _ => None,
            }
        }
        _ => None,
    };
    
    // Get END
    let end: Option<u64> = match info.get("END") {
        Some(Some(InfoValue::Integer(n))) => Some(*n as u64),
        _ => None,
    };
    
    // Parse genotypes
    let mut genotypes = HashMap::new();
    let samples = record.samples();
    
    for (i, sample_name) in sample_names.iter().enumerate() {
        if let Some(sample) = samples.get_index(i) {
            if let Some(Some(gt_value)) = sample.get("GT") {
                if let SampleValue::Genotype(gt) = gt_value {
                    // Build genotype string from alleles using the GenotypeIter trait
                    let mut alleles: Vec<String> = Vec::new();
                    let mut phased = false;
                    let mut first = true;
                    for result in GenotypeIter::iter(&gt) {
                        if let Ok((pos, phasing)) = result {
                            if !first && phasing == Phasing::Phased {
                                phased = true;
                            }
                            match pos {
                                Some(idx) => alleles.push(idx.to_string()),
                                None => alleles.push(".".to_string()),
                            }
                            first = false;
                        }
                    }
                    let sep = if phased { "|" } else { "/" };
                    let gt_str = alleles.join(sep);
                    genotypes.insert(sample_name.clone(), Genotype::from_str(&gt_str));
                }
            }
        }
    }
    
    Ok(Some(SvRecord {
        chrom,
        pos,
        id,
        ref_allele,
        alt_allele,
        sv_type,
        sv_len,
        end,
        genotypes,
    }))
}

/// Parse an SV VCF file (plain text format, supports .vcf and .vcf.gz)
fn parse_sv_vcf_text<P: AsRef<Path>>(path: P) -> Result<(Vec<String>, Vec<SvRecord>), SvVcfError> {
    let reader = open_vcf_reader(path)?;

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
