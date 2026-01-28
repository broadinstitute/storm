//! TRGT VCF parser
//!
//! Parses TRGT VCF files extracting TRID, AL (allele lengths), and GT fields.
//! This is a placeholder for criterion 2.

use std::collections::HashMap;
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
    /// Tandem repeat ID (from INFO TRID)
    pub trid: String,
    /// Sample allele data (sample_id -> TrgtGenotype)
    pub genotypes: HashMap<String, TrgtGenotype>,
}

/// TRGT genotype with allele lengths
#[derive(Debug, Clone)]
pub struct TrgtGenotype {
    /// Genotype allele indices
    pub gt: Vec<Option<u8>>,
    /// Allele lengths in bp
    pub allele_lengths: Vec<i64>,
}

/// Parse a TRGT VCF file and return records
pub fn parse_trgt_vcf<P: AsRef<Path>>(_path: P) -> Result<(Vec<String>, Vec<TrgtRecord>), TrgtVcfError> {
    // Placeholder - will be implemented for criterion 2
    Ok((Vec::new(), Vec::new()))
}
