//! VCF parsing modules for STORM
//!
//! This module handles parsing of:
//! - Structural Variant (SV) VCFs
//! - TRGT VCFs for tandem repeat genotypes

pub mod sv;
pub mod trgt;

pub use sv::{SvRecord, SvType, parse_sv_vcf};
pub use trgt::{TrgtRecord, parse_trgt_vcf};
