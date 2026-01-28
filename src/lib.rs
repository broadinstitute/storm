//! STORM: Structural & Tandem-Repeat Optimized Regression Models
//!
//! A Rust crate with Python front-end for association testing of
//! structural variants and tandem repeats using long-read data.

#[cfg(feature = "python")]
use pyo3::prelude::*;

// Core modules
pub mod vcf;
pub mod catalog;
pub mod mapping;
pub mod testunit;
pub mod resolver;
pub mod cache;
pub mod explain;
pub mod plan;
pub mod encoding;
pub mod glm;

// Re-exports
pub use vcf::sv::{parse_sv_vcf, SvRecord, SvType, SvVcfError, Genotype};
pub use vcf::trgt::{parse_trgt_vcf, TrgtRecord, TrgtGenotype, TrgtVcfError};
pub use catalog::{
    Catalog, CatalogEntry, CatalogError,
    parse_trexplorer_bed, BedRecord, BedError,
    parse_trexplorer_json, JsonRecord, JsonError,
};
pub use mapping::{map_svs_to_catalog, map_svs_with_overlaps, SvMapping, MappingStats, compute_mapping_stats};
pub use testunit::{TestUnit, TestUnitType, TestUnitBuilder, DataSource};
pub use resolver::{Resolver, ResolvedGenotype, PresenceSource, AlleleSource, LocusSvAlleles};
pub use cache::{
    ArrowCache, ArrowCacheError, Provenance, write_parquet, read_parquet, ParquetError,
    parquet_cache::{write_cache_to_dir, read_cache_from_dir},
};
pub use explain::{explain_genotype, explain_locus, format_genotype_compact};
pub use plan::{Plan, PlanError, Rule, RuleCondition, Encoding, Model};
pub use encoding::{encode_s, encode_m, encode_d, encode_tail, encode_categorical, encode_binary, CategoricalBin, apply_encoding};
pub use glm::{GlmError, AssociationResult, linear_regression, logistic_regression, binomi_rare_test, firth_logistic_regression, categorical_regression, run_association};

// Example library function
#[cfg(feature = "python")]
#[pyfunction]
pub fn add(a: i32, b: i32) -> i32 {
    a + b
}

/// Returns the version of the storm package as a string.
///
/// This function returns the version that was set in Cargo.toml at compile time.
#[cfg(feature = "python")]
#[pyfunction]
fn _version() -> PyResult<String> {
    Ok(env!("CARGO_PKG_VERSION").to_string())
}

#[cfg(feature = "python")]
#[pymodule]
fn storm(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(add, m)?)?;
    m.add_function(wrap_pyfunction!(_version, m)?)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    #[test]
    #[cfg(feature = "python")]
    fn test_add() {
        assert_eq!(super::add(2, 2), 4);
        assert_eq!(super::add(-1, 1), 0);
        assert_eq!(super::add(0, 0), 0);
    }

    #[test]
    #[cfg(feature = "python")]
    fn test_version() {
        let version = super::_version().unwrap();
        assert!(!version.is_empty());
        // Version should be in format x.y.z
        assert!(version.matches(r"^\d+\.\d+\.\d+$"));
    }
}
