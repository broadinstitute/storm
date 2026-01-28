//! Results module for association testing output
//!
//! Handles writing results to Parquet with all required metadata.

use std::collections::HashMap;
use std::path::Path;
use std::sync::Arc;
use arrow::array::{
    ArrayRef, StringBuilder, Float64Builder, UInt64Builder, Int64Builder,
};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use thiserror::Error;

use crate::glm::AssociationResult;
use crate::cache::parquet_cache::write_parquet;

/// Errors from results operations
#[derive(Error, Debug)]
pub enum ResultsError {
    #[error("Arrow error: {0}")]
    Arrow(#[from] arrow::error::ArrowError),
    #[error("Parquet error: {0}")]
    Parquet(#[from] crate::cache::ParquetError),
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
}

/// Build Arrow RecordBatch from association results
pub fn build_results_batch(results: &[AssociationResult]) -> Result<RecordBatch, ResultsError> {
    let mut unit_id_builder = StringBuilder::new();
    let mut beta_builder = Float64Builder::new();
    let mut se_builder = Float64Builder::new();
    let mut statistic_builder = Float64Builder::new();
    let mut p_value_builder = Float64Builder::new();
    let mut n_samples_builder = UInt64Builder::new();
    let mut n_carriers_builder = UInt64Builder::new();
    let mut call_rate_builder = Float64Builder::new();
    let mut model_builder = StringBuilder::new();
    let mut encoding_builder = StringBuilder::new();

    for result in results {
        unit_id_builder.append_value(&result.unit_id);
        beta_builder.append_value(result.beta);
        se_builder.append_value(result.se);
        statistic_builder.append_value(result.statistic);
        p_value_builder.append_value(result.p_value);
        n_samples_builder.append_value(result.n_samples as u64);
        n_carriers_builder.append_value(result.n_carriers as u64);
        call_rate_builder.append_value(result.call_rate);
        model_builder.append_value(&result.model);
        encoding_builder.append_value(&result.encoding);
    }

    let schema = Schema::new(vec![
        Field::new("unit_id", DataType::Utf8, false),
        Field::new("beta", DataType::Float64, false),
        Field::new("se", DataType::Float64, false),
        Field::new("statistic", DataType::Float64, false),
        Field::new("p_value", DataType::Float64, false),
        Field::new("n_samples", DataType::UInt64, false),
        Field::new("n_carriers", DataType::UInt64, false),
        Field::new("call_rate", DataType::Float64, false),
        Field::new("model", DataType::Utf8, false),
        Field::new("encoding", DataType::Utf8, false),
    ]);

    Ok(RecordBatch::try_new(
        Arc::new(schema),
        vec![
            Arc::new(unit_id_builder.finish()) as ArrayRef,
            Arc::new(beta_builder.finish()) as ArrayRef,
            Arc::new(se_builder.finish()) as ArrayRef,
            Arc::new(statistic_builder.finish()) as ArrayRef,
            Arc::new(p_value_builder.finish()) as ArrayRef,
            Arc::new(n_samples_builder.finish()) as ArrayRef,
            Arc::new(n_carriers_builder.finish()) as ArrayRef,
            Arc::new(call_rate_builder.finish()) as ArrayRef,
            Arc::new(model_builder.finish()) as ArrayRef,
            Arc::new(encoding_builder.finish()) as ArrayRef,
        ],
    )?)
}

/// Write results to Parquet file
pub fn write_results_parquet<P: AsRef<Path>>(
    path: P,
    results: &[AssociationResult],
) -> Result<(), ResultsError> {
    let batch = build_results_batch(results)?;
    write_parquet(path, &batch)?;
    Ok(())
}

/// Summary statistics for results
#[derive(Debug, Clone)]
pub struct ResultsSummary {
    /// Total tests run
    pub n_tests: usize,
    /// Tests with p < 0.05
    pub n_nominal: usize,
    /// Tests with p < 5e-8 (genome-wide)
    pub n_genome_wide: usize,
    /// Minimum p-value
    pub min_p: f64,
    /// Average call rate
    pub mean_call_rate: f64,
    /// Tests by model
    pub by_model: HashMap<String, usize>,
    /// Tests by encoding
    pub by_encoding: HashMap<String, usize>,
}

impl ResultsSummary {
    /// Compute summary from results
    pub fn from_results(results: &[AssociationResult]) -> Self {
        let n_tests = results.len();
        let n_nominal = results.iter().filter(|r| r.p_value < 0.05).count();
        let n_genome_wide = results.iter().filter(|r| r.p_value < 5e-8).count();
        
        let min_p = results
            .iter()
            .map(|r| r.p_value)
            .fold(f64::INFINITY, f64::min);
        
        let mean_call_rate = if n_tests > 0 {
            results.iter().map(|r| r.call_rate).sum::<f64>() / n_tests as f64
        } else {
            0.0
        };
        
        let mut by_model = HashMap::new();
        let mut by_encoding = HashMap::new();
        
        for result in results {
            *by_model.entry(result.model.clone()).or_insert(0) += 1;
            *by_encoding.entry(result.encoding.clone()).or_insert(0) += 1;
        }
        
        ResultsSummary {
            n_tests,
            n_nominal,
            n_genome_wide,
            min_p,
            mean_call_rate,
            by_model,
            by_encoding,
        }
    }
}

/// Rare-variant ladder for progressive testing
#[derive(Debug, Clone)]
pub struct RareVariantLadder {
    /// Thresholds for carrier count
    pub thresholds: Vec<usize>,
    /// Models to use at each threshold
    pub models: Vec<crate::plan::Model>,
}

impl Default for RareVariantLadder {
    fn default() -> Self {
        RareVariantLadder {
            thresholds: vec![5, 10, 20, 50],
            models: vec![
                crate::plan::Model::BinomiRare,  // Very rare: exact test
                crate::plan::Model::Firth,       // Rare: Firth bias correction
                crate::plan::Model::Logistic,    // Uncommon: standard logistic
                crate::plan::Model::Linear,      // Common: linear
            ],
        }
    }
}

impl RareVariantLadder {
    /// Create a new ladder with custom thresholds
    pub fn new(thresholds: Vec<usize>, models: Vec<crate::plan::Model>) -> Self {
        RareVariantLadder { thresholds, models }
    }

    /// Select model based on carrier count
    pub fn select_model(&self, carrier_count: usize) -> crate::plan::Model {
        for (i, &threshold) in self.thresholds.iter().enumerate() {
            if carrier_count < threshold {
                return self.models.get(i).cloned().unwrap_or(crate::plan::Model::Linear);
            }
        }
        self.models.last().cloned().unwrap_or(crate::plan::Model::Linear)
    }
}

/// Covariate data for association testing
#[derive(Debug, Clone, Default)]
pub struct Covariates {
    /// Sample IDs in order
    pub sample_ids: Vec<String>,
    /// Covariate values (name -> values aligned with sample_ids)
    pub values: HashMap<String, Vec<f64>>,
    /// Principal components (PC1, PC2, ... aligned with sample_ids)
    pub pcs: Vec<Vec<f64>>,
}

impl Covariates {
    /// Create new empty covariates
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a covariate
    pub fn add(&mut self, name: &str, values: Vec<f64>) {
        self.values.insert(name.to_string(), values);
    }

    /// Add principal components
    pub fn add_pcs(&mut self, pcs: Vec<Vec<f64>>) {
        self.pcs = pcs;
    }

    /// Get all covariates as a matrix (for GLM)
    pub fn as_matrix(&self, names: &[String], num_pcs: usize) -> Vec<Vec<f64>> {
        let n = self.sample_ids.len();
        let mut result = Vec::new();

        // Add named covariates
        for name in names {
            if let Some(vals) = self.values.get(name) {
                result.push(vals.clone());
            }
        }

        // Add PCs
        for (i, pc) in self.pcs.iter().enumerate() {
            if i >= num_pcs {
                break;
            }
            result.push(pc.clone());
        }

        result
    }

    /// Get number of samples
    pub fn n_samples(&self) -> usize {
        self.sample_ids.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_build_results_batch() {
        let results = vec![
            AssociationResult {
                unit_id: "unit1".to_string(),
                beta: 0.5,
                se: 0.1,
                statistic: 5.0,
                p_value: 0.001,
                n_samples: 100,
                n_carriers: 20,
                call_rate: 0.95,
                model: "Linear".to_string(),
                encoding: "S".to_string(),
                metadata: HashMap::new(),
            },
            AssociationResult {
                unit_id: "unit2".to_string(),
                beta: -0.3,
                se: 0.15,
                statistic: -2.0,
                p_value: 0.046,
                n_samples: 100,
                n_carriers: 15,
                call_rate: 0.98,
                model: "Logistic".to_string(),
                encoding: "M".to_string(),
                metadata: HashMap::new(),
            },
        ];

        let batch = build_results_batch(&results).unwrap();
        assert_eq!(batch.num_rows(), 2);
        assert_eq!(batch.num_columns(), 10);
    }

    #[test]
    fn test_results_summary() {
        let results = vec![
            AssociationResult {
                unit_id: "unit1".to_string(),
                beta: 0.5,
                se: 0.1,
                statistic: 5.0,
                p_value: 0.001,
                n_samples: 100,
                n_carriers: 20,
                call_rate: 0.95,
                model: "Linear".to_string(),
                encoding: "S".to_string(),
                metadata: HashMap::new(),
            },
            AssociationResult {
                unit_id: "unit2".to_string(),
                beta: -0.3,
                se: 0.15,
                statistic: -2.0,
                p_value: 0.10,
                n_samples: 100,
                n_carriers: 15,
                call_rate: 0.98,
                model: "Logistic".to_string(),
                encoding: "M".to_string(),
                metadata: HashMap::new(),
            },
        ];

        let summary = ResultsSummary::from_results(&results);
        
        assert_eq!(summary.n_tests, 2);
        assert_eq!(summary.n_nominal, 1);
        assert_eq!(summary.n_genome_wide, 0);
        assert_eq!(summary.min_p, 0.001);
        assert!((summary.mean_call_rate - 0.965).abs() < 0.001);
    }

    #[test]
    fn test_rare_variant_ladder() {
        let ladder = RareVariantLadder::default();
        
        // Very rare
        assert_eq!(ladder.select_model(3), crate::plan::Model::BinomiRare);
        
        // Rare
        assert_eq!(ladder.select_model(7), crate::plan::Model::Firth);
        
        // Uncommon
        assert_eq!(ladder.select_model(15), crate::plan::Model::Logistic);
        
        // Common
        assert_eq!(ladder.select_model(100), crate::plan::Model::Linear);
    }

    #[test]
    fn test_covariates() {
        let mut cov = Covariates::new();
        cov.sample_ids = vec!["S1".to_string(), "S2".to_string(), "S3".to_string()];
        cov.add("age", vec![25.0, 30.0, 35.0]);
        cov.add("sex", vec![0.0, 1.0, 0.0]);
        cov.add_pcs(vec![
            vec![0.1, 0.2, 0.3],
            vec![0.4, 0.5, 0.6],
        ]);

        let matrix = cov.as_matrix(&["age".to_string()], 1);
        assert_eq!(matrix.len(), 2); // age + 1 PC
        assert_eq!(matrix[0], vec![25.0, 30.0, 35.0]);
        assert_eq!(matrix[1], vec![0.1, 0.2, 0.3]);
    }
}
