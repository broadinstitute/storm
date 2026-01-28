//! Cache module for Arrow/Parquet storage
//!
//! Provides caching of:
//! - Test units table
//! - Genotypes table
//! - Catalog table
//! - Features table
//! - Provenance metadata

mod arrow_cache;
pub mod parquet_cache;

pub use arrow_cache::{ArrowCache, ArrowCacheError};
pub use parquet_cache::{write_parquet, read_parquet, ParquetError, write_cache_to_dir, read_cache_from_dir};

use std::collections::HashMap;
use std::sync::Arc;
use arrow::array::{
    ArrayRef, StringArray,
    StringBuilder, Int64Builder, UInt64Builder, BooleanBuilder,
};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use chrono::Utc;

use crate::testunit::TestUnit;
use crate::resolver::ResolvedGenotype;
use crate::catalog::CatalogEntry;

/// Build Arrow RecordBatch for test units
pub fn build_test_units_batch(units: &[TestUnit]) -> Result<RecordBatch, arrow::error::ArrowError> {
    let mut id_builder = StringBuilder::new();
    let mut chrom_builder = StringBuilder::new();
    let mut start_builder = UInt64Builder::new();
    let mut end_builder = UInt64Builder::new();
    let mut unit_type_builder = StringBuilder::new();
    let mut source_builder = StringBuilder::new();
    let mut catalog_id_builder = StringBuilder::new();
    let mut motif_builder = StringBuilder::new();

    for unit in units {
        id_builder.append_value(&unit.id);
        chrom_builder.append_value(&unit.chrom);
        start_builder.append_value(unit.start);
        end_builder.append_value(unit.end);
        unit_type_builder.append_value(format!("{:?}", unit.unit_type));
        source_builder.append_value(format!("{:?}", unit.source));
        catalog_id_builder.append_option(unit.catalog_id.as_deref());
        motif_builder.append_option(unit.motif.as_deref());
    }

    let schema = Schema::new(vec![
        Field::new("id", DataType::Utf8, false),
        Field::new("chrom", DataType::Utf8, false),
        Field::new("start", DataType::UInt64, false),
        Field::new("end", DataType::UInt64, false),
        Field::new("unit_type", DataType::Utf8, false),
        Field::new("source", DataType::Utf8, false),
        Field::new("catalog_id", DataType::Utf8, true),
        Field::new("motif", DataType::Utf8, true),
    ]);

    RecordBatch::try_new(
        Arc::new(schema),
        vec![
            Arc::new(id_builder.finish()) as ArrayRef,
            Arc::new(chrom_builder.finish()) as ArrayRef,
            Arc::new(start_builder.finish()) as ArrayRef,
            Arc::new(end_builder.finish()) as ArrayRef,
            Arc::new(unit_type_builder.finish()) as ArrayRef,
            Arc::new(source_builder.finish()) as ArrayRef,
            Arc::new(catalog_id_builder.finish()) as ArrayRef,
            Arc::new(motif_builder.finish()) as ArrayRef,
        ],
    )
}

/// Build Arrow RecordBatch for genotypes
pub fn build_genotypes_batch(
    genotypes: &[(String, Vec<ResolvedGenotype>)], // (unit_id, genotypes)
) -> Result<RecordBatch, arrow::error::ArrowError> {
    let mut unit_id_builder = StringBuilder::new();
    let mut sample_id_builder = StringBuilder::new();
    let mut is_present_builder = BooleanBuilder::new();
    let mut presence_source_builder = StringBuilder::new();
    let mut allele1_builder = Int64Builder::new();
    let mut allele2_builder = Int64Builder::new();
    let mut allele_source_builder = StringBuilder::new();

    for (unit_id, gts) in genotypes {
        for gt in gts {
            unit_id_builder.append_value(unit_id);
            sample_id_builder.append_value(&gt.sample_id);
            is_present_builder.append_value(gt.is_present);
            presence_source_builder.append_value(format!("{:?}", gt.presence_source));
            allele1_builder.append_option(gt.allele1);
            allele2_builder.append_option(gt.allele2);
            allele_source_builder.append_value(format!("{:?}", gt.allele_source));
        }
    }

    let schema = Schema::new(vec![
        Field::new("unit_id", DataType::Utf8, false),
        Field::new("sample_id", DataType::Utf8, false),
        Field::new("is_present", DataType::Boolean, false),
        Field::new("presence_source", DataType::Utf8, false),
        Field::new("allele1", DataType::Int64, true),
        Field::new("allele2", DataType::Int64, true),
        Field::new("allele_source", DataType::Utf8, false),
    ]);

    RecordBatch::try_new(
        Arc::new(schema),
        vec![
            Arc::new(unit_id_builder.finish()) as ArrayRef,
            Arc::new(sample_id_builder.finish()) as ArrayRef,
            Arc::new(is_present_builder.finish()) as ArrayRef,
            Arc::new(presence_source_builder.finish()) as ArrayRef,
            Arc::new(allele1_builder.finish()) as ArrayRef,
            Arc::new(allele2_builder.finish()) as ArrayRef,
            Arc::new(allele_source_builder.finish()) as ArrayRef,
        ],
    )
}

/// Build Arrow RecordBatch for catalog entries
pub fn build_catalog_batch(entries: &[&CatalogEntry]) -> Result<RecordBatch, arrow::error::ArrowError> {
    let mut id_builder = StringBuilder::new();
    let mut chrom_builder = StringBuilder::new();
    let mut start_builder = UInt64Builder::new();
    let mut end_builder = UInt64Builder::new();
    let mut motif_builder = StringBuilder::new();
    let mut repeat_count_builder = UInt64Builder::new();
    let mut gene_builder = StringBuilder::new();
    let mut disease_builder = StringBuilder::new();
    let mut normal_max_builder = UInt64Builder::new();
    let mut pathogenic_min_builder = UInt64Builder::new();

    for entry in entries {
        id_builder.append_value(&entry.id);
        chrom_builder.append_value(&entry.chrom);
        start_builder.append_value(entry.start);
        end_builder.append_value(entry.end);
        motif_builder.append_value(&entry.motif);
        repeat_count_builder.append_value(entry.repeat_count as u64);
        gene_builder.append_option(entry.gene.as_deref());
        disease_builder.append_option(entry.disease.as_deref());
        normal_max_builder.append_option(entry.normal_max.map(|x| x as u64));
        pathogenic_min_builder.append_option(entry.pathogenic_min.map(|x| x as u64));
    }

    let schema = Schema::new(vec![
        Field::new("id", DataType::Utf8, false),
        Field::new("chrom", DataType::Utf8, false),
        Field::new("start", DataType::UInt64, false),
        Field::new("end", DataType::UInt64, false),
        Field::new("motif", DataType::Utf8, false),
        Field::new("repeat_count", DataType::UInt64, false),
        Field::new("gene", DataType::Utf8, true),
        Field::new("disease", DataType::Utf8, true),
        Field::new("normal_max", DataType::UInt64, true),
        Field::new("pathogenic_min", DataType::UInt64, true),
    ]);

    RecordBatch::try_new(
        Arc::new(schema),
        vec![
            Arc::new(id_builder.finish()) as ArrayRef,
            Arc::new(chrom_builder.finish()) as ArrayRef,
            Arc::new(start_builder.finish()) as ArrayRef,
            Arc::new(end_builder.finish()) as ArrayRef,
            Arc::new(motif_builder.finish()) as ArrayRef,
            Arc::new(repeat_count_builder.finish()) as ArrayRef,
            Arc::new(gene_builder.finish()) as ArrayRef,
            Arc::new(disease_builder.finish()) as ArrayRef,
            Arc::new(normal_max_builder.finish()) as ArrayRef,
            Arc::new(pathogenic_min_builder.finish()) as ArrayRef,
        ],
    )
}

/// Build Arrow RecordBatch for features (computed values for modeling)
pub fn build_features_batch(
    features: &[(String, String, HashMap<String, f64>)], // (unit_id, sample_id, feature_map)
) -> Result<RecordBatch, arrow::error::ArrowError> {
    // Collect all feature names
    let mut all_features: Vec<String> = features
        .iter()
        .flat_map(|(_, _, fm)| fm.keys().cloned())
        .collect::<std::collections::HashSet<_>>()
        .into_iter()
        .collect();
    all_features.sort();

    let mut unit_id_builder = StringBuilder::new();
    let mut sample_id_builder = StringBuilder::new();
    let mut feature_builders: Vec<arrow::array::Float64Builder> = all_features
        .iter()
        .map(|_| arrow::array::Float64Builder::new())
        .collect();

    for (unit_id, sample_id, fm) in features {
        unit_id_builder.append_value(unit_id);
        sample_id_builder.append_value(sample_id);
        for (i, feat_name) in all_features.iter().enumerate() {
            feature_builders[i].append_option(fm.get(feat_name).copied());
        }
    }

    let mut fields = vec![
        Field::new("unit_id", DataType::Utf8, false),
        Field::new("sample_id", DataType::Utf8, false),
    ];
    for feat_name in &all_features {
        fields.push(Field::new(feat_name, DataType::Float64, true));
    }

    let mut arrays: Vec<ArrayRef> = vec![
        Arc::new(unit_id_builder.finish()),
        Arc::new(sample_id_builder.finish()),
    ];
    for mut builder in feature_builders {
        arrays.push(Arc::new(builder.finish()));
    }

    RecordBatch::try_new(Arc::new(Schema::new(fields)), arrays)
}

/// Provenance metadata for cache
#[derive(Debug, Clone)]
pub struct Provenance {
    /// Storm version
    pub storm_version: String,
    /// Creation timestamp
    pub created_at: String,
    /// Input file paths
    pub input_files: Vec<String>,
    /// Number of samples
    pub num_samples: usize,
    /// Number of test units
    pub num_test_units: usize,
    /// Additional metadata
    pub metadata: HashMap<String, String>,
}

impl Default for Provenance {
    fn default() -> Self {
        Provenance {
            storm_version: env!("CARGO_PKG_VERSION").to_string(),
            created_at: Utc::now().to_rfc3339(),
            input_files: Vec::new(),
            num_samples: 0,
            num_test_units: 0,
            metadata: HashMap::new(),
        }
    }
}

impl Provenance {
    /// Create a new provenance record
    pub fn new() -> Self {
        Self::default()
    }

    /// Convert to Arrow RecordBatch
    pub fn to_batch(&self) -> Result<RecordBatch, arrow::error::ArrowError> {
        let schema = Schema::new(vec![
            Field::new("key", DataType::Utf8, false),
            Field::new("value", DataType::Utf8, false),
        ]);

        let mut keys = vec![
            "storm_version",
            "created_at",
            "input_files",
            "num_samples",
            "num_test_units",
        ];
        let mut values = vec![
            self.storm_version.clone(),
            self.created_at.clone(),
            self.input_files.join(","),
            self.num_samples.to_string(),
            self.num_test_units.to_string(),
        ];

        for (k, v) in &self.metadata {
            keys.push(k);
            values.push(v.clone());
        }

        let keys_array = StringArray::from(keys.iter().map(|s| s.to_string()).collect::<Vec<_>>());
        let values_array = StringArray::from(values);

        RecordBatch::try_new(
            Arc::new(schema),
            vec![
                Arc::new(keys_array) as ArrayRef,
                Arc::new(values_array) as ArrayRef,
            ],
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_build_test_units_batch() {
        let units = vec![
            TestUnit::from_sv("sv1", "chr1", 1000, 1500, "DEL"),
            TestUnit::from_true_repeat("TR001", "chr1", 10000, 10050, "CAG", "TR001"),
        ];

        let batch = build_test_units_batch(&units).expect("Failed to build batch");
        assert_eq!(batch.num_rows(), 2);
        assert_eq!(batch.num_columns(), 8);
    }

    #[test]
    fn test_build_genotypes_batch() {
        use crate::resolver::{ResolvedGenotype, PresenceSource, AlleleSource};

        let gts = vec![(
            "unit1".to_string(),
            vec![
                ResolvedGenotype {
                    sample_id: "SAMPLE1".to_string(),
                    is_present: true,
                    presence_source: PresenceSource::SvVcf,
                    allele1: Some(30),
                    allele2: Some(45),
                    allele_source: AlleleSource::SvProxy,
                    raw_gt: Some("0/1".to_string()),
                },
            ],
        )];

        let batch = build_genotypes_batch(&gts).expect("Failed to build batch");
        assert_eq!(batch.num_rows(), 1);
        assert_eq!(batch.num_columns(), 7);
    }

    #[test]
    fn test_build_catalog_batch() {
        let catalog = crate::catalog::Catalog::from_bed("fixtures/trexplorer.bed")
            .expect("Failed to load catalog");
        let entries: Vec<&CatalogEntry> = catalog.entries.values().collect();

        let batch = build_catalog_batch(&entries).expect("Failed to build batch");
        assert_eq!(batch.num_rows(), 4);
        assert_eq!(batch.num_columns(), 10);
    }

    #[test]
    fn test_build_features_batch() {
        let mut features = Vec::new();
        let mut fm = HashMap::new();
        fm.insert("S".to_string(), 75.0);
        fm.insert("M".to_string(), 45.0);
        fm.insert("D".to_string(), 15.0);
        features.push(("unit1".to_string(), "SAMPLE1".to_string(), fm));

        let batch = build_features_batch(&features).expect("Failed to build batch");
        assert_eq!(batch.num_rows(), 1);
        assert_eq!(batch.num_columns(), 5); // unit_id, sample_id, D, M, S
    }

    #[test]
    fn test_provenance() {
        let mut prov = Provenance::new();
        prov.input_files = vec!["file1.vcf".to_string(), "file2.vcf".to_string()];
        prov.num_samples = 100;
        prov.num_test_units = 500;

        let batch = prov.to_batch().expect("Failed to build provenance batch");
        assert_eq!(batch.num_rows(), 5);
        assert_eq!(batch.num_columns(), 2);
    }
}
