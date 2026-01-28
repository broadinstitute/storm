//! Parquet file I/O for cache

use std::fs::File;
use std::path::Path;
use arrow::record_batch::RecordBatch;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::arrow::ArrowWriter;
use parquet::file::properties::WriterProperties;
use thiserror::Error;

/// Errors that can occur with Parquet I/O
#[derive(Error, Debug)]
pub enum ParquetError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
    #[error("Parquet error: {0}")]
    Parquet(#[from] parquet::errors::ParquetError),
    #[error("Arrow error: {0}")]
    Arrow(#[from] arrow::error::ArrowError),
}

/// Write a RecordBatch to a Parquet file
pub fn write_parquet<P: AsRef<Path>>(
    path: P,
    batch: &RecordBatch,
) -> Result<(), ParquetError> {
    let file = File::create(path)?;
    
    let props = WriterProperties::builder()
        .set_compression(parquet::basic::Compression::SNAPPY)
        .build();

    let mut writer = ArrowWriter::try_new(file, batch.schema(), Some(props))?;
    writer.write(batch)?;
    writer.close()?;
    
    Ok(())
}

/// Read a RecordBatch from a Parquet file
pub fn read_parquet<P: AsRef<Path>>(path: P) -> Result<RecordBatch, ParquetError> {
    let file = File::open(path.as_ref())?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(file)?;
    let schema = builder.schema().clone();
    let mut reader = builder.build()?;
    
    // Read all batches and concatenate (for simplicity, assume one batch)
    if let Some(batch) = reader.next() {
        Ok(batch?)
    } else {
        // Return an empty batch with the correct schema
        Ok(RecordBatch::new_empty(schema))
    }
}

/// Write full cache to a directory
pub fn write_cache_to_dir<P: AsRef<Path>>(
    dir: P,
    cache: &super::ArrowCache,
) -> Result<(), ParquetError> {
    let dir = dir.as_ref();
    std::fs::create_dir_all(dir)?;

    if let Some(batch) = &cache.test_units {
        write_parquet(dir.join("test_units.parquet"), batch)?;
    }
    if let Some(batch) = &cache.genotypes {
        write_parquet(dir.join("genotypes.parquet"), batch)?;
    }
    if let Some(batch) = &cache.catalog {
        write_parquet(dir.join("catalog.parquet"), batch)?;
    }
    if let Some(batch) = &cache.features {
        write_parquet(dir.join("features.parquet"), batch)?;
    }
    if let Some(batch) = &cache.provenance {
        write_parquet(dir.join("provenance.parquet"), batch)?;
    }

    Ok(())
}

/// Read cache from a directory
pub fn read_cache_from_dir<P: AsRef<Path>>(dir: P) -> Result<super::ArrowCache, ParquetError> {
    let dir = dir.as_ref();
    let mut cache = super::ArrowCache::new();

    let tu_path = dir.join("test_units.parquet");
    if tu_path.exists() {
        cache.test_units = Some(read_parquet(tu_path)?);
    }

    let gt_path = dir.join("genotypes.parquet");
    if gt_path.exists() {
        cache.genotypes = Some(read_parquet(gt_path)?);
    }

    let cat_path = dir.join("catalog.parquet");
    if cat_path.exists() {
        cache.catalog = Some(read_parquet(cat_path)?);
    }

    let feat_path = dir.join("features.parquet");
    if feat_path.exists() {
        cache.features = Some(read_parquet(feat_path)?);
    }

    let prov_path = dir.join("provenance.parquet");
    if prov_path.exists() {
        cache.provenance = Some(read_parquet(prov_path)?);
    }

    Ok(cache)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testunit::TestUnit;
    use crate::cache::build_test_units_batch;
    use std::collections::HashMap;

    #[test]
    fn test_write_read_parquet() {
        let units = vec![
            TestUnit::from_sv("sv1", "chr1", 1000, 1500, "DEL"),
            TestUnit::from_sv("sv2", "chr2", 2000, 2500, "INS"),
        ];

        let batch = build_test_units_batch(&units).expect("Failed to build batch");
        
        let temp_dir = std::env::temp_dir();
        let path = temp_dir.join("storm_test_units.parquet");
        
        // Write
        write_parquet(&path, &batch).expect("Failed to write parquet");
        
        // Read back
        let read_batch = read_parquet(&path).expect("Failed to read parquet");
        
        assert_eq!(read_batch.num_rows(), 2);
        assert_eq!(read_batch.num_columns(), 8);

        // Cleanup
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn test_write_cache_to_dir() {
        use crate::cache::{ArrowCache, Provenance};
        use crate::resolver::{ResolvedGenotype, PresenceSource, AlleleSource};

        let mut cache = ArrowCache::new();

        // Populate cache
        let units = vec![TestUnit::from_sv("sv1", "chr1", 1000, 1500, "DEL")];
        cache.set_test_units(&units).unwrap();

        let gts = vec![(
            "sv1".to_string(),
            vec![ResolvedGenotype {
                sample_id: "S1".to_string(),
                is_present: true,
                presence_source: PresenceSource::SvVcf,
                allele1: Some(30),
                allele2: Some(30),
                allele_source: AlleleSource::Reference,
                raw_gt: None,
            }],
        )];
        cache.set_genotypes(&gts).unwrap();

        let catalog = crate::catalog::Catalog::from_bed("fixtures/trexplorer.bed").unwrap();
        let entries: Vec<&crate::catalog::CatalogEntry> = catalog.entries.values().collect();
        cache.set_catalog(&entries).unwrap();

        let mut fm = HashMap::new();
        fm.insert("S".to_string(), 60.0);
        let features: Vec<(String, String, HashMap<String, f64>)> = vec![
            ("sv1".to_string(), "S1".to_string(), fm),
        ];
        cache.set_features(&features).unwrap();

        let prov = Provenance::new();
        cache.set_provenance(&prov).unwrap();

        // Write to directory
        let temp_dir = std::env::temp_dir().join("storm_cache_test");
        write_cache_to_dir(&temp_dir, &cache).expect("Failed to write cache");

        // Read back
        let read_cache = read_cache_from_dir(&temp_dir).expect("Failed to read cache");
        assert!(read_cache.test_units.is_some());
        assert!(read_cache.genotypes.is_some());
        assert!(read_cache.catalog.is_some());
        assert!(read_cache.provenance.is_some());

        // Cleanup
        std::fs::remove_dir_all(&temp_dir).ok();
    }
}
