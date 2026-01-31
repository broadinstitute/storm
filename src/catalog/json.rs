//! TRExplorer JSON file parser
//! Supports plain .json and gzip .json.gz (standard gzip, not BGZF).

use std::fs::File;
use std::io::{BufReader, Read};
use std::path::Path;
use flate2::read::GzDecoder;
use serde::de::{self, Deserializer};
use serde::Deserialize;
use serde_json::Value;
use thiserror::Error;

/// Errors that can occur during JSON parsing
#[derive(Error, Debug)]
pub enum JsonError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
    #[error("JSON parse error: {0}")]
    Json(#[from] serde_json::Error),
}

/// A parsed JSON record from TRExplorer.
///
/// This supports:
/// - The small fixture format used in this repo (lowercase keys like `id`, `chrom`, `start`, ...).
/// - The TRExplorer "EH with annotations" export format (capitalized keys like `LocusId`,
///   `ReferenceRegion`, `Motif`, `NumRepeatsInReference`, etc.).
#[derive(Debug, Clone)]
pub struct JsonRecord {
    /// Repeat ID
    pub id: String,
    /// Chromosome
    pub chrom: String,
    /// Start position (0-based)
    pub start: u64,
    /// End position
    pub end: u64,
    /// Repeat motif
    pub motif: String,
    /// Repeat count
    pub repeat_count: u32,
    /// Strand
    pub strand: String,
    /// Repeat type
    pub repeat_type: String,
    /// Pathogenicity status
    pub pathogenicity: String,
    /// Associated gene
    pub gene: Option<String>,
    /// Associated disease
    pub disease: Option<String>,
    /// Normal maximum repeat count
    pub normal_max: Option<u32>,
    /// Pathogenic minimum repeat count
    pub pathogenic_min: Option<u32>,
}

impl<'de> Deserialize<'de> for JsonRecord {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let v = Value::deserialize(deserializer)?;
        let obj = v
            .as_object()
            .ok_or_else(|| de::Error::custom("expected JSON object"))?;

        let get_str = |k: &str| obj.get(k).and_then(|v| v.as_str()).map(|s| s.to_string());
        let get_u64 = |k: &str| obj.get(k).and_then(|v| v.as_u64());

        let id = get_str("id")
            .or_else(|| get_str("LocusId"))
            .ok_or_else(|| de::Error::missing_field("id/LocusId"))?;

        // Prefer explicit (chrom, start, end); otherwise parse ReferenceRegion "chr1:10000-10108".
        let (chrom, start, end) = if let (Some(chrom), Some(start), Some(end)) =
            (get_str("chrom"), get_u64("start"), get_u64("end"))
        {
            (chrom, start, end)
        } else if let Some(rr) = get_str("ReferenceRegion") {
            let (chrom_part, range_part) = rr
                .split_once(':')
                .ok_or_else(|| de::Error::custom(format!("invalid ReferenceRegion: {rr}")))?;
            let (start_s, end_s) = range_part
                .split_once('-')
                .ok_or_else(|| de::Error::custom(format!("invalid ReferenceRegion: {rr}")))?;
            let start: u64 = start_s
                .parse()
                .map_err(|_| de::Error::custom(format!("invalid ReferenceRegion start: {rr}")))?;
            let end: u64 = end_s
                .parse()
                .map_err(|_| de::Error::custom(format!("invalid ReferenceRegion end: {rr}")))?;
            (chrom_part.to_string(), start, end)
        } else {
            return Err(de::Error::custom(
                "missing (chrom,start,end) or ReferenceRegion",
            ));
        };

        let motif = get_str("motif")
            .or_else(|| get_str("Motif"))
            .unwrap_or_default();

        let repeat_count = get_u64("repeat_count")
            .or_else(|| get_u64("NumRepeatsInReference"))
            .unwrap_or(0)
            .min(u32::MAX as u64) as u32;

        let strand = get_str("strand")
            .or_else(|| get_str("Strand"))
            .unwrap_or_else(|| "+".to_string());

        let repeat_type = get_str("repeat_type")
            .or_else(|| get_str("RepeatType"))
            .or_else(|| get_str("VariantType"))
            .unwrap_or_else(|| "unknown".to_string());

        let pathogenicity = get_str("pathogenicity").unwrap_or_else(|| "unknown".to_string());

        let gene = get_str("gene")
            .or_else(|| get_str("GencodeGeneName"))
            .or_else(|| get_str("ManeGeneName"))
            .or_else(|| get_str("RefseqGeneName"));

        let disease = get_str("disease")
            .or_else(|| get_str("Disease"))
            .or_else(|| {
                // Some exports use longer field names (e.g. "...Disease...").
                obj.iter()
                    .find_map(|(k, v)| (k.contains("Disease")).then_some(v))
                    .and_then(|v| v.as_str())
                    .map(|s| s.to_string())
            });

        let normal_max = obj
            .get("normal_max")
            .and_then(|v| v.as_u64())
            .or_else(|| obj.get("NormalMax").and_then(|v| v.as_u64()))
            .map(|v| v.min(u32::MAX as u64) as u32);

        let pathogenic_min = obj
            .get("pathogenic_min")
            .and_then(|v| v.as_u64())
            .or_else(|| obj.get("PathogenicMin").and_then(|v| v.as_u64()))
            .map(|v| v.min(u32::MAX as u64) as u32);

        Ok(JsonRecord {
            id,
            chrom,
            start,
            end,
            motif,
            repeat_count,
            strand,
            repeat_type,
            pathogenicity,
            gene,
            disease,
            normal_max,
            pathogenic_min,
        })
    }
}

fn open_json_reader<P: AsRef<Path>>(path: P) -> Result<Box<dyn Read>, JsonError> {
    let path = path.as_ref();
    let file = File::open(path)?;
    let reader: Box<dyn Read> = if path.extension().map(|e| e == "gz").unwrap_or(false) {
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };
    Ok(reader)
}

/// Parse a TRExplorer JSON file (plain .json or .json.gz)
pub fn parse_trexplorer_json<P: AsRef<Path>>(path: P) -> Result<Vec<JsonRecord>, JsonError> {
    let reader = open_json_reader(path)?;
    let records: Vec<JsonRecord> = serde_json::from_reader(reader)?;
    Ok(records)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_trexplorer_json() {
        let records = parse_trexplorer_json("fixtures/trexplorer.json").expect("Failed to parse JSON");
        
        assert_eq!(records.len(), 4);

        let rec = &records[0];
        assert_eq!(rec.id, "TR001");
        assert_eq!(rec.chrom, "chr1");
        assert_eq!(rec.start, 10000);
        assert_eq!(rec.end, 10050);
        assert_eq!(rec.gene, Some("HTT".to_string()));
        assert_eq!(rec.disease, None);
        assert_eq!(rec.normal_max, Some(26));
        assert_eq!(rec.pathogenic_min, Some(36));

        let rec = &records[2];
        assert_eq!(rec.id, "TR003");
        assert_eq!(rec.gene, Some("FMR1".to_string()));
        assert_eq!(rec.disease, Some("Fragile X".to_string()));
    }
}
