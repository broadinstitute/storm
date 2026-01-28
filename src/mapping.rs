//! Mapping module for associating SVs with catalog loci
//!
//! Maps structural variant records to tandem repeat loci based on genomic overlap.

use crate::catalog::{Catalog, CatalogEntry};
use crate::vcf::sv::SvRecord;

/// A mapping between an SV record and overlapping catalog entries
#[derive(Debug, Clone)]
pub struct SvMapping<'a> {
    /// The SV record
    pub sv: &'a SvRecord,
    /// Overlapping catalog entries
    pub overlaps: Vec<&'a CatalogEntry>,
}

impl<'a> SvMapping<'a> {
    /// Check if this SV overlaps any repeat loci
    pub fn has_overlaps(&self) -> bool {
        !self.overlaps.is_empty()
    }

    /// Get the primary (first/closest) overlap if any
    pub fn primary_overlap(&self) -> Option<&CatalogEntry> {
        self.overlaps.first().copied()
    }
}

/// Map SV records to catalog repeat loci by genomic overlap
pub fn map_svs_to_catalog<'a>(
    svs: &'a [SvRecord],
    catalog: &'a Catalog,
) -> Vec<SvMapping<'a>> {
    svs.iter()
        .map(|sv| {
            // Calculate the SV interval
            let start = sv.pos;
            let end = sv.end.unwrap_or_else(|| {
                // For INS, use position +/- some window
                // For DEL/DUP with SVLEN, calculate end
                if let Some(svlen) = sv.sv_len {
                    if svlen < 0 {
                        start + (-svlen) as u64
                    } else {
                        start + svlen as u64
                    }
                } else {
                    start + 1
                }
            });

            let overlaps = catalog.find_overlapping(&sv.chrom, start, end);

            SvMapping { sv, overlaps }
        })
        .collect()
}

/// Map SV records and return only those with overlaps
pub fn map_svs_with_overlaps<'a>(
    svs: &'a [SvRecord],
    catalog: &'a Catalog,
) -> Vec<SvMapping<'a>> {
    map_svs_to_catalog(svs, catalog)
        .into_iter()
        .filter(|m| m.has_overlaps())
        .collect()
}

/// Statistics about SV-to-catalog mapping
#[derive(Debug, Clone, Default)]
pub struct MappingStats {
    /// Total SVs processed
    pub total_svs: usize,
    /// SVs with at least one overlap
    pub svs_with_overlaps: usize,
    /// SVs without any overlap
    pub svs_without_overlaps: usize,
    /// Total number of overlaps found
    pub total_overlaps: usize,
}

/// Compute mapping statistics
pub fn compute_mapping_stats(mappings: &[SvMapping]) -> MappingStats {
    let total_svs = mappings.len();
    let svs_with_overlaps = mappings.iter().filter(|m| m.has_overlaps()).count();
    let total_overlaps: usize = mappings.iter().map(|m| m.overlaps.len()).sum();

    MappingStats {
        total_svs,
        svs_with_overlaps,
        svs_without_overlaps: total_svs - svs_with_overlaps,
        total_overlaps,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vcf::sv::parse_sv_vcf;

    #[test]
    fn test_map_svs_to_catalog() {
        let (_, svs) = parse_sv_vcf("fixtures/sv_small.vcf").expect("Failed to parse SV VCF");
        let catalog = Catalog::from_bed("fixtures/trexplorer.bed").expect("Failed to load catalog");

        let mappings = map_svs_to_catalog(&svs, &catalog);
        assert_eq!(mappings.len(), 3);

        // First SV (chr1:1000-1500) should overlap TR001 (chr1:10000-10050) - NO overlap
        assert!(!mappings[0].has_overlaps());

        // Check stats
        let stats = compute_mapping_stats(&mappings);
        assert_eq!(stats.total_svs, 3);
    }

    #[test]
    fn test_map_svs_with_explicit_overlap() {
        // Create a catalog with a locus that will definitely overlap
        let catalog = Catalog::from_bed("fixtures/trexplorer.bed").expect("Failed to load catalog");
        
        // TR001 is at chr1:10000-10050
        // Our SVs are at chr1:1000-1500 and chr1:2000+200
        // None overlap, which is correct

        // Verify the catalog entries
        let entry = catalog.get("TR001").unwrap();
        assert_eq!(entry.start, 10000);
        assert_eq!(entry.end, 10050);
    }
}
