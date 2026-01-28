//! Genotype Resolver module
//!
//! Handles resolution of genotypes from multiple sources:
//! - SV presence (carrier status)
//! - SV-derived allele lengths (from INS/DEL SVLEN)
//! - TRGT true allele lengths
//!
//! The resolver separates "presence" (is the variant present) from 
//! "allele values" (what are the measured lengths).

use std::collections::HashMap;
use crate::vcf::sv::{SvRecord, SvType, Genotype};
use crate::vcf::trgt::TrgtRecord;
use crate::catalog::CatalogEntry;

/// Source of presence information
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PresenceSource {
    /// Presence from SV VCF genotype
    SvVcf,
    /// Presence from TRGT VCF (non-reference allele)
    TrgtVcf,
    /// Unknown/missing
    Unknown,
}

/// Source of allele length information
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AlleleSource {
    /// Allele length from SV SVLEN (proxy)
    SvProxy,
    /// Allele length from TRGT AL field (true)
    TrgtTrue,
    /// Reference length from catalog
    Reference,
    /// Unknown/missing
    Unknown,
}

/// A resolved genotype for a sample at a locus
#[derive(Debug, Clone)]
pub struct ResolvedGenotype {
    /// Sample ID
    pub sample_id: String,
    /// Is the variant present (any non-ref allele)
    pub is_present: bool,
    /// Source of presence call
    pub presence_source: PresenceSource,
    /// First allele length (bp)
    pub allele1: Option<i64>,
    /// Second allele length (bp)
    pub allele2: Option<i64>,
    /// Source of allele length information
    pub allele_source: AlleleSource,
    /// Raw genotype string
    pub raw_gt: Option<String>,
}

impl ResolvedGenotype {
    /// Create a missing genotype
    pub fn missing(sample_id: &str) -> Self {
        ResolvedGenotype {
            sample_id: sample_id.to_string(),
            is_present: false,
            presence_source: PresenceSource::Unknown,
            allele1: None,
            allele2: None,
            allele_source: AlleleSource::Unknown,
            raw_gt: None,
        }
    }

    /// Sum of allele lengths (L1 + L2) for S encoding
    pub fn sum_lengths(&self) -> Option<i64> {
        match (self.allele1, self.allele2) {
            (Some(a1), Some(a2)) => Some(a1 + a2),
            _ => None,
        }
    }

    /// Maximum allele length for M encoding
    pub fn max_length(&self) -> Option<i64> {
        match (self.allele1, self.allele2) {
            (Some(a1), Some(a2)) => Some(a1.max(a2)),
            (Some(a), None) | (None, Some(a)) => Some(a),
            _ => None,
        }
    }

    /// Absolute difference of allele lengths for D encoding
    pub fn diff_lengths(&self) -> Option<i64> {
        match (self.allele1, self.allele2) {
            (Some(a1), Some(a2)) => Some((a1 - a2).abs()),
            _ => None,
        }
    }
}

/// SV alleles grouped by locus for proxy resolution
#[derive(Debug, Clone, Default)]
pub struct LocusSvAlleles {
    /// List of (SVLEN, is_insertion) from overlapping SVs
    pub alleles: Vec<(i64, bool)>,
}

impl LocusSvAlleles {
    /// Add an SV's allele contribution
    pub fn add_sv(&mut self, sv: &SvRecord, gt: &Genotype) {
        // Only count if sample has non-ref allele
        let has_alt = gt.alleles.iter().any(|a| matches!(a, Some(x) if *x > 0));
        if !has_alt {
            return;
        }

        if let Some(svlen) = sv.sv_len {
            let is_insertion = matches!(sv.sv_type, SvType::Ins);
            self.alleles.push((svlen, is_insertion));
        }
    }

    /// Reconstruct diploid lengths from INS/DEL alleles
    /// Returns (shorter_allele, longer_allele) relative to reference
    pub fn reconstruct_diploid(&self, ref_length: i64) -> Option<(i64, i64)> {
        if self.alleles.is_empty() {
            return Some((ref_length, ref_length));
        }

        // Sum all contributions
        let mut total_delta: i64 = 0;
        for (svlen, is_ins) in &self.alleles {
            if *is_ins {
                total_delta += svlen.abs();
            } else {
                total_delta -= svlen.abs();
            }
        }

        // For simplicity, assume one allele changed and one is reference
        // More sophisticated logic would use genotype phasing
        let alt_length = ref_length + total_delta;
        
        if alt_length >= ref_length {
            Some((ref_length, alt_length))
        } else {
            Some((alt_length, ref_length))
        }
    }
}

/// The genotype resolver
pub struct Resolver {
    /// Reference lengths from catalog (locus_id -> ref_length)
    ref_lengths: HashMap<String, i64>,
}

impl Resolver {
    /// Create a new resolver with catalog reference lengths
    pub fn new() -> Self {
        Resolver {
            ref_lengths: HashMap::new(),
        }
    }

    /// Add reference lengths from catalog entries
    pub fn add_catalog_refs(&mut self, entries: &[&CatalogEntry]) {
        for entry in entries {
            // Reference length = motif_len * repeat_count
            let ref_len = entry.motif.len() as i64 * entry.repeat_count as i64;
            self.ref_lengths.insert(entry.id.clone(), ref_len);
        }
    }

    /// Set a reference length directly
    pub fn set_ref_length(&mut self, locus_id: &str, length: i64) {
        self.ref_lengths.insert(locus_id.to_string(), length);
    }

    /// Resolve genotypes for a locus from SV data only (proxy)
    pub fn resolve_from_svs(
        &self,
        locus_id: &str,
        svs: &[&SvRecord],
        sample_names: &[String],
    ) -> Vec<ResolvedGenotype> {
        let ref_len = self.ref_lengths.get(locus_id).copied().unwrap_or(0);

        sample_names
            .iter()
            .map(|sample| {
                let mut locus_alleles = LocusSvAlleles::default();
                let mut has_any_alt = false;
                let mut raw_gts: Vec<String> = Vec::new();

                for sv in svs {
                    if let Some(gt) = sv.genotypes.get(sample) {
                        let has_alt = gt.alleles.iter().any(|a| matches!(a, Some(x) if *x > 0));
                        if has_alt {
                            has_any_alt = true;
                            locus_alleles.add_sv(sv, gt);
                        }
                        let gt_str = if gt.phased {
                            gt.alleles.iter().map(|a| match a {
                                Some(x) => x.to_string(),
                                None => ".".to_string(),
                            }).collect::<Vec<_>>().join("|")
                        } else {
                            gt.alleles.iter().map(|a| match a {
                                Some(x) => x.to_string(),
                                None => ".".to_string(),
                            }).collect::<Vec<_>>().join("/")
                        };
                        raw_gts.push(gt_str);
                    }
                }

                let (allele1, allele2) = locus_alleles
                    .reconstruct_diploid(ref_len)
                    .unwrap_or((ref_len, ref_len));

                ResolvedGenotype {
                    sample_id: sample.clone(),
                    is_present: has_any_alt,
                    presence_source: if has_any_alt {
                        PresenceSource::SvVcf
                    } else {
                        PresenceSource::Unknown
                    },
                    allele1: Some(allele1),
                    allele2: Some(allele2),
                    allele_source: if has_any_alt {
                        AlleleSource::SvProxy
                    } else {
                        AlleleSource::Reference
                    },
                    raw_gt: if raw_gts.is_empty() {
                        None
                    } else {
                        Some(raw_gts.join(";"))
                    },
                }
            })
            .collect()
    }

    /// Resolve genotypes from TRGT data (true repeats)
    pub fn resolve_from_trgt(
        &self,
        trgt_record: &TrgtRecord,
        sample_names: &[String],
    ) -> Vec<ResolvedGenotype> {
        sample_names
            .iter()
            .map(|sample| {
                if let Some(trgt_gt) = trgt_record.genotypes.get(sample) {
                    let has_alt = trgt_gt.gt.iter().any(|a| matches!(a, Some(x) if *x > 0));
                    let (allele1, allele2) = trgt_gt.diploid_lengths().unwrap_or((0, 0));

                    let gt_str = if trgt_gt.phased {
                        trgt_gt.gt.iter().map(|a| match a {
                            Some(x) => x.to_string(),
                            None => ".".to_string(),
                        }).collect::<Vec<_>>().join("|")
                    } else {
                        trgt_gt.gt.iter().map(|a| match a {
                            Some(x) => x.to_string(),
                            None => ".".to_string(),
                        }).collect::<Vec<_>>().join("/")
                    };

                    ResolvedGenotype {
                        sample_id: sample.clone(),
                        is_present: has_alt,
                        presence_source: PresenceSource::TrgtVcf,
                        allele1: Some(allele1),
                        allele2: Some(allele2),
                        allele_source: AlleleSource::TrgtTrue,
                        raw_gt: Some(gt_str),
                    }
                } else {
                    ResolvedGenotype::missing(sample)
                }
            })
            .collect()
    }

    /// Resolve genotypes with TRGT overriding proxy alleles
    /// Presence comes from SV, allele lengths from TRGT when available
    pub fn resolve_merged(
        &self,
        locus_id: &str,
        svs: &[&SvRecord],
        trgt_record: Option<&TrgtRecord>,
        sample_names: &[String],
    ) -> Vec<ResolvedGenotype> {
        let sv_resolved = self.resolve_from_svs(locus_id, svs, sample_names);

        if let Some(trgt) = trgt_record {
            let trgt_resolved = self.resolve_from_trgt(trgt, sample_names);

            // Create lookup for TRGT genotypes
            let trgt_map: HashMap<&str, &ResolvedGenotype> = trgt_resolved
                .iter()
                .map(|g| (g.sample_id.as_str(), g))
                .collect();

            // Merge: use SV presence, but TRGT allele lengths when available
            sv_resolved
                .into_iter()
                .map(|mut sv_gt| {
                    if let Some(trgt_gt) = trgt_map.get(sv_gt.sample_id.as_str()) {
                        // Override allele lengths with TRGT data
                        if trgt_gt.allele1.is_some() && trgt_gt.allele2.is_some() {
                            sv_gt.allele1 = trgt_gt.allele1;
                            sv_gt.allele2 = trgt_gt.allele2;
                            sv_gt.allele_source = AlleleSource::TrgtTrue;
                        }
                    }
                    sv_gt
                })
                .collect()
        } else {
            sv_resolved
        }
    }
}

impl Default for Resolver {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vcf::sv::parse_sv_vcf;
    use crate::vcf::trgt::parse_trgt_vcf;

    #[test]
    fn test_resolved_genotype_encodings() {
        let gt = ResolvedGenotype {
            sample_id: "SAMPLE1".to_string(),
            is_present: true,
            presence_source: PresenceSource::SvVcf,
            allele1: Some(30),
            allele2: Some(45),
            allele_source: AlleleSource::SvProxy,
            raw_gt: Some("0/1".to_string()),
        };

        assert_eq!(gt.sum_lengths(), Some(75));  // S encoding
        assert_eq!(gt.max_length(), Some(45));   // M encoding
        assert_eq!(gt.diff_lengths(), Some(15)); // D encoding
    }

    #[test]
    fn test_locus_sv_alleles_reconstruct() {
        let mut alleles = LocusSvAlleles::default();
        
        // Simulate a 50bp insertion
        alleles.alleles.push((50, true));
        
        let (a1, a2) = alleles.reconstruct_diploid(100).unwrap();
        assert_eq!((a1, a2), (100, 150));

        // Test deletion
        let mut alleles2 = LocusSvAlleles::default();
        alleles2.alleles.push((-30, false));
        
        let (a1, a2) = alleles2.reconstruct_diploid(100).unwrap();
        assert_eq!((a1, a2), (70, 100));
    }

    #[test]
    fn test_resolve_from_trgt() {
        let (samples, records) = parse_trgt_vcf("fixtures/trgt_small.vcf").unwrap();
        
        let resolver = Resolver::new();
        let resolved = resolver.resolve_from_trgt(&records[0], &samples);
        
        assert_eq!(resolved.len(), 2);
        
        let gt1 = &resolved[0];
        assert_eq!(gt1.sample_id, "SAMPLE1");
        assert!(gt1.is_present);
        assert_eq!(gt1.presence_source, PresenceSource::TrgtVcf);
        assert_eq!(gt1.allele1, Some(9));
        assert_eq!(gt1.allele2, Some(12));
        assert_eq!(gt1.allele_source, AlleleSource::TrgtTrue);
    }

    #[test]
    fn test_presence_and_allele_sources_tracked() {
        let (samples, records) = parse_trgt_vcf("fixtures/trgt_small.vcf").unwrap();
        
        let resolver = Resolver::new();
        let resolved = resolver.resolve_from_trgt(&records[0], &samples);
        
        // Check that sources are properly recorded
        for gt in &resolved {
            assert_eq!(gt.presence_source, PresenceSource::TrgtVcf);
            assert_eq!(gt.allele_source, AlleleSource::TrgtTrue);
        }
    }

    #[test]
    fn test_merged_resolution_trgt_overrides_proxy() {
        let (sv_samples, _svs) = parse_sv_vcf("fixtures/sv_small.vcf").unwrap();
        let (_, trgt_records) = parse_trgt_vcf("fixtures/trgt_small.vcf").unwrap();
        
        let mut resolver = Resolver::new();
        resolver.set_ref_length("test_locus", 50);
        
        // Empty SV list, but with TRGT data
        let sv_refs: Vec<&SvRecord> = Vec::new();
        let resolved = resolver.resolve_merged(
            "test_locus",
            &sv_refs,
            Some(&trgt_records[0]),
            &sv_samples,
        );
        
        // Should have TRGT allele lengths
        let gt1 = &resolved[0];
        assert_eq!(gt1.allele_source, AlleleSource::TrgtTrue);
        assert_eq!(gt1.allele1, Some(9));
        assert_eq!(gt1.allele2, Some(12));
    }
}
