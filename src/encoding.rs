//! Encoding module for genotype transformations
//!
//! Implements various encodings for repeat lengths:
//! - S: Sum of allele lengths (L1 + L2)
//! - M: Maximum allele length
//! - D: Absolute difference |L1 - L2|
//! - Tail: Indicator for exceeding a threshold
//! - Categorical: Binned categories

use crate::resolver::ResolvedGenotype;

/// Encode genotype as sum of allele lengths (S encoding)
/// Returns L1 + L2
pub fn encode_s(gt: &ResolvedGenotype) -> Option<f64> {
    match (gt.allele1, gt.allele2) {
        (Some(a1), Some(a2)) => Some((a1 + a2) as f64),
        _ => None,
    }
}

/// Encode genotype as maximum allele length (M encoding)
/// Returns max(L1, L2)
pub fn encode_m(gt: &ResolvedGenotype) -> Option<f64> {
    match (gt.allele1, gt.allele2) {
        (Some(a1), Some(a2)) => Some(a1.max(a2) as f64),
        (Some(a), None) | (None, Some(a)) => Some(a as f64),
        _ => None,
    }
}

/// Encode genotype as absolute difference (D encoding)
/// Returns |L1 - L2|
pub fn encode_d(gt: &ResolvedGenotype) -> Option<f64> {
    match (gt.allele1, gt.allele2) {
        (Some(a1), Some(a2)) => Some((a1 - a2).abs() as f64),
        _ => None,
    }
}

/// Encode genotype as tail indicator (Tail encoding)
/// Returns 1.0 if max(L1, L2) >= threshold, 0.0 otherwise
pub fn encode_tail(gt: &ResolvedGenotype, threshold: i64) -> Option<f64> {
    match (gt.allele1, gt.allele2) {
        (Some(a1), Some(a2)) => {
            let max_allele = a1.max(a2);
            Some(if max_allele >= threshold { 1.0 } else { 0.0 })
        }
        _ => None,
    }
}

/// Categorical bin definition
#[derive(Debug, Clone)]
pub struct CategoricalBin {
    /// Bin label
    pub label: String,
    /// Lower bound (inclusive)
    pub lower: i64,
    /// Upper bound (exclusive, None means infinity)
    pub upper: Option<i64>,
}

impl CategoricalBin {
    /// Create a new bin
    pub fn new(label: &str, lower: i64, upper: Option<i64>) -> Self {
        CategoricalBin {
            label: label.to_string(),
            lower,
            upper,
        }
    }

    /// Check if value falls in this bin
    pub fn contains(&self, value: i64) -> bool {
        if value < self.lower {
            return false;
        }
        if let Some(upper) = self.upper {
            value < upper
        } else {
            true
        }
    }
}

/// Encode genotype into categorical bins
/// Returns the bin index (0-based) for the maximum allele
pub fn encode_categorical(gt: &ResolvedGenotype, bins: &[CategoricalBin]) -> Option<usize> {
    let max_allele = match (gt.allele1, gt.allele2) {
        (Some(a1), Some(a2)) => a1.max(a2),
        (Some(a), None) | (None, Some(a)) => a,
        _ => return None,
    };

    for (i, bin) in bins.iter().enumerate() {
        if bin.contains(max_allele) {
            return Some(i);
        }
    }

    None
}

/// Encode genotype as binary presence (1.0 if present, 0.0 if absent)
pub fn encode_binary(gt: &ResolvedGenotype) -> f64 {
    if gt.is_present { 1.0 } else { 0.0 }
}

/// Apply encoding to a list of genotypes
pub fn apply_encoding(
    genotypes: &[ResolvedGenotype],
    encoding: &crate::plan::Encoding,
    threshold: Option<i64>,
    bins: Option<&[CategoricalBin]>,
) -> Vec<Option<f64>> {
    genotypes
        .iter()
        .map(|gt| match encoding {
            crate::plan::Encoding::S => encode_s(gt),
            crate::plan::Encoding::M => encode_m(gt),
            crate::plan::Encoding::D => encode_d(gt),
            crate::plan::Encoding::Tail => encode_tail(gt, threshold.unwrap_or(0)),
            crate::plan::Encoding::Categorical => {
                bins.and_then(|b| encode_categorical(gt, b).map(|i| i as f64))
            }
            crate::plan::Encoding::Binary => Some(encode_binary(gt)),
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::resolver::{PresenceSource, AlleleSource};

    fn make_gt(a1: i64, a2: i64, present: bool) -> ResolvedGenotype {
        ResolvedGenotype {
            sample_id: "test".to_string(),
            is_present: present,
            presence_source: PresenceSource::TrgtVcf,
            allele1: Some(a1),
            allele2: Some(a2),
            allele_source: AlleleSource::TrgtTrue,
            raw_gt: None,
        }
    }

    #[test]
    fn test_encode_s() {
        let gt = make_gt(30, 45, true);
        assert_eq!(encode_s(&gt), Some(75.0));

        let gt2 = make_gt(10, 10, false);
        assert_eq!(encode_s(&gt2), Some(20.0));
    }

    #[test]
    fn test_encode_m() {
        let gt = make_gt(30, 45, true);
        assert_eq!(encode_m(&gt), Some(45.0));

        let gt2 = make_gt(100, 50, true);
        assert_eq!(encode_m(&gt2), Some(100.0));
    }

    #[test]
    fn test_encode_d() {
        let gt = make_gt(30, 45, true);
        assert_eq!(encode_d(&gt), Some(15.0));

        let gt2 = make_gt(45, 30, true);
        assert_eq!(encode_d(&gt2), Some(15.0)); // Absolute value
    }

    #[test]
    fn test_encode_tail() {
        let gt = make_gt(30, 45, true);
        
        assert_eq!(encode_tail(&gt, 40), Some(1.0)); // 45 >= 40
        assert_eq!(encode_tail(&gt, 50), Some(0.0)); // 45 < 50
        assert_eq!(encode_tail(&gt, 45), Some(1.0)); // 45 >= 45 (inclusive)
    }

    #[test]
    fn test_encode_categorical() {
        let bins = vec![
            CategoricalBin::new("normal", 0, Some(36)),
            CategoricalBin::new("intermediate", 36, Some(40)),
            CategoricalBin::new("pathogenic", 40, None),
        ];

        let gt_normal = make_gt(20, 25, false);
        assert_eq!(encode_categorical(&gt_normal, &bins), Some(0));

        let gt_intermediate = make_gt(35, 38, true);
        assert_eq!(encode_categorical(&gt_intermediate, &bins), Some(1));

        let gt_pathogenic = make_gt(35, 45, true);
        assert_eq!(encode_categorical(&gt_pathogenic, &bins), Some(2));
    }

    #[test]
    fn test_encode_binary() {
        let gt_present = make_gt(30, 45, true);
        assert_eq!(encode_binary(&gt_present), 1.0);

        let gt_absent = make_gt(30, 30, false);
        assert_eq!(encode_binary(&gt_absent), 0.0);
    }

    #[test]
    fn test_apply_encoding() {
        let gts = vec![
            make_gt(30, 45, true),
            make_gt(20, 20, false),
            make_gt(50, 60, true),
        ];

        let s_encoded = apply_encoding(&gts, &crate::plan::Encoding::S, None, None);
        assert_eq!(s_encoded, vec![Some(75.0), Some(40.0), Some(110.0)]);

        let m_encoded = apply_encoding(&gts, &crate::plan::Encoding::M, None, None);
        assert_eq!(m_encoded, vec![Some(45.0), Some(20.0), Some(60.0)]);

        let binary_encoded = apply_encoding(&gts, &crate::plan::Encoding::Binary, None, None);
        assert_eq!(binary_encoded, vec![Some(1.0), Some(0.0), Some(1.0)]);
    }
}
