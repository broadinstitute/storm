//! Explain module for printing resolved genotype details
//!
//! Provides detailed human-readable output about genotype resolution.

use std::io::Write;
use crate::resolver::ResolvedGenotype;
use crate::testunit::TestUnit;

/// Print detailed explanation of a resolved genotype
pub fn explain_genotype<W: Write>(
    writer: &mut W,
    unit: &TestUnit,
    genotype: &ResolvedGenotype,
) -> std::io::Result<()> {
    writeln!(writer, "=== Genotype Explanation ===")?;
    writeln!(writer)?;
    
    // Test unit info
    writeln!(writer, "Test Unit: {}", unit.id)?;
    writeln!(writer, "  Location: {}:{}-{}", unit.chrom, unit.start, unit.end)?;
    writeln!(writer, "  Type: {:?}", unit.unit_type)?;
    writeln!(writer, "  Data Source: {:?}", unit.source)?;
    if let Some(ref catalog_id) = unit.catalog_id {
        writeln!(writer, "  Catalog ID: {}", catalog_id)?;
    }
    if let Some(ref motif) = unit.motif {
        writeln!(writer, "  Motif: {}", motif)?;
    }
    writeln!(writer)?;

    // Sample genotype info
    writeln!(writer, "Sample: {}", genotype.sample_id)?;
    writeln!(writer, "  Presence: {}", if genotype.is_present { "PRESENT" } else { "ABSENT" })?;
    writeln!(writer, "  Presence Source: {:?}", genotype.presence_source)?;
    writeln!(writer)?;

    // Allele info
    writeln!(writer, "Alleles:")?;
    match (genotype.allele1, genotype.allele2) {
        (Some(a1), Some(a2)) => {
            writeln!(writer, "  Allele 1: {} bp", a1)?;
            writeln!(writer, "  Allele 2: {} bp", a2)?;
        }
        (Some(a1), None) => {
            writeln!(writer, "  Allele 1: {} bp", a1)?;
            writeln!(writer, "  Allele 2: MISSING")?;
        }
        (None, Some(a2)) => {
            writeln!(writer, "  Allele 1: MISSING")?;
            writeln!(writer, "  Allele 2: {} bp", a2)?;
        }
        (None, None) => {
            writeln!(writer, "  Both alleles: MISSING")?;
        }
    }
    writeln!(writer, "  Allele Source: {:?}", genotype.allele_source)?;
    writeln!(writer)?;

    // Computed encodings
    writeln!(writer, "Computed Encodings:")?;
    if let Some(sum) = genotype.sum_lengths() {
        writeln!(writer, "  S (sum): {} bp", sum)?;
    }
    if let Some(max) = genotype.max_length() {
        writeln!(writer, "  M (max): {} bp", max)?;
    }
    if let Some(diff) = genotype.diff_lengths() {
        writeln!(writer, "  D (diff): {} bp", diff)?;
    }
    writeln!(writer)?;

    // Raw genotype if available
    if let Some(ref raw) = genotype.raw_gt {
        writeln!(writer, "Raw Genotype: {}", raw)?;
    }

    Ok(())
}

/// Print explanation for multiple genotypes at a locus
pub fn explain_locus<W: Write>(
    writer: &mut W,
    unit: &TestUnit,
    genotypes: &[ResolvedGenotype],
) -> std::io::Result<()> {
    writeln!(writer, "=== Locus Summary ===")?;
    writeln!(writer)?;
    writeln!(writer, "Test Unit: {}", unit.id)?;
    writeln!(writer, "Location: {}:{}-{}", unit.chrom, unit.start, unit.end)?;
    writeln!(writer, "Type: {:?}", unit.unit_type)?;
    if let Some(ref motif) = unit.motif {
        writeln!(writer, "Motif: {}", motif)?;
    }
    writeln!(writer)?;

    // Summary statistics
    let present_count = genotypes.iter().filter(|g| g.is_present).count();
    let total_count = genotypes.len();
    let call_rate = if total_count > 0 {
        (genotypes.iter().filter(|g| g.allele1.is_some() && g.allele2.is_some()).count() as f64)
            / (total_count as f64)
    } else {
        0.0
    };

    writeln!(writer, "Summary Statistics:")?;
    writeln!(writer, "  Total samples: {}", total_count)?;
    writeln!(writer, "  Present (carriers): {}", present_count)?;
    writeln!(writer, "  Absent (non-carriers): {}", total_count - present_count)?;
    writeln!(writer, "  Call rate: {:.2}%", call_rate * 100.0)?;
    writeln!(writer)?;

    // Per-sample table
    writeln!(writer, "Sample Details:")?;
    writeln!(writer, "{:<15} {:>8} {:>10} {:>10} {:>12}", 
        "Sample", "Present", "Allele1", "Allele2", "Source")?;
    writeln!(writer, "{}", "-".repeat(60))?;

    for gt in genotypes {
        let present_str = if gt.is_present { "Yes" } else { "No" };
        let a1_str = gt.allele1.map(|a| format!("{}", a)).unwrap_or_else(|| ".".to_string());
        let a2_str = gt.allele2.map(|a| format!("{}", a)).unwrap_or_else(|| ".".to_string());
        let source_str = format!("{:?}", gt.allele_source);
        
        writeln!(writer, "{:<15} {:>8} {:>10} {:>10} {:>12}",
            gt.sample_id, present_str, a1_str, a2_str, source_str)?;
    }

    Ok(())
}

/// Format a genotype as a compact string
pub fn format_genotype_compact(genotype: &ResolvedGenotype) -> String {
    let a1 = genotype.allele1.map(|a| a.to_string()).unwrap_or_else(|| ".".to_string());
    let a2 = genotype.allele2.map(|a| a.to_string()).unwrap_or_else(|| ".".to_string());
    let present = if genotype.is_present { "+" } else { "-" };
    format!("{}[{}/{}]", present, a1, a2)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::resolver::{PresenceSource, AlleleSource};

    #[test]
    fn test_explain_genotype() {
        let unit = TestUnit::from_true_repeat("TR001", "chr1", 10000, 10050, "CAG", "TR001");
        let gt = ResolvedGenotype {
            sample_id: "SAMPLE1".to_string(),
            is_present: true,
            presence_source: PresenceSource::TrgtVcf,
            allele1: Some(30),
            allele2: Some(45),
            allele_source: AlleleSource::TrgtTrue,
            raw_gt: Some("0/1".to_string()),
        };

        let mut output = Vec::new();
        explain_genotype(&mut output, &unit, &gt).expect("Failed to write");
        
        let output_str = String::from_utf8(output).expect("Invalid UTF-8");
        assert!(output_str.contains("Genotype Explanation"));
        assert!(output_str.contains("SAMPLE1"));
        assert!(output_str.contains("PRESENT"));
        assert!(output_str.contains("30 bp"));
        assert!(output_str.contains("45 bp"));
        assert!(output_str.contains("S (sum): 75 bp"));
        assert!(output_str.contains("M (max): 45 bp"));
        assert!(output_str.contains("D (diff): 15 bp"));
    }

    #[test]
    fn test_explain_locus() {
        let unit = TestUnit::from_true_repeat("TR001", "chr1", 10000, 10050, "CAG", "TR001");
        let gts = vec![
            ResolvedGenotype {
                sample_id: "SAMPLE1".to_string(),
                is_present: true,
                presence_source: PresenceSource::TrgtVcf,
                allele1: Some(30),
                allele2: Some(45),
                allele_source: AlleleSource::TrgtTrue,
                raw_gt: None,
            },
            ResolvedGenotype {
                sample_id: "SAMPLE2".to_string(),
                is_present: false,
                presence_source: PresenceSource::TrgtVcf,
                allele1: Some(30),
                allele2: Some(30),
                allele_source: AlleleSource::Reference,
                raw_gt: None,
            },
        ];

        let mut output = Vec::new();
        explain_locus(&mut output, &unit, &gts).expect("Failed to write");
        
        let output_str = String::from_utf8(output).expect("Invalid UTF-8");
        assert!(output_str.contains("Locus Summary"));
        assert!(output_str.contains("Total samples: 2"));
        assert!(output_str.contains("Present (carriers): 1"));
        assert!(output_str.contains("SAMPLE1"));
        assert!(output_str.contains("SAMPLE2"));
    }

    #[test]
    fn test_format_genotype_compact() {
        let gt = ResolvedGenotype {
            sample_id: "S1".to_string(),
            is_present: true,
            presence_source: PresenceSource::SvVcf,
            allele1: Some(30),
            allele2: Some(45),
            allele_source: AlleleSource::SvProxy,
            raw_gt: None,
        };

        assert_eq!(format_genotype_compact(&gt), "+[30/45]");

        let gt_absent = ResolvedGenotype {
            sample_id: "S2".to_string(),
            is_present: false,
            presence_source: PresenceSource::Unknown,
            allele1: None,
            allele2: None,
            allele_source: AlleleSource::Unknown,
            raw_gt: None,
        };

        assert_eq!(format_genotype_compact(&gt_absent), "-[./.]");
    }
}
