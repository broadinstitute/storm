//! STORM CLI - Structural & Tandem-Repeat Optimized Regression Models
//!
//! Command-line interface for the STORM association testing framework.

use std::path::PathBuf;
use clap::{Parser, Subcommand};

use storm::{
    build_cache, verify_cache, explain_genotype, explain_locus,
    read_cache_from_dir,
};

/// STORM: Structural & Tandem-Repeat Optimized Regression Models
#[derive(Parser)]
#[command(name = "storm")]
#[command(version = env!("CARGO_PKG_VERSION"))]
#[command(about = "Association testing for structural variants and tandem repeats")]
#[command(long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Option<Commands>,
}

#[derive(Subcommand)]
enum Commands {
    /// Build or verify a STORM cache
    Cache {
        #[command(subcommand)]
        action: CacheAction,
    },
    /// Explain genotypes for a test unit
    Explain {
        /// Test unit ID to explain
        test_unit_id: String,

        /// Path to cache directory
        #[arg(long, default_value = "storm_cache")]
        cache_dir: PathBuf,

        /// Specific sample ID to explain (optional)
        #[arg(long)]
        sample: Option<String>,
    },
}

#[derive(Subcommand)]
enum CacheAction {
    /// Build a cache from input files
    Build {
        /// Path to structural variant VCF file (required)
        #[arg(long, required = true)]
        sv_vcf: PathBuf,

        /// Path to TRGT VCF file (optional)
        #[arg(long)]
        trgt_vcf: Option<PathBuf>,

        /// Path to TRExplorer BED catalog file (optional)
        #[arg(long)]
        catalog_bed: Option<PathBuf>,

        /// Path to TRExplorer JSON catalog file (optional)
        #[arg(long)]
        catalog_json: Option<PathBuf>,

        /// Output directory for cache files
        #[arg(long, default_value = "storm_cache")]
        output_dir: PathBuf,
    },
    /// Verify an existing cache
    Verify {
        /// Path to cache directory
        #[arg(long, default_value = "storm_cache")]
        cache_dir: PathBuf,
    },
}

fn main() {
    let cli = Cli::parse();

    let result = match cli.command {
        Some(Commands::Cache { action }) => match action {
            CacheAction::Build {
                sv_vcf,
                trgt_vcf,
                catalog_bed,
                catalog_json,
                output_dir,
            } => run_cache_build(sv_vcf, trgt_vcf, catalog_bed, catalog_json, output_dir),
            CacheAction::Verify { cache_dir } => run_cache_verify(cache_dir),
        },
        Some(Commands::Explain {
            test_unit_id,
            cache_dir,
            sample,
        }) => run_explain(test_unit_id, cache_dir, sample),
        None => {
            // No subcommand provided, print help
            println!("STORM v{}", env!("CARGO_PKG_VERSION"));
            println!("Use --help for usage information");
            Ok(())
        }
    };

    if let Err(e) = result {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }
}

fn run_cache_build(
    sv_vcf: PathBuf,
    trgt_vcf: Option<PathBuf>,
    catalog_bed: Option<PathBuf>,
    catalog_json: Option<PathBuf>,
    output_dir: PathBuf,
) -> Result<(), Box<dyn std::error::Error>> {
    println!("Building STORM cache...");
    println!("  SV VCF: {}", sv_vcf.display());
    if let Some(ref p) = trgt_vcf {
        println!("  TRGT VCF: {}", p.display());
    }
    if let Some(ref p) = catalog_bed {
        println!("  Catalog BED: {}", p.display());
    }
    if let Some(ref p) = catalog_json {
        println!("  Catalog JSON: {}", p.display());
    }
    println!("  Output: {}", output_dir.display());

    // Validate inputs exist
    if !sv_vcf.exists() {
        return Err(format!("SV VCF file not found: {}", sv_vcf.display()).into());
    }
    if let Some(ref p) = trgt_vcf {
        if !p.exists() {
            return Err(format!("TRGT VCF file not found: {}", p.display()).into());
        }
    }
    if let Some(ref p) = catalog_bed {
        if !p.exists() {
            return Err(format!("Catalog BED file not found: {}", p.display()).into());
        }
    }
    if let Some(ref p) = catalog_json {
        if !p.exists() {
            return Err(format!("Catalog JSON file not found: {}", p.display()).into());
        }
    }

    // Run the build pipeline
    let stats = build_cache(
        sv_vcf.to_str().unwrap(),
        trgt_vcf.as_ref().map(|p| p.to_str().unwrap()),
        catalog_bed.as_ref().map(|p| p.to_str().unwrap()),
        catalog_json.as_ref().map(|p| p.to_str().unwrap()),
        output_dir.to_str().unwrap(),
    )?;

    println!("\nCache built successfully!");
    println!("  Test units: {}", stats.num_test_units);
    println!("  Samples: {}", stats.num_samples);
    println!("  Genotypes: {}", stats.num_genotypes);
    println!("  Output files written to: {}", output_dir.display());

    Ok(())
}

fn run_cache_verify(cache_dir: PathBuf) -> Result<(), Box<dyn std::error::Error>> {
    println!("Verifying STORM cache at {}...", cache_dir.display());

    if !cache_dir.exists() {
        return Err(format!("Cache directory not found: {}", cache_dir.display()).into());
    }

    let result = verify_cache(cache_dir.to_str().unwrap())?;

    if result.is_valid {
        println!("\n✓ Cache is valid!");
        println!("  Test units: {}", result.num_test_units);
        println!("  Genotypes: {}", result.num_genotypes);
        println!("  Catalog entries: {}", result.num_catalog_entries);
        println!("  Features: {}", result.num_features);
    } else {
        println!("\n✗ Cache validation failed:");
        for error in &result.errors {
            println!("  - {}", error);
        }
        return Err("Cache validation failed".into());
    }

    Ok(())
}

fn run_explain(
    test_unit_id: String,
    cache_dir: PathBuf,
    sample: Option<String>,
) -> Result<(), Box<dyn std::error::Error>> {
    if !cache_dir.exists() {
        return Err(format!("Cache directory not found: {}", cache_dir.display()).into());
    }

    let cache = read_cache_from_dir(&cache_dir)?;

    let explanation = if let Some(sample_id) = sample {
        explain_genotype(&cache, &test_unit_id, &sample_id)?
    } else {
        explain_locus(&cache, &test_unit_id)?
    };

    println!("{}", explanation);

    Ok(())
}
