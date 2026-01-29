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
        /// Path to structural variant VCF or BCF file (required)
        #[arg(long, required = true)]
        sv_vcf: PathBuf,

        /// Paths to TRGT VCF files (optional, supports .vcf and .vcf.gz)
        /// Can be specified multiple times: --trgt-vcf file1.vcf --trgt-vcf file2.vcf.gz
        #[arg(long)]
        trgt_vcf: Vec<PathBuf>,
        
        /// Directory containing TRGT VCF files (will glob for *.vcf and *.vcf.gz)
        #[arg(long)]
        trgt_dir: Option<PathBuf>,
        
        /// File containing list of TRGT VCF paths (one path per line)
        #[arg(long)]
        trgt_list: Option<PathBuf>,

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
                trgt_dir,
                trgt_list,
                catalog_bed,
                catalog_json,
                output_dir,
            } => run_cache_build(sv_vcf, trgt_vcf, trgt_dir, trgt_list, catalog_bed, catalog_json, output_dir),
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
    mut trgt_vcf: Vec<PathBuf>,
    trgt_dir: Option<PathBuf>,
    trgt_list: Option<PathBuf>,
    catalog_bed: Option<PathBuf>,
    catalog_json: Option<PathBuf>,
    output_dir: PathBuf,
) -> Result<(), Box<dyn std::error::Error>> {
    // Expand trgt_dir to individual files
    if let Some(dir) = trgt_dir {
        if !dir.exists() {
            return Err(format!("TRGT directory not found: {}", dir.display()).into());
        }
        // Glob for VCF files in the directory
        for entry in std::fs::read_dir(&dir)? {
            let entry = entry?;
            let path = entry.path();
            if let Some(ext) = path.extension() {
                if ext == "vcf" || ext == "gz" {
                    // Check if it's a .vcf or .vcf.gz file
                    let name = path.file_name().unwrap().to_str().unwrap();
                    if name.ends_with(".vcf") || name.ends_with(".vcf.gz") {
                        trgt_vcf.push(path);
                    }
                }
            }
        }
    }
    
    // Read paths from trgt_list file
    if let Some(list_path) = trgt_list {
        if !list_path.exists() {
            return Err(format!("TRGT list file not found: {}", list_path.display()).into());
        }
        let contents = std::fs::read_to_string(&list_path)?;
        for line in contents.lines() {
            let line = line.trim();
            if !line.is_empty() && !line.starts_with('#') {
                trgt_vcf.push(PathBuf::from(line));
            }
        }
    }
    
    println!("Building STORM cache...");
    println!("  SV VCF: {}", sv_vcf.display());
    if !trgt_vcf.is_empty() {
        println!("  TRGT VCF files: {} file(s)", trgt_vcf.len());
        for p in trgt_vcf.iter().take(5) {
            println!("    - {}", p.display());
        }
        if trgt_vcf.len() > 5 {
            println!("    ... and {} more", trgt_vcf.len() - 5);
        }
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
    for p in &trgt_vcf {
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

    // Convert TRGT paths for build_cache
    let trgt_paths: Vec<&str> = trgt_vcf.iter().map(|p| p.to_str().unwrap()).collect();
    let trgt_paths_opt: Option<&[&str]> = if trgt_paths.is_empty() {
        None
    } else {
        Some(&trgt_paths)
    };

    // Run the build pipeline
    let stats = build_cache(
        sv_vcf.to_str().unwrap(),
        trgt_paths_opt,
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
