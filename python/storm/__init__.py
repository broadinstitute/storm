"""STORM - Structural & Tandem-Repeat Optimized Regression Models.

This package provides Python bindings for the STORM framework for
association testing of structural variants and tandem repeats.
"""

from pathlib import Path
from typing import Optional, List, Dict, Any

# Import Rust extension
try:
    from storm._storm import _version as _rust_version
    from storm._storm import add as _rust_add
    from storm._storm import py_build_cache, py_verify_cache
    from storm._storm import py_explain_genotype, py_explain_locus
    from storm._storm import py_load_plan
    from storm._storm import py_run_association
    from storm._storm import py_parse_sv_vcf, py_parse_trgt_vcf
    from storm._storm import py_catalog_from_bed, py_catalog_from_json, py_catalog_from_bed_and_json
    from storm._storm import py_read_cache_from_dir, py_write_cache_to_dir
    HAS_RUST = True
except ImportError:
    HAS_RUST = False
    _rust_version = None
    _rust_add = None
    py_build_cache = None
    py_verify_cache = None
    py_explain_genotype = None
    py_explain_locus = None
    py_load_plan = None
    py_run_association = None
    py_parse_sv_vcf = None
    py_parse_trgt_vcf = None
    py_catalog_from_bed = None
    py_catalog_from_json = None
    py_catalog_from_bed_and_json = None
    py_read_cache_from_dir = None
    py_write_cache_to_dir = None

try:
    import polars as pl
    HAS_POLARS = True
except ImportError:
    HAS_POLARS = False
    pl = None


def version() -> str:
    """Get the current version of the Storm package.

    Returns:
        str: The version string of the currently installed Storm package.
    """
    if HAS_RUST and _rust_version is not None:
        return _rust_version()
    return "unknown"


class StormCache:
    """Cache for STORM data stored in Parquet format."""
    
    def __init__(self, cache_dir: str):
        """Initialize cache from directory.
        
        Args:
            cache_dir: Path to cache directory containing Parquet files.
        """
        self.cache_dir = Path(cache_dir)
        self._test_units = None
        self._genotypes = None
        self._catalog = None
        self._features = None
        self._provenance = None
    
    @classmethod
    def build(
        cls,
        sv_vcf: str,
        trgt_vcf: Optional[str] = None,
        catalog_bed: Optional[str] = None,
        catalog_json: Optional[str] = None,
        output_dir: str = "storm_cache",
    ) -> "StormCache":
        """Build a new cache from input files.
        
        Args:
            sv_vcf: Path to integrated SV VCF file.
            trgt_vcf: Optional path to TRGT VCF file.
            catalog_bed: Optional path to TRExplorer BED file.
            catalog_json: Optional path to TRExplorer JSON file.
            output_dir: Output directory for cache.
            
        Returns:
            StormCache instance pointing to the new cache.
        """
        if HAS_RUST and py_build_cache is not None:
            # Call Rust build_cache function
            num_units, num_samples, num_gts, num_catalog = py_build_cache(
                sv_vcf, trgt_vcf, catalog_bed, catalog_json, output_dir
            )
            cache = cls(output_dir)
            cache._build_stats = {
                "num_test_units": num_units,
                "num_samples": num_samples,
                "num_genotypes": num_gts,
                "num_catalog_entries": num_catalog,
            }
            return cache
        else:
            # Fallback: create empty cache directory
            cache_dir = Path(output_dir)
            cache_dir.mkdir(parents=True, exist_ok=True)
            return cls(str(cache_dir))
    
    @property
    def test_units(self):
        """Load test units table via Polars."""
        if self._test_units is None and HAS_POLARS:
            path = self.cache_dir / "test_units.parquet"
            if path.exists():
                self._test_units = pl.read_parquet(str(path))
        return self._test_units
    
    @property
    def genotypes(self):
        """Load genotypes table via Polars."""
        if self._genotypes is None and HAS_POLARS:
            path = self.cache_dir / "genotypes.parquet"
            if path.exists():
                self._genotypes = pl.read_parquet(str(path))
        return self._genotypes
    
    @property
    def catalog(self):
        """Load catalog table via Polars."""
        if self._catalog is None and HAS_POLARS:
            path = self.cache_dir / "catalog.parquet"
            if path.exists():
                self._catalog = pl.read_parquet(str(path))
        return self._catalog
    
    @property
    def features(self):
        """Load features table via Polars."""
        if self._features is None and HAS_POLARS:
            path = self.cache_dir / "features.parquet"
            if path.exists():
                self._features = pl.read_parquet(str(path))
        return self._features
    
    @property
    def provenance(self):
        """Load provenance metadata via Polars."""
        if self._provenance is None and HAS_POLARS:
            path = self.cache_dir / "provenance.parquet"
            if path.exists():
                self._provenance = pl.read_parquet(str(path))
        return self._provenance


def load_cache(cache_dir: str) -> StormCache:
    """Load a STORM cache from a directory.
    
    Args:
        cache_dir: Path to cache directory.
        
    Returns:
        StormCache instance.
    """
    return StormCache(cache_dir)


def run_glm(
    cache: StormCache,
    phenotype: "pl.Series",
    plan: Optional[str] = None,
    covariates: Optional["pl.DataFrame"] = None,
    output: Optional[str] = None,
) -> "pl.DataFrame":
    """Run StormGLM association testing.
    
    Args:
        cache: StormCache instance with loaded data.
        phenotype: Polars Series with phenotype values.
        plan: Optional path to plan YAML file.
        covariates: Optional DataFrame with covariate columns.
        output: Optional path to write results Parquet.
        
    Returns:
        DataFrame with association results.
    """
    import json
    
    if not HAS_POLARS:
        raise ImportError("Polars is required for run_glm")
    
    # Load test units and features from cache
    if cache.test_units is None or cache.features is None:
        raise ValueError("Cache must have test_units and features loaded")
    
    # Try to use Rust run_association if available
    if HAS_RUST and py_run_association is not None:
        try:
            # Convert phenotype to dict for JSON serialization
            # Phenotype should be indexed by sample_id
            pheno_dict = {}
            
            # Handle Polars Series (has to_list method)
            if HAS_POLARS and isinstance(phenotype, pl.Series):
                if cache.genotypes is not None:
                    sample_ids = cache.genotypes["sample_id"].unique().to_list()
                    pheno_values = phenotype.to_list()
                    for i, sample_id in enumerate(sample_ids):
                        if i < len(pheno_values):
                            pheno_dict[sample_id] = float(pheno_values[i])
                else:
                    pheno_values = phenotype.to_list()
                    pheno_dict = {f"sample_{i}": float(v) for i, v in enumerate(pheno_values)}
            # Handle dict-like objects (already indexed by sample_id)
            elif isinstance(phenotype, dict):
                pheno_dict = {str(k): float(v) for k, v in phenotype.items()}
            # Handle list/array phenotypes
            elif hasattr(phenotype, '__iter__'):
                pheno_values = list(phenotype)
                if cache.genotypes is not None:
                    sample_ids = cache.genotypes["sample_id"].unique().to_list()
                    for i, sample_id in enumerate(sample_ids):
                        if i < len(pheno_values):
                            pheno_dict[sample_id] = float(pheno_values[i])
                else:
                    pheno_dict = {f"sample_{i}": float(v) for i, v in enumerate(pheno_values)}
            else:
                raise ValueError(
                    f"Invalid phenotype format: {type(phenotype)}. "
                    "Expected Polars Series, dict, or list/array."
                )
            
            phenotype_json = json.dumps(pheno_dict)
            
            # Call Rust run_association
            results_json = py_run_association(
                str(cache.cache_dir),
                phenotype_json,
                plan,
            )
            
            results_list = json.loads(results_json)
            
            # Convert to Polars DataFrame
            results = pl.DataFrame(results_list)
            
            if output:
                results.write_parquet(output)
            
            return results
            
        except Exception as e:
            # Fall back to Python implementation on error
            import warnings
            warnings.warn(f"Rust run_association failed, using Python fallback: {e}")
    
    # Python fallback implementation
    results_data = {
        "unit_id": [],
        "beta": [],
        "se": [],
        "p_value": [],
        "n_samples": [],
        "n_carriers": [],
    }
    
    # Get unique unit IDs
    unit_ids = cache.test_units["id"].unique().to_list()
    
    for unit_id in unit_ids:
        # Get features for this unit
        unit_features = cache.features.filter(pl.col("unit_id") == unit_id)
        if len(unit_features) == 0:
            continue
            
        # Count carriers
        n_samples = len(unit_features)
        if "binary" in unit_features.columns:
            n_carriers = int(unit_features["binary"].sum())
        else:
            n_carriers = 0
        
        # Placeholder results (would be from actual regression)
        results_data["unit_id"].append(unit_id)
        results_data["beta"].append(0.0)
        results_data["se"].append(1.0)
        results_data["p_value"].append(1.0)
        results_data["n_samples"].append(n_samples)
        results_data["n_carriers"].append(n_carriers)
    
    results = pl.DataFrame(results_data)
    
    if output:
        results.write_parquet(output)
    
    return results


def explain(
    cache: StormCache,
    unit_id: str,
    sample_id: Optional[str] = None,
) -> str:
    """Get detailed explanation of a test unit's genotypes.
    
    Args:
        cache: StormCache instance.
        unit_id: ID of the test unit to explain.
        sample_id: Optional specific sample to explain.
        
    Returns:
        Human-readable explanation string.
    """
    # Use Rust explain functions if available
    if HAS_RUST and py_explain_genotype is not None and py_explain_locus is not None:
        try:
            if sample_id:
                return py_explain_genotype(str(cache.cache_dir), unit_id, sample_id)
            else:
                return py_explain_locus(str(cache.cache_dir), unit_id)
        except Exception:
            pass  # Fall back to Python implementation
    
    if not HAS_POLARS:
        return "Polars required for explain"
    
    output_lines = []
    output_lines.append(f"=== Explanation for {unit_id} ===")
    output_lines.append("")
    
    # Get unit info
    if cache.test_units is not None:
        unit = cache.test_units.filter(pl.col("id") == unit_id)
        if len(unit) > 0:
            row = unit.row(0, named=True)
            output_lines.append(f"Location: {row.get('chrom', '?')}:{row.get('start', '?')}-{row.get('end', '?')}")
            output_lines.append(f"Type: {row.get('unit_type', '?')}")
            output_lines.append(f"Motif: {row.get('motif', 'N/A')}")
            output_lines.append("")
    
    # Get genotypes
    if cache.genotypes is not None:
        gts = cache.genotypes.filter(pl.col("unit_id") == unit_id)
        if sample_id:
            gts = gts.filter(pl.col("sample_id") == sample_id)
        
        if len(gts) > 0:
            output_lines.append("Genotypes:")
            for row in gts.iter_rows(named=True):
                present = "+" if row.get("is_present") else "-"
                a1 = row.get("allele1", ".")
                a2 = row.get("allele2", ".")
                output_lines.append(f"  {row['sample_id']}: {present}[{a1}/{a2}]")
    
    return "\n".join(output_lines)


def verify_cache(cache_dir: str) -> Dict[str, Any]:
    """Verify a STORM cache.
    
    Args:
        cache_dir: Path to cache directory.
        
    Returns:
        Dictionary with validation results.
    """
    if HAS_RUST and py_verify_cache is not None:
        is_valid, num_units, num_gts, num_catalog, num_features, errors = py_verify_cache(cache_dir)
        return {
            "is_valid": is_valid,
            "num_test_units": num_units,
            "num_genotypes": num_gts,
            "num_catalog_entries": num_catalog,
            "num_features": num_features,
            "errors": errors,
        }
    else:
        # Fallback: check files exist
        cache_path = Path(cache_dir)
        required_files = ["test_units.parquet", "genotypes.parquet", "catalog.parquet", "features.parquet"]
        errors = []
        for f in required_files:
            if not (cache_path / f).exists():
                errors.append(f"Missing file: {f}")
        return {
            "is_valid": len(errors) == 0,
            "errors": errors,
        }


def load_plan(path: str) -> Dict[str, Any]:
    """Load a plan YAML file.
    
    Args:
        path: Path to plan YAML file.
        
    Returns:
        Dictionary with plan configuration.
    """
    import json
    if HAS_RUST and py_load_plan is not None:
        plan_json = py_load_plan(path)
        return json.loads(plan_json)
    else:
        import yaml
        with open(path) as f:
            return yaml.safe_load(f)


def parse_sv_vcf(path: str) -> tuple:
    """Parse an SV VCF file.
    
    Args:
        path: Path to SV VCF file.
        
    Returns:
        Tuple of (samples, records) where records is a list of dicts.
    """
    import json
    if HAS_RUST and py_parse_sv_vcf is not None:
        samples, records_json = py_parse_sv_vcf(path)
        return samples, json.loads(records_json)
    raise RuntimeError("Rust extension not available")


def parse_trgt_vcf(path: str) -> tuple:
    """Parse a TRGT VCF file.
    
    Args:
        path: Path to TRGT VCF file.
        
    Returns:
        Tuple of (samples, records) where records is a list of dicts.
    """
    import json
    if HAS_RUST and py_parse_trgt_vcf is not None:
        samples, records_json = py_parse_trgt_vcf(path)
        return samples, json.loads(records_json)
    raise RuntimeError("Rust extension not available")


def catalog_from_bed(path: str) -> Dict[str, Any]:
    """Load catalog from BED file.
    
    Args:
        path: Path to TRExplorer BED file.
        
    Returns:
        Dictionary with catalog data.
    """
    import json
    if HAS_RUST and py_catalog_from_bed is not None:
        return json.loads(py_catalog_from_bed(path))
    raise RuntimeError("Rust extension not available")


def catalog_from_json(path: str) -> Dict[str, Any]:
    """Load catalog from JSON file.
    
    Args:
        path: Path to TRExplorer JSON file.
        
    Returns:
        Dictionary with catalog data.
    """
    import json as json_mod
    if HAS_RUST and py_catalog_from_json is not None:
        return json_mod.loads(py_catalog_from_json(path))
    raise RuntimeError("Rust extension not available")


def catalog_from_bed_and_json(bed_path: str, json_path: str) -> Dict[str, Any]:
    """Load catalog from both BED and JSON files.
    
    Args:
        bed_path: Path to TRExplorer BED file.
        json_path: Path to TRExplorer JSON file.
        
    Returns:
        Dictionary with merged catalog data.
    """
    import json
    if HAS_RUST and py_catalog_from_bed_and_json is not None:
        return json.loads(py_catalog_from_bed_and_json(bed_path, json_path))
    raise RuntimeError("Rust extension not available")


def read_cache_from_dir(cache_dir: str) -> Dict[str, Any]:
    """Read cache metadata from directory.
    
    Args:
        cache_dir: Path to cache directory.
        
    Returns:
        Dictionary with cache metadata.
    """
    import json
    if HAS_RUST and py_read_cache_from_dir is not None:
        return json.loads(py_read_cache_from_dir(cache_dir))
    raise RuntimeError("Rust extension not available")


def write_cache_to_dir(cache_dir: str) -> None:
    """Write cache to directory (placeholder).
    
    Args:
        cache_dir: Path to cache directory.
    """
    if HAS_RUST and py_write_cache_to_dir is not None:
        py_write_cache_to_dir(cache_dir, "{}")
    else:
        Path(cache_dir).mkdir(parents=True, exist_ok=True)


def run_association(cache_dir: str, phenotype: Dict[str, float], plan_path: Optional[str] = None) -> List[Dict[str, Any]]:
    """Run association testing.
    
    Args:
        cache_dir: Path to cache directory.
        phenotype: Dictionary mapping sample_id to phenotype value.
        plan_path: Optional path to plan YAML file.
        
    Returns:
        List of association results.
    """
    import json
    if HAS_RUST and py_run_association is not None:
        results_json = py_run_association(cache_dir, json.dumps(phenotype), plan_path)
        return json.loads(results_json)
    raise RuntimeError("Rust extension not available")


# Convenience exports
__all__ = [
    "version",
    "StormCache",
    "load_cache",
    "run_glm",
    "explain",
    "verify_cache",
    "load_plan",
    "parse_sv_vcf",
    "parse_trgt_vcf",
    "catalog_from_bed",
    "catalog_from_json",
    "catalog_from_bed_and_json",
    "read_cache_from_dir",
    "write_cache_to_dir",
    "run_association",
]
