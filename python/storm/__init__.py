"""STORM - Structural & Tandem-Repeat Optimized Regression Models.

This package provides Python bindings for the STORM framework for
association testing of structural variants and tandem repeats.
"""

from pathlib import Path
from typing import Optional, List, Dict, Any

try:
    import storm.storm as _storm  # noqa: F403
    HAS_RUST = True
except ImportError:
    HAS_RUST = False
    _storm = None

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
    if HAS_RUST:
        return _storm._version()
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
        cache_dir = Path(output_dir)
        cache_dir.mkdir(parents=True, exist_ok=True)
        
        # TODO: Call Rust build_cache function when available
        # For now, return empty cache
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
    if not HAS_POLARS:
        raise ImportError("Polars is required for run_glm")
    
    # TODO: Implement full GLM pipeline
    # For now, return empty results frame
    results = pl.DataFrame({
        "unit_id": [],
        "beta": [],
        "se": [],
        "p_value": [],
        "n_samples": [],
        "n_carriers": [],
    })
    
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


# Convenience exports
__all__ = [
    "version",
    "StormCache",
    "load_cache",
    "run_glm",
    "explain",
]
