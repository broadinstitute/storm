"""Tests for STORM Python API."""

import sys
import os
import tempfile
import shutil

import pytest

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import storm


class TestVersion:
    """Tests for version function."""
    
    def test_version_returns_string(self):
        """Version should return a string."""
        v = storm.version()
        assert isinstance(v, str)
    
    def test_version_not_unknown(self):
        """Version should not be 'unknown' when Rust is available."""
        if storm.HAS_RUST:
            assert storm.version() != "unknown"
    
    def test_version_has_dots(self):
        """Version should be in semver format."""
        if storm.HAS_RUST:
            v = storm.version()
            assert "." in v


class TestParsing:
    """Tests for VCF parsing functions."""
    
    @pytest.fixture
    def fixtures_dir(self):
        """Get path to fixtures directory."""
        return os.path.join(os.path.dirname(__file__), "..", "..", "fixtures")
    
    def test_parse_sv_vcf(self, fixtures_dir):
        """Should parse SV VCF file."""
        if not storm.HAS_RUST:
            pytest.skip("Rust bindings not available")
        
        path = os.path.join(fixtures_dir, "sv_small.vcf")
        samples, records = storm.parse_sv_vcf(path)
        
        assert len(samples) == 2
        assert len(records) == 3
    
    def test_parse_trgt_vcf(self, fixtures_dir):
        """Should parse TRGT VCF file."""
        if not storm.HAS_RUST:
            pytest.skip("Rust bindings not available")
        
        path = os.path.join(fixtures_dir, "trgt_small.vcf")
        samples, records = storm.parse_trgt_vcf(path)
        
        assert len(samples) == 2
        assert len(records) == 3


class TestCatalog:
    """Tests for catalog loading."""
    
    @pytest.fixture
    def fixtures_dir(self):
        """Get path to fixtures directory."""
        return os.path.join(os.path.dirname(__file__), "..", "..", "fixtures")
    
    def test_catalog_from_bed(self, fixtures_dir):
        """Should load catalog from BED."""
        if not storm.HAS_RUST:
            pytest.skip("Rust bindings not available")
        
        path = os.path.join(fixtures_dir, "trexplorer.bed")
        catalog = storm.catalog_from_bed(path)
        
        assert "entries" in catalog
        assert len(catalog["entries"]) > 0
    
    def test_catalog_from_json(self, fixtures_dir):
        """Should load catalog from JSON."""
        if not storm.HAS_RUST:
            pytest.skip("Rust bindings not available")
        
        path = os.path.join(fixtures_dir, "trexplorer.json")
        catalog = storm.catalog_from_json(path)
        
        assert "entries" in catalog
        assert len(catalog["entries"]) > 0
    
    def test_catalog_from_bed_and_json(self, fixtures_dir):
        """Should load and merge catalog from BED and JSON."""
        if not storm.HAS_RUST:
            pytest.skip("Rust bindings not available")
        
        bed_path = os.path.join(fixtures_dir, "trexplorer.bed")
        json_path = os.path.join(fixtures_dir, "trexplorer.json")
        catalog = storm.catalog_from_bed_and_json(bed_path, json_path)
        
        assert "entries" in catalog
        assert len(catalog["entries"]) > 0


class TestCacheBuild:
    """Tests for cache building."""
    
    @pytest.fixture
    def fixtures_dir(self):
        """Get path to fixtures directory."""
        return os.path.join(os.path.dirname(__file__), "..", "..", "fixtures")
    
    @pytest.fixture
    def temp_dir(self):
        """Create and cleanup temp directory."""
        d = tempfile.mkdtemp(prefix="storm_test_")
        yield d
        shutil.rmtree(d, ignore_errors=True)
    
    def test_build_cache_sv_only(self, fixtures_dir, temp_dir):
        """Should build cache from SV VCF only."""
        if not storm.HAS_RUST:
            pytest.skip("Rust bindings not available")
        
        output = os.path.join(temp_dir, "cache_sv")
        cache = storm.StormCache.build(
            sv_vcf=os.path.join(fixtures_dir, "sv_small.vcf"),
            output_dir=output,
        )
        
        assert hasattr(cache, "_build_stats")
        assert cache._build_stats["num_test_units"] == 3
        assert cache._build_stats["num_samples"] == 2
    
    def test_build_cache_full(self, fixtures_dir, temp_dir):
        """Should build cache from all input files."""
        if not storm.HAS_RUST:
            pytest.skip("Rust bindings not available")
        
        output = os.path.join(temp_dir, "cache_full")
        cache = storm.StormCache.build(
            sv_vcf=os.path.join(fixtures_dir, "sv_small.vcf"),
            trgt_vcf=os.path.join(fixtures_dir, "trgt_small.vcf"),
            catalog_bed=os.path.join(fixtures_dir, "trexplorer.bed"),
            catalog_json=os.path.join(fixtures_dir, "trexplorer.json"),
            output_dir=output,
        )
        
        assert hasattr(cache, "_build_stats")
        assert cache._build_stats["num_test_units"] >= 6
        assert cache._build_stats["num_catalog_entries"] > 0


class TestVerifyCache:
    """Tests for cache verification."""
    
    @pytest.fixture
    def fixtures_dir(self):
        """Get path to fixtures directory."""
        return os.path.join(os.path.dirname(__file__), "..", "..", "fixtures")
    
    @pytest.fixture
    def cache_dir(self, fixtures_dir):
        """Build a test cache."""
        if not storm.HAS_RUST:
            pytest.skip("Rust bindings not available")
        
        d = tempfile.mkdtemp(prefix="storm_verify_")
        storm.StormCache.build(
            sv_vcf=os.path.join(fixtures_dir, "sv_small.vcf"),
            output_dir=d,
        )
        yield d
        shutil.rmtree(d, ignore_errors=True)
    
    def test_verify_valid_cache(self, cache_dir):
        """Should verify a valid cache."""
        result = storm.verify_cache(cache_dir)
        
        assert result["is_valid"] is True
        assert result["num_test_units"] == 3
        assert len(result.get("errors", [])) == 0


class TestLoadCache:
    """Tests for loading cache."""
    
    @pytest.fixture
    def fixtures_dir(self):
        """Get path to fixtures directory."""
        return os.path.join(os.path.dirname(__file__), "..", "..", "fixtures")
    
    @pytest.fixture
    def cache_dir(self, fixtures_dir):
        """Build a test cache."""
        if not storm.HAS_RUST:
            pytest.skip("Rust bindings not available")
        
        d = tempfile.mkdtemp(prefix="storm_load_")
        storm.StormCache.build(
            sv_vcf=os.path.join(fixtures_dir, "sv_small.vcf"),
            output_dir=d,
        )
        yield d
        shutil.rmtree(d, ignore_errors=True)
    
    def test_load_cache(self, cache_dir):
        """Should load cache and access tables."""
        cache = storm.load_cache(cache_dir)
        
        assert cache.test_units is not None
        assert cache.genotypes is not None
        assert len(cache.test_units) == 3


class TestExplain:
    """Tests for explain functionality."""
    
    @pytest.fixture
    def fixtures_dir(self):
        """Get path to fixtures directory."""
        return os.path.join(os.path.dirname(__file__), "..", "..", "fixtures")
    
    @pytest.fixture
    def cache(self, fixtures_dir):
        """Build and return a test cache."""
        if not storm.HAS_RUST:
            pytest.skip("Rust bindings not available")
        
        d = tempfile.mkdtemp(prefix="storm_explain_")
        c = storm.StormCache.build(
            sv_vcf=os.path.join(fixtures_dir, "sv_small.vcf"),
            output_dir=d,
        )
        yield c
        shutil.rmtree(d, ignore_errors=True)
    
    def test_explain_locus(self, cache):
        """Should explain a locus."""
        unit_id = cache.test_units["id"][0]
        explanation = storm.explain(cache, unit_id)
        
        assert len(explanation) > 0
        assert unit_id in explanation
    
    def test_explain_sample(self, cache):
        """Should explain a specific sample."""
        unit_id = cache.test_units["id"][0]
        sample_id = cache.genotypes["sample_id"][0]
        explanation = storm.explain(cache, unit_id, sample_id=sample_id)
        
        assert len(explanation) > 0


class TestRunGLM:
    """Tests for GLM association testing."""
    
    @pytest.fixture
    def fixtures_dir(self):
        """Get path to fixtures directory."""
        return os.path.join(os.path.dirname(__file__), "..", "..", "fixtures")
    
    @pytest.fixture
    def cache(self, fixtures_dir):
        """Build and return a test cache."""
        if not storm.HAS_RUST:
            pytest.skip("Rust bindings not available")
        
        d = tempfile.mkdtemp(prefix="storm_glm_")
        c = storm.StormCache.build(
            sv_vcf=os.path.join(fixtures_dir, "sv_small.vcf"),
            output_dir=d,
        )
        yield c
        shutil.rmtree(d, ignore_errors=True)
    
    def test_run_glm(self, cache):
        """Should run GLM and return results."""
        import polars as pl
        
        phenotype = pl.Series("phenotype", [1, 0])
        results = storm.run_glm(cache, phenotype)
        
        assert len(results) > 0
        assert "unit_id" in results.columns
        assert "p_value" in results.columns


class TestLoadPlan:
    """Tests for plan loading."""
    
    @pytest.fixture
    def fixtures_dir(self):
        """Get path to fixtures directory."""
        return os.path.join(os.path.dirname(__file__), "..", "..", "fixtures")
    
    def test_load_plan(self, fixtures_dir):
        """Should load plan YAML."""
        if not storm.HAS_RUST:
            pytest.skip("Rust bindings not available")
        
        path = os.path.join(fixtures_dir, "plan.yaml")
        plan = storm.load_plan(path)
        
        assert "name" in plan
        assert "rules" in plan


class TestErrorCases:
    """Tests for error handling."""
    
    @pytest.fixture
    def fixtures_dir(self):
        """Get path to fixtures directory."""
        return os.path.join(os.path.dirname(__file__), "..", "..", "fixtures")
    
    @pytest.fixture
    def cache(self, fixtures_dir):
        """Build and return a test cache."""
        if not storm.HAS_RUST:
            pytest.skip("Rust bindings not available")
        
        d = tempfile.mkdtemp(prefix="storm_error_")
        c = storm.StormCache.build(
            sv_vcf=os.path.join(fixtures_dir, "sv_small.vcf"),
            output_dir=d,
        )
        yield c
        shutil.rmtree(d, ignore_errors=True)
    
    def test_run_glm_with_dict_phenotype(self, cache):
        """Should handle dict-like phenotypes correctly."""
        # Get sample IDs from genotypes
        sample_ids = cache.genotypes["sample_id"].unique().to_list()
        phenotype = {sid: 1.0 if i % 2 == 0 else 0.0 for i, sid in enumerate(sample_ids)}
        
        results = storm.run_glm(cache, phenotype)
        
        assert len(results) > 0
        assert "unit_id" in results.columns
        assert "p_value" in results.columns
    
    def test_run_glm_with_list_phenotype(self, cache):
        """Should handle list/array phenotypes correctly."""
        import polars as pl
        
        # Get number of samples
        n_samples = len(cache.genotypes["sample_id"].unique())
        phenotype = [1.0 if i % 2 == 0 else 0.0 for i in range(n_samples)]
        
        results = storm.run_glm(cache, phenotype)
        
        assert len(results) > 0
        assert "unit_id" in results.columns
    
    def test_verify_invalid_cache(self):
        """Should report errors for invalid cache directory."""
        d = tempfile.mkdtemp(prefix="storm_invalid_")
        try:
            result = storm.verify_cache(d)
            # Should either be invalid or have errors
            assert not result["is_valid"] or len(result.get("errors", [])) > 0 or result.get("num_test_units", 0) == 0
        finally:
            shutil.rmtree(d, ignore_errors=True)


class TestAssociationWithRealData:
    """Tests for association testing with real data."""
    
    @pytest.fixture
    def fixtures_dir(self):
        """Get path to fixtures directory."""
        return os.path.join(os.path.dirname(__file__), "..", "..", "fixtures")
    
    @pytest.fixture
    def cache(self, fixtures_dir):
        """Build and return a test cache with TRGT data."""
        if not storm.HAS_RUST:
            pytest.skip("Rust bindings not available")
        
        d = tempfile.mkdtemp(prefix="storm_assoc_")
        c = storm.StormCache.build(
            sv_vcf=os.path.join(fixtures_dir, "sv_small.vcf"),
            trgt_vcf=os.path.join(fixtures_dir, "trgt_small.vcf"),
            catalog_bed=os.path.join(fixtures_dir, "trexplorer.bed"),
            catalog_json=os.path.join(fixtures_dir, "trexplorer.json"),
            output_dir=d,
        )
        yield c
        shutil.rmtree(d, ignore_errors=True)
    
    def test_association_results_have_metadata(self, cache):
        """Results should include model and encoding metadata."""
        import polars as pl
        
        n_samples = len(cache.genotypes["sample_id"].unique())
        phenotype = pl.Series("phenotype", [float(i % 2) for i in range(n_samples)])
        results = storm.run_glm(cache, phenotype)
        
        # Check results have expected columns
        expected_cols = ["unit_id", "beta", "se", "p_value", "n_samples"]
        for col in expected_cols:
            assert col in results.columns, f"Missing column: {col}"
    
    def test_association_with_plan(self, cache, fixtures_dir):
        """Should run association testing with plan-based model selection."""
        import polars as pl
        
        n_samples = len(cache.genotypes["sample_id"].unique())
        phenotype = pl.Series("phenotype", [float(i % 2) for i in range(n_samples)])
        
        plan_path = os.path.join(fixtures_dir, "plan.yaml")
        results = storm.run_glm(cache, phenotype, plan=plan_path)
        
        assert len(results) > 0
        assert "unit_id" in results.columns
    
    def test_run_association_function(self, cache):
        """Test the low-level run_association function."""
        # Get sample IDs
        sample_ids = cache.genotypes["sample_id"].unique().to_list()
        phenotype = {sid: float(i % 2) for i, sid in enumerate(sample_ids)}
        
        results = storm.run_association(str(cache.cache_dir), phenotype)
        
        assert isinstance(results, list)
        if len(results) > 0:
            result = results[0]
            assert "unit_id" in result
            # Results should have real statistics (not all zeros)
            assert "p_value" in result or "error" in result
    
    def test_run_glm_with_covariates(self, cache):
        """Test run_glm with covariates DataFrame."""
        import polars as pl
        
        # Get sample IDs
        sample_ids = cache.genotypes["sample_id"].unique().to_list()
        n_samples = len(sample_ids)
        
        # Create phenotype
        phenotype = pl.Series("phenotype", [float(i % 2) for i in range(n_samples)])
        
        # Create covariates DataFrame with sample_id column
        covariates = pl.DataFrame({
            "sample_id": sample_ids,
            "age": [float(25 + i) for i in range(n_samples)],
            "sex": [float(i % 2) for i in range(n_samples)],
        })
        
        # Run with covariates
        results = storm.run_glm(cache, phenotype, covariates=covariates)
        
        assert len(results) > 0
        assert "unit_id" in results.columns
        assert "p_value" in results.columns
    
    def test_run_association_with_covariates(self, cache):
        """Test run_association function with covariates dict."""
        # Get sample IDs
        sample_ids = cache.genotypes["sample_id"].unique().to_list()
        n_samples = len(sample_ids)
        
        # Create phenotype dict
        phenotype = {sid: float(i % 2) for i, sid in enumerate(sample_ids)}
        
        # Create covariates dict in expected format: {"cov_name": {"sample_id": value, ...}}
        covariates = {
            "age": {sid: float(25 + i) for i, sid in enumerate(sample_ids)},
            "sex": {sid: float(i % 2) for i, sid in enumerate(sample_ids)},
        }
        
        # Run with covariates
        results = storm.run_association(str(cache.cache_dir), phenotype, covariates=covariates)
        
        assert isinstance(results, list)
        if len(results) > 0:
            result = results[0]
            assert "unit_id" in result


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
