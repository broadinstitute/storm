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


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
