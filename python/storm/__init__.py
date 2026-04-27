"""Storm - A high-performance data processing framework.

This package provides Python bindings for the Storm data processing framework.
Storm is designed for high-performance data processing and analysis tasks.

Pure-Python modules (e.g. ``storm.tr_sidecar``) work without the compiled
extension; :func:`version` falls back when ``storm.storm`` is not built.

Public helpers are loaded lazily via :pep:`562` so ``import storm`` does not
import the native extension first (avoids partial package initialization on
some hosts).
"""

from __future__ import annotations

import importlib

__all__ = [
    "annotate_svs",
    "build_tr_sidecar",
    "build_feature_inventory",
    "annotate_saige_predictors",
    "build_long_predictor_tables",
    "build_feature_vcf_row_lookup",
    "align_matrix_cols_to_manifest",
    "build_dense_saige_marker_matrix",
    "dense_marker_matrix_from_long",
    "saige_dosage_vcf_metadata",
    "export_saige_dosage_vcf",
    "export_saige_stratum_vcfs",
    "saige_synthetic_covar_col_list",
    "build_synthetic_saige_pheno_covar_table",
    "export_saige_phenotype_covariate_tsv",
    "write_unrelated_identity_sparse_grm_for_saige",
    "print_feature_inventory_stats",
    "print_long_predictor_stats",
    "print_tr_annotation_summary",
    "build_predictor_feature_qc",
    "export_long_predictor_tables",
    "version",
]


def __getattr__(name: str):
    if name not in __all__:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    mod = importlib.import_module("._api", __name__)
    return getattr(mod, name)


def __dir__() -> list[str]:
    return list(__all__)
