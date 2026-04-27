"""Storm - A high-performance data processing framework.

This package provides Python bindings for the Storm data processing framework.
Storm is designed for high-performance data processing and analysis tasks.

Pure-Python modules (e.g. ``storm.tr_sidecar``) work without the compiled
extension; :func:`version` falls back when ``storm.storm`` is not built.

The implementation lives in :mod:`storm._api`, which avoids importing the
native ``storm.storm`` extension at import time (only :func:`version` loads it
when needed).
"""

from __future__ import annotations

import importlib

# Wheel / environment sanity check (present only in packages that bind API below).
__STORM_PACKAGE_API__ = "dynamic-v1"

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

_api_mod = importlib.import_module(__name__ + "._api")
for _public in __all__:
    globals()[_public] = getattr(_api_mod, _public)
del _api_mod, _public
