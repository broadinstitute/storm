"""Storm - A high-performance data processing framework.

This package provides Python bindings for the Storm data processing framework.
Storm is designed for high-performance data processing and analysis tasks.

Pure-Python modules (e.g. ``storm.tr_sidecar``) work without the compiled
extension; ``version()`` falls back when ``storm.storm`` is not built.

:func:`annotate_svs` and :func:`build_tr_sidecar` use lazy imports so
``import storm`` stays lightweight.
"""

from __future__ import annotations

import importlib.metadata

try:
    import storm.storm as _storm_ext  # type: ignore[import-not-found]
except ImportError:
    _storm_ext = None

__all__ = [
    "annotate_svs",
    "build_tr_sidecar",
    "build_feature_inventory",
    "annotate_saige_predictors",
    "build_long_predictor_tables",
    "print_feature_inventory_stats",
    "print_long_predictor_stats",
    "print_tr_annotation_summary",
    "version",
]


def annotate_svs(mt, *, tr_sidecar_sites=None, keep_star=False, left_aligned=False):
    """Normalize multi-allelic SVs and optionally join TR ``sites.parquet``.

    See :mod:`storm.sv_annotate`.
    """
    from storm.sv_annotate import annotate_svs as _annotate_svs

    return _annotate_svs(
        mt,
        tr_sidecar_sites=tr_sidecar_sites,
        keep_star=keep_star,
        left_aligned=left_aligned,
    )


def build_tr_sidecar(
    vcf_path,
    catalog_path,
    out_dir,
    *,
    with_entries=False,
    batch_size=50_000,
):
    """Build TR sidecar Parquet files from a VCF/BCF and TR catalog."""
    from storm.tr_sidecar import build_tr_sidecar as _build_tr_sidecar

    return _build_tr_sidecar(
        vcf_path,
        catalog_path,
        out_dir,
        with_entries=with_entries,
        batch_size=batch_size,
    )


def build_feature_inventory(
    mt,
    *,
    tr_field="tr",
    standard_label="standard_sv",
    tr_quant_label="tr_quantitative",
):
    """Build a feature-level SAIGE inventory from an annotated MatrixTable."""
    from storm.saige_prep import build_feature_inventory as _build_feature_inventory

    return _build_feature_inventory(
        mt,
        tr_field=tr_field,
        standard_label=standard_label,
        tr_quant_label=tr_quant_label,
    )


def annotate_saige_predictors(
    mt,
    *,
    tr_field="tr",
    standard_label="standard_sv",
    tr_quant_label="tr_quantitative",
):
    """Annotate entry-level SAIGE predictors on a MatrixTable."""
    from storm.saige_prep import annotate_saige_predictors as _annotate_saige_predictors

    return _annotate_saige_predictors(
        mt,
        tr_field=tr_field,
        standard_label=standard_label,
        tr_quant_label=tr_quant_label,
    )


def build_long_predictor_tables(
    mt,
    *,
    sample_field="s",
    tr_field="tr",
    standard_label="standard_sv",
    tr_quant_label="tr_quantitative",
):
    """Build long-form sample-feature predictor tables for SAIGE staging."""
    from storm.saige_prep import build_long_predictor_tables as _build_long_predictor_tables

    return _build_long_predictor_tables(
        mt,
        sample_field=sample_field,
        tr_field=tr_field,
        standard_label=standard_label,
        tr_quant_label=tr_quant_label,
    )


def print_feature_inventory_stats(
    feature_inventory,
    *,
    tr_quant_label="tr_quantitative",
    top_n=20,
):
    """Print common feature-inventory summary tables for SAIGE prep."""
    from storm.saige_prep import print_feature_inventory_stats as _print_feature_inventory_stats

    return _print_feature_inventory_stats(
        feature_inventory,
        tr_quant_label=tr_quant_label,
        top_n=top_n,
    )


def print_long_predictor_stats(
    standard_long,
    tr_quant_long,
    *,
    top_n=20,
):
    """Print common long-predictor summary tables for SAIGE prep."""
    from storm.saige_prep import print_long_predictor_stats as _print_long_predictor_stats

    return _print_long_predictor_stats(
        standard_long,
        tr_quant_long,
        top_n=top_n,
    )


def print_tr_annotation_summary(
    mt,
    *,
    tr_field="tr",
    run_detailed_breakdowns=False,
    top_n=20,
):
    """Print TR annotation summary with optional detailed breakdowns."""
    from storm.saige_prep import print_tr_annotation_summary as _print_tr_annotation_summary

    return _print_tr_annotation_summary(
        mt,
        tr_field=tr_field,
        run_detailed_breakdowns=run_detailed_breakdowns,
        top_n=top_n,
    )


def version() -> str:
    """Get the current version of the Storm package.

    Returns:
        str: The version string of the currently installed Storm package.
    """
    if _storm_ext is not None:
        return str(_storm_ext._version())
    try:
        return importlib.metadata.version("storm")
    except importlib.metadata.PackageNotFoundError:
        return "0.0.0"
