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
    "print_feature_inventory_stats",
    "print_long_predictor_stats",
    "print_tr_annotation_summary",
    "build_predictor_feature_qc",
    "export_long_predictor_tables",
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


def build_predictor_feature_qc(
    predictor_long,
    *,
    feature_id_field="feature_id",
    predictor_field="predictor",
):
    """Build feature-level QC summary from a long predictor table."""
    from storm.saige_prep import build_predictor_feature_qc as _build_predictor_feature_qc

    return _build_predictor_feature_qc(
        predictor_long,
        feature_id_field=feature_id_field,
        predictor_field=predictor_field,
    )


def export_long_predictor_tables(
    standard_long,
    tr_quant_long,
    *,
    out_dir,
    prefix="saige",
    include_qc=True,
    include_sample_manifest=True,
    carriers_only=False,
    output_format="tsv",
    sample_manifest=None,
    sample_id_field="sample_id",
):
    """Export SAIGE staging tables and optional QC/manifest artifacts."""
    from storm.saige_prep import export_long_predictor_tables as _export_long_predictor_tables

    return _export_long_predictor_tables(
        standard_long,
        tr_quant_long,
        out_dir=out_dir,
        prefix=prefix,
        include_qc=include_qc,
        include_sample_manifest=include_sample_manifest,
        carriers_only=carriers_only,
        output_format=output_format,
        sample_manifest=sample_manifest,
        sample_id_field=sample_id_field,
    )


def build_feature_vcf_row_lookup(mt, *, feature_field="allele_id"):
    """Table keyed by feature id with locus/alleles/rsid for SAIGE VCF densification."""
    from storm.saige_prep import build_feature_vcf_row_lookup as _fn

    return _fn(mt, feature_field=feature_field)


def align_matrix_cols_to_manifest(mt, manifest_ht, *, sample_id_field="sample_id", col_key="s"):
    """Reorder MatrixTable columns to match a sample manifest table."""
    from storm.saige_prep import align_matrix_cols_to_manifest as _fn

    return _fn(mt, manifest_ht, sample_id_field=sample_id_field, col_key=col_key)


def build_dense_saige_marker_matrix(
    mt,
    *,
    stratum,
    tr_field="tr",
    standard_label="standard_sv",
    tr_quant_label="tr_quantitative",
    sample_field="s",
    fill_value=0.0,
):
    """Dense DS MatrixTable for one SAIGE stratum."""
    from storm.saige_prep import build_dense_saige_marker_matrix as _fn

    return _fn(
        mt,
        stratum=stratum,
        tr_field=tr_field,
        standard_label=standard_label,
        tr_quant_label=tr_quant_label,
        sample_field=sample_field,
        fill_value=fill_value,
    )


def dense_marker_matrix_from_long(
    predictor_long,
    row_lookup,
    *,
    feature_id_field="feature_id",
    sample_id_field="sample_id",
    predictor_field="predictor",
    fill_value=0.0,
):
    """Dense DS matrix from a sparse long predictor table."""
    from storm.saige_prep import dense_marker_matrix_from_long as _fn

    return _fn(
        predictor_long,
        row_lookup,
        feature_id_field=feature_id_field,
        sample_id_field=sample_id_field,
        predictor_field=predictor_field,
        fill_value=fill_value,
    )


def saige_dosage_vcf_metadata():
    """VCF header metadata dict for FORMAT/DS (``hail.export_vcf``)."""
    from storm.saige_prep import saige_dosage_vcf_metadata as _fn

    return _fn()


def export_saige_dosage_vcf(mt, output_path, *, parallel=None, tabix=False):
    """Write ``DS``-only entries to ``.vcf.bgz`` for SAIGE Step 2."""
    from storm.saige_prep import export_saige_dosage_vcf as _fn

    return _fn(mt, output_path, parallel=parallel, tabix=tabix)


def export_saige_stratum_vcfs(
    mt,
    manifest_ht,
    *,
    out_dir,
    prefix="saige",
    sample_id_field="sample_id",
    tr_field="tr",
    standard_label="standard_sv",
    tr_quant_label="tr_quantitative",
    parallel=None,
    tabix=False,
):
    """Export both strata as dosage VCFs with manifest column order."""
    from storm.saige_prep import export_saige_stratum_vcfs as _fn

    return _fn(
        mt,
        manifest_ht,
        out_dir=out_dir,
        prefix=prefix,
        sample_id_field=sample_id_field,
        tr_field=tr_field,
        standard_label=standard_label,
        tr_quant_label=tr_quant_label,
        parallel=parallel,
        tabix=tabix,
    )


def saige_synthetic_covar_col_list(*, n_pcs=5):
    """Covariate column list string for SAIGE ``--covarColList`` (synthetic cohort)."""
    from storm.saige_prep import saige_synthetic_covar_col_list as _fn

    return _fn(n_pcs=n_pcs)


def build_synthetic_saige_pheno_covar_table(manifest_ht, *, sample_id_field="sample_id", iid_output_field="IID", n_pcs=5):
    """Toy phenotype + covariates aligned to a sample manifest (SAIGE Step 1)."""
    from storm.saige_prep import build_synthetic_saige_pheno_covar_table as _fn

    return _fn(
        manifest_ht,
        sample_id_field=sample_id_field,
        iid_output_field=iid_output_field,
        n_pcs=n_pcs,
    )


def export_saige_phenotype_covariate_tsv(pheno_covar_ht, path, *, sort_by_iid="IID"):
    """Export tab-delimited phenotype/covariate file for SAIGE Step 1."""
    from storm.saige_prep import export_saige_phenotype_covariate_tsv as _fn

    return _fn(pheno_covar_ht, path, sort_by_iid=sort_by_iid)


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
