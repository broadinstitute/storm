"""Hail helpers for preparing Storm SV/TR data for SAIGE workflows.

This module packages notebook-only logic into reusable functions so users can:

1. classify variants into SAIGE feature strata
2. project row-level repeat annotations into sample-level predictors
3. export long-form staging tables for downstream marker-file conversion
4. densify sparse long tables and write SAIGE Step 2 VCF inputs (FORMAT ``DS``)
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from storm.sv_annotate import _require_hail

if TYPE_CHECKING:
    import hail as hl


def _vcf_export_row_value_fields(mt: "hl.MatrixTable") -> dict[str, "hl.Expression"]:
    """Non-key row fields to keep for :func:`hail.export_vcf`.

    ``locus`` and ``alleles`` are MatrixTable row keys and must not be passed to
    :meth:`MatrixTable.select_rows`; they are preserved automatically.
    """
    hl = _require_hail()
    out: dict[str, hl.Expression] = {}
    if "rsid" in mt.row:
        out["rsid"] = mt.rsid
    else:
        out["rsid"] = hl.missing(hl.tstr)
    for name in ("qual", "filters", "info"):
        if name in mt.row:
            out[name] = getattr(mt, name)
    return out


def build_feature_inventory(
    mt: "hl.MatrixTable",
    *,
    tr_field: str = "tr",
    standard_label: str = "standard_sv",
    tr_quant_label: str = "tr_quantitative",
) -> "hl.Table":
    """Create a feature-level SAIGE inventory table from an annotated MatrixTable.

    Expects ``mt`` to contain split-allele rows (e.g. from :func:`storm.annotate_svs`)
    and optional row struct ``tr`` with Storm sidecar fields.
    """
    hl = _require_hail()

    rows = mt.rows()
    tr = getattr(rows, tr_field)
    fi0 = rows.select(
        feature_id=rows.allele_id,
        svtype=hl.if_else(hl.is_defined(tr), tr.svtype, ""),
        svlen=rows.svlen,
        tr_locus_id=hl.if_else(hl.is_defined(tr), tr.tr_locus_id, ""),
        rule_applicable=hl.if_else(hl.is_defined(tr), tr.rule_applicable, hl.missing(hl.tbool)),
        motif=hl.if_else(hl.is_defined(tr), tr.motif, ""),
        repeat_units_estimate=hl.if_else(
            hl.is_defined(tr),
            tr.repeat_units_estimate,
            hl.missing(hl.tfloat64),
        ),
    )
    fi1 = fi0.annotate(
        has_repeat_units=hl.is_defined(fi0.repeat_units_estimate)
        & hl.is_finite(fi0.repeat_units_estimate),
    )
    return fi1.annotate(
        feature_class=hl.if_else(fi1.has_repeat_units, tr_quant_label, standard_label),
    )


def annotate_saige_predictors(
    mt: "hl.MatrixTable",
    *,
    tr_field: str = "tr",
    standard_label: str = "standard_sv",
    tr_quant_label: str = "tr_quantitative",
) -> "hl.MatrixTable":
    """Annotate entry-level predictors for SAIGE staging.

    Entry annotations added:
    - ``predictor_standard``: dosage for ``standard_sv`` rows
    - ``predictor_tr_quant``: ``dosage * repeat_units_estimate`` for ``tr_quantitative`` rows
    """
    hl = _require_hail()

    mt1 = mt.annotate_rows(
        has_repeat_units=hl.is_defined(getattr(mt, tr_field))
        & hl.is_defined(getattr(mt, tr_field).repeat_units_estimate)
        & hl.is_finite(getattr(mt, tr_field).repeat_units_estimate),
    )
    mt1 = mt1.annotate_rows(
        feature_class=hl.if_else(mt1.has_repeat_units, tr_quant_label, standard_label),
    )
    tr = getattr(mt1, tr_field)

    dosage = hl.if_else(hl.is_defined(mt1.dosage), mt1.dosage, mt1.GT.n_alt_alleles())
    return mt1.annotate_entries(
        predictor_standard=hl.or_missing(
            mt1.feature_class == standard_label,
            hl.float64(dosage),
        ),
        predictor_tr_quant=hl.or_missing(
            mt1.feature_class == tr_quant_label,
            hl.float64(dosage) * tr.repeat_units_estimate,
        ),
    )


def build_long_predictor_tables(
    mt: "hl.MatrixTable",
    *,
    sample_field: str = "s",
    tr_field: str = "tr",
    standard_label: str = "standard_sv",
    tr_quant_label: str = "tr_quantitative",
) -> dict[str, "hl.Table"]:
    """Build long-form sample-feature predictor tables for SAIGE staging.

    Returns a dict with:
    - ``standard_long``: sample-level dosage predictors
    - ``tr_quant_long``: dosage-weighted repeat predictors and TR metadata
    """
    hl = _require_hail()
    mt1 = annotate_saige_predictors(
        mt,
        tr_field=tr_field,
        standard_label=standard_label,
        tr_quant_label=tr_quant_label,
    )
    et = mt1.entries()
    et_std = et.filter(hl.is_defined(et.predictor_standard))
    standard_long = et_std.select(
        sample_id=hl.str(getattr(et_std, sample_field)),
        feature_id=et_std.allele_id,
        feature_class=et_std.feature_class,
        predictor=et_std.predictor_standard,
    )

    et_tr = et.filter(hl.is_defined(et.predictor_tr_quant))
    tr = getattr(et_tr, tr_field)
    tr_quant_long = et_tr.select(
        sample_id=hl.str(getattr(et_tr, sample_field)),
        feature_id=et_tr.allele_id,
        feature_class=et_tr.feature_class,
        predictor=et_tr.predictor_tr_quant,
        repeat_units_estimate=tr.repeat_units_estimate,
        tr_locus_id=tr.tr_locus_id,
        rule_applicable=tr.rule_applicable,
    )
    return {"standard_long": standard_long, "tr_quant_long": tr_quant_long}


def print_feature_inventory_stats(
    feature_inventory: "hl.Table",
    *,
    tr_quant_label: str = "tr_quantitative",
    top_n: int = 20,
) -> None:
    """Print common summary tables for SAIGE feature inventory."""
    hl = _require_hail()
    print("Feature strata counts")
    feature_inventory.group_by(feature_class=feature_inventory.feature_class).aggregate(
        n=hl.agg.count()
    ).order_by(hl.desc("n")).show(10)

    print("TR quantitative features: catalog metadata snapshot")
    fi_tr = feature_inventory.filter(feature_inventory.feature_class == tr_quant_label)
    fi_tr.group_by(
        has_catalog_locus=fi_tr.tr_locus_id != "",
        rule_applicable=fi_tr.rule_applicable,
    ).aggregate(n=hl.agg.count()).order_by(hl.desc("n")).show(top_n)


def print_long_predictor_stats(
    standard_long: "hl.Table",
    tr_quant_long: "hl.Table",
    *,
    top_n: int = 20,
) -> None:
    """Print common summary outputs for long-form SAIGE predictor tables."""
    hl = _require_hail()
    print("Entry-level staging counts")
    print("standard_sv non-missing predictors:", standard_long.count())
    print("tr_quantitative non-missing predictors:", tr_quant_long.count())

    print(f"\nPer-feature sample counts (tr_quantitative top {top_n})")
    tr_quant_long.group_by(feature_id=tr_quant_long.feature_id).aggregate(
        n_samples=hl.agg.count(),
        n_nonzero=hl.agg.count_where(tr_quant_long.predictor != 0),
    ).order_by(hl.desc("n_nonzero"), hl.desc("n_samples")).show(top_n)


def build_predictor_feature_qc(
    predictor_long: "hl.Table",
    *,
    feature_id_field: str = "feature_id",
    predictor_field: str = "predictor",
) -> "hl.Table":
    """Aggregate feature-level QC metrics from a long-form predictor table."""
    hl = _require_hail()
    fid = getattr(predictor_long, feature_id_field)
    pred = getattr(predictor_long, predictor_field)
    return predictor_long.group_by(feature_id=fid).aggregate(
        n_samples=hl.agg.count(),
        n_nonzero=hl.agg.count_where(pred != 0),
        predictor_mean=hl.agg.mean(pred),
        predictor_stdev=hl.agg.stats(pred).stdev,
        predictor_min=hl.agg.min(pred),
        predictor_max=hl.agg.max(pred),
    )


def export_long_predictor_tables(
    standard_long: "hl.Table",
    tr_quant_long: "hl.Table",
    *,
    out_dir: str | Path,
    prefix: str = "saige",
    include_qc: bool = True,
    include_sample_manifest: bool = True,
    carriers_only: bool = False,
    output_format: str = "tsv",
    sample_manifest: "hl.Table | None" = None,
    sample_id_field: str = "sample_id",
) -> dict[str, str]:
    """Export SAIGE staging tables and optional QC/manifest artifacts.

    Returns a dict of emitted artifact paths.
    """
    hl = _require_hail()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    if output_format not in {"tsv", "parquet", "ht"}:
        raise ValueError("output_format must be one of: 'tsv', 'parquet', 'ht'")

    std_to_write = standard_long
    tr_to_write = tr_quant_long
    if carriers_only:
        std_to_write = std_to_write.filter(std_to_write.predictor != 0)
        tr_to_write = tr_to_write.filter(tr_to_write.predictor != 0)
        # Sparse export needs sample universe for correct zero-filling downstream.
        include_sample_manifest = True

    if output_format == "tsv":
        standard_path = str(out / f"{prefix}.standard_long.tsv.bgz")
        tr_quant_path = str(out / f"{prefix}.tr_quant_long.tsv.bgz")
        std_to_write.export(standard_path)
        tr_to_write.export(tr_quant_path)
    elif output_format == "parquet":
        standard_path = str(out / f"{prefix}.standard_long.parquet")
        tr_quant_path = str(out / f"{prefix}.tr_quant_long.parquet")
        try:
            std_to_write.to_spark(flatten=True).write.mode("overwrite").parquet(standard_path)
            tr_to_write.to_spark(flatten=True).write.mode("overwrite").parquet(tr_quant_path)
        except NotImplementedError:
            # LocalBackend does not support to_spark(); stage as native Hail tables.
            output_format = "ht"
            standard_path = str(out / f"{prefix}.standard_long.ht")
            tr_quant_path = str(out / f"{prefix}.tr_quant_long.ht")
            std_to_write.write(standard_path, overwrite=True)
            tr_to_write.write(tr_quant_path, overwrite=True)
    else:
        standard_path = str(out / f"{prefix}.standard_long.ht")
        tr_quant_path = str(out / f"{prefix}.tr_quant_long.ht")
        std_to_write.write(standard_path, overwrite=True)
        tr_to_write.write(tr_quant_path, overwrite=True)

    emitted = {
        "standard_long": standard_path,
        "tr_quant_long": tr_quant_path,
    }

    if include_qc:
        std_qc = build_predictor_feature_qc(std_to_write).order_by(
            hl.desc("n_nonzero"), hl.desc("n_samples")
        )
        tr_qc = build_predictor_feature_qc(tr_to_write).order_by(
            hl.desc("n_nonzero"), hl.desc("n_samples")
        )
        std_qc_path = str(out / f"{prefix}.standard_feature_qc.tsv.bgz")
        tr_qc_path = str(out / f"{prefix}.tr_quant_feature_qc.tsv.bgz")
        std_qc.export(std_qc_path)
        tr_qc.export(tr_qc_path)
        emitted["standard_feature_qc"] = std_qc_path
        emitted["tr_quant_feature_qc"] = tr_qc_path

    if include_sample_manifest:
        if sample_manifest is None:
            # Fallback (potentially expensive): derive from long tables.
            s_std = standard_long.select("sample_id")
            s_tr = tr_quant_long.select("sample_id")
            sample_manifest_ht = s_std.union(s_tr).distinct().order_by("sample_id")
        else:
            sample_manifest_ht = sample_manifest.select(
                sample_id=hl.str(getattr(sample_manifest, sample_id_field))
            ).distinct().order_by("sample_id")
        sample_manifest_path = str(out / f"{prefix}.sample_manifest.tsv")
        sample_manifest_ht.export(sample_manifest_path)
        emitted["sample_manifest"] = sample_manifest_path

    emitted["carriers_only"] = str(bool(carriers_only))
    emitted["output_format"] = output_format

    return emitted


def build_feature_vcf_row_lookup(
    mt: "hl.MatrixTable",
    *,
    feature_field: str = "allele_id",
) -> "hl.Table":
    """Map ``feature_field`` (default ``allele_id``) to VCF row fields for densification.

    Use with sparse long predictor tables and :func:`dense_marker_matrix_from_long`.
    """
    hl = _require_hail()
    rows = mt.rows()
    fid = getattr(rows, feature_field)
    sel = rows.select(
        feature_id=fid,
        locus=rows.locus,
        alleles=rows.alleles,
    )
    if "rsid" in rows:
        sel = sel.annotate(rsid=rows.rsid)
    else:
        sel = sel.annotate(rsid=hl.missing(hl.tstr))
    return sel.key_by("feature_id")


def align_matrix_cols_to_manifest(
    mt: "hl.MatrixTable",
    manifest_ht: "hl.Table",
    *,
    sample_id_field: str = "sample_id",
    col_key: str = "s",
) -> "hl.MatrixTable":
    """Reorder columns to match ``manifest_ht`` (ordered by ``sample_id_field``).

    Notes
    -----
    Uses ``MatrixTable.choose_cols``, which requires the full column order in memory
    on the driver. Fine for biobank-scale *sample* counts in typical Hail setups, but
    avoid calling in tight loops on huge clusters without revisiting.
    """
    if col_key not in mt.col_key:
        raise ValueError(f"expected column key field {col_key!r}, found {list(mt.col_key)}")
    sid_expr = getattr(manifest_ht, sample_id_field)
    # Hail Table.distinct() requires a key; manifests from mt.cols().select(...) are often unkeyed.
    ordered_tbl = (
        manifest_ht.select(_manifest_sid=sid_expr).key_by("_manifest_sid").distinct()
    )
    ordered = ordered_tbl.order_by("_manifest_sid")["_manifest_sid"].collect()
    present = getattr(mt, col_key).collect()
    idx_map = {sid: i for i, sid in enumerate(present)}
    missing = [sid for sid in ordered if sid not in idx_map]
    if missing:
        raise ValueError(
            "manifest contains sample_ids not present in matrix columns "
            f"(showing up to 5): {missing[:5]!r}"
        )
    return mt.choose_cols([idx_map[sid] for sid in ordered])


def build_dense_saige_marker_matrix(
    mt: "hl.MatrixTable",
    *,
    stratum: str,
    tr_field: str = "tr",
    standard_label: str = "standard_sv",
    tr_quant_label: str = "tr_quantitative",
    sample_field: str = "s",
    fill_value: float = 0.0,
) -> "hl.MatrixTable":
    """Dense per-stratum matrix with entries ``DS`` for SAIGE ``--vcfField=DS``.

    Keeps only rows whose ``feature_class`` matches ``stratum``, then fills missing
    predictors with ``fill_value`` (non-carriers and filtered genotypes).
    """
    hl = _require_hail()
    if stratum not in {standard_label, tr_quant_label}:
        raise ValueError(f"stratum must be {standard_label!r} or {tr_quant_label!r}, got {stratum!r}")

    mt1 = annotate_saige_predictors(
        mt,
        tr_field=tr_field,
        standard_label=standard_label,
        tr_quant_label=tr_quant_label,
    )
    mt1 = mt1.filter_rows(mt1.feature_class == stratum)
    pred = (
        mt1.predictor_standard if stratum == standard_label else mt1.predictor_tr_quant
    )
    mt1 = mt1.select_entries(DS=hl.float64(hl.or_else(pred, hl.float64(fill_value))))
    s = getattr(mt1, sample_field)
    mt1 = mt1.key_cols_by(s=hl.str(s))
    return mt1.select_rows(**_vcf_export_row_value_fields(mt1))


def dense_marker_matrix_from_long(
    predictor_long: "hl.Table",
    row_lookup: "hl.Table",
    *,
    feature_id_field: str = "feature_id",
    sample_id_field: str = "sample_id",
    predictor_field: str = "predictor",
    fill_value: float = 0.0,
) -> "hl.MatrixTable":
    """Build a dense ``DS`` matrix from a sparse long table + feature row lookup.

    ``row_lookup`` must be keyed by ``feature_id`` (see :func:`build_feature_vcf_row_lookup`).
    """
    hl = _require_hail()
    pl = predictor_long
    fid = getattr(pl, feature_id_field)
    sid = getattr(pl, sample_id_field)
    pred = getattr(pl, predictor_field)
    rk = row_lookup[fid]
    joined = pl.filter(hl.is_defined(rk.locus)).annotate(
        locus=rk.locus,
        alleles=rk.alleles,
        rsid=rk.rsid,
    )
    mt = joined.select(
        locus=joined.locus,
        alleles=joined.alleles,
        rsid=joined.rsid,
        sample_id=sid,
        feature_id=fid,
        predictor=pred,
    ).to_matrix_table(
        row_key=["feature_id"],
        col_key=["sample_id"],
        row_fields=["locus", "alleles", "rsid"],
    )
    mt = mt.annotate_entries(DS=hl.float64(hl.or_else(mt.predictor, hl.float64(fill_value))))
    mt = mt.drop("predictor")
    mt = mt.key_rows_by(locus=mt.locus, alleles=mt.alleles)
    mt = mt.key_cols_by(s=hl.str(mt.sample_id))
    drop_row: list[str] = [n for n in ("feature_id", "qual", "filters", "info") if n in mt.row]
    if drop_row:
        mt = mt.drop(*drop_row)
    return mt


def saige_dosage_vcf_metadata() -> dict[str, dict[str, dict[str, str]]]:
    """VCF header fragment for FORMAT/DS (pass to :func:`hail.export_vcf` ``metadata``)."""
    return {
        "format": {
            "DS": {
                "Description": "Storm SAIGE predictor (dosage or dosage-weighted TR proxy)",
                "Number": "1",
                "Type": "Float",
            },
        },
    }


def export_saige_dosage_vcf(
    mt: "hl.MatrixTable",
    output_path: str | Path,
    *,
    parallel: str | None = None,
    tabix: bool = False,
) -> str:
    """Write a MatrixTable with entry field ``DS`` to a blocked gzip VCF.

    SAIGE Step 2: pass ``--vcfFile=...``, ``--vcfFileIndex=...``, ``--vcfField=DS``.
    """
    hl = _require_hail()
    path = str(output_path)
    if not path.endswith(".bgz"):
        raise ValueError("output_path should end with .vcf.bgz for blocked gzip (tabix / SAIGE)")
    hl.export_vcf(
        mt,
        path,
        metadata=saige_dosage_vcf_metadata(),
        parallel=parallel,
        tabix=tabix,
    )
    return path


def export_saige_stratum_vcfs(
    mt: "hl.MatrixTable",
    manifest_ht: "hl.Table",
    *,
    out_dir: str | Path,
    prefix: str = "saige",
    sample_id_field: str = "sample_id",
    tr_field: str = "tr",
    standard_label: str = "standard_sv",
    tr_quant_label: str = "tr_quantitative",
    parallel: str | None = None,
    tabix: bool = False,
) -> dict[str, str]:
    """Export both SAIGE strata as dosage VCFs with identical column order.

    Column order follows ``manifest_ht`` ordered by ``sample_id_field``.
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    m_std = build_dense_saige_marker_matrix(
        mt,
        stratum=standard_label,
        tr_field=tr_field,
        standard_label=standard_label,
        tr_quant_label=tr_quant_label,
    )
    m_tr = build_dense_saige_marker_matrix(
        mt,
        stratum=tr_quant_label,
        tr_field=tr_field,
        standard_label=standard_label,
        tr_quant_label=tr_quant_label,
    )
    m_std = align_matrix_cols_to_manifest(m_std, manifest_ht, sample_id_field=sample_id_field)
    m_tr = align_matrix_cols_to_manifest(m_tr, manifest_ht, sample_id_field=sample_id_field)
    p_std = str(out / f"{prefix}.{standard_label}.ds.vcf.bgz")
    p_tr = str(out / f"{prefix}.{tr_quant_label}.ds.vcf.bgz")
    export_saige_dosage_vcf(m_std, p_std, parallel=parallel, tabix=tabix)
    export_saige_dosage_vcf(m_tr, p_tr, parallel=parallel, tabix=tabix)
    return {
        standard_label: p_std,
        tr_quant_label: p_tr,
    }


def saige_synthetic_covar_col_list(*, n_pcs: int = 5) -> str:
    """Comma-separated covariate column names for ``--covarColList`` (synthetic table)."""
    if n_pcs < 0 or n_pcs > 50:
        raise ValueError("n_pcs must be between 0 and 50")
    cols = ["sex", "age", "seq_depth_log"]
    cols += [f"PC{i}" for i in range(1, n_pcs + 1)]
    return ",".join(cols)


def build_synthetic_saige_pheno_covar_table(
    manifest_ht: "hl.Table",
    *,
    sample_id_field: str = "sample_id",
    iid_output_field: str = "IID",
    n_pcs: int = 5,
) -> "hl.Table":
    """Build a **placeholder** phenotype + covariate table keyed to a sample manifest.

    Values are **deterministic** from ``sample_id`` (via Hail ``hash``), so the same
    cohort always gets the same synthetic fields. Replace with real EHR / QC tables
    for production.

    Columns:

    - ``IID``: copy of sample id (use with ``--sampleIDColinphenoFile=IID``).
    - ``pheno_binary``: 0/1 (roughly balanced).
    - ``pheno_quant``: float (unrelated toy outcome).
    - ``sex``, ``age``, ``seq_depth_log``, ``PC1`` … ``PC{n_pcs}``: toy covariates.

    SAIGE Step 1 examples: ``--phenoCol=pheno_binary`` or ``pheno_quant``,
    ``--covarColList`` from :func:`saige_synthetic_covar_col_list`.
    """
    if n_pcs < 0 or n_pcs > 50:
        raise ValueError("n_pcs must be between 0 and 50")
    hl = _require_hail()

    sid = getattr(manifest_ht, sample_id_field)
    t0 = manifest_ht.select(_sid=sid)
    t1 = t0.annotate(
        _h1=hl.abs(hl.int64(hl.hash(t0._sid))),
        _h2=hl.abs(hl.int64(hl.hash(hl.str(t0._sid) + hl.literal("|storm_syn|")))),
    )
    sex = hl.int32(t1._h1 % 2)
    age = hl.int32(30 + (t1._h1 % 50))
    seq_depth_log = hl.log10(10.0 + hl.float64(t1._h2 % 90))
    pheno_binary = hl.int32(((t1._h1 + t1._h2) % 100) < 42)
    pheno_quant = hl.float64(t1._h1 % 100) * 0.1 + hl.float64(t1._h2 % 300) * 0.01

    pc_fields: dict[str, hl.Expression] = {}
    for k in range(1, n_pcs + 1):
        pc_fields[f"PC{k}"] = hl.float64(
            (t1._h1 // hl.int64(7 * k) + t1._h2 * hl.int64(k + 1)) % 2000
        ) / 1000.0 - 1.0

    return t1.select(
        **{
            iid_output_field: t1._sid,
            "pheno_binary": pheno_binary,
            "pheno_quant": pheno_quant,
            "sex": sex,
            "age": age,
            "seq_depth_log": seq_depth_log,
            **pc_fields,
        }
    )


def export_saige_phenotype_covariate_tsv(
    pheno_covar_ht: "hl.Table",
    path: str | Path,
    *,
    sort_by_iid: str | None = "IID",
) -> str:
    """Write a tab-delimited phenotype/covariate file for SAIGE Step 1.

    Parameters
    ----------
    sort_by_iid
        If set, column name to ``order_by`` before export (stable cohort ordering).
        Pass ``None`` to preserve table order.
    """
    _require_hail()
    out = str(path)
    ht = pheno_covar_ht
    if sort_by_iid is not None:
        if sort_by_iid not in ht.row:
            raise ValueError(f"sort_by_iid={sort_by_iid!r} not in table row fields")
        ht = ht.order_by(sort_by_iid)
    ht.export(out, delimiter="\t")
    return out


def print_tr_annotation_summary(
    mt: "hl.MatrixTable",
    *,
    tr_field: str = "tr",
    run_detailed_breakdowns: bool = False,
    top_n: int = 20,
) -> dict[str, int]:
    """Print a fast TR annotation summary, with optional detailed breakdowns.

    Returns a dict-like struct of core counts for programmatic reuse.
    """
    hl = _require_hail()
    rows = mt.rows()
    tr = getattr(rows, tr_field)
    rows = rows.select(tr=tr)
    tr = rows.tr

    core = rows.aggregate(
        hl.struct(
            rows_total=hl.agg.count(),
            rows_with_tr_struct=hl.agg.count_where(hl.is_defined(tr)),
            rows_with_catalog_locus=hl.agg.count_where(hl.is_defined(tr) & (tr.tr_locus_id != "")),
            rows_match_method_best_overlap=hl.agg.count_where(
                hl.is_defined(tr) & (tr.match_method == "best_overlap")
            ),
            rows_match_method_no_overlap=hl.agg.count_where(
                hl.is_defined(tr) & (tr.match_method == "no_overlap")
            ),
            rows_with_motif=hl.agg.count_where(hl.is_defined(tr) & (tr.motif != "")),
            rows_with_motif_matches_catalog_true=hl.agg.count_where(
                hl.is_defined(tr)
                & hl.is_defined(tr.motif_matches_catalog)
                & tr.motif_matches_catalog
            ),
            rows_with_motif_matches_catalog_false=hl.agg.count_where(
                hl.is_defined(tr)
                & hl.is_defined(tr.motif_matches_catalog)
                & (~tr.motif_matches_catalog)
            ),
            rows_with_repeat_units_estimate=hl.agg.count_where(
                hl.is_defined(tr)
                & hl.is_defined(tr.repeat_units_estimate)
                & hl.is_finite(tr.repeat_units_estimate)
            ),
            rows_rule_applicable_true=hl.agg.count_where(
                hl.is_defined(tr) & hl.is_defined(tr.rule_applicable) & tr.rule_applicable
            ),
            rows_rule_applicable_false=hl.agg.count_where(
                hl.is_defined(tr) & hl.is_defined(tr.rule_applicable) & (~tr.rule_applicable)
            ),
        )
    )

    print("Core counts (single-pass aggregate)")
    for k, v in core.items():
        print(f"{k}: {v}")
    if core.rows_total:
        print("\nSelected percentages")
        print(f"catalog_locus_pct: {100.0 * core.rows_with_catalog_locus / core.rows_total:.3f}%")
        print(
            f"repeat_units_estimate_pct: {100.0 * core.rows_with_repeat_units_estimate / core.rows_total:.3f}%"
        )
        print(
            f"rule_applicable_true_pct: {100.0 * core.rows_rule_applicable_true / core.rows_total:.3f}%"
        )

    if run_detailed_breakdowns:
        print("\nBy match_method")
        rows.group_by(
            match_method=hl.if_else(hl.is_defined(tr), tr.match_method, "tr_missing")
        ).aggregate(n=hl.agg.count()).order_by(hl.desc("n")).show(top_n)

        print("\nBy svtype")
        rows.group_by(
            svtype=hl.if_else(hl.is_defined(tr), tr.svtype, "tr_missing")
        ).aggregate(n=hl.agg.count()).order_by(hl.desc("n")).show(top_n)

        print("\nBy motif presence x motif_matches_catalog")
        rows.group_by(
            has_motif=hl.if_else(hl.is_defined(tr), tr.motif != "", False),
            motif_matches_catalog=hl.if_else(
                hl.is_defined(tr), tr.motif_matches_catalog, hl.missing(hl.tbool)
            ),
        ).aggregate(n=hl.agg.count()).order_by(hl.desc("n")).show(top_n)

        print("\nBy rule_applicable")
        rows.group_by(
            rule_applicable=hl.if_else(
                hl.is_defined(tr), tr.rule_applicable, hl.missing(hl.tbool)
            )
        ).aggregate(n=hl.agg.count()).order_by(hl.desc("n")).show(10)

        print("\nTop genes among matched catalog loci")
        rows.filter(hl.is_defined(tr) & (tr.tr_locus_id != "")).group_by(
            gene=tr.gene
        ).aggregate(n=hl.agg.count()).order_by(hl.desc("n")).show(top_n)

    return dict(core.items())
