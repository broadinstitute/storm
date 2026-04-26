"""Hail helpers for preparing Storm SV/TR data for SAIGE workflows.

This module packages notebook-only logic into reusable functions so users can:

1. classify variants into SAIGE feature strata
2. project row-level repeat annotations into sample-level predictors
3. export long-form staging tables for downstream marker-file conversion
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from storm.sv_annotate import _require_hail

if TYPE_CHECKING:
    import hail as hl


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
