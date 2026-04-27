"""Frozen sidecar schema version identifiers and Parquet column documentation.

Site table (one row per variant) — ``sites.parquet``:

- ``schema_version`` (string): format of this row; bump when columns change.
- ``ruleset_version`` (string): TR rules / catalog generation version.
- ``variant_id`` (string): ``INFO/ID`` when set (must match Hail ``allele_id``), else the
  VCF third-column ID, else synthetic ``chrom:pos:ref:alt``.
- ``contig`` (string): chromosome with ``chr`` prefix when possible.
- ``pos_start`` / ``pos_end`` (int32): 0-based half-open span (pysam coordinates).
- ``ref`` (string), ``alt`` (string): first alternate allele only.
- ``svtype``, ``svlen``, ``end`` (optional): from INFO when present.
- ``motif``, ``motif_size``, ``detected``, ``has_tr_flag``: TR annotations from INFO.
- Catalog fields: ``tr_locus_id``, ``gene``, ``motif_primary``, ``motifs_accepted``,
  ``rule_group``, ``inheritance``, threshold strings, flags.
- ``catalog_overlap_bp`` (int32): bases overlapping chosen catalog interval.
- ``match_method`` (string): e.g. ``best_overlap``.
- ``motif_matches_catalog`` (bool): whether INFO motif matches catalog (nullable).
- ``repeat_units_estimate`` (float64): ``|SVLEN|/MOTIF_SIZE`` when computable (nullable).
- ``rule_applicable`` (bool): whether automated rule application is meaningful for v1.

Entry table (optional) — ``entries.parquet``:

- ``schema_version`` (use :data:`ENTRY_SCHEMA_VERSION`), ``ruleset_version``, ``variant_id``, ``sample``
- ``gt_str`` (string): diploid GT summary.
- ``is_carrier`` (bool): any non-reference allele.
- ``num_alt`` (int32): count of non-ref alleles in GT.
- ``format_json`` (string): JSON object of requested per-sample FORMAT fields (expansion callers),
  e.g. ``CN`` / repeat intervals; ``"{}"`` when none were requested or missing.
"""

from __future__ import annotations

import pyarrow as pa

# Per-row version written in ``sites.parquet`` rows.
SCHEMA_VERSION = "1"
# Per-row version written in ``entries.parquet`` rows (bumped when entry columns change).
ENTRY_SCHEMA_VERSION = "2"
RULESET_VERSION = "1"

SITE_ARROW_SCHEMA = pa.schema(
    [
        ("schema_version", pa.string()),
        ("ruleset_version", pa.string()),
        ("variant_id", pa.string()),
        ("contig", pa.string()),
        ("pos_start", pa.int32()),
        ("pos_end", pa.int32()),
        ("ref", pa.string()),
        ("alt", pa.string()),
        ("svtype", pa.string()),
        ("svlen", pa.int32()),
        ("end", pa.int32()),
        ("motif", pa.string()),
        ("motif_size", pa.int32()),
        ("detected", pa.string()),
        ("has_tr_flag", pa.bool_()),
        ("tr_locus_id", pa.string()),
        ("gene", pa.string()),
        ("motif_primary", pa.string()),
        ("motifs_accepted", pa.string()),
        ("rule_group", pa.string()),
        ("inheritance", pa.string()),
        ("benign_max_repeats", pa.string()),
        ("pathogenic_min_repeats", pa.string()),
        ("motif_change_required", pa.bool_()),
        ("contraction_pathogenic", pa.bool_()),
        ("prioritize_in_sv_scan", pa.bool_()),
        ("sv_scan_default_applicable", pa.bool_()),
        ("catalog_notes", pa.string()),
        ("catalog_overlap_bp", pa.int32()),
        ("match_method", pa.string()),
        ("motif_matches_catalog", pa.bool_()),
        ("repeat_units_estimate", pa.float64()),
        ("rule_applicable", pa.bool_()),
    ]
)

ENTRY_ARROW_SCHEMA = pa.schema(
    [
        ("schema_version", pa.string()),
        ("ruleset_version", pa.string()),
        ("variant_id", pa.string()),
        ("sample", pa.string()),
        ("gt_str", pa.string()),
        ("is_carrier", pa.bool_()),
        ("num_alt", pa.int32()),
        ("format_json", pa.string()),
    ]
)
