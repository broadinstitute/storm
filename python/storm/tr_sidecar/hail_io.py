"""Load TR sidecar ``sites.parquet`` into Hail without reading the whole file into memory.

Terra / Hail 0.2.134 often lacks ``hl.import_parquet``; the fallback must not call
``pq.read_table(...).to_pylist()`` on whole-genome site tables (driver OOM).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any


def sites_table_row_type(hl: Any) -> Any:
    """Row type for :func:`read_tr_sidecar_sites_table` (matches ``sites.parquet`` columns)."""
    return hl.tstruct(
        schema_version=hl.tstr,
        ruleset_version=hl.tstr,
        variant_id=hl.tstr,
        contig=hl.tstr,
        pos_start=hl.tint32,
        pos_end=hl.tint32,
        ref=hl.tstr,
        alt=hl.tstr,
        svtype=hl.tstr,
        svlen=hl.tint32,
        end=hl.tint32,
        motif=hl.tstr,
        motif_size=hl.tint32,
        detected=hl.tstr,
        has_tr_flag=hl.tbool,
        tr_locus_id=hl.tstr,
        gene=hl.tstr,
        motif_primary=hl.tstr,
        motifs_accepted=hl.tstr,
        rule_group=hl.tstr,
        inheritance=hl.tstr,
        benign_max_repeats=hl.tstr,
        pathogenic_min_repeats=hl.tstr,
        motif_change_required=hl.tbool,
        contraction_pathogenic=hl.tbool,
        prioritize_in_sv_scan=hl.tbool,
        sv_scan_default_applicable=hl.tbool,
        catalog_notes=hl.tstr,
        catalog_overlap_bp=hl.tint32,
        match_method=hl.tstr,
        motif_matches_catalog=hl.tbool,
        repeat_units_estimate=hl.tfloat64,
        rule_applicable=hl.tbool,
    )


def read_tr_sidecar_sites_table(
    hl: Any,
    sites_path: str | Path,
    *,
    parquet_chunk_rows: int = 200_000,
) -> Any:
    """Return an unkeyed Hail Table of sidecar sites (key with ``variant_id`` yourself).

    Uses ``hl.import_parquet`` when available; otherwise streams Parquet row groups
    in chunks of at most ``parquet_chunk_rows`` and unions Hail tables.
    """
    path = Path(sites_path)
    if not path.is_file():
        raise FileNotFoundError(path)

    row_type = sites_table_row_type(hl)

    if hasattr(hl, "import_parquet"):
        return hl.import_parquet(str(path))

    import pyarrow.parquet as pq

    pf = pq.ParquetFile(str(path))
    chunks: list[Any] = []
    for batch in pf.iter_batches(batch_size=max(1, int(parquet_chunk_rows))):
        rows = batch.to_pylist()
        if not rows:
            continue
        chunks.append(hl.Table.parallelize(rows, schema=row_type))

    if not chunks:
        return hl.Table.parallelize([], schema=row_type)
    if len(chunks) == 1:
        return chunks[0]
    acc = chunks[0]
    for t in chunks[1:]:
        acc = acc.union(t)
    return acc
