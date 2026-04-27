"""Hail: split multi-allelic SVs, per-allele IDs, dosage/phased, optional TR sidecar join.

Install Hail per https://hail.is/docs/0.2/getting_started.html (often conda).
Optional pip extra: ``pip install 'storm[hail]'`` (may not match all Hail installs).
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    import hail as hl


def _require_hail() -> Any:
    try:
        import hail as hl
    except ImportError as e:
        raise ImportError(
            "annotate_svs requires hail. Install hail in your environment "
            "(see https://hail.is/docs/0.2/getting_started.html)."
        ) from e
    return hl


def annotate_svs(
    mt: "hl.MatrixTable",
    *,
    tr_sidecar_sites: str | Path | None = None,
    keep_star: bool = False,
    left_aligned: bool = False,
) -> "hl.MatrixTable":
    """Prepare an SV MatrixTable: min-rep, parse per-ALT IDs, split multi, join TR sidecar.

    Row fields added (after split): ``allele_id``, ``allele_ids``, ``ref_len``, ``alt_len``,
    ``svlen`` (uses ``INFO/SVLEN`` when defined, else ``len(ALT) - len(REF)``), ``is_indel``,
    ``id_parse_ok`` (from pre-split ID vs ALT count).

    Entry fields: ``dosage`` (``GT.n_alt_alleles()``), ``phased``.

    If ``tr_sidecar_sites`` is set, joins ``sites.parquet`` on ``allele_id`` = ``variant_id``
    into row struct ``tr``.

    Parameters
    ----------
    mt
        From ``hl.import_vcf(..., reference_genome='GRCh38')`` (or compatible).
    tr_sidecar_sites
        Path to ``sites.parquet`` from :func:`storm.tr_sidecar.build_tr_sidecar`.
    keep_star
        Passed to :func:`hail.split_multi_hts`.
    left_aligned
        Passed to :func:`hail.split_multi_hts`.
    """
    hl = _require_hail()

    mr0 = hl.min_rep(mt.locus, mt.alleles)
    mt1 = mt.key_rows_by(locus=mr0.locus, alleles=mr0.alleles)

    # Must stay in sync with pysam ``variant_id`` in tr_sidecar (``INFO/ID`` first,
    # then the VCF column-3 id available as ``rsid`` in Hail).
    id_info = mt1.info.ID
    id_from_info = hl.or_missing(
        hl.is_defined(id_info)
        & (hl.str(id_info) != "")
        & (hl.str(id_info) != "."),
        hl.str(id_info),
    )
    id_from_rsid = hl.or_missing(
        hl.is_defined(mt1.rsid)
        & (hl.str(mt1.rsid) != "")
        & (hl.str(mt1.rsid) != "."),
        hl.str(mt1.rsid),
    )
    raw_id = hl.coalesce(id_from_info, id_from_rsid, "")
    looks_json = raw_id.matches(r"^\s*\[")
    id_clean = hl.if_else(
        looks_json,
        raw_id.replace(r"^\s*\[", "")
        .replace(r"\]\s*$", "")
        .replace('"', "")
        .replace("'", ""),
        raw_id,
    )

    n_alt = hl.len(mt1.alleles) - 1
    comma_items = id_clean.split(",").map(lambda s: s.strip())
    id_groups = hl.or_missing(
        (n_alt > 0) & (hl.len(comma_items) == n_alt),
        comma_items.map(lambda x: hl.array([x])),
    )

    mt1 = mt1.annotate_rows(
        id_groups=id_groups,
        id_parse_ok=(n_alt > 0) & (hl.len(comma_items) == n_alt),
    )

    mt2 = hl.split_multi_hts(mt1, keep_star=keep_star, left_aligned=left_aligned)
    mr1 = hl.min_rep(mt2.locus, mt2.alleles)
    mt2 = mt2.key_rows_by(locus=mr1.locus, alleles=mr1.alleles)

    has_group = hl.is_defined(mt2.id_groups) & (mt2.a_index - 1 < hl.len(mt2.id_groups))
    ref_len = hl.len(mt2.alleles[0])
    alt_len = hl.len(mt2.alleles[1])
    svlen_seq = alt_len - ref_len

    # INFO/SVLEN can be scalar or array (e.g. Number=A in multiallelic records).
    # Use a dtype-aware expression so split rows pick the appropriate ALT value.
    svlen_from_info = None
    try:
        svlen_dtype = mt2.info.SVLEN.dtype
        if svlen_dtype in (hl.tint32, hl.tint64):
            svlen_from_info = hl.int32(mt2.info.SVLEN)
        elif svlen_dtype in (hl.tarray(hl.tint32), hl.tarray(hl.tint64)):
            svlen_from_info = hl.or_missing(
                hl.is_defined(mt2.info.SVLEN) & (mt2.a_index - 1 < hl.len(mt2.info.SVLEN)),
                hl.int32(mt2.info.SVLEN[mt2.a_index - 1]),
            )
    except Exception:
        svlen_from_info = None

    svlen = (
        hl.if_else(
            hl.is_defined(svlen_from_info),
            svlen_from_info,
            hl.int32(svlen_seq),
        )
        if svlen_from_info is not None
        else hl.int32(svlen_seq)
    )

    mt2 = mt2.annotate_rows(
        allele_ids=hl.or_missing(has_group, mt2.id_groups[mt2.a_index - 1]),
        allele_id=hl.or_missing(has_group, mt2.id_groups[mt2.a_index - 1][0]),
        ref_len=ref_len,
        alt_len=alt_len,
        svlen=svlen,
        is_indel=ref_len != alt_len,
    )

    gt_defined = hl.is_defined(mt2.GT)
    mt2 = mt2.annotate_entries(
        dosage=hl.or_missing(gt_defined, mt2.GT.n_alt_alleles()),
        phased=hl.or_missing(gt_defined, mt2.GT.phased),
    )

    if tr_sidecar_sites is not None:
        sites_path = str(tr_sidecar_sites)
        if sites_path.endswith(".ht"):
            sites = hl.read_table(sites_path)
        else:
            from storm.tr_sidecar.hail_io import read_tr_sidecar_sites_table

            sites = read_tr_sidecar_sites_table(hl, sites_path)
        sites = sites.key_by("variant_id")
        mt2 = mt2.annotate_rows(tr=sites[mt2.allele_id])

    return mt2
