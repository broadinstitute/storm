from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from .catalog import TrLocus

import pysam
import pyarrow as pa
import pyarrow.parquet as pq

from . import schema
from .catalog import load_catalog, norm_chrom
from .match import best_catalog_match, compute_rule_applicable, motif_matches


@dataclass
class SidecarResult:
    out_dir: Path
    sites_path: Path
    entries_path: Path | None
    manifest_path: Path
    n_variants: int


def _info_scalar(rec: pysam.VariantRecord, key: str) -> Any:
    if key not in rec.info:
        return None
    v = rec.info[key]
    if isinstance(v, tuple):
        return v[0] if len(v) == 1 else v
    return v


def _info_first_atomic(rec: pysam.VariantRecord, key: str) -> Any:
    """First INFO value when the field is scalar or multi-value (e.g. Number=A ``SVLEN``)."""
    if key not in rec.info:
        return None
    v = rec.info[key]
    if isinstance(v, tuple):
        return v[0] if len(v) else None
    if isinstance(v, list):
        return v[0] if len(v) else None
    return v


def _info_flag(rec: pysam.VariantRecord, key: str) -> bool:
    return key in rec.info


def _first_alt(rec: pysam.VariantRecord) -> str:
    if not rec.alts:
        return ""
    return str(rec.alts[0])


def variant_id_for_record(rec: pysam.VariantRecord) -> str:
    """ID used in ``sites.parquet``; must match Hail ``allele_id`` from :func:`storm.annotate_svs`.

    ``annotate_svs`` takes ``INFO/ID`` (when set), else the VCF's third-column ID
    (``mt.rsid`` in Hail). A bare ``pysam`` ``record.id`` alone can disagree with
    ``INFO/ID``; always prefer the INFO field to avoid join misses.
    """
    info_id = _info_scalar(rec, "ID")
    if info_id is not None:
        s = str(info_id).strip()
        if s and s not in (".", ""):
            # Hail: comma-list aligned with alts; first token is the first alternate.
            return s.split(",")[0].strip()
    rid = rec.id
    if rid and str(rid) not in (".", ""):
        s = str(rid).strip()
        return s.split(",")[0].strip()
    alt = _first_alt(rec)
    return f"{rec.chrom}:{rec.pos}:{rec.ref}:{alt}"


def _variant_span(rec: pysam.VariantRecord) -> tuple[int, int]:
    """0-based half-open span for catalog overlap.

    Prefer ``record.stop`` (htslib applies INFO END when the header declares it).
    Otherwise use INFO ``END``, then approximate with ``|SVLEN|`` for 1 bp symbolic rows.
    """
    vs, ve = rec.start, rec.stop
    if ve <= vs:
        ve = vs + 1
    end_v = _info_scalar(rec, "END")
    if end_v is not None:
        try:
            end_i = int(end_v)
            ve = max(ve, end_i)
        except (TypeError, ValueError):
            pass
    if ve <= vs + 1:
        svlen_v = _info_scalar(rec, "SVLEN")
        if svlen_v is not None:
            try:
                sl = abs(int(svlen_v))
                if sl > 0:
                    ve = max(ve, vs + sl)
            except (TypeError, ValueError):
                pass
    return vs, ve


def _repeat_units_estimate(rec: pysam.VariantRecord) -> float | None:
    svlen = _info_first_atomic(rec, "SVLEN")
    ms = _info_first_atomic(rec, "MOTIF_SIZE")
    if svlen is None or ms is None:
        return None
    try:
        sl = int(svlen)
        m = int(ms)
    except (TypeError, ValueError):
        return None
    if m <= 0:
        return None
    return abs(float(sl)) / float(m)


def _site_row(
    rec: pysam.VariantRecord,
    loc: TrLocus | None,
    overlap: int,
    match_method: str,
    motif_match: bool | None,
) -> dict[str, Any]:
    vid = variant_id_for_record(rec)
    contig = norm_chrom(rec.chrom)
    vs, ve = _variant_span(rec)

    svtype_v = _info_scalar(rec, "SVTYPE")
    svtype = str(svtype_v) if svtype_v is not None else ""

    svlen_v = _info_scalar(rec, "SVLEN")
    svlen: int | None
    try:
        svlen = int(svlen_v) if svlen_v is not None else None
    except (TypeError, ValueError):
        svlen = None

    end_v = _info_scalar(rec, "END")
    try:
        end_i = int(end_v) if end_v is not None else None
    except (TypeError, ValueError):
        end_i = None

    motif_v = _info_scalar(rec, "MOTIF")
    motif = str(motif_v) if motif_v is not None else ""

    ms_v = _info_scalar(rec, "MOTIF_SIZE")
    try:
        motif_size = int(ms_v) if ms_v is not None else None
    except (TypeError, ValueError):
        motif_size = None

    det_v = _info_scalar(rec, "DETECTED")
    detected = str(det_v) if det_v is not None else ""

    has_tr = _info_flag(rec, "TR")

    if loc is None:
        return {
            "schema_version": schema.SCHEMA_VERSION,
            "ruleset_version": schema.RULESET_VERSION,
            "variant_id": vid,
            "contig": contig,
            "pos_start": vs,
            "pos_end": ve,
            "ref": rec.ref,
            "alt": _first_alt(rec),
            "svtype": svtype,
            "svlen": svlen if svlen is not None else None,
            "end": end_i if end_i is not None else None,
            "motif": motif,
            "motif_size": motif_size if motif_size is not None else None,
            "detected": detected,
            "has_tr_flag": has_tr,
            "tr_locus_id": "",
            "gene": "",
            "motif_primary": "",
            "motifs_accepted": "",
            "rule_group": "",
            "inheritance": "",
            "benign_max_repeats": "",
            "pathogenic_min_repeats": "",
            "motif_change_required": False,
            "contraction_pathogenic": False,
            "prioritize_in_sv_scan": False,
            "sv_scan_default_applicable": False,
            "catalog_notes": "",
            "catalog_overlap_bp": overlap,
            "match_method": match_method,
            "motif_matches_catalog": motif_match,
            "repeat_units_estimate": _repeat_units_estimate(rec),
            "rule_applicable": False,
        }

    ru = _repeat_units_estimate(rec)
    mm = motif_match
    if motif and loc.motif_primary:
        mm = motif_matches(motif, loc.motif_primary, loc.motifs_accepted)

    return {
        "schema_version": schema.SCHEMA_VERSION,
        "ruleset_version": schema.RULESET_VERSION,
        "variant_id": vid,
        "contig": contig,
        "pos_start": vs,
        "pos_end": ve,
        "ref": rec.ref,
        "alt": _first_alt(rec),
        "svtype": svtype,
        "svlen": svlen if svlen is not None else None,
        "end": end_i if end_i is not None else None,
        "motif": motif,
        "motif_size": motif_size if motif_size is not None else None,
        "detected": detected,
        "has_tr_flag": has_tr,
        "tr_locus_id": loc.tr_locus_id,
        "gene": loc.gene,
        "motif_primary": loc.motif_primary,
        "motifs_accepted": loc.motifs_accepted,
        "rule_group": loc.rule_group,
        "inheritance": loc.inheritance,
        "benign_max_repeats": loc.benign_max_repeats,
        "pathogenic_min_repeats": loc.pathogenic_min_repeats,
        "motif_change_required": loc.motif_change_required,
        "contraction_pathogenic": loc.contraction_pathogenic,
        "prioritize_in_sv_scan": loc.prioritize_in_sv_scan,
        "sv_scan_default_applicable": loc.sv_scan_default_applicable,
        "catalog_notes": loc.notes,
        "catalog_overlap_bp": overlap,
        "match_method": match_method,
        "motif_matches_catalog": mm,
        "repeat_units_estimate": ru,
        "rule_applicable": compute_rule_applicable(loc, mm),
    }


def _gt_summary(rec: pysam.VariantRecord, sample: str) -> tuple[str, bool, int]:
    try:
        gt = rec.samples[sample]["GT"]
    except (KeyError, TypeError):
        return "./.", False, 0
    if gt is None:
        return "./.", False, 0
    gtl = list(gt)
    gt_str = "/".join("." if x is None else str(int(x)) for x in gtl)
    alts = sum(1 for x in gtl if x is not None and int(x) > 0)
    is_carrier = alts > 0
    return gt_str, is_carrier, alts


def build_tr_sidecar(
    vcf_path: str | Path,
    catalog_path: str | Path,
    out_dir: str | Path,
    *,
    with_entries: bool = False,
    batch_size: int = 50_000,
) -> SidecarResult:
    """Stream ``vcf_path``, join ``catalog_path``, write Parquet sidecar under ``out_dir``."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    catalog_by_chrom = load_catalog(catalog_path)

    sites_path = out / "sites.parquet"
    entries_path = out / "entries.parquet" if with_entries else None
    manifest_path = out / "manifest.json"

    sites_buffer: list[dict[str, Any]] = []
    entries_buffer: list[dict[str, Any]] = []
    n_variants = 0

    writer: pq.ParquetWriter | None = None
    entry_writer: pq.ParquetWriter | None = None

    try:
        with pysam.VariantFile(str(vcf_path)) as vf:
            samples = list(vf.header.samples)
            for rec in vf:
                n_variants += 1
                chrom = norm_chrom(rec.chrom)
                locs = catalog_by_chrom.get(chrom, [])
                vs, ve = _variant_span(rec)
                loc, ov, mmeth = best_catalog_match(chrom, vs, ve, locs)
                motif_v = _info_scalar(rec, "MOTIF")
                motif_s = str(motif_v) if motif_v is not None else ""
                mm = (
                    motif_matches(motif_s, loc.motif_primary, loc.motifs_accepted)
                    if loc and motif_s
                    else None
                )
                row = _site_row(rec, loc, ov, mmeth, mm)
                sites_buffer.append(row)

                if with_entries and samples:
                    vid = row["variant_id"]
                    for s in samples:
                        gt_str, is_carrier, num_alt = _gt_summary(rec, s)
                        entries_buffer.append(
                            {
                                "schema_version": schema.SCHEMA_VERSION,
                                "ruleset_version": schema.RULESET_VERSION,
                                "variant_id": vid,
                                "sample": s,
                                "gt_str": gt_str,
                                "is_carrier": is_carrier,
                                "num_alt": num_alt,
                            }
                        )

                if len(sites_buffer) >= batch_size:
                    table = pa.Table.from_pylist(
                        sites_buffer, schema=schema.SITE_ARROW_SCHEMA
                    )
                    if writer is None:
                        writer = pq.ParquetWriter(
                            str(sites_path), schema.SITE_ARROW_SCHEMA
                        )
                    writer.write_table(table)
                    sites_buffer.clear()

                    if with_entries and entry_writer is None and entries_path:
                        entry_writer = pq.ParquetWriter(
                            str(entries_path), schema=schema.ENTRY_ARROW_SCHEMA
                        )
                    if entry_writer is not None and entries_buffer:
                        et = pa.Table.from_pylist(
                            entries_buffer, schema=schema.ENTRY_ARROW_SCHEMA
                        )
                        entry_writer.write_table(et)
                        entries_buffer.clear()

        if sites_buffer:
            table = pa.Table.from_pylist(
                sites_buffer, schema=schema.SITE_ARROW_SCHEMA
            )
            if writer is None:
                writer = pq.ParquetWriter(str(sites_path), schema.SITE_ARROW_SCHEMA)
            writer.write_table(table)
            sites_buffer.clear()

        if with_entries and entries_buffer:
            if entry_writer is None and entries_path:
                entry_writer = pq.ParquetWriter(
                    str(entries_path), schema=schema.ENTRY_ARROW_SCHEMA
                )
            et = pa.Table.from_pylist(
                entries_buffer, schema=schema.ENTRY_ARROW_SCHEMA
            )
            assert entry_writer is not None
            entry_writer.write_table(et)
            entries_buffer.clear()

        if writer is None:
            empty = pa.Table.from_pylist([], schema=schema.SITE_ARROW_SCHEMA)
            pq.write_table(empty, str(sites_path))

        if with_entries and entries_path and entry_writer is None:
            empty_e = pa.Table.from_pylist([], schema=schema.ENTRY_ARROW_SCHEMA)
            pq.write_table(empty_e, str(entries_path))
    finally:
        if writer is not None:
            writer.close()
        if entry_writer is not None:
            entry_writer.close()

    catalog_sha = _sha256_file(Path(catalog_path))
    manifest = {
        "schema_version": schema.SCHEMA_VERSION,
        "ruleset_version": schema.RULESET_VERSION,
        "sites_path": str(sites_path.resolve()),
        "entries_path": str(entries_path.resolve()) if entries_path else None,
        "vcf_path": str(Path(vcf_path).resolve()),
        "catalog_path": str(Path(catalog_path).resolve()),
        "catalog_sha256": catalog_sha,
        "n_variants_processed": n_variants,
        "with_entries": with_entries,
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    return SidecarResult(
        out_dir=out,
        sites_path=sites_path,
        entries_path=entries_path,
        manifest_path=manifest_path,
        n_variants=n_variants,
    )


def _sha256_file(path: Path) -> str:
    import hashlib

    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(65536), b""):
            h.update(block)
    return h.hexdigest()
