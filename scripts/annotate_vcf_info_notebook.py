"""
Notebook-friendly INFO overlay for VCF/BCF files.

This module mirrors the behavior of `scripts/annotate_vcf_info.py` while exposing
Python-callable functions that are easier to run in Jupyter cells.
"""

from __future__ import annotations

import logging
import os
import time
from dataclasses import dataclass
from typing import Any, Iterator, Literal, TypeVar

import pysam

T = TypeVar("T")
OnInfoConflict = Literal["prefer_primary", "prefer_secondary", "error"]

logger = logging.getLogger("annotate_vcf_info_notebook")


@dataclass(frozen=True)
class AnnotationStats:
    output_path: str
    secondary_loaded_ids: int
    secondary_duplicate_ids: int
    secondary_skipped_no_id: int
    primary_rows_written: int
    primary_rows_merged: int
    primary_skipped_no_id: int
    secondary_unmatched_ids: int
    elapsed_seconds: float
    stream_write_seconds: float


def setup_notebook_logging(level: str = "INFO") -> None:
    """Configure logger output once per notebook session."""
    if logging.getLogger().handlers:
        return
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="%(asctime)s [%(levelname)s] %(message)s",
    )


def id_key(rec: pysam.VariantRecord) -> str | None:
    rid = rec.id
    if rid is None or rid == "" or rid == ".":
        return None
    return str(rid)


def _info_value_equal(a: Any, b: Any) -> bool:
    if a is b:
        return True
    try:
        return bool(a == b)
    except Exception:
        return False


def _vcf_meta_quote(s: str) -> str:
    return s.replace("\\", "\\\\").replace('"', '\\"')


def augment_header(dst: pysam.VariantHeader, src: pysam.VariantHeader) -> list[str]:
    notes: list[str] = []

    for key in src.info:
        if key in dst.info:
            continue
        m = src.info[key]
        dst.info.add(m.name, m.number, m.type, m.description)
        notes.append(f"added INFO/{key} from secondary header")

    for key in src.formats:
        if key in dst.formats:
            continue
        m = src.formats[key]
        dst.formats.add(m.name, m.number, m.type, m.description)
        notes.append(f"added FORMAT/{key} from secondary header")

    for key in src.filters:
        if key in dst.filters:
            continue
        desc = _vcf_meta_quote(src.filters[key].description)
        dst.add_line(f'##FILTER=<ID={key},Description="{desc}">')
        notes.append(f"added FILTER/{key} from secondary header")

    for contig in src.contigs:
        if contig in dst.contigs:
            continue
        try:
            ln = src.get_reference_length(contig)
        except ValueError:
            ln = None
        if ln is not None:
            dst.add_line(f"##contig=<ID={contig},length={ln}>")
        else:
            dst.add_line(f"##contig=<ID={contig}>")
        notes.append(f"added contig {contig} from secondary header")

    return notes


def merge_headers(
    primary: pysam.VariantHeader, secondary: pysam.VariantHeader
) -> tuple[pysam.VariantHeader, list[str]]:
    out = primary.copy()
    notes = augment_header(out, secondary)
    return out, notes


def assert_same_samples(h1: pysam.VariantHeader, h2: pysam.VariantHeader) -> None:
    s1, s2 = list(h1.samples), list(h2.samples)
    if s1 != s2:
        raise ValueError(
            "Sample columns differ; annotate requires identical sample names and order.\n"
            f"  primary:   {s1}\n  secondary: {s2}"
        )


def _try_tqdm_notebook() -> Any:
    try:
        from tqdm.auto import tqdm  # type: ignore[import-not-found]

        return tqdm
    except ImportError:
        return None


def _progress_iterator(
    it: Iterator[T],
    *,
    desc: str,
    tqdm_fn: Any,
    log_every: int,
    total: int | None = None,
) -> Iterator[T]:
    if tqdm_fn is not None:
        yield from tqdm_fn(it, desc=desc, total=total, unit="var", mininterval=0.3)
        return

    n = 0
    t0 = time.perf_counter()
    for x in it:
        n += 1
        if log_every > 0 and n % log_every == 0:
            dt = time.perf_counter() - t0
            rate = n / dt if dt > 0 else 0.0
            if total is not None:
                logger.info(
                    "%s: %d / %d (%.1f%%) - %.0f var/s",
                    desc,
                    n,
                    total,
                    100.0 * n / total,
                    rate,
                )
            else:
                logger.info("%s: %d - %.0f var/s", desc, n, rate)
        yield x
    dt = time.perf_counter() - t0
    if n:
        logger.info("%s: finished %d in %.1fs (%.0f var/s)", desc, n, dt, n / dt if dt > 0 else 0.0)


def _apply_secondary_info_dict_in_place(
    rec: pysam.VariantRecord, secondary: dict[str, Any], on_conflict: OnInfoConflict
) -> None:
    locus = f"ID {rec.id!r} @ {rec.chrom}:{rec.pos}"
    for key, sval in secondary.items():
        if key not in rec.info:
            rec.info[key] = sval
            continue
        if _info_value_equal(rec.info[key], sval):
            continue
        if on_conflict == "error":
            raise ValueError(f"INFO conflict at {locus} on tag {key!r}: {rec.info[key]!r} vs {sval!r}")
        if on_conflict == "prefer_primary":
            continue
        if on_conflict == "prefer_secondary":
            rec.info[key] = sval
            continue
        raise ValueError(f"unknown on_conflict mode: {on_conflict!r}")


def _materialize_info(rec: pysam.VariantRecord) -> dict[str, Any]:
    return {k: rec.info[k] for k in rec.info}


def annotate_vcf_info(
    primary_vcf: str,
    secondary_vcf: str,
    output: str,
    *,
    on_info_conflict: OnInfoConflict = "prefer_primary",
    allow_extra_secondary_rows: bool = False,
    allow_different_samples: bool = False,
    progress: bool = True,
    progress_every: int = 2000,
    log_level: str = "INFO",
) -> AnnotationStats:
    """
    Overlay INFO fields from `secondary_vcf` onto `primary_vcf`, joined by ID.

    Returns summary stats for quick inspection in notebook cells.
    """
    setup_notebook_logging(log_level)
    tqdm_fn = _try_tqdm_notebook() if progress else None
    log_every = max(0, int(progress_every))
    run_start = time.perf_counter()

    logger.info("Primary:   %s", primary_vcf)
    logger.info("Secondary: %s", secondary_vcf)
    logger.info("Output:    %s", output)

    read_start = time.perf_counter()
    with pysam.VariantFile(primary_vcf) as vf1, pysam.VariantFile(secondary_vcf) as vf2:
        if not allow_different_samples:
            assert_same_samples(vf1.header, vf2.header)

        merged_h, header_notes = merge_headers(vf1.header, vf2.header)
        if header_notes:
            logger.info("Merged %d new header line(s) from secondary", len(header_notes))

        overlay: dict[str, dict[str, Any]] = {}
        secondary_duplicate_ids = 0
        secondary_skipped_no_id = 0

        for rec in _progress_iterator(
            vf2.fetch(),
            desc="Read secondary",
            tqdm_fn=tqdm_fn,
            log_every=log_every,
        ):
            key = id_key(rec)
            if key is None:
                secondary_skipped_no_id += 1
                continue
            if key in overlay:
                secondary_duplicate_ids += 1
            overlay[key] = _materialize_info(rec)

        secondary_loaded_ids = len(overlay)
        logger.info("Loaded %d ID(s) from secondary (%.1fs)", secondary_loaded_ids, time.perf_counter() - read_start)

    unseen: set[str] = set(overlay)
    stream_start = time.perf_counter()
    primary_rows_merged = 0
    primary_rows_written = 0
    primary_skipped_no_id = 0

    with pysam.VariantFile(primary_vcf) as vf1, pysam.VariantFile(output, "w", header=merged_h) as out:
        for rec in _progress_iterator(
            vf1.fetch(),
            desc="Write",
            tqdm_fn=tqdm_fn,
            log_every=log_every,
        ):
            primary_rows_written += 1
            rec.translate(merged_h)
            key = id_key(rec)
            if key is None:
                primary_skipped_no_id += 1
            elif key in overlay:
                _apply_secondary_info_dict_in_place(rec, overlay[key], on_info_conflict)
                primary_rows_merged += 1
                unseen.discard(key)
            out.write(rec)

    if unseen and not allow_extra_secondary_rows:
        try:
            os.remove(output)
        except OSError as exc:
            logger.warning("Could not remove output %s: %s", output, exc)
        example = next(iter(unseen))
        raise ValueError(
            f"{len(unseen)} ID(s) in secondary not found in primary (example: {example!r}). "
            "Set allow_extra_secondary_rows=True to ignore this."
        )

    if unseen:
        logger.warning("Ignored %d ID(s) only in secondary", len(unseen))

    elapsed_seconds = time.perf_counter() - run_start
    stream_write_seconds = time.perf_counter() - stream_start
    logger.info(
        "Wrote %d row(s); merged secondary INFO for %d of them. Elapsed %.1fs.",
        primary_rows_written,
        primary_rows_merged,
        elapsed_seconds,
    )

    return AnnotationStats(
        output_path=output,
        secondary_loaded_ids=secondary_loaded_ids,
        secondary_duplicate_ids=secondary_duplicate_ids,
        secondary_skipped_no_id=secondary_skipped_no_id,
        primary_rows_written=primary_rows_written,
        primary_rows_merged=primary_rows_merged,
        primary_skipped_no_id=primary_skipped_no_id,
        secondary_unmatched_ids=len(unseen),
        elapsed_seconds=elapsed_seconds,
        stream_write_seconds=stream_write_seconds,
    )
