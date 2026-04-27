from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path


def norm_chrom(contig: str) -> str:
    c = str(contig).strip()
    if not c:
        return c
    return c if c.startswith("chr") else f"chr{c}"


@dataclass(frozen=True)
class TrLocus:
    tr_locus_id: str
    chrom: str
    start: int
    end: int
    gene: str
    motif_primary: str
    motifs_accepted: str
    rule_group: str
    inheritance: str
    benign_max_repeats: str
    pathogenic_min_repeats: str
    motif_change_required: bool
    contraction_pathogenic: bool
    prioritize_in_sv_scan: bool
    sv_scan_default_applicable: bool
    notes: str


def load_catalog(path: str | Path) -> dict[str, list[TrLocus]]:
    """Load TSV into ``chrom -> [TrLocus, ...]`` (chrom normalized with ``chr``)."""
    p = Path(path)
    by_chrom: dict[str, list[TrLocus]] = {}
    with p.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            chrom = norm_chrom(row["chrom"])
            loc = TrLocus(
                tr_locus_id=row["tr_locus_id"].strip(),
                chrom=chrom,
                start=int(row["start"]),
                end=int(row["end"]),
                gene=row.get("gene", "").strip(),
                motif_primary=row.get("motif_primary", "").strip(),
                motifs_accepted=row.get("motifs_accepted", "").strip(),
                rule_group=row.get("rule_group", "").strip(),
                inheritance=row.get("inheritance", "").strip(),
                benign_max_repeats=row.get("benign_max_repeats", "").strip(),
                pathogenic_min_repeats=row.get("pathogenic_min_repeats", "").strip(),
                motif_change_required=_parse_bool(
                    row.get("motif_change_required", "0")
                ),
                contraction_pathogenic=_parse_bool(
                    row.get("contraction_pathogenic", "0")
                ),
                prioritize_in_sv_scan=_parse_bool(
                    row.get("prioritize_in_sv_scan", "1")
                ),
                sv_scan_default_applicable=_parse_bool(
                    row.get("sv_scan_default_applicable", "1")
                ),
                notes=row.get("notes", "").strip(),
            )
            by_chrom.setdefault(chrom, []).append(loc)
    return by_chrom


def _parse_bool(s: str) -> bool:
    return str(s).strip() in ("1", "true", "True", "yes", "YES")
