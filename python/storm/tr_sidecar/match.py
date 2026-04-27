from __future__ import annotations

from .catalog import TrLocus


def overlap_bp(a0: int, a1: int, b0: int, b1: int) -> int:
    """Overlap width between half-open intervals [a0,a1) and [b0,b1)."""
    return max(0, min(a1, b1) - max(a0, b0))


def best_catalog_match(
    contig: str,
    var_start: int,
    var_end: int,
    candidates: list[TrLocus],
) -> tuple[TrLocus | None, int, str]:
    """Return (best_locus, overlap_bp, match_method)."""
    if not candidates:
        return None, 0, "no_catalog_on_chrom"
    best: TrLocus | None = None
    best_ov = -1
    for loc in candidates:
        ov = overlap_bp(var_start, var_end, loc.start, loc.end)
        if ov > best_ov:
            best_ov = ov
            best = loc
        elif ov == best_ov and ov > 0 and best is not None:
            # deterministic tie-break
            t = (loc.chrom, loc.start, loc.end, loc.tr_locus_id)
            u = (best.chrom, best.start, best.end, best.tr_locus_id)
            if t < u:
                best = loc
    if best_ov <= 0:
        return None, 0, "no_overlap"
    return best, best_ov, "best_overlap"


def motif_matches(
    vcf_motif: str | None,
    motif_primary: str,
    motifs_accepted: str,
) -> bool | None:
    """True/False if ``vcf_motif`` present; None if no VCF motif to compare."""
    if not vcf_motif:
        return None
    vm = vcf_motif.strip().upper()
    if not vm:
        return None
    primary = motif_primary.strip().upper()
    if primary and vm == primary:
        return True
    accepted = [x.strip().upper() for x in motifs_accepted.split(",") if x.strip()]
    if vm in accepted:
        return True
    return False


def compute_rule_applicable(
    loc: TrLocus | None,
    motif_match: bool | None,
) -> bool:
    if loc is None:
        return False
    if not loc.sv_scan_default_applicable:
        return False
    if loc.motif_change_required:
        # v1: cannot verify motif change from INFO alone
        return False
    if motif_match is False:
        return False
    return True
