"""SV + TRGT Plotly visualization (loaded by ``notebooks/visualize.ipynb``)."""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional

import pysam
import plotly.graph_objects as go
import plotly.io as pio

# Jupyter Classic / Lab: try notebook_connected first; fall back to browser if needed.
pio.renderers.default = "notebook_connected"


def find_data_dir() -> Path:
    cwd = Path.cwd()
    for p in (cwd, cwd.parent):
        if (p / "1000920.phase2.bcf").exists():
            return p
    return cwd


def parse_locus(spec: str) -> tuple[str, int, int]:
    """Parse ``chrom:start-end`` into (chrom, start, end) with 1-based inclusive coordinates.

    Commas in numbers are allowed (e.g. ``chr6:31,803,187-32,050,925``).
    Works for any contig name without spaces.
    """
    s = spec.strip().replace(",", "").replace(" ", "")
    m = re.match(r"^([^:]+):(\d+)-(\d+)$", s)
    if not m:
        raise ValueError(
            f"Bad locus {spec!r}; expected chrom:start-end (e.g. chr6:31803187-32050925)"
        )
    chrom, a, b = m.group(1), int(m.group(2)), int(m.group(3))
    if b < a:
        raise ValueError(f"Locus end < start: {spec!r}")
    return chrom, a, b


# Vertical bands (plot y) for haplotype tracks
Y = {
    "trgt_h1": 3.2,
    "trgt_h2": 2.65,
    "sv_h1": 2.05,
    "sv_h2": 1.5,
    "ucsc_simple": 0.76,
    "ucsc_rmsk": 0.52,
    "ucsc_segdup": 0.28,
}
TRACK_H = 0.24
TRACK_UCSC = 0.17

# TRGT and SV variant strokes share this width; INS rects use the same thickness in data units.
VARIANT_LINE_WIDTH_PX = 10
# y-axis view + nominal inner plot height (must stay consistent with ``update_layout`` below).
FIG_Y_RANGE = (0.05, 3.38)
FIG_INNER_PLOT_H_PX = 268.0  # ~ height 500 − top/bottom margins; maps stroke px ↔ y data for INS

# INS has no reference span; drawn width is max(insertion length, INS_PAD_MIN_BP) from VCF POS.
INS_PAD_MIN_BP = 20

# Per-base sequence strips: cap points for very long TR alleles.
INSSEQ_MAX_BP = 20000
INSSEQ_STRIP_H = 0.06
INSSEQ_GAP_BELOW_MAIN = 0.02

# IGV-style nucleotide colors.
BASE_COLOR = {
    "A": "#00A000",
    "C": "#0000E0",
    "G": "#F09000",
    "T": "#E00000",
}
BASE_COLOR_OTHER = "#888888"


def base_color(ch: str) -> str:
    return BASE_COLOR.get((ch or "?").upper(), BASE_COLOR_OTHER)


def insseq_group_id(*, chrom: str, pos: int, hap: int, kind: str, vid: str) -> str:
    return f"{chrom}|{pos}|H{hap + 1}|{kind}|{vid}"


def _variant_stroke_height_data() -> float:
    """Vertical span in y data units matching ``VARIANT_LINE_WIDTH_PX`` on the default figure."""
    y0, y1 = FIG_Y_RANGE
    return (VARIANT_LINE_WIDTH_PX / FIG_INNER_PLOT_H_PX) * (y1 - y0)


@dataclass
class SVRow:
    pos: int  # 1-based VCF POS (pysam rec.pos is 0-based; we convert at parse)
    svtype: str
    svlen: int
    gt: tuple[Optional[int], Optional[int]]
    phased: bool
    vid: str
    qual: float
    filter: str
    info: str
    # Per-haplotype inserted sequence (concrete ALT only); None if not INS / not carried / symbolic.
    ins_seq: tuple[Optional[str], Optional[str]] = (None, None)


@dataclass
class TRGTRow:
    pos: int  # 1-based first ref base (from pysam rec.pos + 1)
    end: int
    trid: str
    gt: tuple[Optional[int], Optional[int]]
    phased: bool
    al: tuple[int, ...]
    mc: tuple[int, ...]
    info: str
    ref: str  # VCF REF (full repeat + padding base)
    vcf_alts: tuple[str, ...]  # VCF ALT consensus strings (same order as alleles 1..n)


def gt_is_called(gt: Any) -> bool:
    if gt is None:
        return False
    if not isinstance(gt, tuple):
        return False
    return all(a is not None for a in gt)


def gt_is_homozygous_reference(gt: tuple[Optional[int], ...]) -> bool:
    """True if every allele is reference (0/0, 0|0, haploid 0)."""
    return all((a or 0) == 0 for a in gt)


def _is_symbolic_alt(alt: str) -> bool:
    return len(alt) >= 2 and alt[0] == "<" and alt[-1] == ">"


def _info_svlen_abs(rec: Any) -> Optional[int]:
    raw = rec.info.get("SVLEN")
    if isinstance(raw, tuple):
        raw = raw[0] if raw else None
    if raw is None:
        return None
    return abs(int(raw))


def _deletion_bp_symbolic(rec: Any) -> int:
    """Deleted reference length for symbolic DEL (SVLEN, END, or ref span)."""
    sl = _info_svlen_abs(rec)
    if sl is not None:
        return sl
    raw_end = rec.info.get("END")
    if isinstance(raw_end, tuple):
        raw_end = raw_end[0] if raw_end else None
    if raw_end is not None:
        return max(1, int(raw_end) - int(rec.pos))
    return max(1, int(rec.stop) - int(rec.pos))


def insertion_bp_from_ref_alt(ref: str, alt: str, rec: Any) -> int:
    """Inserted bases from REF/ALT; symbolic INS/DUP uses |SVLEN| when present."""
    if not alt:
        return 0
    if _is_symbolic_alt(alt):
        u = alt.upper()
        if "INS" in u or "DUP" in u:
            sl = _info_svlen_abs(rec)
            return sl if sl is not None else 0
        return 0
    return max(0, len(alt) - len(ref))


def deletion_bp_from_ref_alt(ref: str, alt: str, rec: Any) -> int:
    """Deleted reference bases from REF/ALT; symbolic DEL uses SVLEN/END/stop."""
    if not alt:
        return 0
    if _is_symbolic_alt(alt):
        if "DEL" in alt.upper():
            return _deletion_bp_symbolic(rec)
        return 0
    return max(0, len(ref) - len(alt))


INFO_VALUE_MAX_CHARS = 30


def _truncate_info_value(s: str, max_chars: int = INFO_VALUE_MAX_CHARS) -> str:
    if len(s) <= max_chars:
        return s
    return s[:max_chars] + "..."


def _html_escape(s: str) -> str:
    # Make INFO safe for HTML rendering in the sidebar.
    return (
        s.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
        .replace("'", "&#39;")
    )


def _info_value_to_str(v: Any) -> str:
    if isinstance(v, (tuple, list)):
        return ",".join(_info_value_to_str(x) for x in v)
    if isinstance(v, (bytes, bytearray)):
        raw = bytes(v).decode("ascii", errors="replace")
        return _html_escape(_truncate_info_value(raw))
    if v is None:
        return "."
    return _html_escape(_truncate_info_value(str(v)))


def record_info_summary(rec: Any) -> str:
    keys = sorted(rec.info.keys())
    if not keys:
        return "&nbsp;&nbsp;."
    return "<br>".join(
        f"&nbsp;&nbsp;{k}={_info_value_to_str(rec.info.get(k))}" for k in keys
    )


def _norm_vcf_seq(x: Any) -> str:
    """VCF REF/ALT as str (pysam may return ``bytes`` in some builds)."""
    if x is None:
        return ""
    if isinstance(x, (bytes, bytearray)):
        return bytes(x).decode("ascii", errors="replace")
    return str(x)


def _sv_hap_ins_seq(
    gt: tuple[Optional[int], ...],
    hap: int,
    ref: str,
    alts: tuple[Any, ...],
    rec: Any,
) -> Optional[str]:
    """Inserted bases on ``hap`` for a VCF record; ``None`` if not a concrete insertion."""
    if hap >= len(gt):
        return None
    ai = gt[hap]
    if ai is None or ai <= 0:
        return None
    if ai - 1 < 0 or ai - 1 >= len(alts):
        return None
    alt = _norm_vcf_seq(alts[ai - 1])
    if insertion_bp_from_ref_alt(ref, alt, rec) <= 0:
        return None
    if _is_symbolic_alt(alt):
        return None
    tail = alt[len(ref) :]
    if not tail:
        return None
    return tail.upper()


def parse_sv_rows(path: Path, chrom: str, start: int, end: int, sample: str) -> list[SVRow]:
    """Load indel SV records from REF/ALT (and symbolic alleles + SVLEN/END).

    Stored ``pos`` is **1-based VCF POS**, matching ``parse_locus`` coordinates.
    """
    rows: list[SVRow] = []
    with pysam.VariantFile(str(path)) as vf:
        if sample not in vf.header.samples:
            raise ValueError(f"Sample {sample!r} not in {path}")
        for rec in vf.fetch(chrom, start - 1, end):
            s = rec.samples[sample]
            gt = s.get("GT")
            if not gt_is_called(gt):
                continue
            gt = tuple(gt)
            if gt_is_homozygous_reference(gt):
                continue
            ref = _norm_vcf_seq(rec.ref)
            alts = tuple(_norm_vcf_seq(a) for a in (rec.alts or ()))
            carried = sorted({a for a in gt if a is not None and a > 0})
            if not carried:
                continue
            ins_lens: list[int] = []
            del_lens: list[int] = []
            for ai in carried:
                if ai - 1 < 0 or ai - 1 >= len(alts):
                    continue
                alt = alts[ai - 1]
                ib = insertion_bp_from_ref_alt(ref, alt, rec)
                db = deletion_bp_from_ref_alt(ref, alt, rec)
                if ib > 0:
                    ins_lens.append(ib)
                if db > 0:
                    del_lens.append(db)
            if not ins_lens and not del_lens:
                continue
            pos = int(rec.pos) + 1  # 1-based VCF POS; aligns with win_start / ref_left
            base_vid = str(rec.id or str(rec.pos))
            qual = float(rec.qual or 0.0)
            flt = ";".join(rec.filter.keys()) if rec.filter else "PASS"
            info_txt = record_info_summary(rec)
            phased = bool(s.phased)
            if del_lens:
                L = max(del_lens)
                suf = "_DEL" if ins_lens else ""
                rows.append(
                    SVRow(
                        pos=pos,
                        svtype="DEL",
                        svlen=L,
                        gt=gt,
                        phased=phased,
                        vid=f"{base_vid}{suf}",
                        qual=qual,
                        filter=flt,
                        info=info_txt,
                    )
                )
            if ins_lens:
                L = max(ins_lens)
                suf = "_INS" if del_lens else ""
                ins0 = _sv_hap_ins_seq(gt, 0, ref, alts, rec)
                ins1 = _sv_hap_ins_seq(gt, 1, ref, alts, rec) if len(gt) > 1 else None
                rows.append(
                    SVRow(
                        pos=pos,
                        svtype="INS",
                        svlen=L,
                        gt=gt,
                        phased=phased,
                        vid=f"{base_vid}{suf}",
                        qual=qual,
                        filter=flt,
                        info=info_txt,
                        ins_seq=(ins0, ins1),
                    )
                )
    return rows


def _fmt_tuple_ints(x: Any) -> tuple[int, ...]:
    if x is None:
        return ()
    if isinstance(x, tuple):
        return tuple(int(v) for v in x)
    return (int(x),)


def parse_trgt_rows(path: Path, chrom: str, start: int, end: int, sample: str) -> list[TRGTRow]:
    """TRGT rows use 1-based ``pos`` and inclusive ``end`` (VCF / INFO END)."""
    rows: list[TRGTRow] = []
    with pysam.VariantFile(str(path)) as vf:
        if sample not in vf.header.samples:
            raise ValueError(f"Sample {sample!r} not in {path}")
        for rec in vf.fetch(chrom, start - 1, end):
            s = rec.samples[sample]
            gt = s.get("GT")
            if not gt_is_called(gt):
                continue
            gt = tuple(gt)
            if gt_is_homozygous_reference(gt):
                continue
            raw_end = rec.info.get("END")
            if isinstance(raw_end, tuple):
                raw_end = raw_end[0] if raw_end else None
            end_i = int(raw_end) if raw_end is not None else int(rec.stop)
            ref_s = _norm_vcf_seq(rec.ref)
            alt_tpl = tuple(_norm_vcf_seq(a) for a in (rec.alts or ()))
            rows.append(
                TRGTRow(
                    pos=int(rec.pos) + 1,
                    end=end_i,
                    trid=str(rec.info.get("TRID", rec.id or "")),
                    gt=gt,
                    phased=bool(s.phased),
                    al=_fmt_tuple_ints(s.get("AL")),
                    mc=_fmt_tuple_ints(s.get("MC")),
                    info=record_info_summary(rec),
                    ref=ref_s,
                    vcf_alts=alt_tpl,
                )
            )
    return rows


def ref_left_linear(win_start: int, r: int) -> float:
    """Left edge of reference base r in linear (no-gap) plot coordinates."""
    return float(r - win_start)


def plot_x_to_ref_bp(x: float, win_start: int, win_end: int) -> float:
    """Map linear plot-x back to a continuous 1-based reference coordinate."""
    ref_w = win_end - win_start + 1
    xf = max(0.0, min(float(ref_w), float(x)))
    return float(win_start) + xf


def _linspace(a: float, b: float, n: int) -> list[float]:
    if n <= 1:
        return [a]
    return [a + (b - a) * i / (n - 1) for i in range(n)]


def display_figure_with_genomic_cursor(
    fig: go.Figure,
    chrom: str,
    win_start: int,
    win_end: int,
    *,
    n_hover_pts: int = 2000,
) -> Any:
    """Jupyter: status line + cursor strip + clickable variant info panel.
    If ipywidgets/FigureWidget unavailable, falls back to ``fig.show()`` and returns ``None``.
    """
    if not any(getattr(t, "name", None) == "_cursor_ref" for t in fig.data):
        x0 = 0.0
        x1 = float(win_end - win_start + 1)
        n = max(200, min(4000, n_hover_pts))
        xs = _linspace(x0, x1, n)
        texts = [
            f"{chrom}:{int(round(plot_x_to_ref_bp(xv, win_start, win_end))):,}"
            for xv in xs
        ]
        traces = list(fig.data)
        vm = [t for t in traces if getattr(t, "name", None) == "_ref_variant_markers"]
        other = [t for t in traces if getattr(t, "name", None) != "_ref_variant_markers"]
        cursor = go.Scatter(
            x=xs,
            y=[1.0] * len(xs),
            mode="markers",
            marker=dict(size=18, opacity=0.004, color="#ffffff"),
            text=texts,
            hovertemplate="%{text}<extra></extra>",
            showlegend=False,
            name="_cursor_ref",
        )
        fig = go.Figure(data=other + [cursor] + vm, layout=fig.layout)
    try:
        from ipywidgets import HTML, HBox, VBox, Layout
        from plotly.graph_objects import FigureWidget
    except ImportError:
        fig.show()
        return None
    status = HTML(
        value='<div style="font-family:ui-monospace,monospace;color:#64748b">'
        "Move mouse over plot…</div>"
    )
    panel = HTML(
        value=(
            '<div style="height:500px;overflow-y:auto;border:1px solid #e2e8f0;'
            'border-radius:6px;padding:8px 10px;background:#f8fafc;'
            'font-family:ui-monospace,monospace;font-size:11px;color:#64748b">'
            'Click a TRGT/SV variant to view details here.</div>'
        ),
        layout=Layout(width="360px"),
    )
    fw = FigureWidget(fig)
    fw.update_layout(hovermode="closest")

    # Use trace *indices* into ``fw.data``: FigureWidget may not apply ``visible``
    # updates to stale trace object references captured at build time.
    insseq_indices_by_group: dict[str, list[int]] = {}
    for i, tr in enumerate(fw.data):
        gid_k: Optional[str] = None
        meta = getattr(tr, "meta", None)
        if isinstance(meta, dict) and meta.get("insseq_group"):
            gid_k = str(meta["insseq_group"])
        else:
            nm = getattr(tr, "name", None) or ""
            if isinstance(nm, str) and nm.startswith("_insseq:"):
                gid_k = nm[len("_insseq:") :]
        if gid_k:
            insseq_indices_by_group.setdefault(gid_k, []).append(i)

    def _insseq_group_is_visible(gid: str) -> bool:
        for j in insseq_indices_by_group.get(gid, []):
            if fw.data[j].visible:
                return True
        return False

    def _cd_cell_detail(cell: Any) -> str:
        if isinstance(cell, str):
            return cell
        try:
            if len(cell) >= 1:  # type: ignore[arg-type]
                return str(cell[0])
        except Exception:
            pass
        return str(cell)

    def _cd_cell_group(cell: Any) -> str:
        if isinstance(cell, str):
            return ""
        try:
            if len(cell) >= 2 and cell[1] is not None:  # type: ignore[arg-type]
                return str(cell[1]).strip()
        except Exception:
            pass
        return ""

    def on_hover(trace: Any, points: Any, state: Any) -> None:
        if not getattr(points, "point_inds", None):
            return
        name = getattr(trace, "name", None)
        if name not in ("_cursor_ref", "_ref_variant_markers"):
            return
        idx = points.point_inds[0]
        tx = trace.text[idx] if trace.text is not None else ""
        status.value = (
            f'<div style="font-family:ui-monospace,monospace;font-size:14px">'
            f"<b>{tx}</b></div>"
        )

    def on_unhover(_trace: Any, _points: Any, _state: Any) -> None:
        status.value = '<div style="font-family:ui-monospace,monospace;color:#64748b">—</div>'

    def _point_text(trace: Any, idx: int) -> str:
        def _from_attr(attr: str) -> str:
            v = getattr(trace, attr, None)
            if v is None:
                return ""
            if isinstance(v, str):
                return v
            try:
                if 0 <= idx < len(v):
                    if attr == "customdata":
                        return _cd_cell_detail(v[idx])
                    return str(v[idx])
            except Exception:
                return ""
            return ""

        # Prefer customdata for sidebar; fall back to hover/text.
        for attr in ("customdata", "hovertext", "text"):
            out = _from_attr(attr)
            if out:
                return out
        return ""

    def _point_customdata(trace: Any, idx: int) -> str:
        v = getattr(trace, "customdata", None)
        if v is None:
            return ""
        try:
            if 0 <= idx < len(v):
                out = _cd_cell_detail(v[idx])
                if out.strip():
                    return out
        except Exception:
            return ""
        return ""

    def _point_insseq_group(trace: Any, idx: int) -> str:
        v = getattr(trace, "customdata", None)
        if v is None:
            return ""
        try:
            if 0 <= idx < len(v):
                return _cd_cell_group(v[idx])
        except Exception:
            return ""
        return ""

    def _trace_insseq_group_scan(trace: Any) -> str:
        """Filled polygons sometimes report a point index whose customdata row lacks the group; scan the trace."""
        v = getattr(trace, "customdata", None)
        if v is None:
            return ""
        try:
            for item in v:
                g = _cd_cell_group(item)
                if g:
                    return g
        except Exception:
            return ""
        return ""

    def _first_nonempty_customdata(trace: Any) -> str:
        v = getattr(trace, "customdata", None)
        if v is None:
            return ""
        if isinstance(v, str):
            return v if v.strip() else ""
        try:
            for item in v:
                s = _cd_cell_detail(item)
                if s.strip():
                    return s
        except Exception:
            return ""
        return ""

    def on_variant_click(trace: Any, points: Any, state: Any) -> None:
        if not getattr(points, "point_inds", None):
            return
        idx = points.point_inds[0]
        gid = (
            _point_insseq_group(trace, idx) or _trace_insseq_group_scan(trace)
        ).strip()
        detail = (
            _point_customdata(trace, idx)
            or _first_nonempty_customdata(trace)
            or _point_text(trace, idx)
        )
        if detail.strip():
            panel.value = (
                '<div style="height:500px;overflow-y:auto;border:1px solid #e2e8f0;'
                'border-radius:6px;padding:8px 10px;background:#f8fafc;'
                'font-family:ui-monospace,monospace;font-size:11px;line-height:1.35">'
                f"{detail}</div>"
            )

        if not gid:
            return
        # Per-variant toggle: multiple strips can stay visible at once.
        if _insseq_group_is_visible(gid):
            for j in insseq_indices_by_group.get(gid, []):
                fw.data[j].visible = False
        else:
            for j in insseq_indices_by_group.get(gid, []):
                fw.data[j].visible = True

    for tr in fw.data:
        nm = getattr(tr, "name", None) or ""
        if nm in ("_cursor_ref", "_ref_variant_markers"):
            tr.on_hover(on_hover)
            if hasattr(tr, "on_unhover"):
                tr.on_unhover(on_unhover)
        if nm.startswith("SV ") or nm.startswith("TRGT "):
            tr.on_click(on_variant_click)
    return VBox([status, HBox([fw, panel], layout=Layout(width="100%"))])


@dataclass
class UCSCInterval:
    """1-based inclusive coordinates on the reference assembly."""

    start: int
    end: int
    hover: str


def _ucsc_mysql_str(v: Any) -> str:
    """pymysql returns some columns (e.g. simpleRepeat.sequence) as bytes."""
    if v is None:
        return ""
    if isinstance(v, (bytes, bytearray)):
        return bytes(v).decode("ascii", errors="replace")
    return str(v)


def fetch_ucsc_repeat_tracks(
    genome: str,
    chrom: str,
    win_start: int,
    win_end: int,
) -> tuple[list[UCSCInterval], list[UCSCInterval], list[UCSCInterval], Optional[str]]:
    """Load simpleRepeat, RepeatMasker (rmsk), and segmental dups overlapping the window.

    Uses UCSC public MySQL (``genome-mysql.soe.ucsc.edu``). Coordinates in the DB are
    0-based half-open; returned intervals are 1-based inclusive to match the rest of this notebook.
    """
    import pymysql

    ws, we = win_start, win_end
    simple_rows: list[UCSCInterval] = []
    rmsk_rows: list[UCSCInterval] = []
    seg_rows: list[UCSCInterval] = []
    try:
        conn = pymysql.connect(
            host="genome-mysql.soe.ucsc.edu",
            user="genome",
            database=genome,
            connect_timeout=12,
            read_timeout=90,
        )
        try:
            with conn.cursor() as cur:
                cur.execute(
                    """
                    SELECT chromStart, chromEnd, name, period, copyNum, sequence
                    FROM simpleRepeat
                    WHERE chrom=%s AND chromEnd > %s AND chromStart < %s
                    """,
                    (chrom, ws - 1, we),
                )
                for cs, ce, name, period, copy_num, seq in cur.fetchall():
                    s1, e1 = cs + 1, ce
                    name_s = _ucsc_mysql_str(name)
                    seq_s = _ucsc_mysql_str(seq)
                    if len(seq_s) > 48:
                        seq_s = seq_s[:48] + "…"
                    hover = (
                        f"<b>simpleRepeat</b><br>{chrom}:{s1}-{e1}<br>"
                        f"{name_s} period={period} copyNum={copy_num}"
                        + (f"<br>{seq_s}" if seq_s else "")
                    )
                    simple_rows.append(UCSCInterval(s1, e1, hover))

                cur.execute(
                    """
                    SELECT genoStart, genoEnd, repName, repClass, repFamily, strand
                    FROM rmsk
                    WHERE genoName=%s AND genoEnd > %s AND genoStart < %s
                    """,
                    (chrom, ws - 1, we),
                )
                for gs, ge, rep_name, rep_class, rep_family, strand in cur.fetchall():
                    s1, e1 = gs + 1, ge
                    hover = (
                        f"<b>RepeatMasker</b><br>{chrom}:{s1}-{e1}<br>"
                        f"{_ucsc_mysql_str(rep_name)} / {_ucsc_mysql_str(rep_class)} / "
                        f"{_ucsc_mysql_str(rep_family)}<br>strand {_ucsc_mysql_str(strand)}"
                    )
                    rmsk_rows.append(UCSCInterval(s1, e1, hover))

                cur.execute(
                    """
                    SELECT chromStart, chromEnd, name, otherChrom, otherStart, otherEnd
                    FROM genomicSuperDups
                    WHERE chrom=%s AND chromEnd > %s AND chromStart < %s
                    """,
                    (chrom, ws - 1, we),
                )
                for cs, ce, name, och, os, oe in cur.fetchall():
                    s1, e1 = cs + 1, ce
                    pair_s = _ucsc_mysql_str(name)
                    och_s = _ucsc_mysql_str(och)
                    other = (
                        f"{och_s}:{int(os) + 1}-{oe}" if och_s and os is not None and oe is not None else och_s
                    )
                    hover = (
                        f"<b>SegDup</b><br>{chrom}:{s1}-{e1}<br>"
                        f"pair: {pair_s}<br>other: {other}"
                    )
                    seg_rows.append(UCSCInterval(s1, e1, hover))
        finally:
            conn.close()
    except Exception as e:
        return [], [], [], str(e)

    return simple_rows, rmsk_rows, seg_rows, None


def allele_on_haplotype(gt: tuple[Optional[int], ...], hap: int) -> Optional[int]:
    if len(gt) <= hap:
        return None
    return gt[hap]


def sv_affects_hap(r: SVRow, hap: int) -> bool:
    """Haplotype carries a non-reference allele for this SV record."""
    a0, a1 = r.gt[0], r.gt[1] if len(r.gt) > 1 else None
    if a1 is None:
        return (a0 or 0) > 0
    ah = allele_on_haplotype(r.gt, hap)
    return ah is not None and ah > 0


def _variant_header_fields(*, vid: str, pos: int, vtype: str, length: int) -> str:
    return (
        f"<b>ID:</b> {vid}<br>"
        f"<b>POS:</b> {pos}<br>"
        f"<b>TYPE:</b> {vtype}<br>"
        f"<b>Length:</b> {length}"
    )


def variant_sidebar_detail(*, vid: str, pos: int, vtype: str, length: int, info: str) -> str:
    return _variant_header_fields(vid=vid, pos=pos, vtype=vtype, length=length) + (
        f"<br><b>INFO</b><br>{info}"
    )


def trgt_detail(r: TRGTRow) -> str:
    return variant_sidebar_detail(
        vid=r.trid,
        pos=r.pos,
        vtype="TRGT",
        length=trgt_ref_len_bp(r),
        info=r.info,
    )


def trgt_hover(r: TRGTRow, hap: int, chrom: str) -> str:
    # Hover stays compact; full detail goes to the click sidebar.
    return _variant_header_fields(
        vid=r.trid, pos=r.pos, vtype="TRGT", length=trgt_ref_len_bp(r)
    )


def trgt_ref_len_bp(r: TRGTRow) -> int:
    """Reference repeat span (bp) with 1-based inclusive ``pos``..``end``."""
    return max(1, int(r.end) - int(r.pos) + 1)


def trgt_max_expansion_bp(r: TRGTRow) -> int:
    """Max (AL − ref span) over alternate alleles in GT; 0 if no expansion."""
    ref_bp = trgt_ref_len_bp(r)
    m = 0
    for ai in r.gt:
        if ai is None or ai <= 0 or ai >= len(r.al):
            continue
        m = max(m, int(r.al[ai]) - ref_bp)
    return max(0, m)


def trgt_hap_meets_min_allele_delta(
    r: TRGTRow, hap: int, min_allele_bp: int
) -> bool:
    """True if |alt - ref| (AL vs ref span, bp) is at least min_allele_bp."""
    if min_allele_bp <= 0:
        return True
    ai = allele_on_haplotype(r.gt, hap)
    if ai is None or ai >= len(r.al):
        return False
    ref_bp = trgt_ref_len_bp(r)
    return abs(int(r.al[ai]) - ref_bp) >= min_allele_bp


def trgt_allele_seq_for_hap(r: TRGTRow, hap: int) -> Optional[str]:
    """Consensus sequence for the carried ALT on ``hap``; strips TRGT padding base when it matches REF[0]."""
    ai = allele_on_haplotype(r.gt, hap)
    if ai is None or ai <= 0:
        return None
    if ai - 1 < 0 or ai - 1 >= len(r.vcf_alts):
        return None
    seq = r.vcf_alts[ai - 1].upper()
    if r.ref and seq and seq[0] == r.ref[0]:
        seq = seq[1:]
    return seq if seq else None


def sv_ins_seq_for_hap(r: SVRow, hap: int) -> Optional[str]:
    if r.svtype != "INS":
        return None
    if hap >= len(r.ins_seq):
        return None
    return r.ins_seq[hap]


def _clip_insseq(seq: str) -> str:
    if len(seq) <= INSSEQ_MAX_BP:
        return seq
    return seq[:INSSEQ_MAX_BP]


def _add_insseq_bar_trace(
    fig: go.Figure,
    seq: str,
    ref_left_pos: float,
    y_lo: float,
    y_hi: float,
    group_id: str,
) -> None:
    """One Bar trace: one column per bp at true genomic width (centers at half-integer offsets)."""
    seq = _clip_insseq(seq)
    if not seq:
        return
    n = len(seq)
    h = y_hi - y_lo
    xs = [ref_left_pos + i + 0.5 for i in range(n)]
    colors = [base_color(seq[i]) for i in range(n)]
    fig.add_trace(
        go.Bar(
            x=xs,
            y=[h] * n,
            base=[y_lo] * n,
            width=[1.0] * n,
            marker=dict(color=colors, line=dict(width=0)),
            hoverinfo="skip",
            showlegend=False,
            visible=False,
            name=f"_insseq:{group_id}",
            meta=dict(insseq_group=group_id),
            xaxis="x",
            yaxis="y",
        )
    )


def trgt_sites_visible_count(rows: list[TRGTRow], min_allele_bp: int) -> int:
    """Sites with at least one haplotype passing the allele-delta filter."""
    n = 0
    for r in rows:
        if any(trgt_hap_meets_min_allele_delta(r, h, min_allele_bp) for h in (0, 1)):
            n += 1
    return n


def sv_detail(r: SVRow, hap: int) -> str:
    return variant_sidebar_detail(
        vid=r.vid,
        pos=r.pos,
        vtype=r.svtype,
        length=r.svlen,
        info=r.info,
    )


def sv_hover(r: SVRow, hap: int) -> str:
    # Hover stays compact; full detail goes to the click sidebar.
    return _variant_header_fields(vid=r.vid, pos=r.pos, vtype=r.svtype, length=r.svlen)


def nice_reference_ticks(ws: int, we: int, max_ticks: int = 9) -> list[int]:
    """Pick readable 1-based coordinates on a 1/2/5×10ⁿ step (not window-start + fixed stride)."""
    span = we - ws + 1
    if span <= 1:
        return [ws, we]
    nt = min(max_ticks, max(2, span))
    raw = span / (nt - 1)
    exp = int(math.floor(math.log10(raw)))
    f = raw / (10**exp)
    if f <= 1.5:
        nf = 1
    elif f <= 3:
        nf = 2
    elif f <= 7:
        nf = 5
    else:
        nf = 10
    step = max(1, int(nf * 10**exp))
    t0 = math.ceil(ws / step) * step
    ticks: list[int] = []
    t = t0
    while t <= we:
        ticks.append(int(t))
        t += step
    if not ticks:
        return [ws, we]
    # If a "nice" grid tick sits just inside a window end, Plotly stacks two labels (overlap).
    merge_bp = max(5000, step // 5)
    if ticks[0] != ws:
        if ticks[0] - ws < merge_bp:
            ticks[0] = ws
        else:
            ticks.insert(0, ws)
    if ticks[-1] != we:
        if we - ticks[-1] < merge_bp:
            ticks[-1] = we
        else:
            ticks.append(we)
    return ticks


def _expand_multiseg_line(
    segments: list[tuple[float, float]],
    y_mid: float,
    hovers: list[str],
    min_span_data: float,
    *,
    anchor_left: bool = False,
) -> tuple[list[Any], list[Any], list[str]]:
    """Rebuild Scatter line arrays; widen short spans to min_span_data (plot x units).
    If anchor_left, keep the left edge at xa (VCF start); else center the widened span.
    """
    lx: list[Any] = []
    ly: list[Any] = []
    lt: list[str] = []
    ms = max(0.0, float(min_span_data))
    for (xa, xb), h in zip(segments, hovers):
        w = xb - xa
        if w <= 0:
            continue
        w2 = max(w, ms) if ms > 0 else w
        if anchor_left:
            xa2 = xa
            xb2 = xa + w2
        else:
            cx = 0.5 * (xa + xb)
            xa2 = cx - 0.5 * w2
            xb2 = cx + 0.5 * w2
        lx.extend([xa2, xb2, None])
        ly.extend([y_mid, y_mid, None])
        lt.extend([h, h, ""])
    return lx, ly, lt


def _expand_multiseg_rects(
    segments: list[tuple[float, float]],
    y_lo: float,
    y_hi: float,
    hovers: list[str],
    min_span_data: float,
    *,
    anchor_left: bool = True,
) -> tuple[list[Any], list[Any], list[str]]:
    """Closed rectangles (for INS): same width rules as ``_expand_multiseg_line``.
    Each rectangle is a closed loop separated by None for Plotly ``fill='toself'``.
    """
    lx: list[Any] = []
    ly: list[Any] = []
    lt: list[str] = []
    ms = max(0.0, float(min_span_data))
    for (xa, xb), h in zip(segments, hovers):
        w = xb - xa
        if w <= 0:
            continue
        w2 = max(w, ms) if ms > 0 else w
        if anchor_left:
            xa2 = xa
            xb2 = xa + w2
        else:
            cx = 0.5 * (xa + xb)
            xa2 = cx - 0.5 * w2
            xb2 = cx + 0.5 * w2
        lx.extend([xa2, xb2, xb2, xa2, xa2, None])
        ly.extend([y_lo, y_lo, y_hi, y_hi, y_lo, None])
        lt.extend([h, h, h, h, h, ""])
    return lx, ly, lt


def _min_span_plot_units_for_window(
    win_start: int,
    win_end: int,
    min_px: float,
    *,
    layout_width_px: float = 700.0,
    margin_l: float = 52.0,
    margin_r: float = 132.0,
) -> float:
    """Linear x-axis span so a segment is ~min_px wide at default zoom."""
    if min_px <= 0:
        return 0.0
    x0 = ref_left_linear(win_start, win_start)
    x1 = ref_left_linear(win_start, win_end + 1)
    if x1 <= x0:
        return 0.0
    inner = max(120.0, float(layout_width_px) - margin_l - margin_r)
    return float(min_px) * (x1 - x0) / inner


def build_figure(
    sv_path: Path,
    trgt_path: Path,
    chrom: str,
    win_start: int,
    win_end: int,
    sample: Optional[str] = None,
    *,
    genome: str = "hg38",
    ucsc_tracks: bool = True,
    min_allele_bp: int = 5,
    min_segment_px: float = 6.0,
) -> go.Figure:
    with pysam.VariantFile(str(sv_path)) as vf:
        samples = list(vf.header.samples)
    if not samples:
        raise ValueError("No samples in SV VCF")
    sample = sample or samples[0]

    sv_rows = parse_sv_rows(sv_path, chrom, win_start, win_end, sample)
    trgt_rows = parse_trgt_rows(trgt_path, chrom, win_start, win_end, sample)

    ref_left = lambda r: ref_left_linear(win_start, r)
    ref_span = float(win_end - win_start + 1)
    ms_plot = _min_span_plot_units_for_window(
        win_start, win_end, min_segment_px
    )

    ucsc_err: Optional[str] = None
    n_ucsc_simple = n_ucsc_rmsk = n_ucsc_seg = 0
    simple_u: list[UCSCInterval] = []
    rmsk_u: list[UCSCInterval] = []
    seg_u: list[UCSCInterval] = []
    if ucsc_tracks:
        simple_u, rmsk_u, seg_u, ucsc_err = fetch_ucsc_repeat_tracks(
            genome, chrom, win_start, win_end
        )
        n_ucsc_simple, n_ucsc_rmsk, n_ucsc_seg = len(simple_u), len(rmsk_u), len(seg_u)

    fig = go.Figure()

    # Reference axis line (light)
    x0, x1 = ref_left(win_start), ref_left(win_end + 1)
    fig.add_trace(
        go.Scatter(
            x=[x0, x1],
            y=[1.0, 1.0],
            mode="lines",
            line=dict(color="#bbb", width=2),
            name="reference",
            showlegend=False,
            hoverinfo="skip",
        )
    )

    def band_y(key: str) -> tuple[float, float]:
        yc = Y[key]
        return yc - TRACK_H / 2, yc + TRACK_H / 2

    # --- TRGT: ref span bars per haplotype ---
    colors_trgt = ("#3b82f6", "#60a5fa")
    variants_legend_title = True
    for hap in (0, 1):
        key = f"trgt_h{hap + 1}"
        y_lo, y_hi = band_y(key)
        y_mid = 0.5 * (y_lo + y_hi)
        segs: list[tuple[float, float]] = []
        hovs: list[str] = []
        dets: list[str] = []
        gids_tr: list[str] = []
        n_tr = 0
        strip_y_hi = y_lo - INSSEQ_GAP_BELOW_MAIN
        strip_y_lo = strip_y_hi - INSSEQ_STRIP_H
        for r in trgt_rows:
            if allele_on_haplotype(r.gt, hap) is None:
                continue
            if not trgt_hap_meets_min_allele_delta(r, hap, min_allele_bp):
                continue
            xa = ref_left(r.pos)
            xb = ref_left(r.end + 1)
            if xb <= xa:
                continue
            n_tr += 1
            segs.append((xa, xb))
            hovs.append(trgt_hover(r, hap, chrom))
            dets.append(trgt_detail(r))
            raw_sq = trgt_allele_seq_for_hap(r, hap)
            gid = ""
            if raw_sq:
                gid = insseq_group_id(
                    chrom=chrom, pos=r.pos, hap=hap, kind="TRGT", vid=r.trid
                )
                _add_insseq_bar_trace(
                    fig,
                    raw_sq,
                    float(ref_left(r.pos)),
                    strip_y_lo,
                    strip_y_hi,
                    gid,
                )
            gids_tr.append(gid)
        lx, ly, lt = _expand_multiseg_line(
            segs, y_mid, hovs, ms_plot, anchor_left=True
        )
        cd_trgt: list[list[str]] = []
        for (xa, xb), d, gid in zip(segs, dets, gids_tr):
            if xb <= xa:
                continue
            cd_trgt.extend([[d, gid], [d, gid], ["", ""]])
        leg_trgt: dict[str, Any] = {"legendgroup": "variants", "legendrank": 1}
        if variants_legend_title:
            leg_trgt["legendgrouptitle"] = dict(text="Variants")
            variants_legend_title = False
        if segs:
            fig.add_trace(
                go.Scatter(
                    x=lx,
                    y=ly,
                    mode="lines",
                    line=dict(width=VARIANT_LINE_WIDTH_PX, color=colors_trgt[hap], simplify=False),
                    name=f"TRGT H{hap + 1} ({n_tr})",
                    hoverinfo="text",
                    text=lt,
                    customdata=cd_trgt,
                    **leg_trgt,
                )
            )
        else:
            # No segments (e.g. filtered): still show TRGT in legend
            fig.add_trace(
                go.Scatter(
                    x=[],
                    y=[],
                    mode="lines",
                    line=dict(width=VARIANT_LINE_WIDTH_PX, color=colors_trgt[hap], simplify=False),
                    name=f"TRGT H{hap + 1} ({n_tr})",
                    hoverinfo="skip",
                    showlegend=True,
                    **leg_trgt,
                )
            )

    # --- SV: DEL and INS per haplotype ---
    col_del = "#dc2626"
    col_ins = "#16a34a"
    for hap in (0, 1):
        key = f"sv_h{hap + 1}"
        y_lo, y_hi = band_y(key)
        y_mid = 0.5 * (y_lo + y_hi)

        n_sv_del = sum(
            1 for r in sv_rows if r.svtype == "DEL" and sv_affects_hap(r, hap)
        )
        segs_d: list[tuple[float, float]] = []
        hovs_d: list[str] = []
        dets_d: list[str] = []
        for r in sv_rows:
            if r.svtype != "DEL" or not sv_affects_hap(r, hap):
                continue
            xa = ref_left(r.pos)
            xb = ref_left(r.pos + r.svlen)
            segs_d.append((xa, xb))
            hovs_d.append(sv_hover(r, hap))
            dets_d.append(sv_detail(r, hap))
        lx_d, ly_d, lt_d = _expand_multiseg_line(
            segs_d, y_mid, hovs_d, ms_plot, anchor_left=True
        )
        cd_sv_del: list[list[str]] = []
        for (xa, xb), d in zip(segs_d, dets_d):
            if xb <= xa:
                continue
            cd_sv_del.extend([[d, ""], [d, ""], ["", ""]])
        ins_y_half = 0.5 * _variant_stroke_height_data()
        if segs_d:
            leg_sv: dict[str, Any] = {"legendgroup": "variants", "legendrank": 1}
            if variants_legend_title:
                leg_sv["legendgrouptitle"] = dict(text="Variants")
                variants_legend_title = False
            fig.add_trace(
                go.Scatter(
                    x=lx_d,
                    y=ly_d,
                    mode="lines",
                    line=dict(width=VARIANT_LINE_WIDTH_PX, color=col_del, simplify=False),
                    name=f"SV DEL H{hap + 1} ({n_sv_del})",
                    hoverinfo="text",
                    text=lt_d,
                    customdata=cd_sv_del,
                    **leg_sv,
                )
            )

        ins_rows = [
            r for r in sv_rows if r.svtype == "INS" and sv_affects_hap(r, hap)
        ]
        n_sv_ins = len(ins_rows)
        if ins_rows:
            leg_si: dict[str, Any] = {"legendgroup": "variants", "legendrank": 1}
            if variants_legend_title:
                leg_si["legendgrouptitle"] = dict(text="Variants")
                variants_legend_title = False
            show_leg = True
            for r in ins_rows:
                xa = ref_left(r.pos)
                span_bp = max(INS_PAD_MIN_BP, int(r.svlen))
                xb = xa + float(span_bp)
                hov = sv_hover(r, hap)
                det = sv_detail(r, hap)
                iy_lo, iy_hi = y_mid - ins_y_half, y_mid + ins_y_half
                rx_i, ry_i, _ = _expand_multiseg_rects(
                    [(xa, xb)], iy_lo, iy_hi, [hov], ms_plot, anchor_left=True
                )
                raw_ins = sv_ins_seq_for_hap(r, hap)
                # Symbolic <INS> / missing ALT sequence: still draw bp-wide strip as N so
                # zoom + click-toggle match TRGT (true bases when concrete REF/ALT).
                seq_strip = raw_ins if raw_ins else None
                if not seq_strip and int(r.svlen) > 0:
                    seq_strip = "N" * min(int(r.svlen), INSSEQ_MAX_BP)
                gid_ins = ""
                if seq_strip:
                    gid_ins = insseq_group_id(
                        chrom=chrom, pos=r.pos, hap=hap, kind="SV_INS", vid=r.vid
                    )
                    sy_hi = iy_lo - INSSEQ_GAP_BELOW_MAIN
                    sy_lo = sy_hi - INSSEQ_STRIP_H
                    _add_insseq_bar_trace(
                        fig, seq_strip, float(ref_left(r.pos)), sy_lo, sy_hi, gid_ins
                    )
                # Ensure every clickable point carries full sidebar detail.
                # Plotly includes a trailing separator point; keep customdata non-empty there too.
                cd_ins = [[det, gid_ins]] * 6
                fig.add_trace(
                    go.Scatter(
                        x=rx_i,
                        y=ry_i,
                        mode="lines",
                        hoveron="fills",
                        fill="toself",
                        fillcolor="rgba(22, 163, 74, 0.8)",
                        # Non-zero line width helps Plotly register clicks on filled polygons
                        # (width=0 INS rects often ignored by FigureWidget click handlers).
                        line=dict(width=2, color="rgba(22, 163, 74, 0.45)"),
                        name=f"SV INS H{hap + 1} ({n_sv_ins})",
                        hoverinfo="text",
                        text=hov,
                        hovertext=hov,
                        hovertemplate="%{hovertext}<extra></extra>",
                        customdata=cd_ins,
                        showlegend=show_leg,
                        **leg_si,
                    )
                )
                show_leg = False

    # --- UCSC annotation tracks (after variants so legend groups stay ordered) ---
    if ucsc_tracks:
        ucsc_specs = [
            ("ucsc_segdup", seg_u, "#db2777", "SegDup (UCSC)"),
            ("ucsc_rmsk", rmsk_u, "#ca8a04", "RepeatMasker"),
            ("ucsc_simple", simple_u, "#0d9488", "simpleRepeat"),
        ]
        ucsc_legend_title = True
        for ykey, rows, color, legname in ucsc_specs:
            yc = Y[ykey]
            y_lo = yc - TRACK_UCSC / 2
            y_hi = yc + TRACK_UCSC / 2
            y_mid = 0.5 * (y_lo + y_hi)
            segs_uc: list[tuple[float, float]] = []
            hovs_uc: list[str] = []
            n_uc = 0
            for u in rows:
                xa = ref_left(u.start)
                xb = ref_left(u.end + 1)
                if xb <= xa:
                    continue
                n_uc += 1
                segs_uc.append((xa, xb))
                hovs_uc.append(u.hover)
            lx, ly, lt = _expand_multiseg_line(
                segs_uc, y_mid, hovs_uc, ms_plot, anchor_left=True
            )
            if segs_uc:
                leg_uc: dict[str, Any] = {"legendgroup": "ucsc", "legendrank": 2}
                if ucsc_legend_title:
                    leg_uc["legendgrouptitle"] = dict(text="UCSC")
                    ucsc_legend_title = False
                fig.add_trace(
                    go.Scatter(
                        x=lx,
                        y=ly,
                        mode="lines",
                        line=dict(width=9, color=color, simplify=False),
                        opacity=0.6,
                        name=f"{legname} ({n_uc})",
                        hoverinfo="text",
                        text=lt,
                        **leg_uc,
                    )
                )

    # Reference line: one circular marker per distinct anchor locus (SV + TRGT only).
    # Hover should be locus-only (no variant details); variant tracks keep their own tooltips.
    ref_marker_pos: set[int] = set()
    for hap in (0, 1):
        for r in trgt_rows:
            if allele_on_haplotype(r.gt, hap) is None:
                continue
            if not trgt_hap_meets_min_allele_delta(r, hap, min_allele_bp):
                continue
            xa_t = ref_left(r.pos)
            xb_t = ref_left(r.end + 1)
            if xb_t <= xa_t:
                continue
            ref_marker_pos.add(r.pos)
    for hap in (0, 1):
        for r in sv_rows:
            if r.svtype not in ("DEL", "INS") or not sv_affects_hap(r, hap):
                continue
            ref_marker_pos.add(r.pos)
    if ref_marker_pos:
        sorted_pos = sorted(ref_marker_pos)
        xs_vm = [ref_left(p) for p in sorted_pos]
        texts_vm = [f"{chrom}:{p:,}" for p in sorted_pos]
        fig.add_trace(
            go.Scatter(
                x=xs_vm,
                y=[1.0] * len(xs_vm),
                mode="markers",
                marker=dict(
                    size=10,
                    color="#475569",
                    symbol="circle",
                    line=dict(width=1, color="#ffffff"),
                ),
                text=texts_vm,
                hovertemplate="%{text}<extra></extra>",
                showlegend=False,
                name="_ref_variant_markers",
            )
        )

    # X tick marks at reference coordinates (~8 ticks)
    tick_refs = nice_reference_ticks(win_start, win_end, max_ticks=9)
    tick_x = [ref_left(r) for r in tick_refs]
    tick_text = [f"{r:,}" for r in tick_refs]

    n_trgt_vis = trgt_sites_visible_count(trgt_rows, min_allele_bp)
    status_bits = [
        f"SV: {len(sv_rows)} sites",
        f"TRGT: {n_trgt_vis} sites",
        "0/0 excluded",
        (
            f"TRGT min |alt - ref|: {min_allele_bp} bp"
            if min_allele_bp > 0
            else "TRGT min |alt - ref|: off (0)"
        ),
    ]
    if min_allele_bp > 0 and len(trgt_rows) > n_trgt_vis:
        status_bits.append(
            f"{len(trgt_rows) - n_trgt_vis} TRGT sites below threshold"
        )
    status_line = " | ".join(status_bits)

    fig.update_layout(
        barmode="overlay",
        title=dict(
            text=f"<b>{sample}</b> ({chrom}:{win_start:,}-{win_end:,})",
            x=0.02,
            xanchor="left",
            y=0.99,
            yanchor="top",
            pad=dict(t=2, b=0),
        ),
        height=500,
        margin=dict(l=52, r=24, t=64, b=168),
        xaxis=dict(
            type="linear",
            title=dict(text="Reference position", standoff=10),
            tickvals=tick_x,
            ticktext=tick_text,
            zeroline=False,
            range=[ref_left(win_start) - 0.02 * ref_span, ref_left(win_end + 1) + 0.02 * ref_span],
        ),
        yaxis=dict(
            range=[FIG_Y_RANGE[0], FIG_Y_RANGE[1]],
            fixedrange=True,  # zoom/pan/scroll only affect x (genomic) axis
            tickmode="array",
            tickvals=[
                Y["ucsc_segdup"],
                Y["ucsc_rmsk"],
                Y["ucsc_simple"],
                1.0,
                Y["sv_h2"],
                Y["sv_h1"],
                Y["trgt_h2"],
                Y["trgt_h1"],
            ],
            ticktext=[
                "SegDup",
                "RMask",
                "simRep",
                "ref",
                "SV H2",
                "SV H1",
                "TRGT H2",
                "TRGT H1",
            ],
            title="",
        ),
        showlegend=False,
        legend=dict(
            orientation="v",
            yanchor="top",
            y=1,
            x=1.02,
            xanchor="left",
            xref="paper",
            yref="paper",
            font=dict(size=10),
            bgcolor="rgba(255,255,255,0.94)",
            bordercolor="#e2e8f0",
            borderwidth=1,
            tracegroupgap=24,
        ),
        hoverlabel=dict(align="left"),
    )

    # yref="paper" is the *plot domain* (0=bottom of axes box). yanchor="bottom" at y=0
    # draws text upward into the chart; use negative y + top anchor to sit in the margin
    # below ticks and the x-axis title.
    fig.add_annotation(
        x=-0.15,
        y=-0.50,
        xref="paper",
        yref="paper",
        xanchor="left",
        yanchor="top",
        showarrow=False,
        align="left",
        text=status_line,
        font=dict(size=11, color="#64748b"),
    )

    return fig
