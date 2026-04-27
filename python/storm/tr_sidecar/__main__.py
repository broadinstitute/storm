from __future__ import annotations

import argparse
import sys

from .builder import build_tr_sidecar


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description="Build TR annotation sidecar (Parquet) from an SV VCF/BCF and catalog TSV."
    )
    p.add_argument("--vcf", required=True, help="Input VCF/BCF path")
    p.add_argument("--catalog", required=True, help="tr_loci_hg38.tsv (or compatible)")
    p.add_argument("--out", required=True, help="Output directory")
    p.add_argument(
        "--with-entries",
        action="store_true",
        help="Also write entries.parquet (per-sample GT + optional FORMAT JSON)",
    )
    p.add_argument(
        "--entry-format-fields",
        type=str,
        default="",
        help=(
            "Comma-separated FORMAT names to store per (variant,sample) in entries.parquet "
            "as JSON (e.g. CN for copy number; caller-specific repeat fields). "
            "Only applies with --with-entries."
        ),
    )
    p.add_argument(
        "--batch-size",
        type=int,
        default=50_000,
        help="Records to buffer before writing a Parquet row group",
    )
    args = p.parse_args(argv)
    fmt = tuple(x.strip() for x in args.entry_format_fields.split(",") if x.strip())
    res = build_tr_sidecar(
        args.vcf,
        args.catalog,
        args.out,
        with_entries=args.with_entries,
        batch_size=args.batch_size,
        entry_format_fields=fmt,
    )
    print(
        f"Wrote {res.sites_path} ({res.n_variants} variants)",
        file=sys.stderr,
    )
    if res.entries_path:
        print(f"Wrote {res.entries_path}", file=sys.stderr)
    print(f"Manifest {res.manifest_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
