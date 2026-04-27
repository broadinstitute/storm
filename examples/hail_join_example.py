#!/usr/bin/env python3
"""
Hail workflow: import SV VCF, run ``storm.annotate_svs`` (split multi, allele IDs,
optional TR sidecar join on ``allele_id`` = ``variant_id``).

Requires Hail **0.2.134** (Terra) or compatible 0.2.x; see
https://hail.is/docs/0.2/getting_started.html. Often installed via conda or your
cluster image rather than ``pip install 'storm[hail]'`` (which pins 0.2.134).

Usage::

    hailctl dataproc submit ... hail_join_example.py -- \\
        --vcf gs://.../union.bcf \\
        --sites gs://.../sidecar/sites.parquet

Notebook::

    import hail as hl
    import storm

    mt0 = hl.import_vcf(vcf_path, force_bgz=True, reference_genome="GRCh38")
    mt1 = storm.annotate_svs(mt0, tr_sidecar_sites="gs://.../sites.parquet")

Legacy (pre-split join on ``rsid`` only) — use only if your VCF is already biallelic
and IDs match ``variant_id`` in the sidecar without splitting::

    mt = hl.import_vcf(...)
    mt = mt.annotate_rows(variant_id=mt.rsid)
    sites = hl.import_parquet(sites_path).key_by("variant_id")
    mt = mt.annotate_rows(tr=sites[mt.variant_id])

Whole-genome ``sites.parquet``: :func:`storm.annotate_svs` uses chunked Parquet reads
when ``hl.import_parquet`` is missing (e.g. Terra 0.2.134), so the driver does not
load the entire sites table into memory.
"""

from __future__ import annotations

import argparse


def run(
    vcf: str,
    sites_parquet: str | None,
    reference_genome: str = "GRCh38",
) -> None:
    import hail as hl

    import storm

    hl.init(default_reference=reference_genome)

    mt0 = hl.import_vcf(vcf, reference_genome=reference_genome, force_bgz=True)
    mt1 = storm.annotate_svs(
        mt0,
        tr_sidecar_sites=sites_parquet,
    )

    mt1.describe()
    # hits = mt1.filter_rows(hl.is_defined(mt1.tr))


def main() -> None:
    p = argparse.ArgumentParser(description="import_vcf + storm.annotate_svs")
    p.add_argument("--vcf", required=True)
    p.add_argument(
        "--sites",
        default=None,
        help="Optional sites.parquet from storm-tr-sidecar (join on allele_id)",
    )
    p.add_argument("--reference-genome", default="GRCh38")
    args = p.parse_args()
    run(args.vcf, args.sites, reference_genome=args.reference_genome)


if __name__ == "__main__":
    main()
