from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pyarrow.parquet as pq

# Run from repo root: python -m pytest tests/test_tr_sidecar.py
# or: cd python && python -m pytest ../tests/test_tr_sidecar.py
import sys

ROOT = Path(__file__).resolve().parents[1]
PY = ROOT / "python"
if str(PY) not in sys.path:
    sys.path.insert(0, str(PY))

from storm.tr_sidecar import build_tr_sidecar  # noqa: E402


MINI_CATALOG = """tr_locus_id	chrom	start	end	gene	motif_primary	motifs_accepted	rule_group	inheritance	benign_max_repeats	pathogenic_min_repeats	motif_change_required	contraction_pathogenic	prioritize_in_sv_scan	sv_scan_default_applicable	notes
TEST_chr22	chr22	24990000	25100000	TESTGENE	T		T	test	AD	10	20	0	0	1	1	test_locus
"""


def _write_vcf(path: Path, pos: int = 25_003_401, chrom: str = "chr22") -> None:
    # 1-based POS matches common VCF; pysam start = POS-1
    lines = [
        "##fileformat=VCFv4.2",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        f"##contig=<ID={chrom},length=50818468>",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1",
        f"{chrom}\t{pos}\ttr-test-1\tN\t<INS>\t.\t.\tSVTYPE=INS;SVLEN=22;MOTIF=T;MOTIF_SIZE=1;END=25003422\tGT\t0/1",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


class TrSidecarTest(unittest.TestCase):
    def test_build_sidecar_match_and_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            cat = td_path / "cat.tsv"
            cat.write_text(MINI_CATALOG, encoding="utf-8")
            vcf = td_path / "t.vcf"
            _write_vcf(vcf)
            out = td_path / "out"
            res = build_tr_sidecar(vcf, cat, out, with_entries=True)
            self.assertEqual(res.n_variants, 1)
            st = pq.read_table(res.sites_path)
            self.assertEqual(st.num_rows, 1)
            row = st.to_pylist()[0]
            self.assertEqual(row["variant_id"], "tr-test-1")
            self.assertEqual(row["tr_locus_id"], "TEST_chr22")
            self.assertEqual(row["gene"], "TESTGENE")
            self.assertEqual(row["catalog_overlap_bp"], 22)
            self.assertTrue(row["motif_matches_catalog"])
            self.assertTrue(row["rule_applicable"])
            self.assertEqual(row["repeat_units_estimate"], 22.0)

            et = pq.read_table(res.entries_path)
            self.assertEqual(et.num_rows, 1)
            er = et.to_pylist()[0]
            self.assertEqual(er["sample"], "s1")
            self.assertTrue(er["is_carrier"])

            man = json.loads(res.manifest_path.read_text(encoding="utf-8"))
            self.assertEqual(man["n_variants_processed"], 1)
            self.assertTrue(man["with_entries"])
            self.assertIn("catalog_sha256", man)

    def test_repeat_units_estimate_uses_first_svlen_when_tuple(self) -> None:
        """Number=A SVLEN is a tuple in pysam; estimate must not be dropped."""
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            cat = td_path / "cat.tsv"
            cat.write_text(MINI_CATALOG, encoding="utf-8")
            vcf = td_path / "t.vcf"
            lines = [
                "##fileformat=VCFv4.2",
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                "##contig=<ID=chr22,length=50818468>",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1",
                "chr22\t25003401\ttr-test-2\tN\t<DEL>,<INS>\t.\t.\t"
                "SVTYPE=DEL;SVLEN=22,30;MOTIF=T;MOTIF_SIZE=1;END=25003422\tGT\t0/1",
                "",
            ]
            vcf.write_text("\n".join(lines), encoding="utf-8")
            out = td_path / "out"
            res = build_tr_sidecar(vcf, cat, out, with_entries=False)
            row = pq.read_table(res.sites_path).to_pylist()[0]
            self.assertEqual(row["repeat_units_estimate"], 22.0)

    def test_variant_id_prefers_info_id_over_column_id(self) -> None:
        """``annotate_svs`` join key uses INFO/ID; sidecar must match, not the short VCF ID column."""
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            cat = td_path / "cat.tsv"
            cat.write_text(MINI_CATALOG, encoding="utf-8")
            vcf = td_path / "t.vcf"
            # Column 3 is a short id; long descriptive id is in INFO/ID (like many union VCFs).
            lines = [
                "##fileformat=VCFv4.2",
                '##INFO=<ID=ID,Number=1,Type=String,Description="Variant ID">',
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                f"##contig=<ID=chr22,length=50818468>",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1",
                "chr22\t25003401\t22_12345\tN\t<INS>\t.\t.\t"
                "ID=chr22-25003401-INS-22;SVTYPE=INS;SVLEN=22;MOTIF=T;MOTIF_SIZE=1;END=25003422"
                "\tGT\t0/1",
                "",
            ]
            vcf.write_text("\n".join(lines), encoding="utf-8")
            out = td_path / "out"
            res = build_tr_sidecar(vcf, cat, out, with_entries=False)
            row = pq.read_table(res.sites_path).to_pylist()[0]
            self.assertEqual(row["variant_id"], "chr22-25003401-INS-22")

    def test_no_overlap(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            cat = td_path / "cat.tsv"
            cat.write_text(MINI_CATALOG, encoding="utf-8")
            vcf = td_path / "t.vcf"
            _write_vcf(vcf, pos=1_000_000)
            out = td_path / "out"
            res = build_tr_sidecar(vcf, cat, out, with_entries=False)
            row = pq.read_table(res.sites_path).to_pylist()[0]
            self.assertEqual(row["tr_locus_id"], "")
            self.assertEqual(row["match_method"], "no_overlap")


if __name__ == "__main__":
    unittest.main()
