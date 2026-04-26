"""Tests for storm.saige_prep helpers (requires hail)."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PY = ROOT / "python"
if str(PY) not in sys.path:
    sys.path.insert(0, str(PY))

import pyarrow as pa  # noqa: E402
import pyarrow.parquet as pq  # noqa: E402

try:
    import hail as hl

    _HAS_HAIL = True
except ImportError:
    hl = None  # type: ignore[assignment]
    _HAS_HAIL = False

import storm  # noqa: E402


def _write_vcf(path: Path) -> None:
    lines = [
        "##fileformat=VCFv4.2",
        "##INFO=<ID=ID,Number=1,Type=String,Description=\"Per-ALT ids\">",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "##contig=<ID=chr1,length=248956422>",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2",
        "chr1\t1000\trsid_a,rsid_b\tG\tA,C\t.\t.\tID=id_a,id_b\tGT\t0/1\t1/1",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def _write_sidecar(path: Path) -> None:
    table = pa.Table.from_pylist(
        [
            {
                "variant_id": "id_a",
                "svtype": "INS",
                "tr_locus_id": "L1",
                "rule_applicable": True,
                "motif": "T",
                "repeat_units_estimate": 5.0,
            },
            {
                "variant_id": "id_b",
                "svtype": "DEL",
                "tr_locus_id": "",
                "rule_applicable": False,
                "motif": "",
                "repeat_units_estimate": None,
            },
        ]
    )
    pq.write_table(table, str(path))


class SaigePrepTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not _HAS_HAIL:
            raise unittest.SkipTest("hail not installed")
        try:
            hl.init_local(quiet=True, default_reference="GRCh38")
        except Exception as e:  # pragma: no cover
            raise unittest.SkipTest(f"Hail init failed: {e}") from e

    @classmethod
    def tearDownClass(cls) -> None:
        if _HAS_HAIL:
            hl.stop()

    def test_inventory_and_long_predictors(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            vcf = td_path / "m.vcf"
            sites = td_path / "sites.parquet"
            _write_vcf(vcf)
            _write_sidecar(sites)

            mt0 = hl.import_vcf(str(vcf), reference_genome="GRCh38", force=True)
            mt1 = storm.annotate_svs(mt0, tr_sidecar_sites=sites)

            fi = storm.build_feature_inventory(mt1)
            self.assertEqual(fi.count(), 2)
            rows = {r.feature_id: r for r in fi.collect()}
            self.assertEqual(rows["id_a"].feature_class, "tr_quantitative")
            self.assertEqual(rows["id_b"].feature_class, "standard_sv")

            tables = storm.build_long_predictor_tables(mt1)
            std = tables["standard_long"].collect()
            trq = tables["tr_quant_long"].collect()

            # id_b contributes standard predictors for both samples
            self.assertEqual(len(std), 2)
            # id_a contributes tr_quant predictors for both samples
            self.assertEqual(len(trq), 2)

            by_sample = {r.sample_id: r.predictor for r in trq}
            # s1 GT=0/1 -> dosage 1 * repeat_units_estimate 5.0
            self.assertAlmostEqual(by_sample["s1"], 5.0)
            # s2 GT=1/1 -> dosage 2 * repeat_units_estimate 5.0
            self.assertAlmostEqual(by_sample["s2"], 10.0)

    def test_dense_marker_matrix_and_long_replay(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            vcf = td_path / "m.vcf"
            sites = td_path / "sites.parquet"
            _write_vcf(vcf)
            _write_sidecar(sites)

            mt0 = hl.import_vcf(str(vcf), reference_genome="GRCh38", force=True)
            mt1 = storm.annotate_svs(mt0, tr_sidecar_sites=sites)

            manifest = mt1.cols().select(sample_id=hl.str(mt1.s)).order_by("sample_id")

            m_std = storm.build_dense_saige_marker_matrix(
                mt1, stratum="standard_sv", standard_label="standard_sv", tr_quant_label="tr_quantitative"
            )
            m_tr = storm.build_dense_saige_marker_matrix(
                mt1, stratum="tr_quantitative", standard_label="standard_sv", tr_quant_label="tr_quantitative"
            )
            self.assertEqual(m_std.count_rows(), 1)
            self.assertEqual(m_tr.count_rows(), 1)
            self.assertEqual(m_std.count_cols(), 2)
            self.assertEqual(m_tr.count_cols(), 2)

            m_std_a = storm.align_matrix_cols_to_manifest(m_std, manifest)
            entries = m_std_a.entries().collect()
            self.assertEqual(len(entries), 2)
            by_s = {e.s: e.DS for e in entries}
            self.assertAlmostEqual(by_s["s1"], 1.0)
            self.assertAlmostEqual(by_s["s2"], 2.0)

            lookup = storm.build_feature_vcf_row_lookup(mt1)
            long_tables = storm.build_long_predictor_tables(mt1)
            m_from_long = storm.dense_marker_matrix_from_long(
                long_tables["standard_long"], lookup, fill_value=0.0
            )
            self.assertEqual(m_from_long.count_rows(), 1)
            self.assertEqual(m_from_long.count_cols(), 2)

            out_vcf = td_path / "out.ds.vcf.bgz"
            storm.export_saige_dosage_vcf(m_std_a, str(out_vcf))
            self.assertTrue(out_vcf.is_file())
            back = hl.import_vcf(str(out_vcf), reference_genome="GRCh38", force=True)
            self.assertIn("DS", back.entry)
            self.assertEqual(back.count_rows(), 1)
            row = back.entries().collect()
            self.assertEqual(len(row), 2)
            by_s2 = {r.s: r.DS for r in row}
            self.assertAlmostEqual(by_s2["s1"], 1.0)
            self.assertAlmostEqual(by_s2["s2"], 2.0)

    def test_export_both_strata(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            vcf = td_path / "m.vcf"
            sites = td_path / "sites.parquet"
            _write_vcf(vcf)
            _write_sidecar(sites)
            mt0 = hl.import_vcf(str(vcf), reference_genome="GRCh38", force=True)
            mt1 = storm.annotate_svs(mt0, tr_sidecar_sites=sites)
            manifest = mt1.cols().select(sample_id=hl.str(mt1.s)).order_by("sample_id")
            paths = storm.export_saige_stratum_vcfs(mt1, manifest, out_dir=str(td_path), prefix="t")
            self.assertIn("standard_sv", paths)
            self.assertIn("tr_quantitative", paths)
            self.assertTrue(Path(paths["standard_sv"]).is_file())
            self.assertTrue(Path(paths["tr_quantitative"]).is_file())


if __name__ == "__main__":
    unittest.main()
