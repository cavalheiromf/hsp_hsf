import json
import math
import tempfile
import unittest
import csv
from pathlib import Path

from scripts.validate_salmon_quant import QC_FIELDS, collect_qc, validate_quant, write_qc


def write_quant(
    root: Path,
    percent_mapped: float = 75.0,
    *,
    meta_overrides: dict | None = None,
    lib_payload: object | None = None,
    quant_body: str = "transcript:t1\t100\t80\t1000000\t10\n",
    elapsed_seconds: object = 12,
) -> None:
    (root / "aux_info").mkdir(parents=True)
    (root / "quant.sf").write_text(
        "Name\tLength\tEffectiveLength\tTPM\tNumReads\n" + quant_body,
        encoding="utf-8",
    )
    mapped = int(percent_mapped) if isinstance(percent_mapped, (int, float)) and math.isfinite(percent_mapped) else 75
    meta = {
        "salmon_version": "1.11.4", "index_seq_hash": "seq", "index_name_hash": "name",
        "num_processed": 100, "num_mapped": mapped,
        "percent_mapped": percent_mapped, "library_types": ["IU"],
    }
    if meta_overrides:
        meta.update(meta_overrides)
    (root / "aux_info/meta_info.json").write_text(json.dumps(meta), encoding="utf-8")
    if lib_payload is None:
        lib_payload = {"expected_format": "IU", "compatible_fragment_ratio": 0.98}
    (root / "lib_format_counts.json").write_text(json.dumps(lib_payload), encoding="utf-8")
    (root / "run_metrics.json").write_text(json.dumps({"elapsed_seconds": elapsed_seconds}), encoding="utf-8")


class SalmonValidationTest(unittest.TestCase):
    def test_valid_quantification_returns_pass_metrics(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            write_quant(root, 75.0)
            metrics = validate_quant(root)
            self.assertEqual(metrics["status"], "pass")
            self.assertEqual(metrics["mapping_flag"], "pass")
            self.assertEqual(metrics["num_processed"], 100)
            self.assertEqual(metrics["expected_format"], "IU")
            self.assertEqual(metrics["compatible_fragment_ratio"], 0.98)

    def test_mapping_below_70_warns(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            write_quant(root, 65.0)
            self.assertEqual(validate_quant(root)["mapping_flag"], "warning")

    def test_mapping_below_50_blocks(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            write_quant(root, 49.0)
            self.assertEqual(validate_quant(root)["mapping_flag"], "block")

    def test_missing_quant_sf_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "aux_info").mkdir()
            with self.assertRaisesRegex(ValueError, "quant.sf"):
                validate_quant(root)

    def test_rejects_nonfinite_out_of_range_or_inconsistent_percent_mapped(self):
        for percent, mapped in ((float("nan"), 75), (float("inf"), 75), (-0.1, 0), (100.1, 100), (80.0, 75)):
            with self.subTest(percent=percent, mapped=mapped), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                write_quant(root, percent, meta_overrides={"num_mapped": mapped})
                with self.assertRaisesRegex(ValueError, "percent_mapped"):
                    validate_quant(root)

    def test_requires_complete_typed_meta_info_fields(self):
        invalid_overrides = (
            {"library_types": []}, {"library_types": ["IU", 1]},
            {"index_seq_hash": ""}, {"index_name_hash": 3},
        )
        for overrides in invalid_overrides:
            with self.subTest(overrides=overrides), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                write_quant(root, meta_overrides=overrides)
                with self.assertRaises(ValueError):
                    validate_quant(root)

    def test_requires_valid_library_format_metrics(self):
        invalid_payloads = (
            [], {"expected_format": ""}, {"expected_format": "IU"},
            {"expected_format": "IU", "compatible_fragment_ratio": float("nan")},
            {"expected_format": "IU", "compatible_fragment_ratio": 1.1},
        )
        for payload in invalid_payloads:
            with self.subTest(payload=payload), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                write_quant(root, lib_payload=payload)
                with self.assertRaises(ValueError):
                    validate_quant(root)

    def test_rejects_malformed_or_nonfinite_quantification_rows(self):
        invalid_rows = (
            "transcript:t1\t100\t80\t1\n",
            "transcript:t1\t0\t80\t1\t1\n",
            "transcript:t1\t100\t0\t1\t1\n",
            "transcript:t1\t100\t80\t-1\t1\n",
            "transcript:t1\t100\t80\t1\t-1\n",
            "transcript:t1\t100\t80\tnan\t1\n",
        )
        for body in invalid_rows:
            with self.subTest(body=body), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                write_quant(root, quant_body=body)
                with self.assertRaisesRegex(ValueError, "quant.sf"):
                    validate_quant(root)

    def test_requires_nonnegative_integer_elapsed_seconds(self):
        for elapsed in (-1, 1.5, "12"):
            with self.subTest(elapsed=elapsed), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                write_quant(root, elapsed_seconds=elapsed)
                with self.assertRaisesRegex(ValueError, "elapsed_seconds"):
                    validate_quant(root)

    def test_collect_qc_filters_canaries_and_writes_aligned_tsv(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = root / "manifest.tsv"
            manifest.write_text(
                "sample_id\tspecies\tcanary\nS1\tsp_a\ttrue\nS2\tsp_b\tfalse\n",
                encoding="utf-8",
            )
            write_quant(root / "quant" / "sp_a" / "S1")
            write_quant(root / "quant" / "sp_b" / "S2")
            rows = collect_qc(manifest, root / "quant", "canary")
            output = root / "qc.tsv"
            write_qc(rows, output)
            with output.open(encoding="utf-8") as handle:
                written = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([row["sample_id"] for row in written], ["S1"])
            self.assertEqual(list(written[0]), QC_FIELDS)
            self.assertEqual(written[0]["expected_format"], "IU")


if __name__ == "__main__":
    unittest.main()
