import json
import tempfile
import unittest
from pathlib import Path

from scripts.validate_salmon_quant import validate_quant


def write_quant(root: Path, percent_mapped: float = 75.0) -> None:
    (root / "aux_info").mkdir(parents=True)
    (root / "quant.sf").write_text(
        "Name\tLength\tEffectiveLength\tTPM\tNumReads\ntranscript:t1\t100\t80\t1000000\t10\n",
        encoding="utf-8",
    )
    (root / "aux_info/meta_info.json").write_text(json.dumps({
        "salmon_version": "1.11.4", "index_seq_hash": "seq", "index_name_hash": "name",
        "num_processed": 100, "num_mapped": int(percent_mapped),
        "percent_mapped": percent_mapped, "library_types": ["IU"],
    }), encoding="utf-8")
    (root / "lib_format_counts.json").write_text(json.dumps({"expected_format": "IU"}), encoding="utf-8")
    (root / "run_metrics.json").write_text(json.dumps({"elapsed_seconds": 12}), encoding="utf-8")


class SalmonValidationTest(unittest.TestCase):
    def test_valid_quantification_returns_pass_metrics(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            write_quant(root, 75.0)
            metrics = validate_quant(root)
            self.assertEqual(metrics["status"], "pass")
            self.assertEqual(metrics["mapping_flag"], "pass")
            self.assertEqual(metrics["num_processed"], 100)

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


if __name__ == "__main__":
    unittest.main()
