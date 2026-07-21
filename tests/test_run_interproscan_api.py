import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "scripts" / "run_interproscan_api.py"
SPEC = importlib.util.spec_from_file_location("run_interproscan_api", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class InterProScanApiHelpersTest(unittest.TestCase):
    def test_fasta_round_trip_and_chunking(self):
        with tempfile.TemporaryDirectory() as directory:
            fasta = Path(directory) / "candidates.fa"
            fasta.write_text(">p1 description\nAAAA\n>p2\nCC\n>p3\nGGG\n", encoding="utf-8")
            records = MODULE.parse_fasta(fasta)

        self.assertEqual(records[0], ("p1 description", "AAAA"))
        self.assertEqual([len(item) for item in MODULE.chunks(records, 2)], [2, 1])
        self.assertEqual(MODULE.format_fasta(records).count(">"), 3)

    def test_empty_sequence_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            fasta = Path(directory) / "empty.fa"
            fasta.write_text(">p1\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                MODULE.parse_fasta(fasta)


if __name__ == "__main__":
    unittest.main()
