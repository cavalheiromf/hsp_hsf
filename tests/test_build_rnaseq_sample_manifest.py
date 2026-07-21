import csv
import tempfile
import unittest
from pathlib import Path

from scripts.build_rnaseq_sample_manifest import build_rows, select_row, write_manifest


class RnaSeqManifestTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        (self.root / "data/fastq").mkdir(parents=True)
        for species in ("triticum_aestivum", "glycine_max"):
            index = self.root / "results/rnaseq" / species / "index"
            index.mkdir(parents=True)
            (index / "versionInfo.json").write_text("{}\n", encoding="utf-8")

        self.csv_path = self.root / "metadata.csv"
        with self.csv_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=["Bioproject", "SRA", "Species", "Tissue", "Condition", "Layout"],
            )
            writer.writeheader()
            writer.writerow({
                "Bioproject": "PRJ1", "SRA": "SRR33083344", "Species": "Triticum aestivum",
                "Tissue": "Leaf", "Condition": "Alana_elevated-CO2+drought_biol_rep1",
                "Layout": "Paired",
            })
            writer.writerow({
                "Bioproject": "PRJ2", "SRA": "SRR26553587", "Species": "Glycine max",
                "Tissue": "Leaf", "Condition": "Controlled water treatment 2",
                "Layout": "Paired",
            })
        for sample in ("SRR33083344", "SRR26553587"):
            for mate in (1, 2):
                (self.root / f"data/fastq/{sample}_{mate}.fastq.gz").write_bytes(b"gzip-fixture")

    def tearDown(self):
        self.temp.cleanup()

    def test_builds_normalized_rows_and_marks_canaries(self):
        rows = build_rows(self.csv_path, self.root, expected_samples=2)
        self.assertEqual([row["species"] for row in rows], ["triticum_aestivum", "glycine_max"])
        self.assertEqual([row["replicate"] for row in rows], ["1", "2"])
        self.assertTrue(all(row["canary"] == "true" for row in rows))
        self.assertEqual(rows[0]["r1"], "data/fastq/SRR33083344_1.fastq.gz")
        self.assertEqual(rows[0]["salmon_index"], "results/rnaseq/triticum_aestivum/index")

    def test_selects_canary_by_zero_based_array_index(self):
        manifest = self.root / "manifest.tsv"
        write_manifest(build_rows(self.csv_path, self.root, expected_samples=2), manifest)
        self.assertEqual(select_row(manifest, "canary", 1)["sample_id"], "SRR26553587")
        self.assertEqual(select_row(manifest, "all", 0)["sample_id"], "SRR33083344")

    def test_writes_lf_terminated_tsv(self):
        manifest = self.root / "manifest.tsv"
        write_manifest(build_rows(self.csv_path, self.root, expected_samples=2), manifest)
        self.assertNotIn(b"\r\n", manifest.read_bytes())

    def test_rejects_negative_array_index(self):
        manifest = self.root / "manifest.tsv"
        write_manifest(build_rows(self.csv_path, self.root, expected_samples=2), manifest)
        with self.assertRaisesRegex(ValueError, "Array index -1 is outside scope all"):
            select_row(manifest, "all", -1)

    def test_rejects_missing_fastq_mate(self):
        (self.root / "data/fastq/SRR26553587_2.fastq.gz").unlink()
        with self.assertRaisesRegex(ValueError, "Missing R2"):
            build_rows(self.csv_path, self.root, expected_samples=2)

    def test_rejects_duplicate_sample(self):
        text = self.csv_path.read_text(encoding="utf-8")
        self.csv_path.write_text(text + text.splitlines()[1] + "\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "Duplicate sample"):
            build_rows(self.csv_path, self.root, expected_samples=2)


if __name__ == "__main__":
    unittest.main()
