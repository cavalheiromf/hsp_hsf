import gzip
import tempfile
import unittest
from pathlib import Path

from scripts.build_tx2gene import extract_mapping, write_mapping


class Tx2GeneTest(unittest.TestCase):
    def test_preserves_transcript_prefix_and_removes_gene_prefix(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            gff = root / "annotation.gff3.gz"
            with gzip.open(gff, "wt", encoding="utf-8") as handle:
                handle.write("##gff-version 3\n")
                handle.write("1\ttest\tmRNA\t1\t10\t.\t+\t.\tID=transcript:t1;Parent=gene:g1\n")
                handle.write("1\ttest\tmRNA\t20\t30\t.\t+\t.\tID=transcript:t2;Parent=gene:g1\n")
            fasta = root / "transcriptome.fa"
            fasta.write_text(">transcript:t1 CDS=1-9\nAAA\n>transcript:t2 CDS=1-9\nCCC\n", encoding="utf-8")
            self.assertEqual(extract_mapping(gff, fasta), [("transcript:t1", "g1"), ("transcript:t2", "g1")])

    def test_rejects_unmapped_indexed_transcript(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            gff = root / "annotation.gff3"
            gff.write_text(
                "1\ttest\tmRNA\t1\t10\t.\t+\t.\tID=transcript:t1;Parent=gene:g1\n",
                encoding="utf-8",
            )
            fasta = root / "transcriptome.fa"
            fasta.write_text(">transcript:t1\nAAA\n>transcript:t2\nCCC\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Indexed transcripts missing from GFF3: transcript:t2"):
                extract_mapping(gff, fasta)

    def test_maps_non_mrna_transcript_features(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            gff = root / "annotation.gff3"
            gff.write_text(
                "1\ttest\ttRNA\t1\t10\t.\t+\t.\tID=transcript:t1;Parent=gene:g1\n",
                encoding="utf-8",
            )
            fasta = root / "transcriptome.fa"
            fasta.write_text(">transcript:t1\nAAA\n", encoding="utf-8")
            self.assertEqual(extract_mapping(gff, fasta), [("transcript:t1", "g1")])

    def test_writes_lf_terminated_tsv(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "tx2gene.tsv"
            write_mapping([("transcript:t1", "g1")], output)
            self.assertEqual(
                output.read_bytes(),
                b"transcript_id\tgene_id\ntranscript:t1\tg1\n",
            )

    def test_rejects_transcript_with_multiple_parent_genes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            gff = root / "annotation.gff3"
            gff.write_text(
                "1\ttest\tmRNA\t1\t10\t.\t+\t.\tID=transcript:t1;Parent=gene:g1,gene:g2\n",
                encoding="utf-8",
            )
            fasta = root / "transcriptome.fa"
            fasta.write_text(">transcript:t1\nAAA\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "multiple parent genes"):
                extract_mapping(gff, fasta)


if __name__ == "__main__":
    unittest.main()
