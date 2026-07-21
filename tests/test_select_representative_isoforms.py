import tempfile
import unittest
from pathlib import Path

from scripts.select_representative_isoforms import (
    canonical_transcripts,
    parse_fasta,
    select_representatives,
)


class RepresentativeIsoformTests(unittest.TestCase):
    def test_canonical_then_longest_then_stable_id(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta = root / "proteins.fa"
            gff3 = root / "annotation.gff3"
            fasta.write_text(
                ">p1 gene:g1 transcript:t1\nAAAA\n"
                ">p2 gene:g1 transcript:t2\nAAAAAAAA\n"
                ">p4 gene:g2 transcript:t4\nAAAAAA\n"
                ">p3 gene:g2 transcript:t3\nAAAAAA\n"
                ">p5 gene:g3 transcript:t5\nAAA\n"
                ">p6 gene:g3 transcript:t6\nAAAAA\n",
                encoding="utf-8",
            )
            gff3.write_text(
                "##gff-version 3\n"
                "1\ttest\tmRNA\t1\t10\t.\t+\t.\tID=transcript:t1;Parent=gene:g1;tag=Ensembl_canonical\n",
                encoding="utf-8",
            )

            proteins = parse_fasta(fasta)
            selected = list(select_representatives(proteins, canonical_transcripts(gff3)))
            observed = {protein.gene_id: (protein.protein_id, rule) for protein, rule, _ in selected}

            self.assertEqual(observed["g1"], ("p1", "canonical"))
            self.assertEqual(observed["g2"], ("p3", "length_tie_id"))
            self.assertEqual(observed["g3"], ("p6", "longest"))

    def test_missing_header_mapping_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            fasta = Path(tmp) / "bad.fa"
            fasta.write_text(">p1\nAAAA\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                parse_fasta(fasta)


if __name__ == "__main__":
    unittest.main()
