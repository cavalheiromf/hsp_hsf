import csv
import gzip
import os
import subprocess
import tempfile
import unittest
from pathlib import Path

from scripts.build_rnaseq_sample_manifest import FIELDS


REPOSITORY = Path(__file__).parents[1]


class SalmonQuantJobTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.project = Path(self.temp.name)
        (self.project / "scripts").symlink_to(REPOSITORY / "scripts", target_is_directory=True)
        fastq = self.project / "data/fastq"
        fastq.mkdir(parents=True)
        for mate in (1, 2):
            with gzip.open(fastq / f"SRR39669459_{mate}.fastq.gz", "wt", encoding="utf-8") as handle:
                handle.write("@read\nACGT\n+\n!!!!\n")
        index = self.project / "results/rnaseq/setaria_viridis/index"
        index.mkdir(parents=True)
        (index / "versionInfo.json").write_text("{}\n", encoding="utf-8")
        self.manifest = self.project / "config/rnaseq_samples.tsv"
        self.manifest.parent.mkdir()
        row = {
            "sample_id": "SRR39669459", "species": "setaria_viridis", "bioproject": "PRJ",
            "tissue": "Leaf", "condition": "control rep2", "replicate": "2", "layout": "paired",
            "r1": "data/fastq/SRR39669459_1.fastq.gz",
            "r2": "data/fastq/SRR39669459_2.fastq.gz",
            "salmon_index": "results/rnaseq/setaria_viridis/index", "canary": "true",
        }
        with self.manifest.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
            writer.writeheader()
            writer.writerow(row)
        self.fake_salmon = self.project / "fake_salmon.py"
        self.fake_salmon.write_text(
            "#!/usr/bin/env python3\n"
            "import json, pathlib, sys\n"
            "out = pathlib.Path(sys.argv[sys.argv.index('-o') + 1])\n"
            "(out / 'aux_info').mkdir(parents=True)\n"
            "(out / 'quant.sf').write_text('Name\\tLength\\tEffectiveLength\\tTPM\\tNumReads\\ntranscript:t1\\t100\\t80\\t1000000\\t80\\n')\n"
            "(out / 'aux_info/meta_info.json').write_text(json.dumps({'salmon_version':'1.11.4','index_seq_hash':'seq','index_name_hash':'name','num_processed':100,'num_mapped':80,'percent_mapped':80.0,'library_types':['IU']}))\n"
            "(out / 'lib_format_counts.json').write_text(json.dumps({'expected_format':'IU','compatible_fragment_ratio':0.98}))\n",
            encoding="utf-8",
        )
        self.fake_salmon.chmod(0o755)
        self.quant_root = self.project / "results/rnaseq/quant"

    def tearDown(self):
        self.temp.cleanup()

    def run_job(self) -> subprocess.CompletedProcess[str]:
        environment = {
            **os.environ,
            "PROJECT_DIR": str(self.project),
            "MANIFEST": str(self.manifest),
            "QUANT_ROOT": str(self.quant_root),
            "QUANT_SCOPE": "canary",
            "SLURM_ARRAY_TASK_ID": "0",
            "SLURM_JOB_ID": "fixture",
            "SLURM_CPUS_PER_TASK": "2",
            "SALMON_BIN": str(self.fake_salmon),
            "MIN_FREE_GB": "0",
        }
        return subprocess.run(
            ["bash", str(REPOSITORY / "jobs/salmon_quant.sbatch")],
            cwd=REPOSITORY,
            env=environment,
            text=True,
            capture_output=True,
        )

    def test_quantifies_and_atomically_finalizes_sample(self):
        result = self.run_job()
        self.assertEqual(result.returncode, 0, result.stderr)
        final = self.quant_root / "setaria_viridis/SRR39669459"
        self.assertTrue((final / "quant.sf").is_file())
        self.assertTrue((final / "run_metrics.json").is_file())
        self.assertFalse(any(self.quant_root.rglob("*.partial.*")))

    def test_valid_completed_sample_is_idempotent(self):
        first = self.run_job()
        self.assertEqual(first.returncode, 0, first.stderr)
        second = self.run_job()
        self.assertEqual(second.returncode, 0, second.stderr)
        self.assertIn("already complete; skipping", second.stdout)


if __name__ == "__main__":
    unittest.main()
