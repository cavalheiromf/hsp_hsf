import csv
import gzip
import os
import subprocess
import tempfile
import time
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
            "import json, os, pathlib, sys, time\n"
            "args = sys.argv[1:]\n"
            "expected = ['quant', '-i', None, '-l', 'A', '-1', None, '-2', None, '--validateMappings', '--seqBias', '--gcBias', '-p', None, '-o', None]\n"
            "if len(args) != len(expected) or any(want is not None and got != want for got, want in zip(args, expected)) or not all(args[index] for index in (2, 6, 8, 13, 15)):\n"
            "    raise SystemExit('unexpected Salmon command: ' + repr(args))\n"
            "log = os.environ.get('FAKE_SALMON_LOG')\n"
            "if log:\n"
            "    with open(log, 'a', encoding='utf-8') as handle: handle.write('invoked\\n')\n"
            "out = pathlib.Path(args[-1])\n"
            "out.mkdir(parents=True)\n"
            "if os.environ.get('FAKE_SALMON_FAIL'):\n"
            "    raise SystemExit(13)\n"
            "time.sleep(float(os.environ.get('FAKE_SALMON_SLEEP', '0')))\n"
            "(out / 'aux_info').mkdir(parents=True)\n"
            "(out / 'quant.sf').write_text('Name\\tLength\\tEffectiveLength\\tTPM\\tNumReads\\ntranscript:t1\\t100\\t80\\t1000000\\t80\\n')\n"
            "meta = {'salmon_version':'1.11.4','index_seq_hash':'seq','index_name_hash':'name','num_processed':100,'num_mapped':80,'percent_mapped':80.0,'library_types':[] if os.environ.get('FAKE_SALMON_INVALID') else ['IU']}\n"
            "(out / 'aux_info/meta_info.json').write_text(json.dumps(meta))\n"
            "(out / 'lib_format_counts.json').write_text(json.dumps({'expected_format':'IU','compatible_fragment_ratio':0.98}))\n",
            encoding="utf-8",
        )
        self.fake_salmon.chmod(0o755)
        self.quant_root = self.project / "results/rnaseq/quant"

    def tearDown(self):
        self.temp.cleanup()

    def job_environment(self, **overrides: str) -> dict[str, str]:
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
        environment.update(overrides)
        return environment

    def run_job(self, **overrides: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            ["bash", str(REPOSITORY / "jobs/salmon_quant.sbatch")],
            cwd=REPOSITORY,
            env=self.job_environment(**overrides),
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

    def test_invalid_final_blocks_before_invoking_salmon(self):
        final = self.quant_root / "setaria_viridis/SRR39669459"
        final.parent.mkdir(parents=True)
        final.write_text("not a quantification directory\n", encoding="utf-8")
        log = self.project / "salmon.log"

        result = self.run_job(FAKE_SALMON_LOG=str(log))

        self.assertEqual(result.returncode, 20, result.stderr)
        self.assertFalse(log.exists())

    def test_old_partial_blocks_before_invoking_salmon(self):
        partial = self.quant_root / "setaria_viridis/SRR39669459.partial.old-job"
        partial.mkdir(parents=True)
        log = self.project / "salmon.log"

        result = self.run_job(FAKE_SALMON_LOG=str(log))

        self.assertEqual(result.returncode, 21, result.stderr)
        self.assertTrue(partial.is_dir())
        self.assertFalse(log.exists())

    def test_salmon_failure_preserves_partial_output(self):
        result = self.run_job(FAKE_SALMON_FAIL="1", SLURM_JOB_ID="failed-salmon")

        self.assertEqual(result.returncode, 13, result.stderr)
        self.assertTrue((self.quant_root / "setaria_viridis/SRR39669459.partial.failed-salmon").is_dir())

    def test_validation_failure_preserves_partial_output(self):
        result = self.run_job(FAKE_SALMON_INVALID="1", SLURM_JOB_ID="failed-validation")

        self.assertNotEqual(result.returncode, 0)
        self.assertTrue((self.quant_root / "setaria_viridis/SRR39669459.partial.failed-validation").is_dir())

    def test_concurrent_jobs_serialize_and_second_skips(self):
        log = self.project / "salmon.log"
        first = subprocess.Popen(
            ["bash", str(REPOSITORY / "jobs/salmon_quant.sbatch")],
            cwd=REPOSITORY,
            env=self.job_environment(SLURM_JOB_ID="first", FAKE_SALMON_LOG=str(log), FAKE_SALMON_SLEEP="1"),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        deadline = time.monotonic() + 5
        while not log.exists() and time.monotonic() < deadline:
            time.sleep(0.01)
        self.assertTrue(log.exists(), "first fake Salmon invocation did not start")
        second = subprocess.Popen(
            ["bash", str(REPOSITORY / "jobs/salmon_quant.sbatch")],
            cwd=REPOSITORY,
            env=self.job_environment(SLURM_JOB_ID="second", FAKE_SALMON_LOG=str(log)),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        first_stdout, first_stderr = first.communicate(timeout=10)
        second_stdout, second_stderr = second.communicate(timeout=10)

        self.assertEqual(first.returncode, 0, first_stderr)
        self.assertEqual(second.returncode, 0, second_stderr)
        self.assertIn("already complete; skipping", second_stdout)
        self.assertEqual(log.read_text(encoding="utf-8").splitlines(), ["invoked"])

    def test_rejects_invalid_resource_parameters_before_salmon(self):
        for overrides in ({"MIN_FREE_GB": "-1"}, {"MIN_FREE_GB": "nope"}, {"SLURM_CPUS_PER_TASK": "0"}, {"SLURM_CPUS_PER_TASK": "nope"}):
            with self.subTest(overrides=overrides):
                log = self.project / f"salmon-{len(overrides)}.log"
                result = self.run_job(**overrides, FAKE_SALMON_LOG=str(log))
                self.assertNotEqual(result.returncode, 0)
                self.assertFalse(log.exists())


if __name__ == "__main__":
    unittest.main()
