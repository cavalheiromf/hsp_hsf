# Salmon RNA-seq Quantification Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Quantify 58 paired-end RNA-seq samples with the correct species-specific Salmon index, validate four canaries before releasing the full array, and create transcript- and gene-level matrices with `tximport`.

**Architecture:** Python utilities normalize sample metadata, derive complete species-specific transcript-to-gene mappings, and validate Salmon outputs. One idempotent Slurm array routes each sample to its index and writes atomically to a per-sample directory; an R/`tximport` stage builds separate matrices for each species.

**Tech Stack:** Python 3.12 standard library, `unittest`, Bash/Slurm, Salmon 1.11.4, R 4.5.1, Bioconductor `tximport`, Ensembl Plants release 63.

## Global Constraints

- Quantify all libraries as paired-end and infer library type with Salmon `-l A`.
- Run Salmon with `--validateMappings`, `--seqBias`, and `--gcBias`.
- Use the four existing transcriptome indices; do not rebuild decoy-aware indices unless canary review requires it.
- Do not trim reads or generate inferential replicates in this stage.
- Use exactly four canaries: SRR33083344, SRR26553587, SRR39669459, and SRR14935704.
- Preserve every per-sample `quant.sf` and Salmon auxiliary metadata.
- Generate `tx2gene` from the same release-63 GFF3 used for each indexed transcriptome.
- Use `tximport(..., countsFromAbundance = "no")`; retain counts, TPM, and effective lengths.
- Build matrices independently by species; never combine features from different species into one matrix.
- Warn below 70% mapping and block automatic release below 50% mapping.
- Stop after validated matrices; differential-expression analysis is out of scope.
- The current workspace exposes an unusable `.git` directory. Run each listed commit after Git metadata is restored; do not replace commits with destructive repository initialization.

---

## File Structure

- `scripts/build_rnaseq_sample_manifest.py`: normalize the curated CSV, validate sample inputs, and select one manifest row for a Slurm array task.
- `scripts/build_tx2gene.py`: extract complete transcript-to-gene relations from release-63 GFF3 and verify them against indexed transcript FASTA identifiers.
- `scripts/validate_salmon_quant.py`: validate one Salmon output directory and consolidate sample QC metrics.
- `jobs/salmon_quant.sbatch`: perform idempotent, atomic per-sample Salmon quantification.
- `scripts/import_salmon_tximport.R`: create transcript- and gene-level matrices independently for each species.
- `tests/test_build_rnaseq_sample_manifest.py`: manifest schema, canary, replicate, routing, and missing-input tests.
- `tests/test_build_tx2gene.py`: GFF3 parsing, prefix preservation, and transcript coverage tests.
- `tests/test_validate_salmon_quant.py`: output-schema and mapping-threshold tests.
- `tests/test_salmon_quant_job.py`: local fake-Salmon integration test for the Slurm script.
- `tests/test_import_salmon_tximport.R`: `tximport` fixture test for counts, TPM, lengths, and gene aggregation.
- `config/rnaseq_samples.tsv`: generated canonical manifest for all 58 samples.
- `results/rnaseq/tx2gene/<species>.tsv`: generated complete mapping for each species.
- `results/rnaseq/quant/<species>/<sample>/`: production Salmon output.
- `results/rnaseq/qc/*.tsv`: canary and full-run QC tables.
- `results/rnaseq/matrices/<species>/`: final per-species matrices and aligned sample metadata.

---

### Task 1: Canonical RNA-seq sample manifest

**Files:**
- Create: `scripts/build_rnaseq_sample_manifest.py`
- Create: `tests/test_build_rnaseq_sample_manifest.py`
- Create on execution: `config/rnaseq_samples.tsv`

**Interfaces:**
- Consumes: curated CSV columns `Bioproject`, `SRA`, `Species`, `Tissue`, `Condition`, and `Layout`; project paths `data/fastq/` and `results/rnaseq/<species>/index/`.
- Produces: `build_rows(metadata_csv: Path, project_root: Path, expected_samples: int | None = 58) -> list[dict[str, str]]`, `write_manifest(rows: list[dict[str, str]], output: Path) -> None`, and `select_row(manifest: Path, scope: str, index: int) -> dict[str, str]`.

- [ ] **Step 1: Write failing manifest tests**

Create `tests/test_build_rnaseq_sample_manifest.py`:

```python
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
```

- [ ] **Step 2: Run the tests and verify the missing module failure**

Run:

```bash
python3 -m unittest tests.test_build_rnaseq_sample_manifest -v
```

Expected: `ModuleNotFoundError: No module named 'scripts.build_rnaseq_sample_manifest'`.

- [ ] **Step 3: Implement the manifest builder and row selector**

Create `scripts/build_rnaseq_sample_manifest.py`:

```python
#!/usr/bin/env python3
"""Build and query the canonical RNA-seq sample manifest."""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path


FIELDS = [
    "sample_id", "species", "bioproject", "tissue", "condition", "replicate",
    "layout", "r1", "r2", "salmon_index", "canary",
]
SPECIES = {
    "Setaria viridis": "setaria_viridis",
    "Sorghum bicolor": "sorghum_bicolor",
    "Triticum aestivum": "triticum_aestivum",
    "Glycine max": "glycine_max",
}
CANARIES = {"SRR33083344", "SRR26553587", "SRR39669459", "SRR14935704"}


def extract_replicate(condition: str) -> str:
    patterns = (
        r"(?:biological\s+)?rep(?:licate)?[_ ]*(\d+)\s*$",
        r"(?:treatment|stress)[ _]*(\d+)\s*$",
        r"(\d+)\s*$",
    )
    for pattern in patterns:
        match = re.search(pattern, condition, flags=re.IGNORECASE)
        if match:
            return match.group(1)
    raise ValueError(f"Cannot derive replicate from condition: {condition}")


def build_rows(
    metadata_csv: Path,
    project_root: Path,
    expected_samples: int | None = 58,
) -> list[dict[str, str]]:
    with metadata_csv.open(newline="", encoding="utf-8-sig") as handle:
        source = list(csv.DictReader(handle))
    rows: list[dict[str, str]] = []
    seen: set[str] = set()
    for item in source:
        sample = item["SRA"].strip()
        if sample in seen:
            raise ValueError(f"Duplicate sample: {sample}")
        seen.add(sample)
        if item["Layout"].strip().lower() != "paired":
            raise ValueError(f"Non-paired layout for {sample}: {item['Layout']}")
        try:
            species = SPECIES[item["Species"].strip()]
        except KeyError as error:
            raise ValueError(f"Unknown species for {sample}: {item['Species']}") from error
        r1 = Path("data/fastq") / f"{sample}_1.fastq.gz"
        r2 = Path("data/fastq") / f"{sample}_2.fastq.gz"
        index = Path("results/rnaseq") / species / "index"
        for label, path in (("R1", r1), ("R2", r2)):
            absolute = project_root / path
            if not absolute.is_file() or absolute.stat().st_size == 0:
                raise ValueError(f"Missing {label} for {sample}: {path}")
        if not (project_root / index / "versionInfo.json").is_file():
            raise ValueError(f"Missing Salmon index for {sample}: {index}")
        condition = item["Condition"].strip()
        rows.append({
            "sample_id": sample,
            "species": species,
            "bioproject": item["Bioproject"].strip(),
            "tissue": item["Tissue"].strip(),
            "condition": condition,
            "replicate": extract_replicate(condition),
            "layout": "paired",
            "r1": str(r1),
            "r2": str(r2),
            "salmon_index": str(index),
            "canary": "true" if sample in CANARIES else "false",
        })
    if expected_samples is not None and len(rows) != expected_samples:
        raise ValueError(f"Expected {expected_samples} samples, found {len(rows)}")
    sample_ids = {row["sample_id"] for row in rows}
    observed_canaries = {row["sample_id"] for row in rows if row["canary"] == "true"}
    required_canaries = CANARIES if expected_samples == 58 else CANARIES & sample_ids
    if observed_canaries != required_canaries:
        raise ValueError("Manifest does not contain exactly the four approved canaries")
    return rows


def write_manifest(rows: list[dict[str, str]], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def select_row(manifest: Path, scope: str, index: int) -> dict[str, str]:
    with manifest.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if scope == "canary":
        rows = [row for row in rows if row["canary"] == "true"]
    elif scope != "all":
        raise ValueError(f"Unknown scope: {scope}")
    try:
        return rows[index]
    except IndexError as error:
        raise ValueError(f"Array index {index} is outside scope {scope} ({len(rows)} rows)") from error


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata-csv", type=Path)
    parser.add_argument("--project-root", type=Path, default=Path("."))
    parser.add_argument("--output", type=Path, default=Path("config/rnaseq_samples.tsv"))
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--scope", choices=("canary", "all"), default="all")
    parser.add_argument("--select-index", type=int)
    args = parser.parse_args()
    if args.select_index is not None:
        if args.manifest is None:
            parser.error("--manifest is required with --select-index")
        print(json.dumps(select_row(args.manifest, args.scope, args.select_index)))
        return
    if args.metadata_csv is None:
        parser.error("--metadata-csv is required when building a manifest")
    rows = build_rows(args.metadata_csv, args.project_root)
    write_manifest(rows, args.output)
    print(f"Wrote {len(rows)} samples to {args.output}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run the focused and complete Python test suites**

Run:

```bash
python3 -m unittest tests.test_build_rnaseq_sample_manifest -v
python3 -m unittest discover -s tests -v
```

Expected: four new manifest tests pass; all previously existing tests remain green.

- [ ] **Step 5: Generate and inspect the production manifest**

Run:

```bash
python3 scripts/build_rnaseq_sample_manifest.py \
  --metadata-csv 'reference/external/SRA list_WD_eCO2 NCBI - Página1.csv' \
  --project-root . \
  --output config/rnaseq_samples.tsv
awk -F '\t' 'NR>1 {species[$2]++; canary[$11]++} END {for (s in species) print s, species[s]; for (c in canary) print "canary=" c, canary[c]}' config/rnaseq_samples.tsv
```

Expected: 58 rows; species totals 18 Setaria, 8 sorghum, 18 wheat, and 14 soybean; `canary=true 4`.

- [ ] **Step 6: Commit the independently testable manifest component**

```bash
git add scripts/build_rnaseq_sample_manifest.py tests/test_build_rnaseq_sample_manifest.py config/rnaseq_samples.tsv
git commit -m "feat: add canonical RNA-seq sample manifest"
```

Expected after Git metadata restoration: one commit containing only the manifest component.

---

### Task 2: Complete species-specific transcript-to-gene mappings

**Files:**
- Create: `scripts/build_tx2gene.py`
- Create: `tests/test_build_tx2gene.py`
- Create on execution: `results/rnaseq/tx2gene/<species>.tsv`

**Interfaces:**
- Consumes: GFF3 `mRNA` features and the first token of each indexed `transcriptome.fa` header.
- Produces: `extract_mapping(gff3: Path, transcriptome: Path) -> list[tuple[str, str]]` where transcript IDs exactly match Salmon feature names and gene IDs omit only the leading `gene:` namespace.

- [ ] **Step 1: Write failing GFF3/FASTA agreement tests**

Create `tests/test_build_tx2gene.py`:

```python
import gzip
import tempfile
import unittest
from pathlib import Path

from scripts.build_tx2gene import extract_mapping


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
```

- [ ] **Step 2: Run the tests and verify the missing module failure**

Run: `python3 -m unittest tests.test_build_tx2gene -v`

Expected: `ModuleNotFoundError: No module named 'scripts.build_tx2gene'`.

- [ ] **Step 3: Implement complete mapping extraction and CLI generation**

Create `scripts/build_tx2gene.py`:

```python
#!/usr/bin/env python3
"""Build complete transcript-to-gene mappings for Salmon/tximport."""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path
from typing import TextIO


def open_text(path: Path) -> TextIO:
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def parse_attributes(text: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value
    return result


def fasta_ids(path: Path) -> list[str]:
    identifiers: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                identifiers.append(line[1:].split()[0])
    if not identifiers:
        raise ValueError(f"No transcript records in {path}")
    if len(identifiers) != len(set(identifiers)):
        raise ValueError(f"Duplicate transcript identifiers in {path}")
    return identifiers


def extract_mapping(gff3: Path, transcriptome: Path) -> list[tuple[str, str]]:
    mapping: dict[str, str] = {}
    with open_text(gff3) as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"Malformed GFF3 row at {gff3}:{line_number}")
            if fields[2] != "mRNA":
                continue
            attributes = parse_attributes(fields[8])
            transcript = attributes.get("ID")
            parents = attributes.get("Parent", "").split(",")
            if not transcript or not parents[0]:
                raise ValueError(f"mRNA lacks ID or Parent at {gff3}:{line_number}")
            if len(parents) != 1:
                raise ValueError(f"Transcript {transcript} has multiple parent genes")
            gene = parents[0].removeprefix("gene:")
            previous = mapping.setdefault(transcript, gene)
            if previous != gene:
                raise ValueError(f"Transcript {transcript} maps to both {previous} and {gene}")

    indexed = fasta_ids(transcriptome)
    missing = [transcript for transcript in indexed if transcript not in mapping]
    extra = sorted(set(mapping) - set(indexed))
    if missing:
        preview = ", ".join(missing[:5])
        raise ValueError(f"Indexed transcripts missing from GFF3: {preview}")
    if extra:
        raise ValueError(f"GFF3 transcripts absent from indexed FASTA: {', '.join(extra[:5])}")
    return [(transcript, mapping[transcript]) for transcript in indexed]


def read_species(path: Path) -> list[str]:
    with path.open(encoding="utf-8") as handle:
        return [row["species_id"] for row in csv.DictReader(handle, delimiter="\t")]


def write_mapping(rows: list[tuple[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(("transcript_id", "gene_id"))
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species-config", type=Path, default=Path("config/species.tsv"))
    parser.add_argument("--reference-root", type=Path, default=Path("data/reference/ensembl_plants_63"))
    parser.add_argument("--rnaseq-root", type=Path, default=Path("results/rnaseq"))
    parser.add_argument("--output-root", type=Path, default=Path("results/rnaseq/tx2gene"))
    args = parser.parse_args()
    for species in read_species(args.species_config):
        rows = extract_mapping(
            args.reference_root / species / "annotation.gff3.gz",
            args.rnaseq_root / species / "transcriptome.fa",
        )
        output = args.output_root / f"{species}.tsv"
        write_mapping(rows, output)
        print(f"{species}: {len(rows)} transcripts -> {output}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run focused and complete tests**

Run:

```bash
python3 -m unittest tests.test_build_tx2gene -v
python3 -m unittest discover -s tests -v
```

Expected: three mapping tests pass and no regression appears.

- [ ] **Step 5: Generate production mappings and compare transcript counts**

Run:

```bash
python3 scripts/build_tx2gene.py
for species in setaria_viridis sorghum_bicolor triticum_aestivum glycine_max; do
  printf '%s\t' "$species"
  awk 'END {print NR-1}' "results/rnaseq/tx2gene/${species}.tsv"
  cat "results/rnaseq/${species}/transcript_count.txt"
done
```

Expected: each mapping row count equals its indexed transcript count: 52,459; 48,559; 146,597; and 89,662, respectively.

- [ ] **Step 6: Commit the mapping component**

```bash
git add scripts/build_tx2gene.py tests/test_build_tx2gene.py results/rnaseq/tx2gene
git commit -m "feat: build complete Salmon transcript-to-gene mappings"
```

---

### Task 3: Salmon output validation and QC consolidation

**Files:**
- Create: `scripts/validate_salmon_quant.py`
- Create: `tests/test_validate_salmon_quant.py`

**Interfaces:**
- Consumes: `quant.sf`, `aux_info/meta_info.json`, and `lib_format_counts.json` from one Salmon output directory.
- Produces: `validate_quant(directory: Path) -> dict[str, str | int | float]` and a CLI that validates `--single` or writes a manifest-aligned QC TSV.

- [ ] **Step 1: Write failing validation and threshold tests**

Create `tests/test_validate_salmon_quant.py`:

```python
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
```

- [ ] **Step 2: Run tests and verify the missing module failure**

Run: `python3 -m unittest tests.test_validate_salmon_quant -v`

Expected: missing-module failure.

- [ ] **Step 3: Implement strict validation and QC collection**

Create `scripts/validate_salmon_quant.py`:

```python
#!/usr/bin/env python3
"""Validate Salmon outputs and consolidate mapping diagnostics."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


QC_FIELDS = [
    "sample_id", "species", "status", "mapping_flag", "num_processed", "num_mapped",
    "percent_mapped", "library_types", "salmon_version", "index_seq_hash", "index_name_hash",
    "elapsed_seconds",
]


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"Missing or empty required file: {path}")


def validate_quant(directory: Path) -> dict[str, str | int | float]:
    quant = directory / "quant.sf"
    meta_path = directory / "aux_info/meta_info.json"
    lib_path = directory / "lib_format_counts.json"
    run_path = directory / "run_metrics.json"
    for path in (quant, meta_path, lib_path, run_path):
        require_file(path)
    with quant.open(encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        first = handle.readline()
    if header != ["Name", "Length", "EffectiveLength", "TPM", "NumReads"] or not first:
        raise ValueError(f"Malformed quant.sf: {quant}")
    meta = json.loads(meta_path.read_text(encoding="utf-8"))
    json.loads(lib_path.read_text(encoding="utf-8"))
    run = json.loads(run_path.read_text(encoding="utf-8"))
    required = ("num_processed", "num_mapped", "percent_mapped", "salmon_version")
    missing = [key for key in required if key not in meta]
    if missing:
        raise ValueError(f"meta_info.json lacks: {', '.join(missing)}")
    processed = int(meta["num_processed"])
    mapped = int(meta["num_mapped"])
    percent = float(meta["percent_mapped"])
    if processed <= 0 or mapped < 0 or mapped > processed:
        raise ValueError(f"Invalid fragment counts in {meta_path}")
    flag = "block" if percent < 50.0 else "warning" if percent < 70.0 else "pass"
    library_types = meta.get("library_types", [])
    return {
        "status": "pass",
        "mapping_flag": flag,
        "num_processed": processed,
        "num_mapped": mapped,
        "percent_mapped": percent,
        "library_types": ",".join(library_types) if isinstance(library_types, list) else str(library_types),
        "salmon_version": str(meta["salmon_version"]),
        "index_seq_hash": str(meta.get("index_seq_hash", "")),
        "index_name_hash": str(meta.get("index_name_hash", "")),
        "elapsed_seconds": int(run["elapsed_seconds"]),
    }


def collect_qc(manifest: Path, quant_root: Path, scope: str) -> list[dict[str, str | int | float]]:
    with manifest.open(encoding="utf-8") as handle:
        samples = list(csv.DictReader(handle, delimiter="\t"))
    if scope == "canary":
        samples = [sample for sample in samples if sample["canary"] == "true"]
    rows = []
    for sample in samples:
        metrics = validate_quant(quant_root / sample["species"] / sample["sample_id"])
        rows.append({"sample_id": sample["sample_id"], "species": sample["species"], **metrics})
    return rows


def write_qc(rows: list[dict[str, str | int | float]], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=QC_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--single", type=Path)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--quant-root", type=Path, default=Path("results/rnaseq/quant"))
    parser.add_argument("--scope", choices=("canary", "all"), default="all")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.single:
        print(json.dumps(validate_quant(args.single), sort_keys=True))
        return
    if not args.manifest or not args.output:
        parser.error("--manifest and --output are required unless --single is used")
    rows = collect_qc(args.manifest, args.quant_root, args.scope)
    write_qc(rows, args.output)
    print(f"Wrote {len(rows)} QC rows to {args.output}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run focused and complete tests**

```bash
python3 -m unittest tests.test_validate_salmon_quant -v
python3 -m unittest discover -s tests -v
```

Expected: four validator tests and all earlier tests pass.

- [ ] **Step 5: Commit the validator**

```bash
git add scripts/validate_salmon_quant.py tests/test_validate_salmon_quant.py
git commit -m "feat: validate Salmon quantification outputs"
```

---

### Task 4: Idempotent Salmon Slurm array

**Files:**
- Create: `jobs/salmon_quant.sbatch`
- Create: `tests/test_salmon_quant_job.py`

**Interfaces:**
- Consumes: JSON emitted by `build_rnaseq_sample_manifest.py --manifest ... --scope ... --select-index ...` and `validate_salmon_quant.py --single`.
- Produces: atomically finalized `results/rnaseq/quant/<species>/<sample>/` directories.

- [ ] **Step 1: Write a failing fake-Salmon integration test**

Create `tests/test_salmon_quant_job.py`:

```python
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
            "(out / 'lib_format_counts.json').write_text(json.dumps({'expected_format':'IU'}))\n",
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
        self.assertTrue((self.quant_root / "setaria_viridis/SRR39669459/quant.sf").is_file())
        self.assertFalse(any(self.quant_root.rglob("*.partial.*")))

    def test_valid_completed_sample_is_idempotent(self):
        first = self.run_job()
        self.assertEqual(first.returncode, 0, first.stderr)
        second = self.run_job()
        self.assertEqual(second.returncode, 0, second.stderr)
        self.assertIn("already complete; skipping", second.stdout)


if __name__ == "__main__":
    unittest.main()
```

- [ ] **Step 2: Run the test and verify the missing job failure**

Run: `python3 -m unittest tests.test_salmon_quant_job -v`

Expected: failure because `jobs/salmon_quant.sbatch` does not exist.

- [ ] **Step 3: Implement the complete Slurm job**

Create `jobs/salmon_quant.sbatch`:

```bash
#!/usr/bin/env bash
#SBATCH --job-name=salmon-quant
#SBATCH --partition=short
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --chdir=/home/mfcaval/github/hsp_hsf
#SBATCH --output=/home/mfcaval/github/hsp_hsf/logs/salmon_quant_%A_%a.out
#SBATCH --error=/home/mfcaval/github/hsp_hsf/logs/salmon_quant_%A_%a.err

set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-/home/mfcaval/github/hsp_hsf}"
MANIFEST="${MANIFEST:-$PROJECT_DIR/config/rnaseq_samples.tsv}"
QUANT_ROOT="${QUANT_ROOT:-$PROJECT_DIR/results/rnaseq/quant}"
QUANT_SCOPE="${QUANT_SCOPE:-canary}"
MIN_FREE_GB="${MIN_FREE_GB:-5}"
task_index="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"
threads="${SLURM_CPUS_PER_TASK:-16}"
job_id="${SLURM_JOB_ID:-manual}"

row_json="$(python3 "$PROJECT_DIR/scripts/build_rnaseq_sample_manifest.py" \
  --manifest "$MANIFEST" --scope "$QUANT_SCOPE" --select-index "$task_index")"
mapfile -t fields < <(python3 -c '
import json, sys
row = json.loads(sys.argv[1])
for key in ("sample_id", "species", "r1", "r2", "salmon_index"):
    print(row[key])
' "$row_json")
sample="${fields[0]}"
species="${fields[1]}"
r1="$PROJECT_DIR/${fields[2]}"
r2="$PROJECT_DIR/${fields[3]}"
index="$PROJECT_DIR/${fields[4]}"
final="$QUANT_ROOT/$species/$sample"
partial="${final}.partial.${job_id}"

echo "Started: $(date --iso-8601=seconds)"
start_epoch="$(date +%s)"
echo "Sample: $sample"
echo "Species: $species"
echo "Scope: $QUANT_SCOPE"
echo "Index: $index"

if [[ -d "$final" ]]; then
    if python3 "$PROJECT_DIR/scripts/validate_salmon_quant.py" --single "$final" >/dev/null; then
        echo "$sample already complete; skipping"
        exit 0
    fi
    echo "Existing final output is invalid and will not be overwritten: $final" >&2
    exit 20
fi
if [[ -e "$partial" ]]; then
    echo "Partial output already exists and requires review: $partial" >&2
    exit 21
fi
for path in "$r1" "$r2" "$index/versionInfo.json"; do
    [[ -s "$path" ]] || { echo "Missing or empty input: $path" >&2; exit 22; }
done
gzip -t "$r1"
gzip -t "$r2"
mkdir -p "$(dirname "$final")"
free_kb="$(df -Pk "$(dirname "$final")" | awk 'NR==2 {print $4}')"
required_kb="$((MIN_FREE_GB * 1024 * 1024))"
if (( free_kb < required_kb )); then
    echo "Insufficient free space: ${free_kb} KiB available, ${required_kb} KiB required" >&2
    exit 23
fi

if [[ -z "${SALMON_BIN:-}" ]]; then
    module load Bio/Salmon/1.11.4
    SALMON_BIN="$(command -v salmon)"
fi
command=(
    "$SALMON_BIN" quant
    -i "$index"
    -l A
    -1 "$r1"
    -2 "$r2"
    --validateMappings
    --seqBias
    --gcBias
    -p "$threads"
    -o "$partial"
)
printf 'Command:'
printf ' %q' "${command[@]}"
printf '\n'
"${command[@]}"
end_epoch="$(date +%s)"
printf '{"elapsed_seconds": %d}\n' "$((end_epoch - start_epoch))" > "$partial/run_metrics.json"
python3 "$PROJECT_DIR/scripts/validate_salmon_quant.py" --single "$partial"
mv "$partial" "$final"
echo "Completed: $(date --iso-8601=seconds)"
```

- [ ] **Step 4: Complete the fake-Salmon fixture and run syntax/integration checks**

Run:

```bash
bash -n jobs/salmon_quant.sbatch
python3 -m unittest tests.test_salmon_quant_job -v
python3 -m unittest discover -s tests -v
```

Expected: Bash syntax succeeds; initial and idempotent rerun tests pass; all Python tests remain green.

- [ ] **Step 5: Commit the Slurm execution layer**

```bash
git add jobs/salmon_quant.sbatch tests/test_salmon_quant_job.py
git commit -m "feat: add idempotent Salmon quantification array"
```

---

### Task 5: Species-specific matrix generation with tximport

**Files:**
- Create: `scripts/import_salmon_tximport.R`
- Create: `tests/test_import_salmon_tximport.R`

**Interfaces:**
- Consumes: manifest rows, per-sample `quant.sf`, and `results/rnaseq/tx2gene/<species>.tsv`.
- Produces: six TSV files per species: transcript counts, transcript TPM, gene counts, gene TPM, gene effective lengths, and aligned sample metadata.

- [ ] **Step 1: Write a failing R fixture test**

Create `tests/test_import_salmon_tximport.R`:

```r
source("scripts/import_salmon_tximport.R")

root <- tempfile("tximport-fixture-")
dir.create(root)
quant_root <- file.path(root, "quant")
output_dir <- file.path(root, "matrices")
tx2gene_path <- file.path(root, "tx2gene.tsv")

write_quant <- function(sample, counts, tpm) {
  directory <- file.path(quant_root, "setaria_viridis", sample)
  dir.create(directory, recursive = TRUE)
  table <- data.frame(
    Name = c("transcript:t1", "transcript:t2"),
    Length = c(100, 200),
    EffectiveLength = c(80, 180),
    TPM = tpm,
    NumReads = counts,
    check.names = FALSE
  )
  write.table(
    table, file.path(directory, "quant.sf"), sep = "\t", quote = FALSE,
    row.names = FALSE
  )
}

write_quant("sample1", c(10, 20), c(400000, 600000))
write_quant("sample2", c(5, 15), c(300000, 700000))
write.table(
  data.frame(
    transcript_id = c("transcript:t1", "transcript:t2"),
    gene_id = c("g1", "g1")
  ),
  tx2gene_path, sep = "\t", quote = FALSE, row.names = FALSE
)
manifest <- data.frame(
  sample_id = c("sample1", "sample2"),
  species = c("setaria_viridis", "setaria_viridis"),
  canary = c("true", "false"),
  stringsAsFactors = FALSE
)

import_species(manifest, quant_root, tx2gene_path, output_dir)

stopifnot(file.exists(file.path(output_dir, "transcript_counts.tsv")))
stopifnot(file.exists(file.path(output_dir, "gene_counts.tsv")))
gene_counts <- read.delim(file.path(output_dir, "gene_counts.tsv"), check.names = FALSE)
stopifnot(identical(gene_counts$feature_id, "g1"))
stopifnot(abs(gene_counts$sample1 - 30) < 1e-8)
stopifnot(abs(gene_counts$sample2 - 20) < 1e-8)
metadata <- read.delim(file.path(output_dir, "sample_metadata.tsv"), check.names = FALSE)
stopifnot(identical(metadata$sample_id, c("sample1", "sample2")))
unlink(root, recursive = TRUE)
```

- [ ] **Step 2: Run the R test and verify the missing script failure**

Run:

```bash
Rscript -e 'stopifnot(requireNamespace("tximport", quietly=TRUE))'
Rscript tests/test_import_salmon_tximport.R
```

Expected: the dependency check passes and the test fails because `scripts/import_salmon_tximport.R` does not exist.

- [ ] **Step 3: Implement per-species tximport aggregation**

Create `scripts/import_salmon_tximport.R`:

```r
#!/usr/bin/env Rscript

write_matrix <- function(matrix, path) {
  output <- data.frame(feature_id = rownames(matrix), matrix, check.names = FALSE)
  write.table(output, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

import_species <- function(sample_rows, quant_root, tx2gene_path, output_dir) {
  if (!requireNamespace("tximport", quietly = TRUE)) {
    stop("Bioconductor package 'tximport' is required")
  }
  files <- file.path(quant_root, sample_rows$species, sample_rows$sample_id, "quant.sf")
  names(files) <- sample_rows$sample_id
  missing <- files[!file.exists(files)]
  if (length(missing) > 0) {
    stop("Missing quant.sf: ", paste(missing, collapse = ", "))
  }
  tx2gene <- read.delim(tx2gene_path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!identical(names(tx2gene), c("transcript_id", "gene_id"))) {
    stop("tx2gene must contain transcript_id and gene_id columns")
  }
  for (file in files) {
    quantified <- read.delim(file, stringsAsFactors = FALSE, check.names = FALSE)$Name
    missing_mapping <- setdiff(quantified, tx2gene$transcript_id)
    extra_mapping <- setdiff(tx2gene$transcript_id, quantified)
    if (length(missing_mapping) > 0 || length(extra_mapping) > 0) {
      stop(
        "Transcript/tx2gene mismatch for ", file,
        ": missing=", length(missing_mapping), ", extra=", length(extra_mapping)
      )
    }
  }
  tx <- tximport::tximport(
    files, type = "salmon", txOut = TRUE, countsFromAbundance = "no"
  )
  gene <- tximport::tximport(
    files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "no"
  )
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write_matrix(tx$counts, file.path(output_dir, "transcript_counts.tsv"))
  write_matrix(tx$abundance, file.path(output_dir, "transcript_tpm.tsv"))
  write_matrix(gene$counts, file.path(output_dir, "gene_counts.tsv"))
  write_matrix(gene$abundance, file.path(output_dir, "gene_tpm.tsv"))
  write_matrix(gene$length, file.path(output_dir, "effective_lengths.tsv"))
  write.table(
    sample_rows,
    file.path(output_dir, "sample_metadata.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  invisible(list(transcript = tx, gene = gene))
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 5) {
    stop("Usage: import_salmon_tximport.R MANIFEST QUANT_ROOT TX2GENE_ROOT OUTPUT_ROOT SCOPE")
  }
  manifest <- read.delim(args[[1]], stringsAsFactors = FALSE, check.names = FALSE)
  quant_root <- args[[2]]
  tx2gene_root <- args[[3]]
  output_root <- args[[4]]
  scope <- args[[5]]
  if (!scope %in% c("canary", "all")) stop("Scope must be canary or all")
  if (scope == "canary") manifest <- manifest[manifest$canary == "true", , drop = FALSE]
  for (species in unique(manifest$species)) {
    rows <- manifest[manifest$species == species, , drop = FALSE]
    import_species(
      rows,
      quant_root,
      file.path(tx2gene_root, paste0(species, ".tsv")),
      file.path(output_root, species)
    )
  }
}

if (sys.nframe() == 0) main()
```

- [ ] **Step 4: Run the R fixture and all Python tests**

```bash
Rscript tests/test_import_salmon_tximport.R
python3 -m unittest discover -s tests -v
```

Expected: the R fixture exits zero and all Python tests pass.

- [ ] **Step 5: Commit matrix generation**

```bash
git add scripts/import_salmon_tximport.R tests/test_import_salmon_tximport.R
git commit -m "feat: aggregate Salmon estimates with tximport"
```

---

### Task 6: Production preflight and four-canary submission

**Files:**
- Verify: `config/rnaseq_samples.tsv`
- Verify: `results/rnaseq/tx2gene/*.tsv`
- Create externally through Slurm: four `results/rnaseq/quant/<species>/<sample>/` directories and four log pairs.

**Interfaces:**
- Consumes: all components from Tasks 1–5.
- Produces: four real, validated Salmon quantifications.

- [ ] **Step 1: Run the complete local validation suite**

```bash
python3 -m unittest discover -s tests -v
Rscript tests/test_import_salmon_tximport.R
bash -n jobs/salmon_quant.sbatch
```

Expected: all commands exit zero.

- [ ] **Step 2: Rebuild deterministic generated inputs**

```bash
python3 scripts/build_rnaseq_sample_manifest.py \
  --metadata-csv 'reference/external/SRA list_WD_eCO2 NCBI - Página1.csv' \
  --project-root . \
  --output config/rnaseq_samples.tsv
python3 scripts/build_tx2gene.py
```

Expected: 58 manifest rows and four mapping files with the exact transcript counts recorded in Task 2.

- [ ] **Step 3: Submit only the four canaries**

```bash
mkdir -p results/rnaseq/qc
canary_job_id="$(sbatch --parsable \
  --array=0-3%4 \
  --export=ALL,QUANT_SCOPE=canary \
  jobs/salmon_quant.sbatch)"
printf '%s\n' "$canary_job_id" | tee results/rnaseq/qc/canary_slurm_job_id.txt
```

Expected: Slurm returns one job ID with four array tasks. Record that ID in the run notes; do not submit the full array.

- [ ] **Step 4: Wait for all four tasks and verify scheduler completion**

```bash
sacct -j "$(cat results/rnaseq/qc/canary_slurm_job_id.txt)" \
  --format=JobID,JobName,State,ExitCode,Elapsed,MaxRSS
```

Expected: all four array tasks are `COMPLETED` with `0:0` exit codes.

- [ ] **Step 5: Build the canary QC report and matrices**

```bash
python3 scripts/validate_salmon_quant.py \
  --manifest config/rnaseq_samples.tsv \
  --quant-root results/rnaseq/quant \
  --scope canary \
  --output results/rnaseq/qc/canary_quantification.tsv
Rscript scripts/import_salmon_tximport.R \
  config/rnaseq_samples.tsv \
  results/rnaseq/quant \
  results/rnaseq/tx2gene \
  results/rnaseq/matrices/canary \
  canary
```

Expected: four QC rows and one-sample transcript/gene matrices under each of the four species directories.

- [ ] **Step 6: Enforce the manual release gate**

Run:

```bash
column -t -s $'\t' results/rnaseq/qc/canary_quantification.tsv
```

Expected: no `mapping_flag=block`; investigate every `warning`. Present the table to the user and obtain explicit approval before Task 7.

- [ ] **Step 7: Commit canary provenance after approval**

```bash
git add results/rnaseq/qc/canary_quantification.tsv results/rnaseq/matrices/canary logs/salmon_quant_*.out logs/salmon_quant_*.err
git commit -m "data: record Salmon canary validation"
```

Do not add FASTQs, indices, or other large sequence data.

---

### Task 7: Full quantification and final per-species matrices

**Files:**
- Create externally through Slurm: remaining 54 per-sample quantification directories and logs.
- Create: `results/rnaseq/qc/quantification.tsv`
- Create: `results/rnaseq/matrices/<species>/*.tsv`

**Interfaces:**
- Consumes: explicit canary approval and the unchanged manifest/index/mapping versions used for canaries.
- Produces: validated full-run QC and final transcript/gene matrices for each species.

- [ ] **Step 1: Submit the complete manifest idempotently**

```bash
full_job_id="$(sbatch --parsable \
  --array=0-57%4 \
  --export=ALL,QUANT_SCOPE=all \
  jobs/salmon_quant.sbatch)"
printf '%s\n' "$full_job_id" | tee results/rnaseq/qc/full_slurm_job_id.txt
```

Expected: four canaries skip as already valid and the other 54 samples quantify. The `%4` concurrency cap limits shared filesystem pressure.

- [ ] **Step 2: Monitor without mutating outputs**

```bash
squeue -j "$(cat results/rnaseq/qc/full_slurm_job_id.txt)"
sacct -j "$(cat results/rnaseq/qc/full_slurm_job_id.txt)" \
  --format=JobID,State,ExitCode,Elapsed,MaxRSS
```

Expected terminal state: every task `COMPLETED` with exit code `0:0`. Resubmit only failed array indices after diagnosing their preserved partial outputs.

- [ ] **Step 3: Consolidate and inspect full-run QC**

```bash
python3 scripts/validate_salmon_quant.py \
  --manifest config/rnaseq_samples.tsv \
  --quant-root results/rnaseq/quant \
  --scope all \
  --output results/rnaseq/qc/quantification.tsv
awk -F '\t' 'NR>1 {flag[$4]++; species[$2]++} END {for (x in flag) print x, flag[x]; for (x in species) print x, species[x]}' \
  results/rnaseq/qc/quantification.tsv
```

Expected: 58 QC rows and zero `block` rows. Review all warnings before matrix publication.

- [ ] **Step 4: Generate final matrices independently for every species**

```bash
Rscript scripts/import_salmon_tximport.R \
  config/rnaseq_samples.tsv \
  results/rnaseq/quant \
  results/rnaseq/tx2gene \
  results/rnaseq/matrices \
  all
```

Expected matrix column counts, excluding `feature_id`: 18 Setaria, 8 sorghum, 18 wheat, and 14 soybean.

- [ ] **Step 5: Validate final feature and sample dimensions**

```bash
for species in setaria_viridis sorghum_bicolor triticum_aestivum glycine_max; do
  echo "$species"
  wc -l "results/rnaseq/matrices/$species/"*.tsv
  head -1 "results/rnaseq/matrices/$species/gene_counts.tsv"
done
```

Expected: transcript matrix rows match each species' `tx2gene` row count plus one header; all six files exist per species; metadata and matrix sample ordering agree.

- [ ] **Step 6: Run the final regression suite**

```bash
python3 -m unittest discover -s tests -v
Rscript tests/test_import_salmon_tximport.R
```

Expected: all tests pass after production data generation.

- [ ] **Step 7: Commit final lightweight outputs and implementation**

```bash
git add config/rnaseq_samples.tsv scripts jobs tests results/rnaseq/tx2gene results/rnaseq/qc results/rnaseq/matrices
git commit -m "feat: complete four-species Salmon quantification"
```

Before committing, run `git status --short` and verify that FASTQs, `.sra` files, transcriptome FASTA files, Salmon indices, and per-sample Salmon working directories remain excluded.
