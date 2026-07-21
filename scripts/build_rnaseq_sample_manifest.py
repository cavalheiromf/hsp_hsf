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
        writer = csv.DictWriter(
            handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def select_row(manifest: Path, scope: str, index: int) -> dict[str, str]:
    with manifest.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if scope == "canary":
        rows = [row for row in rows if row["canary"] == "true"]
    elif scope != "all":
        raise ValueError(f"Unknown scope: {scope}")
    if index < 0:
        raise ValueError(f"Array index {index} is outside scope {scope} ({len(rows)} rows)")
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
