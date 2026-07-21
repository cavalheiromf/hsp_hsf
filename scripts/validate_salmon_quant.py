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
