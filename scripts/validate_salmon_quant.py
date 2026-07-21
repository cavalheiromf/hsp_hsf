#!/usr/bin/env python3
"""Validate Salmon outputs and consolidate mapping diagnostics."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


QC_FIELDS = [
    "sample_id", "species", "status", "mapping_flag", "num_processed", "num_mapped",
    "percent_mapped", "library_types", "salmon_version", "index_seq_hash", "index_name_hash",
    "expected_format", "compatible_fragment_ratio", "elapsed_seconds",
]


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"Missing or empty required file: {path}")


def require_finite_number(value: object, field: str, *, minimum: float | None = None,
                          maximum: float | None = None) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"Invalid {field}")
    number = float(value)
    if (minimum is not None and number < minimum) or (maximum is not None and number > maximum):
        raise ValueError(f"Invalid {field}")
    return number


def validate_quant_rows(quant: Path) -> None:
    with quant.open(encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        if header != ["Name", "Length", "EffectiveLength", "TPM", "NumReads"]:
            raise ValueError(f"Malformed quant.sf: {quant}")
        rows = 0
        for line_number, line in enumerate(handle, start=2):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 5:
                raise ValueError(f"Malformed quant.sf: {quant} line {line_number}")
            try:
                length, effective_length, tpm, reads = (float(value) for value in fields[1:])
            except ValueError as error:
                raise ValueError(f"Malformed quant.sf: {quant} line {line_number}") from error
            if not all(math.isfinite(value) for value in (length, effective_length, tpm, reads)):
                raise ValueError(f"Malformed quant.sf: {quant} line {line_number}")
            if length <= 0 or effective_length <= 0 or tpm < 0 or reads < 0:
                raise ValueError(f"Malformed quant.sf: {quant} line {line_number}")
            rows += 1
    if rows == 0:
        raise ValueError(f"Malformed quant.sf: {quant}")


def validate_quant(directory: Path) -> dict[str, str | int | float]:
    quant = directory / "quant.sf"
    meta_path = directory / "aux_info/meta_info.json"
    lib_path = directory / "lib_format_counts.json"
    run_path = directory / "run_metrics.json"
    for path in (quant, meta_path, lib_path, run_path):
        require_file(path)
    validate_quant_rows(quant)
    meta = json.loads(meta_path.read_text(encoding="utf-8"))
    lib = json.loads(lib_path.read_text(encoding="utf-8"))
    run = json.loads(run_path.read_text(encoding="utf-8"))
    if not isinstance(meta, dict):
        raise ValueError(f"Invalid meta_info.json: {meta_path}")
    required = (
        "num_processed", "num_mapped", "percent_mapped", "salmon_version", "library_types",
        "index_seq_hash", "index_name_hash",
    )
    missing = [key for key in required if key not in meta]
    if missing:
        raise ValueError(f"meta_info.json lacks: {', '.join(missing)}")
    processed = meta["num_processed"]
    mapped = meta["num_mapped"]
    if isinstance(processed, bool) or not isinstance(processed, int) or isinstance(mapped, bool) or not isinstance(mapped, int):
        raise ValueError(f"Invalid fragment counts in {meta_path}")
    if processed <= 0 or mapped < 0 or mapped > processed:
        raise ValueError(f"Invalid fragment counts in {meta_path}")
    percent = require_finite_number(meta["percent_mapped"], "percent_mapped", minimum=0.0, maximum=100.0)
    if abs(percent - (100.0 * mapped / processed)) > 0.1:
        raise ValueError(f"Invalid percent_mapped in {meta_path}")
    library_types = meta["library_types"]
    if not isinstance(library_types, list) or not library_types or not all(isinstance(item, str) and item for item in library_types):
        raise ValueError(f"Invalid library_types in {meta_path}")
    for field in ("salmon_version", "index_seq_hash", "index_name_hash"):
        if not isinstance(meta[field], str) or not meta[field]:
            raise ValueError(f"Invalid {field} in {meta_path}")
    if not isinstance(lib, dict):
        raise ValueError(f"Invalid lib_format_counts.json: {lib_path}")
    expected_format = lib.get("expected_format")
    if not isinstance(expected_format, str) or not expected_format:
        raise ValueError(f"Invalid expected_format in {lib_path}")
    ratio = require_finite_number(lib.get("compatible_fragment_ratio"), "compatible_fragment_ratio", minimum=0.0, maximum=1.0)
    if not isinstance(run, dict):
        raise ValueError(f"Invalid run_metrics.json: {run_path}")
    elapsed_seconds = run.get("elapsed_seconds")
    if isinstance(elapsed_seconds, bool) or not isinstance(elapsed_seconds, int) or elapsed_seconds < 0:
        raise ValueError(f"Invalid elapsed_seconds in {run_path}")
    flag = "block" if percent < 50.0 else "warning" if percent < 70.0 else "pass"
    return {
        "status": "pass",
        "mapping_flag": flag,
        "num_processed": processed,
        "num_mapped": mapped,
        "percent_mapped": percent,
        "library_types": ",".join(library_types) if isinstance(library_types, list) else str(library_types),
        "salmon_version": str(meta["salmon_version"]),
        "index_seq_hash": meta["index_seq_hash"],
        "index_name_hash": meta["index_name_hash"],
        "expected_format": expected_format,
        "compatible_fragment_ratio": ratio,
        "elapsed_seconds": elapsed_seconds,
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
