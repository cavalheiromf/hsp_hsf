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
            attributes = parse_attributes(fields[8])
            transcript = attributes.get("ID")
            parents = attributes.get("Parent", "").split(",")
            if fields[2] != "mRNA" and (not transcript or not parents[0].startswith("gene:")):
                continue
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
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
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
