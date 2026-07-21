#!/usr/bin/env python3
"""Select one canonical-or-longest protein isoform per gene."""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import unquote


HEADER_FIELD = re.compile(r"(?:^|\s)(gene|transcript):([^\s]+)")


@dataclass(frozen=True)
class Protein:
    protein_id: str
    gene_id: str
    transcript_id: str
    sequence: str


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def clean_id(value: str) -> str:
    return value.split(":", 1)[1] if value.startswith(("gene:", "transcript:", "protein:")) else value


def parse_fasta(path: Path) -> list[Protein]:
    proteins: list[Protein] = []
    header: str | None = None
    sequence: list[str] = []

    def emit() -> None:
        if header is None:
            return
        protein_id = header.split()[0]
        fields = {key: clean_id(value) for key, value in HEADER_FIELD.findall(header)}
        if "gene" not in fields or "transcript" not in fields:
            raise ValueError(f"Missing gene/transcript fields in FASTA header: {header}")
        seq = "".join(sequence).replace("*", "").upper()
        if not seq:
            raise ValueError(f"Empty protein sequence: {protein_id}")
        proteins.append(Protein(protein_id, fields["gene"], fields["transcript"], seq))

    with open_text(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                emit()
                header = line[1:]
                sequence = []
            elif header is None:
                raise ValueError("Sequence encountered before first FASTA header")
            else:
                sequence.append(line)
    emit()

    ids = [protein.protein_id for protein in proteins]
    if len(ids) != len(set(ids)):
        raise ValueError("Duplicated protein identifiers in FASTA")
    if not proteins:
        raise ValueError("No proteins found in FASTA")
    return proteins


def parse_attributes(text: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for item in text.rstrip().split(";"):
        if not item or "=" not in item:
            continue
        key, value = item.split("=", 1)
        attributes[key] = unquote(value)
    return attributes


def canonical_transcripts(gff3: Path) -> set[str]:
    canonical: set[str] = set()
    with open_text(gff3) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            feature = fields[2].lower()
            attributes = parse_attributes(fields[8])

            for key in ("canonical_transcript", "canonical_transcript_id"):
                if key in attributes:
                    canonical.update(clean_id(value) for value in attributes[key].split(","))

            if feature in {"mrna", "transcript"}:
                transcript_id = clean_id(attributes.get("ID", ""))
                tags = {tag.lower() for tag in attributes.get("tag", "").split(",")}
                is_canonical = attributes.get("is_canonical", "").lower() in {"1", "true", "yes"}
                tagged_canonical = any(tag == "canonical" or tag.endswith("_canonical") for tag in tags)
                if transcript_id and (tagged_canonical or is_canonical):
                    canonical.add(transcript_id)
    return canonical


def choose(proteins: list[Protein], canonical: set[str]) -> tuple[Protein, str]:
    canonical_options = [protein for protein in proteins if protein.transcript_id in canonical]
    if canonical_options:
        candidates = canonical_options
        base_rule = "canonical"
    else:
        max_length = max(len(protein.sequence) for protein in proteins)
        candidates = [protein for protein in proteins if len(protein.sequence) == max_length]
        base_rule = "longest"
    selected = min(candidates, key=lambda protein: protein.protein_id)
    if base_rule == "longest" and len(candidates) > 1:
        return selected, "length_tie_id"
    return selected, base_rule


def select_representatives(proteins: list[Protein], canonical: set[str]):
    by_gene: dict[str, list[Protein]] = defaultdict(list)
    for protein in proteins:
        by_gene[protein.gene_id].append(protein)
    for gene_id in sorted(by_gene):
        selected, rule = choose(by_gene[gene_id], canonical)
        yield selected, rule, sorted(by_gene[gene_id], key=lambda protein: protein.protein_id)


def write_outputs(proteins: list[Protein], canonical: set[str], fasta_out: Path, mapping_out: Path) -> None:
    fasta_out.parent.mkdir(parents=True, exist_ok=True)
    mapping_out.parent.mkdir(parents=True, exist_ok=True)
    selected_rows = list(select_representatives(proteins, canonical))

    with fasta_out.open("w", encoding="utf-8") as fasta:
        for selected, rule, _ in selected_rows:
            fasta.write(
                f">{selected.protein_id} gene:{selected.gene_id} "
                f"transcript:{selected.transcript_id} selection:{rule}\n"
            )
            for start in range(0, len(selected.sequence), 60):
                fasta.write(selected.sequence[start : start + 60] + "\n")

    fields = [
        "gene_id",
        "selected_transcript_id",
        "selected_protein_id",
        "protein_length",
        "selection_rule",
        "alternative_transcript_ids",
        "alternative_protein_ids",
    ]
    with mapping_out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for selected, rule, alternatives in selected_rows:
            writer.writerow(
                {
                    "gene_id": selected.gene_id,
                    "selected_transcript_id": selected.transcript_id,
                    "selected_protein_id": selected.protein_id,
                    "protein_length": len(selected.sequence),
                    "selection_rule": rule,
                    "alternative_transcript_ids": ";".join(p.transcript_id for p in alternatives),
                    "alternative_protein_ids": ";".join(p.protein_id for p in alternatives),
                }
            )

    if len(selected_rows) != len({selected.gene_id for selected, _, _ in selected_rows}):
        raise AssertionError("Representative output contains duplicated genes")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--proteins", type=Path, required=True)
    parser.add_argument("--gff3", type=Path, required=True)
    parser.add_argument("--fasta-out", type=Path, required=True)
    parser.add_argument("--mapping-out", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    proteins = parse_fasta(args.proteins)
    canonical = canonical_transcripts(args.gff3)
    write_outputs(proteins, canonical, args.fasta_out, args.mapping_out)
    genes = len({protein.gene_id for protein in proteins})
    print(f"Proteins: {len(proteins)}; genes/representatives: {genes}; canonical transcripts: {len(canonical)}")


if __name__ == "__main__":
    try:
        main()
    except ValueError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
