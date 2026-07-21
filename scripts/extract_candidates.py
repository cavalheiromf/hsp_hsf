#!/usr/bin/env python3
"""Merge accepted HMMER hits and extract candidate protein sequences."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from select_representative_isoforms import parse_fasta


DOM_COLUMNS = [
    "target_name", "target_accession", "target_length", "query_name", "query_accession",
    "query_length", "full_evalue", "full_score", "full_bias", "domain_number",
    "domain_count", "domain_cevalue", "domain_ievalue", "domain_score", "domain_bias",
    "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to", "accuracy",
]


def parse_domtblout(path: Path, species: str) -> list[dict[str, str | float]]:
    stem = path.name.removesuffix(".domtblout")
    family, pfam_id = stem.split("__", 1)
    rows: list[dict[str, str | float]] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split(maxsplit=22)
            if len(fields) < 22:
                raise ValueError(f"Malformed domtblout row in {path}")
            record = dict(zip(DOM_COLUMNS, fields[:22]))
            record.update(
                {
                    "species": species,
                    "family": family,
                    "pfam_id": pfam_id,
                    "hmm_coverage": (int(record["hmm_to"]) - int(record["hmm_from"]) + 1) / int(record["query_length"]),
                    "target_coverage": (int(record["ali_to"]) - int(record["ali_from"]) + 1) / int(record["target_length"]),
                    "source_file": str(path),
                }
            )
            rows.append(record)
    return rows


def write_fasta(proteome: Path, candidate_ids: set[str], output: Path) -> None:
    proteins = {protein.protein_id: protein for protein in parse_fasta(proteome)}
    missing = candidate_ids - proteins.keys()
    if missing:
        raise ValueError(f"Candidate IDs absent from proteome: {', '.join(sorted(missing)[:10])}")
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as handle:
        for protein_id in sorted(candidate_ids):
            protein = proteins[protein_id]
            handle.write(f">{protein.protein_id} gene:{protein.gene_id} transcript:{protein.transcript_id}\n")
            for start in range(0, len(protein.sequence), 60):
                handle.write(protein.sequence[start : start + 60] + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--species", required=True)
    parser.add_argument("--proteome", type=Path, required=True)
    parser.add_argument("--hmmer-dir", type=Path, required=True)
    parser.add_argument("--fasta-out", type=Path, required=True)
    parser.add_argument("--hits-out", type=Path, required=True)
    args = parser.parse_args()

    files = sorted(args.hmmer_dir.glob("*.domtblout"))
    if not files:
        raise ValueError(f"No domtblout files in {args.hmmer_dir}")
    rows = [row for path in files for row in parse_domtblout(path, args.species)]
    if not rows:
        raise ValueError(f"No accepted HMMER hits for {args.species}")

    args.hits_out.parent.mkdir(parents=True, exist_ok=True)
    fields = ["species", "family", "pfam_id"] + DOM_COLUMNS + ["hmm_coverage", "target_coverage", "source_file"]
    with args.hits_out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    write_fasta(args.proteome, {str(row["target_name"]) for row in rows}, args.fasta_out)
    print(f"{args.species}: {len(rows)} domain hits; {len({row['target_name'] for row in rows})} candidate proteins")


if __name__ == "__main__":
    main()
