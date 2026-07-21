#!/usr/bin/env python3
"""Create an evidence table for discrepant and newly detected HSFs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from Bio import Align

from select_representative_isoforms import parse_fasta
from validate_historical_hsfs import parse_alignment


FIELDS = [
    "case_type", "species", "historical_id", "gene_id", "protein_id",
    "historical_length", "current_length", "length_delta", "aligned_identity",
    "alignment_gap_columns", "hmm_score", "hmm_evalue", "pfam_domain",
    "interpro_id", "supporting_databases", "interpretation",
]


def align_identity(historical: str, current: str) -> tuple[str, int]:
    """Measure identity on aligned residues, separating indels from substitutions."""
    aligner = Align.PairwiseAligner()
    aligner.match_score = 1
    aligner.mismatch_score = 0
    aligner.open_gap_score = 0
    aligner.extend_gap_score = 0
    alignment = aligner.align(historical, current)[0]
    left, right = alignment[0], alignment[1]
    aligned = [(a, b) for a, b in zip(left, right) if a != "-" and b != "-"]
    matches = sum(a == b for a, b in aligned)
    gaps = sum(a == "-" or b == "-" for a, b in zip(left, right))
    return f"{matches / len(aligned):.4f}" if aligned else "", gaps


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", type=Path, default=Path("results/catalog/hsp_hsf_catalog.tsv"))
    parser.add_argument("--validation", type=Path, default=Path("results/catalog/historical_hsf_validation.tsv"))
    parser.add_argument("--new-candidates", type=Path, default=Path("results/catalog/new_hsf_candidates.tsv"))
    parser.add_argument("--alignment", type=Path, default=Path("reference/external/Alinhamento_HSF.fas"))
    parser.add_argument("--reference-dir", type=Path, default=Path("data/reference/ensembl_plants_63"))
    parser.add_argument("--output", type=Path, default=Path("results/catalog/hsf_discrepancy_review.tsv"))
    args = parser.parse_args()

    catalog = {
        row["protein_id"]: row
        for row in read_rows(args.catalog)
        if row["family"] == "HSF"
    }
    validation = read_rows(args.validation)
    new_candidates = read_rows(args.new_candidates)
    historical = parse_alignment(args.alignment)
    proteins: dict[tuple[str, str], str] = {}
    for species in {row["species"] for row in catalog.values()}:
        proteome = args.reference_dir / species / "proteins_all.fa.gz"
        proteins.update({(species, p.protein_id): p.sequence for p in parse_fasta(proteome)})

    rows: list[dict[str, str]] = []
    for item in validation:
        if item["sequence_status"] != "sequence_changed":
            continue
        current = proteins[(item["species"], item["current_protein_id"])]
        identity, gaps = align_identity(historical[item["historical_id"]], current)
        c = catalog[item["current_protein_id"]]
        rows.append(
            {
                "case_type": "historical_sequence_changed",
                "species": item["species"],
                "historical_id": item["historical_id"],
                "gene_id": item["mapped_gene_id"],
                "protein_id": item["current_protein_id"],
                "historical_length": item["historical_sequence_length"],
                "current_length": item["current_protein_length"],
                "length_delta": str(int(item["current_protein_length"]) - int(item["historical_sequence_length"])),
                "aligned_identity": identity,
                "alignment_gap_columns": str(gaps),
                "hmm_score": c["hmm_score"],
                "hmm_evalue": c["hmm_evalue"],
                "pfam_domain": f"{c['domain_start']}-{c['domain_end']}",
                "interpro_id": c["interpro_id"],
                "supporting_databases": c["interpro_supporting_databases"],
                "interpretation": "indel_only; retain as updated annotation",
            }
        )

    for item in new_candidates:
        c = catalog[item["protein_id"]]
        sequence = proteins[(item["species"], item["protein_id"])]
        rows.append(
            {
                "case_type": "new_release63_candidate",
                "species": item["species"],
                "historical_id": "",
                "gene_id": item["gene_id"],
                "protein_id": item["protein_id"],
                "historical_length": "",
                "current_length": str(len(sequence)),
                "length_delta": "",
                "aligned_identity": "",
                "alignment_gap_columns": "",
                "hmm_score": c["hmm_score"],
                "hmm_evalue": c["hmm_evalue"],
                "pfam_domain": f"{c['domain_start']}-{c['domain_end']}",
                "interpro_id": c["interpro_id"],
                "supporting_databases": c["interpro_supporting_databases"],
                "interpretation": "confirmed PF00447 candidate absent from historical set",
            }
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(sorted(rows, key=lambda row: (row["species"], row["case_type"], row["gene_id"])))
    print(f"Wrote {len(rows)} HSF review cases to {args.output}")


if __name__ == "__main__":
    main()
