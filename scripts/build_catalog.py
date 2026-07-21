#!/usr/bin/env python3
"""Integrate HMMER assignments, representative isoforms, and InterProScan hits."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


CATALOG_FIELDS = [
    "species", "gene_id", "transcript_id", "protein_id", "protein_length",
    "isoform_selection", "family", "pfam_id", "hmm_score", "hmm_evalue",
    "domain_score", "domain_evalue", "domain_start", "domain_end",
    "domain_coverage", "hmm_domain_count", "interpro_signature_id",
    "interpro_id", "interpro_description", "interpro_supporting_databases",
    "classification_status", "reference_release", "assembly_accession",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_interpro(directory: Path) -> dict[str, list[dict[str, str]]]:
    """Index the standard 15-column InterProScan TSV by protein accession."""
    fields = [
        "protein_id", "md5", "length", "database", "signature_id",
        "signature_description", "start", "end", "score", "status", "date",
        "interpro_id", "interpro_description", "go_terms", "pathways",
    ]
    matches: dict[str, list[dict[str, str]]] = defaultdict(list)
    files = sorted(directory.glob("batch_*.tsv"))
    if not files:
        raise ValueError(f"No InterProScan batch TSVs in {directory}")
    for path in files:
        with path.open(encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, start=1):
                if not line.strip() or line.startswith("#"):
                    continue
                values = line.rstrip("\n").split("\t")
                if len(values) != len(fields):
                    raise ValueError(f"Expected 15 columns in {path}:{line_number}")
                record = dict(zip(fields, values))
                matches[record["protein_id"]].append(record)
    return matches


def build_species(
    species: dict[str, str], work_dir: Path, interpro_dir: Path
) -> list[dict[str, str]]:
    species_id = species["species_id"]
    isoforms = {
        row["selected_protein_id"]: row
        for row in read_tsv(work_dir / "isoform_mapping" / f"{species_id}.tsv")
    }
    hits = read_tsv(work_dir / "candidates" / f"{species_id}.hits.tsv")
    interpro = read_interpro(interpro_dir / species_id)

    # A catalog row is a protein-family assignment. For repeated domains of the
    # same family, retain the highest-scoring domain and record the multiplicity.
    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    protein_families: dict[str, set[str]] = defaultdict(set)
    for hit in hits:
        key = (hit["target_name"], hit["family"])
        grouped[key].append(hit)
        protein_families[hit["target_name"]].add(hit["family"])

    rows: list[dict[str, str]] = []
    for (protein_id, family), family_hits in sorted(grouped.items()):
        best = max(family_hits, key=lambda item: float(item["domain_score"]))
        if protein_id not in isoforms:
            raise ValueError(f"Missing isoform metadata for {protein_id}")
        metadata = isoforms[protein_id]
        pfam_id = best["pfam_id"].split(".", 1)[0]
        all_matches = interpro.get(protein_id, [])
        exact = [
            match for match in all_matches
            if match["database"] == "Pfam"
            and match["signature_id"].split(".", 1)[0] == pfam_id
        ]

        if len(protein_families[protein_id]) > 1:
            classification = "multi_family"
        elif exact:
            classification = "confirmed"
        elif any(match["interpro_id"] != "-" for match in all_matches):
            classification = "interpro_related"
        else:
            classification = "hmm_only"

        interpro_ids = sorted({m["interpro_id"] for m in exact if m["interpro_id"] != "-"})
        descriptions = sorted({m["interpro_description"] for m in exact if m["interpro_description"] != "-"})
        databases = sorted({m["database"] for m in all_matches})
        rows.append(
            {
                "species": species_id,
                "gene_id": metadata["gene_id"],
                "transcript_id": metadata["selected_transcript_id"],
                "protein_id": protein_id,
                "protein_length": metadata["protein_length"],
                "isoform_selection": metadata["selection_rule"],
                "family": family,
                "pfam_id": pfam_id,
                "hmm_score": best["full_score"],
                "hmm_evalue": best["full_evalue"],
                "domain_score": best["domain_score"],
                "domain_evalue": best["domain_ievalue"],
                "domain_start": best["ali_from"],
                "domain_end": best["ali_to"],
                "domain_coverage": best["hmm_coverage"],
                "hmm_domain_count": str(len(family_hits)),
                "interpro_signature_id": ";".join(sorted({m["signature_id"] for m in exact})),
                "interpro_id": ";".join(interpro_ids),
                "interpro_description": ";".join(descriptions),
                "interpro_supporting_databases": ";".join(databases),
                "classification_status": classification,
                "reference_release": species["ensembl_release"],
                "assembly_accession": species["assembly_accession"],
            }
        )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species-config", type=Path, default=Path("config/species.tsv"))
    parser.add_argument("--work-dir", type=Path, default=Path("work"))
    parser.add_argument("--interpro-dir", type=Path, default=Path("results/interproscan"))
    parser.add_argument("--output", type=Path, default=Path("results/catalog/hsp_hsf_catalog.tsv"))
    args = parser.parse_args()

    species_rows = read_tsv(args.species_config)
    catalog = [
        row
        for species in species_rows
        for row in build_species(species, args.work_dir, args.interpro_dir)
    ]
    if not catalog:
        raise ValueError("Catalog is empty")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CATALOG_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(catalog)
    print(f"Wrote {len(catalog)} protein-family assignments to {args.output}")


if __name__ == "__main__":
    main()
