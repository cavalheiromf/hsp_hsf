#!/usr/bin/env python3
"""Validate and document the four standardized RNA-seq references."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
from collections import Counter
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fasta_count(path: Path) -> tuple[int, int]:
    records = 0
    residues = 0
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                records += 1
            elif line.strip():
                residues += len(line.strip())
    if not records or not residues:
        raise ValueError(f"Empty FASTA: {path}")
    return records, residues


def gff_counts(path: Path) -> tuple[Counter[str], set[str]]:
    counts: Counter[str] = Counter()
    gene_ids: set[str] = set()
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"Malformed GFF3 row at {path}:{line_number}")
            counts[fields[2]] += 1
            if fields[2] == "gene":
                for attribute in fields[8].split(";"):
                    if attribute.startswith("ID="):
                        gene_id = attribute[3:]
                        gene_ids.add(gene_id.removeprefix("gene:"))
                        break
    if not counts.get("gene") or not counts.get("mRNA"):
        raise ValueError(f"GFF3 has no gene/mRNA features: {path}")
    return counts, gene_ids


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species-config", type=Path, default=Path("config/species.tsv"))
    parser.add_argument("--reference-dir", type=Path, default=Path("data/reference/ensembl_plants_63"))
    parser.add_argument("--mapping-dir", type=Path, default=Path("work/isoform_mapping"))
    parser.add_argument("--output-dir", type=Path, default=Path("results/rnaseq"))
    args = parser.parse_args()

    species_rows = read_tsv(args.species_config)
    manifest_rows: list[dict[str, str | int]] = []
    mapping_rows: list[dict[str, str]] = []
    for species in species_rows:
        species_id = species["species_id"]
        root = args.reference_dir / species_id
        genome = root / "genome.fa.gz"
        gff3 = root / "annotation.gff3.gz"
        proteins = root / "proteins_all.fa.gz"
        for path in (genome, gff3, proteins):
            if not path.is_file() or path.stat().st_size == 0:
                raise ValueError(f"Missing reference file: {path}")

        protein_count, residue_count = fasta_count(proteins)
        counts, gff_gene_ids = gff_counts(gff3)
        selected = read_tsv(args.mapping_dir / f"{species_id}.tsv")
        selected_genes = {row["gene_id"] for row in selected}
        if len(selected_genes) != len(selected):
            raise ValueError(f"Duplicate representative genes for {species_id}")
        if not selected_genes.issubset(gff_gene_ids):
            missing_in_gff = sorted(selected_genes - gff_gene_ids)
            raise ValueError(f"Representative genes absent from GFF3 for {species_id}: {missing_in_gff[:5]}")
        if len(selected) > counts["gene"]:
            raise ValueError(
                f"Representative/GFF3 gene mismatch for {species_id}: "
                f"{len(selected)} vs {counts['gene']}"
            )
        for row in selected:
            mapping_rows.append({"species": species_id, **row})

        manifest_rows.append(
            {
                "species": species_id,
                "scientific_name": species["scientific_name"],
                "release": species["ensembl_release"],
                "assembly": species["assembly"],
                "assembly_accession": species["assembly_accession"],
                "genome": str(genome),
                "annotation": str(gff3),
                "proteins": str(proteins),
                "genes": counts["gene"],
                "mrnas": counts["mRNA"],
                "representatives": len(selected),
                "genes_without_representative_protein": counts["gene"] - len(selected),
                "protein_records": protein_count,
                "protein_residues": residue_count,
                "genome_sha256": sha256(genome),
                "annotation_sha256": sha256(gff3),
                "proteins_sha256": sha256(proteins),
            }
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = args.output_dir / "reference_manifest.tsv"
    mapping_path = args.output_dir / "gene_transcript_protein.tsv"
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(manifest_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)
    with mapping_path.open("w", encoding="utf-8", newline="") as handle:
        fields = list(mapping_rows[0])
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(mapping_rows)
    print(f"Validated {len(manifest_rows)} species and wrote {len(mapping_rows)} representative mappings")


if __name__ == "__main__":
    main()
