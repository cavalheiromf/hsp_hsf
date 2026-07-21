#!/usr/bin/env python3
"""Download and verify frozen Ensembl Plants reference files."""

from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path
from urllib.parse import urlparse


LOCAL_NAMES = {
    "genome_url": "genome.fa.gz",
    "gff3_url": "annotation.gff3.gz",
    "proteome_url": "proteins_all.fa.gz",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=Path("config/species.tsv"))
    parser.add_argument("--outdir", type=Path, default=Path("data/reference/ensembl_plants_63"))
    parser.add_argument("--species", action="append", dest="species_ids")
    parser.add_argument("--retries", type=int, default=3)
    return parser.parse_args()


def read_config(path: Path, selected: set[str] | None) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if selected:
        known = {row["species_id"] for row in rows}
        unknown = selected - known
        if unknown:
            raise ValueError(f"Unknown species: {', '.join(sorted(unknown))}")
        rows = [row for row in rows if row["species_id"] in selected]
    if not rows:
        raise ValueError("No species selected")
    return rows


def checksum_url(file_url: str) -> str:
    return file_url.rsplit("/", 1)[0] + "/CHECKSUMS"


def fetch_text(url: str, retries: int) -> str:
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(url, timeout=120) as response:
                return response.read().decode("utf-8")
        except (OSError, urllib.error.URLError) as error:
            if attempt == retries:
                raise RuntimeError(f"Failed to fetch {url}: {error}") from error
            time.sleep(2**attempt)
    raise AssertionError("unreachable")


def expected_sum(checksums: str, remote_name: str) -> tuple[str, str]:
    for line in checksums.splitlines():
        fields = line.split()
        if len(fields) >= 3 and fields[-1] == remote_name:
            return fields[0], fields[1]
    raise ValueError(f"No official checksum found for {remote_name}")


def local_sum(path: Path) -> tuple[str, str]:
    result = subprocess.run(["sum", str(path)], check=True, capture_output=True, text=True)
    fields = result.stdout.split()
    return fields[0], fields[1]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate(path: Path, expected: tuple[str, str]) -> str:
    observed = local_sum(path)
    if observed != expected:
        raise ValueError(f"Checksum mismatch for {path}: expected {expected}, observed {observed}")
    subprocess.run(["gzip", "-t", str(path)], check=True)
    return sha256(path)


def download(url: str, target: Path, retries: int) -> None:
    partial = target.with_suffix(target.suffix + ".part")
    if partial.exists():
        partial.unlink()
    for attempt in range(1, retries + 1):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "hsp-hsf-reference-pipeline/1.0"})
            with urllib.request.urlopen(request, timeout=120) as response, partial.open("wb") as out:
                shutil.copyfileobj(response, out, length=1024 * 1024)
            partial.replace(target)
            return
        except (OSError, urllib.error.URLError) as error:
            partial.unlink(missing_ok=True)
            if attempt == retries:
                raise RuntimeError(f"Failed to download {url}: {error}") from error
            time.sleep(2**attempt)


def process_species(row: dict[str, str], outdir: Path, retries: int) -> None:
    species_dir = outdir / row["species_id"]
    species_dir.mkdir(parents=True, exist_ok=True)
    manifest_rows: list[dict[str, str | int]] = []

    for url_field, local_name in LOCAL_NAMES.items():
        url = row[url_field]
        remote_name = Path(urlparse(url).path).name
        official = expected_sum(fetch_text(checksum_url(url), retries), remote_name)
        target = species_dir / local_name

        if target.exists():
            digest = validate(target, official)
            status = "existing_valid"
        else:
            print(f"Downloading {row['species_id']} {local_name}", flush=True)
            download(url, target, retries)
            digest = validate(target, official)
            status = "downloaded"

        manifest_rows.append(
            {
                "species_id": row["species_id"],
                "scientific_name": row["scientific_name"],
                "ensembl_release": row["ensembl_release"],
                "assembly": row["assembly"],
                "assembly_accession": row["assembly_accession"],
                "file_type": url_field.removesuffix("_url"),
                "url": url,
                "local_file": local_name,
                "bytes": target.stat().st_size,
                "ensembl_sum": official[0],
                "ensembl_blocks": official[1],
                "sha256": digest,
                "status": status,
            }
        )

    manifest = species_dir / "checksums.tsv"
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(manifest_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)
    print(f"Verified {row['species_id']}: {manifest}")


def main() -> None:
    args = parse_args()
    selected = set(args.species_ids) if args.species_ids else None
    rows = read_config(args.config, selected)
    for row in rows:
        process_species(row, args.outdir, args.retries)


if __name__ == "__main__":
    try:
        main()
    except (RuntimeError, ValueError, subprocess.CalledProcessError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
