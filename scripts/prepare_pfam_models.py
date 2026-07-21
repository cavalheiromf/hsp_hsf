#!/usr/bin/env python3
"""Download Pfam-A and extract the target HSP/HSF profiles."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path


def md5(path: Path) -> str:
    digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_single_row(path: Path) -> dict[str, str]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 1:
        raise ValueError(f"Expected one Pfam release row in {path}")
    return rows[0]


def read_families(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError("No Pfam families configured")
    return rows


def download(url: str, target: Path) -> None:
    partial = target.with_suffix(target.suffix + ".part")
    partial.unlink(missing_ok=True)
    request = urllib.request.Request(url, headers={"User-Agent": "hsp-hsf-reference-pipeline/1.0"})
    with urllib.request.urlopen(request, timeout=180) as response, partial.open("wb") as out:
        shutil.copyfileobj(response, out, length=1024 * 1024)
    partial.replace(target)


def decompress(source: Path, target: Path) -> None:
    partial = target.with_suffix(target.suffix + ".part")
    partial.unlink(missing_ok=True)
    with gzip.open(source, "rb") as compressed, partial.open("wb") as out:
        shutil.copyfileobj(compressed, out, length=1024 * 1024)
    partial.replace(target)


def verify_profile(path: Path, pfam_id: str) -> None:
    text = path.read_text(encoding="utf-8")
    accession = next((line.split()[1].split(".", 1)[0] for line in text.splitlines() if line.startswith("ACC")), None)
    if accession != pfam_id:
        raise ValueError(f"Expected {pfam_id} in {path}, observed {accession}")
    if not any(line.startswith("GA") for line in text.splitlines()):
        raise ValueError(f"Profile {pfam_id} lacks a gathering threshold")


def resolve_accession_versions(database: Path, requested: set[str]) -> dict[str, str]:
    resolved: dict[str, str] = {}
    with database.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("ACC"):
                continue
            accession = line.split()[1]
            base = accession.split(".", 1)[0]
            if base in requested:
                resolved[base] = accession
                if resolved.keys() == requested:
                    break
    missing = requested - resolved.keys()
    if missing:
        raise ValueError(f"Pfam accessions absent from release: {', '.join(sorted(missing))}")
    return resolved


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--release-config", type=Path, default=Path("config/pfam_release.tsv"))
    parser.add_argument("--families", type=Path, default=Path("config/pfam_families.tsv"))
    parser.add_argument("--outdir", type=Path, default=Path("data/pfam"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    release = read_single_row(args.release_config)
    families = read_families(args.families)
    release_dir = args.outdir / release["release"]
    models_dir = release_dir / "models"
    models_dir.mkdir(parents=True, exist_ok=True)

    archive = release_dir / "Pfam-A.hmm.gz"
    database = release_dir / "Pfam-A.hmm"
    if not archive.exists():
        print(f"Downloading Pfam {release['release']} profiles", flush=True)
        download(release["url"], archive)
    observed_md5 = md5(archive)
    if observed_md5 != release["md5"]:
        raise ValueError(f"Pfam archive MD5 mismatch: expected {release['md5']}, observed {observed_md5}")

    if not database.exists():
        print("Decompressing Pfam-A.hmm", flush=True)
        decompress(archive, database)

    hmmfetch = shutil.which("hmmfetch")
    if not hmmfetch:
        raise RuntimeError("hmmfetch not found; load Bio/HMMER3/3.4")

    resolved = resolve_accession_versions(database, {row["pfam_id"] for row in families})
    index = Path(str(database) + ".ssi")
    if not index.exists():
        print("Indexing Pfam-A.hmm for accession lookup", flush=True)
        subprocess.run([hmmfetch, "--index", str(database)], check=True)

    metadata: list[dict[str, str | int]] = []
    for row in families:
        output = models_dir / f"{row['family']}__{row['pfam_id']}.hmm"
        partial = output.with_suffix(output.suffix + ".part")
        partial.unlink(missing_ok=True)
        with partial.open("w", encoding="utf-8") as handle:
            subprocess.run([hmmfetch, str(database), resolved[row["pfam_id"]]], check=True, stdout=handle)
        partial.replace(output)
        verify_profile(output, row["pfam_id"])
        metadata.append(
            {
                **row,
                "pfam_release": release["release"],
                "pfam_accession_version": resolved[row["pfam_id"]],
                "profile_file": output.name,
                "bytes": output.stat().st_size,
                "sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
            }
        )

    metadata_path = release_dir / "metadata.tsv"
    with metadata_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(metadata[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(metadata)
    print(f"Extracted {len(metadata)} profiles: {metadata_path}")


if __name__ == "__main__":
    try:
        main()
    except (OSError, RuntimeError, ValueError, subprocess.CalledProcessError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
