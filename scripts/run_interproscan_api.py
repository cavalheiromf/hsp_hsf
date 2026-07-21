#!/usr/bin/env python3
"""Submit candidate proteins to the EMBL-EBI InterProScan REST service.

The client runs one batch at a time, persists every job identifier, and can
resume without resubmitting completed work.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


BASE_URL = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5"
DEFAULT_APPLICATIONS = (
    "PfamA,SMART,Panther,Gene3d,SuperFamily,PRINTS,"
    "PrositePatterns,PrositeProfiles"
)
TERMINAL_FAILURES = {"ERROR", "FAILURE", "NOT_FOUND"}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def parse_fasta(path: Path) -> list[tuple[str, str]]:
    """Read FASTA records while retaining the complete header."""
    records: list[tuple[str, str]] = []
    header: str | None = None
    sequence: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(sequence)))
                header, sequence = line[1:], []
            elif header is None:
                raise ValueError(f"Sequence before the first FASTA header in {path}")
            else:
                sequence.append(line)
    if header is not None:
        records.append((header, "".join(sequence)))
    if not records or any(not sequence for _, sequence in records):
        raise ValueError(f"Empty FASTA or sequence in {path}")
    return records


def format_fasta(records: list[tuple[str, str]]) -> str:
    blocks = []
    for header, sequence in records:
        wrapped = "\n".join(sequence[i : i + 60] for i in range(0, len(sequence), 60))
        blocks.append(f">{header}\n{wrapped}")
    return "\n".join(blocks) + "\n"


def chunks(records: list[tuple[str, str]], size: int):
    for start in range(0, len(records), size):
        yield records[start : start + size]


def request(
    url: str,
    *,
    data: dict[str, str] | None = None,
    retries: int = 5,
    timeout: int = 120,
) -> bytes:
    encoded = urllib.parse.urlencode(data).encode() if data is not None else None
    req = urllib.request.Request(
        url,
        data=encoded,
        headers={
            "Accept": "text/plain",
            "User-Agent": "hsp-hsf-pilot/1.0 (EMBL-EBI Job Dispatcher client)",
        },
    )
    for attempt in range(retries + 1):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as response:
                return response.read()
        except urllib.error.HTTPError as exc:
            if exc.code not in {429, 500, 502, 503, 504} or attempt == retries:
                detail = exc.read().decode(errors="replace")
                raise RuntimeError(f"HTTP {exc.code} from {url}: {detail}") from exc
        except urllib.error.URLError as exc:
            if attempt == retries:
                raise RuntimeError(f"Request failed for {url}: {exc.reason}") from exc
        time.sleep(min(60, 2**attempt * 5))
    raise AssertionError("unreachable")


def load_manifest(path: Path) -> dict:
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    return {"service": BASE_URL, "created_at": utc_now(), "batches": []}


def save_manifest(path: Path, manifest: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def submit(email: str, title: str, applications: str, fasta: str) -> str:
    payload = {
        "email": email,
        "title": title,
        "stype": "p",
        "appl": applications,
        "goterms": "true",
        "pathways": "true",
        "sequence": fasta,
    }
    return request(f"{BASE_URL}/run", data=payload).decode().strip()


def status(job_id: str) -> str:
    return request(f"{BASE_URL}/status/{job_id}").decode().strip().upper()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="Candidate protein FASTA")
    parser.add_argument("--species", required=True, help="Stable species slug")
    parser.add_argument("--outdir", type=Path, default=Path("results/interproscan"))
    parser.add_argument("--email", default=os.environ.get("INTERPRO_EMAIL"))
    parser.add_argument("--chunk-size", type=int, default=100)
    parser.add_argument("--poll-seconds", type=int, default=30)
    parser.add_argument("--applications", default=DEFAULT_APPLICATIONS)
    args = parser.parse_args()

    if not args.email or "@" not in args.email:
        parser.error("provide --email or set the INTERPRO_EMAIL environment variable")
    if not 1 <= args.chunk_size <= 1000:
        parser.error("--chunk-size must be between 1 and 1000")
    if args.poll_seconds < 10:
        parser.error("--poll-seconds must be at least 10")

    records = parse_fasta(args.input)
    species_dir = args.outdir / args.species
    species_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = species_dir / "jobs.json"
    manifest = load_manifest(manifest_path)
    manifest.update(
        {
            "species": args.species,
            "input": str(args.input),
            "applications": args.applications.split(","),
            "chunk_size": args.chunk_size,
            "sequence_count": len(records),
            "updated_at": utc_now(),
        }
    )

    batches_by_number = {item["batch"]: item for item in manifest["batches"]}
    for batch_number, batch_records in enumerate(chunks(records, args.chunk_size), start=1):
        fasta = format_fasta(batch_records)
        digest = hashlib.sha256(fasta.encode()).hexdigest()
        result_path = species_dir / f"batch_{batch_number:03d}.tsv"
        entry = batches_by_number.get(batch_number)

        if entry and entry.get("sha256") != digest:
            raise RuntimeError(
                f"Batch {batch_number} changed since submission; use a new output directory"
            )
        if entry and entry.get("status") == "FINISHED" and result_path.exists():
            print(f"batch {batch_number}: already finished ({entry['job_id']})", flush=True)
            continue
        if entry and entry.get("status") in TERMINAL_FAILURES:
            raise RuntimeError(
                f"Batch {batch_number} previously ended as {entry['status']}; "
                "inspect jobs.json before resubmitting"
            )

        if entry is None:
            job_id = submit(
                args.email,
                f"hsp_hsf_{args.species}_{batch_number:03d}",
                args.applications,
                fasta,
            )
            entry = {
                "batch": batch_number,
                "sequence_count": len(batch_records),
                "sha256": digest,
                "job_id": job_id,
                "status": "SUBMITTED",
                "submitted_at": utc_now(),
            }
            manifest["batches"].append(entry)
            batches_by_number[batch_number] = entry
            save_manifest(manifest_path, manifest)
            print(f"batch {batch_number}: submitted as {job_id}", flush=True)

        while True:
            current_status = status(entry["job_id"])
            entry["status"] = current_status
            entry["checked_at"] = utc_now()
            manifest["updated_at"] = utc_now()
            save_manifest(manifest_path, manifest)
            print(f"batch {batch_number}: {current_status}", flush=True)
            if current_status == "FINISHED":
                result = request(f"{BASE_URL}/result/{entry['job_id']}/tsv")
                temporary = result_path.with_suffix(".tsv.tmp")
                temporary.write_bytes(result)
                temporary.replace(result_path)
                entry["result"] = str(result_path)
                entry["finished_at"] = utc_now()
                save_manifest(manifest_path, manifest)
                break
            if current_status in TERMINAL_FAILURES:
                raise RuntimeError(
                    f"InterProScan batch {batch_number} ended as {current_status}"
                )
            time.sleep(args.poll_seconds)

    print(f"{args.species}: {len(records)} proteins annotated; manifest: {manifest_path}")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("Interrupted; rerun the same command to resume.", file=sys.stderr)
        raise SystemExit(130)
