#!/usr/bin/env python3
"""Create a reproducibility manifest for a B0 tracking validation run."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_commit(path: Path | None) -> str | None:
    if path is None:
        return None
    try:
        return subprocess.run(
            ["git", "-C", str(path), "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="simulation/reconstruction input file")
    parser.add_argument("--dataset", default="", help="stable dataset/sample identifier")
    parser.add_argument("--material-map", type=Path)
    parser.add_argument("--b0trackers-repo", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--eicrecon-repo", type=Path)
    parser.add_argument("--epic-repo", type=Path)
    parser.add_argument("--primary-pdg", type=int, default=2212)
    parser.add_argument("--primary-status", type=int, default=1)
    parser.add_argument("--events", type=int)
    parser.add_argument("--arg", action="append", default=[], help="important reconstruction argument")
    parser.add_argument("--output", type=Path, default=Path("b0-run-manifest.json"))
    args = parser.parse_args()

    if not args.input.exists():
        parser.error(f"input does not exist: {args.input}")
    if args.material_map and not args.material_map.exists():
        parser.error(f"material map does not exist: {args.material_map}")

    manifest = {
        "input_dataset": args.dataset,
        "input_path": str(args.input.resolve()),
        "input_sha256": sha256(args.input),
        "events": args.events,
        "primary_pdg": args.primary_pdg,
        "primary_status": args.primary_status,
        "detector_config": os.environ.get("DETECTOR_CONFIG", ""),
        "detector_path": os.environ.get("DETECTOR_PATH", ""),
        "b0trackers_commit": git_commit(args.b0trackers_repo),
        "eicrecon_commit": git_commit(args.eicrecon_repo),
        "epic_commit": git_commit(args.epic_repo),
        "material_map_path": str(args.material_map.resolve()) if args.material_map else "",
        "material_map_sha256": sha256(args.material_map) if args.material_map else "",
        "reconstruction_args": args.arg,
    }
    args.output.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
