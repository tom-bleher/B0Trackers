#!/usr/bin/env python3
"""Benchmark B0Trackers overhead and thread scaling with reproducible eicrecon runs."""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import statistics
import subprocess
import time
from pathlib import Path


def _gnu_time() -> str | None:
    candidate = Path("/usr/bin/time")
    if not candidate.exists():
        return None
    try:
        result = subprocess.run(
            [str(candidate), "--version"],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
    except OSError:
        return None
    return str(candidate) if "GNU time" in result.stdout else None


def _run_once(command: list[str], metrics_path: Path) -> dict:
    timer = _gnu_time()
    if timer:
        timed_command = [
            timer,
            "-f",
            '{"wall_s":%e,"user_s":%U,"sys_s":%S,"max_rss_kb":%M}',
            "-o",
            str(metrics_path),
            *command,
        ]
        start = time.perf_counter()
        completed = subprocess.run(timed_command, check=False)
        fallback_wall = time.perf_counter() - start
        if completed.returncode != 0:
            raise RuntimeError(f"command failed with exit code {completed.returncode}: {' '.join(command)}")
        try:
            metrics = json.loads(metrics_path.read_text())
        except (OSError, json.JSONDecodeError):
            metrics = {"wall_s": fallback_wall, "user_s": None, "sys_s": None, "max_rss_kb": None}
        return metrics

    start = time.perf_counter()
    completed = subprocess.run(command, check=False)
    wall = time.perf_counter() - start
    if completed.returncode != 0:
        raise RuntimeError(f"command failed with exit code {completed.returncode}: {' '.join(command)}")
    return {"wall_s": wall, "user_s": None, "sys_s": None, "max_rss_kb": None}


def _aggregate(records: list[dict]) -> dict:
    out: dict[str, float | int | None] = {"runs": len(records)}
    for key in ["wall_s", "user_s", "sys_s", "max_rss_kb", "hist_bytes", "podio_bytes"]:
        values = [float(row[key]) for row in records if row.get(key) is not None]
        out[f"{key}_mean"] = statistics.fmean(values) if values else None
        out[f"{key}_stdev"] = statistics.stdev(values) if len(values) > 1 else 0.0 if values else None
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="simulation input file")
    parser.add_argument("--eicrecon", default="eicrecon")
    parser.add_argument("--threads", default="1,2,4,8", help="comma-separated thread counts")
    parser.add_argument("--runs", type=int, default=3, help="repetitions per mode/thread count")
    parser.add_argument("--events", type=int, help="optional fixed jana:nevents for quicker comparisons")
    parser.add_argument("--arg", action="append", default=[], help="additional eicrecon argument; repeat as needed")
    parser.add_argument("--plugin-arg", default="-Pplugins=B0Trackers")
    parser.add_argument("--output-dir", type=Path, default=Path("b0-runtime"))
    parser.add_argument("--keep-output", action="store_true")
    args = parser.parse_args()

    if args.runs < 1:
        parser.error("--runs must be >= 1")
    threads = [int(item) for item in args.threads.split(",") if item.strip()]
    if not threads or any(value < 1 for value in threads):
        parser.error("--threads must contain positive integers")
    if not args.input.exists():
        parser.error(f"input does not exist: {args.input}")
    if shutil.which(args.eicrecon) is None:
        parser.error(f"cannot find eicrecon executable: {args.eicrecon}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    records: list[dict] = []

    for threads_count in threads:
        for mode in ["baseline", "b0trackers"]:
            for run in range(1, args.runs + 1):
                stem = f"{mode}-t{threads_count}-r{run}"
                podio = args.output_dir / f"{stem}.edm4eic.root"
                hist = args.output_dir / f"{stem}.hists.root"
                time_json = args.output_dir / f"{stem}.time.json"

                command = [
                    args.eicrecon,
                    *args.arg,
                    f"-Pnthreads={threads_count}",
                    f"-Ppodio:output_file={podio}",
                ]
                if args.events is not None:
                    command.append(f"-Pjana:nevents={args.events}")
                if mode == "b0trackers":
                    command.extend([args.plugin_arg, f"-Phistsfile={hist}"])
                command.append(str(args.input))

                print(f"[{mode} threads={threads_count} run={run}] {' '.join(command)}", flush=True)
                metrics = _run_once(command, time_json)
                row = {
                    "mode": mode,
                    "threads": threads_count,
                    "run": run,
                    **metrics,
                    "podio_bytes": podio.stat().st_size if podio.exists() else None,
                    "hist_bytes": hist.stat().st_size if hist.exists() else None,
                }
                records.append(row)

                if not args.keep_output:
                    for path in [podio, hist, time_json]:
                        try:
                            path.unlink()
                        except FileNotFoundError:
                            pass

    csv_path = args.output_dir / "runs.csv"
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)

    summary: dict[str, dict] = {}
    for threads_count in threads:
        for mode in ["baseline", "b0trackers"]:
            subset = [row for row in records if row["threads"] == threads_count and row["mode"] == mode]
            summary[f"{mode}_t{threads_count}"] = _aggregate(subset)

        base = summary[f"baseline_t{threads_count}"]["wall_s_mean"]
        plugin = summary[f"b0trackers_t{threads_count}"]["wall_s_mean"]
        if base and plugin:
            summary[f"overhead_t{threads_count}"] = {
                "wall_ratio": plugin / base,
                "wall_overhead_fraction": plugin / base - 1.0,
            }

    output = {
        "input": str(args.input),
        "runs_per_configuration": args.runs,
        "threads": threads,
        "events": args.events,
        "extra_args": args.arg,
        "summary": summary,
    }
    summary_path = args.output_dir / "summary.json"
    summary_path.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")

    print(f"wrote {csv_path}")
    print(f"wrote {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
