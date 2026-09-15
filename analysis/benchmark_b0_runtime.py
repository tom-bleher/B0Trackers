#!/usr/bin/env python3
"""Benchmark B0Trackers overhead and thread scaling with reproducible eicrecon runs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import random
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


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


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


def _mad(values: list[float]) -> float | None:
    if not values:
        return None
    center = statistics.median(values)
    return statistics.median(abs(value - center) for value in values)


def _aggregate(records: list[dict]) -> dict:
    out: dict[str, float | int | None] = {"runs": len(records)}
    for key in ["wall_s", "user_s", "sys_s", "max_rss_kb", "hist_bytes", "podio_bytes"]:
        values = [float(row[key]) for row in records if row.get(key) is not None]
        out[f"{key}_mean"] = statistics.fmean(values) if values else None
        out[f"{key}_median"] = statistics.median(values) if values else None
        out[f"{key}_mad"] = _mad(values)
        out[f"{key}_stdev"] = statistics.stdev(values) if len(values) > 1 else 0.0 if values else None
    return out


def _command(args, threads_count: int, mode: str, podio: Path, hist: Path) -> list[str]:
    command = [
        args.eicrecon,
        *args.arg,
        f"-Pnthreads={threads_count}",
        f"-Ppodio:output_file={podio}",
        f"-Phistsfile={hist}",
    ]
    if args.events is not None:
        command.append(f"-Pjana:nevents={args.events}")
    if mode == "b0trackers":
        command.append(args.plugin_arg)
        command.extend(args.plugin_extra_arg)
    command.append(str(args.input))
    return command


def _cleanup(paths: list[Path]) -> None:
    for path in paths:
        try:
            path.unlink()
        except FileNotFoundError:
            pass


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="simulation input file")
    parser.add_argument("--eicrecon", default="eicrecon")
    parser.add_argument("--threads", default="1,2,4,8", help="comma-separated thread counts")
    parser.add_argument("--runs", type=int, default=5, help="recorded paired repetitions per thread count")
    parser.add_argument("--warmup", type=int, default=1, help="unrecorded warmup pairs per thread count")
    parser.add_argument("--seed", type=int, default=73021, help="deterministic randomization seed")
    parser.add_argument("--events", type=int, help="optional fixed jana:nevents for quicker comparisons")
    parser.add_argument(
        "--arg", action="append", default=[], help="eicrecon argument applied to baseline and plugin runs"
    )
    parser.add_argument("--plugin-arg", default="-Pplugins=B0Trackers")
    parser.add_argument(
        "--plugin-extra-arg",
        action="append",
        default=[],
        help="argument applied only when B0Trackers is enabled; repeat as needed",
    )
    parser.add_argument("--output-dir", type=Path, default=Path("b0-runtime"))
    parser.add_argument("--keep-output", action="store_true")
    args = parser.parse_args()

    if args.runs < 1:
        parser.error("--runs must be >= 1")
    if args.warmup < 0:
        parser.error("--warmup must be >= 0")
    threads = [int(item) for item in args.threads.split(",") if item.strip()]
    if not threads or any(value < 1 for value in threads):
        parser.error("--threads must contain positive integers")
    if not args.input.exists():
        parser.error(f"input does not exist: {args.input}")
    if shutil.which(args.eicrecon) is None:
        parser.error(f"cannot find eicrecon executable: {args.eicrecon}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)
    records: list[dict] = []

    for threads_count in threads:
        # Warm both paths before recording. Pair order is randomized so neither
        # mode systematically benefits from filesystem/page-cache warming.
        for warmup in range(1, args.warmup + 1):
            modes = ["baseline", "b0trackers"]
            rng.shuffle(modes)
            for mode in modes:
                stem = f"warmup-{mode}-t{threads_count}-r{warmup}"
                podio = args.output_dir / f"{stem}.edm4eic.root"
                hist = args.output_dir / f"{stem}.hists.root"
                time_json = args.output_dir / f"{stem}.time.json"
                command = _command(args, threads_count, mode, podio, hist)
                print(f"[warmup {mode} threads={threads_count}] {' '.join(command)}", flush=True)
                _run_once(command, time_json)
                _cleanup([podio, hist, time_json])

        for run in range(1, args.runs + 1):
            modes = ["baseline", "b0trackers"]
            rng.shuffle(modes)
            for pair_position, mode in enumerate(modes):
                stem = f"{mode}-t{threads_count}-r{run}"
                podio = args.output_dir / f"{stem}.edm4eic.root"
                hist = args.output_dir / f"{stem}.hists.root"
                time_json = args.output_dir / f"{stem}.time.json"
                command = _command(args, threads_count, mode, podio, hist)

                print(f"[{mode} threads={threads_count} run={run}] {' '.join(command)}", flush=True)
                metrics = _run_once(command, time_json)
                row = {
                    "mode": mode,
                    "threads": threads_count,
                    "run": run,
                    "pair_position": pair_position,
                    **metrics,
                    "podio_bytes": podio.stat().st_size if podio.exists() else None,
                    "hist_bytes": hist.stat().st_size if hist.exists() else None,
                }
                records.append(row)

                if not args.keep_output:
                    _cleanup([podio, hist, time_json])

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

        base = summary[f"baseline_t{threads_count}"]["wall_s_median"]
        plugin = summary[f"b0trackers_t{threads_count}"]["wall_s_median"]
        if base and plugin:
            summary[f"overhead_t{threads_count}"] = {
                "wall_ratio_median": plugin / base,
                "wall_overhead_fraction_median": plugin / base - 1.0,
            }

    single_thread = summary.get("b0trackers_t1", {}).get("wall_s_median")
    if single_thread:
        for threads_count in threads:
            wall = summary[f"b0trackers_t{threads_count}"]["wall_s_median"]
            if wall:
                summary[f"b0trackers_scaling_t{threads_count}"] = {
                    "speedup_vs_t1": single_thread / wall,
                    "parallel_efficiency": (single_thread / wall) / threads_count,
                }

    output = {
        "input": str(args.input),
        "input_sha256": _sha256(args.input),
        "runs_per_configuration": args.runs,
        "warmup_pairs": args.warmup,
        "random_seed": args.seed,
        "threads": threads,
        "events": args.events,
        "extra_args": args.arg,
        "plugin_arg": args.plugin_arg,
        "plugin_extra_args": args.plugin_extra_arg,
        "summary": summary,
    }
    summary_path = args.output_dir / "summary.json"
    summary_path.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")

    print(f"wrote {csv_path}")
    print(f"wrote {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
