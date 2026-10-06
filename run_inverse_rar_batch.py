#!/usr/bin/env python
"""Shard the exploratory 175-galaxy inverse RAR diagnostic; no queue submission.

Dry-run prints the complete deterministic assignment without computing disc fields.
Each queued shard writes disjoint galaxy JSON/PNG pairs and a resumable report.
"""
import argparse
import hashlib
import json
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from utils import DEFAULT_DENS_ROOT, get_SBdisk, load_sparc_table

HERE = Path(__file__).resolve().parent
DEFAULT_ROOT = HERE / "outputs/inverse_rar_diagnostic/full175_densknots_p400_q9_ml10_retry20"
DRIVER = HERE / "inverse_rar_diagnostic.py"
PLOTTER = HERE / "plot_inverse_rar_comparison.py"


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def assignments(shards, dens_root=DEFAULT_DENS_ROOT):
    """Largest-first greedy allocation by knot count squared; stable tie breaks."""
    table = load_sparc_table()
    names = list(table.Galaxy)
    if len(names) != 175 or len(set(names)) != 175:
        raise ValueError("expected exactly 175 unique SPARC galaxies")
    items = []
    for name in names:
        rad, _ = get_SBdisk(name, dens_root=dens_root, table=table)
        if not (HERE.parent / "data" / f"{name}_rotmod.dat").is_file():
            # The actual data root is checked below using the diagnostic default.
            pass
        items.append((len(rad) ** 2, str(name), len(rad)))
    bins = [[] for _ in range(shards)]
    costs = [0] * shards
    for cost, name, knots in sorted(items, key=lambda x: (-x[0], x[1])):
        target = min(range(shards), key=lambda i: (costs[i], len(bins[i]), i))
        bins[target].append({"galaxy": name, "knots": knots, "weight": cost})
        costs[target] += cost
    if len({x["galaxy"] for group in bins for x in group}) != 175:
        raise ValueError("shard assignment does not cover the sample once")
    return bins, costs


def existing_ok(path, galaxy, driver_hash, plot_hash):
    if not path.is_file() or not path.with_suffix(".png").is_file():
        return False
    try:
        obj = json.loads(path.read_text())
        provenance = obj["provenance"]
        settings = provenance["args"]
        return (obj["galaxy"] == galaxy
                and provenance["script_sha256"] == driver_hash
                and provenance["plot_script_sha256"] == plot_hash
                and settings["knots_from_dens"] == "True"
                and settings["prob_grid_size"] == "400"
                and settings["ml_samples"] == "9"
                and settings["max_ml"] == "10.0"
                and settings["retry_max_ml"] == "20.0"
                and obj["plot"]["status"] == "created"
                and obj["plot"]["path"] == str(path.with_suffix(".png")))
    except (OSError, ValueError, KeyError, TypeError):
        return False


def save_report(path, report):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(report, indent=2) + "\n")
    temporary.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shard", type=int, help="one-based shard number")
    parser.add_argument("--shards", type=int, default=10)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--audit", action="store_true", help="read-only 175-output audit")
    args = parser.parse_args()
    if args.shards < 1 or args.shards > 175:
        parser.error("shards must be 1..175")
    if args.shard is not None and not 1 <= args.shard <= args.shards:
        parser.error("shard must be in 1..shards")
    if not (args.dry_run or args.audit) and args.shard is None:
        parser.error("--shard is required for computation")
    bins, costs = assignments(args.shards)
    selected = range(args.shards) if args.shard is None else [args.shard - 1]
    if args.dry_run:
        for index in selected:
            print(f"shard {index + 1:02d}: {len(bins[index])} galaxies, weight {costs[index]}; "
                  + ",".join(x["galaxy"] for x in bins[index]))
        return 0
    driver_hash, plot_hash = digest(DRIVER), digest(PLOTTER)
    if args.audit:
        found = {"valid": [], "missing_or_mismatched": []}
        for group in bins:
            for item in group:
                galaxy = item["galaxy"]
                path = args.output_root / "galaxies" / f"{galaxy}_inverse.json"
                key = "valid" if existing_ok(path, galaxy, driver_hash, plot_hash) else "missing_or_mismatched"
                found[key].append(galaxy)
        print(json.dumps({k: {"count": len(v), "galaxies": v} for k, v in found.items()}, indent=2))
        return 0 if len(found["valid"]) == 175 else 1
    index = args.shard - 1
    report_path = args.output_root / f"shard_{args.shard:02d}.json"
    report = {"created_utc": datetime.now(timezone.utc).isoformat(),
              "shard": args.shard, "shards": args.shards,
              "galaxies": [x["galaxy"] for x in bins[index]],
              "estimated_weight": costs[index],
              "driver_sha256": driver_hash, "plotter_sha256": plot_hash,
              "catalogue_sha256": digest(HERE / "SPARC_table.txt"),
              "settings": {"knots": "all .dens radii", "max_ml": 10.0,
                           "retry_max_ml": 20.0, "prob_grid_size": 400,
                           "ml_samples": 9}, "results": {}}
    failures = 0
    for item in bins[index]:
        galaxy = item["galaxy"]
        path = args.output_root / "galaxies" / f"{galaxy}_inverse.json"
        if existing_ok(path, galaxy, driver_hash, plot_hash):
            report["results"][galaxy] = {"status": "existing_valid", "path": str(path)}
            print(f"{galaxy}: existing valid pair", flush=True)
            save_report(report_path, report)
            continue
        if path.exists() or path.with_suffix(".png").exists():
            report["results"][galaxy] = {"status": "collision_needs_review", "path": str(path)}
            failures += 1
            print(f"{galaxy}: existing pair needs review; no overwrite", flush=True)
            save_report(report_path, report)
            continue
        path.parent.mkdir(parents=True, exist_ok=True)
        command = [sys.executable, str(DRIVER), "--galaxy", galaxy,
                   "--output", str(path), "--prob-grid-size", "400",
                   "--ml-samples", "9", "--max-ml", "10",
                   "--retry-max-ml", "20", "--knots-from-dens"]
        started = time.monotonic()
        completed = subprocess.run(command, cwd=HERE, text=True,
                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        elapsed = round(time.monotonic() - started, 2)
        if completed.returncode == 0 and existing_ok(path, galaxy, driver_hash, plot_hash):
            status = json.loads(path.read_text())["status"]
            report["results"][galaxy] = {"status": status, "seconds": elapsed,
                                          "path": str(path)}
            print(f"{galaxy}: {status} ({elapsed}s)", flush=True)
        else:
            failures += 1
            report["results"][galaxy] = {"status": "error", "seconds": elapsed,
                                          "returncode": completed.returncode,
                                          "stdout": completed.stdout[-4000:],
                                          "stderr": completed.stderr[-8000:],
                                          "path": str(path)}
            print(f"{galaxy}: error ({elapsed}s); see {report_path}", flush=True)
        save_report(report_path, report)
    print(f"shard {args.shard:02d}: {len(bins[index])} attempted, {failures} errors", flush=True)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
