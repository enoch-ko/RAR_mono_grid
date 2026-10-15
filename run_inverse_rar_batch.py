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
DEFAULT_ROOT = HERE / "outputs/inverse_rar_diagnostic/comparison_cholesky"
DRIVER = HERE / "inverse_rar_boundary.py"
INVERSE = HERE / "inverse_rar_diagnostic.py"
PROBABILITY = HERE / "probability_recursion.py"
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


def run_settings(args):
    return {"integration": args.integration, "knots_from_dens": "True",
            "prob_grid_size": str(args.prob_grid_size), "ml_samples": str(args.ml_samples),
            "max_ml": "10.0", "retry_max_ml": "20.0", "min_ml": "0.01",
            "length_fraction": "0.3", "gas_ml": "1.0", "bulge_ml": "0.7",
            "gas_ml_dex": "0.04", "bulge_ml_dex": "0.1", "zpadding": str(args.zpadding),
            "boundary_log_tolerance": args.boundary_log_tolerance,
            "boundary_max_refinements": args.boundary_max_refinements,
            "qp_backend": getattr(args, "qp_backend", "slsqp"),
            "osqp_eps": str(getattr(args, "osqp_eps", 1e-9)),
            "osqp_max_iter": str(getattr(args, "osqp_max_iter", 100000)),
            "osqp_time_limit": str(getattr(args, "osqp_time_limit", 10.0))}


def source_hashes():
    return {name: digest(HERE / name) for name in (
        "inverse_rar_diagnostic.py", "inverse_rar_boundary.py", "probability_recursion.py",
        "mono_SPARC_RAR.py", "utils.py", "plot_inverse_rar_comparison.py", "SPARC_table.txt")}


def command_for(galaxy, path, args):
    driver = DRIVER if args.integration == "boundary" else INVERSE
    command = [sys.executable, str(driver), "--galaxy", galaxy, "--output", str(path),
               "--prob-grid-size", str(args.prob_grid_size), "--ml-samples", str(args.ml_samples),
               "--max-ml", "10", "--retry-max-ml", "20", "--knots-from-dens",
               "--zpadding", str(args.zpadding)]
    command += ["--qp-backend", getattr(args, "qp_backend", "slsqp"),
                "--osqp-eps", str(getattr(args, "osqp_eps", 1e-9)),
                "--osqp-max-iter", str(getattr(args, "osqp_max_iter", 100000)),
                "--osqp-time-limit", str(getattr(args, "osqp_time_limit", 10.0))]
    if args.integration == "boundary":
        command += ["--boundary-log-tolerance", str(args.boundary_log_tolerance),
                    "--boundary-max-refinements", str(args.boundary_max_refinements)]
    return command


def existing_ok(path, galaxy, settings, hashes):
    """Matching artifacts, not a scientific convergence endorsement."""
    if not path.is_file() or not path.with_suffix(".png").is_file():
        return False
    try:
        obj = json.loads(path.read_text())
        p, stored = obj["provenance"], obj["provenance"]["args"]
        if not (obj["galaxy"] == galaxy
                and p["script_sha256"] == hashes["inverse_rar_diagnostic.py"]
                and p["driver_sha256"] == hashes["mono_SPARC_RAR.py"]
                and p["plot_script_sha256"] == hashes["plot_inverse_rar_comparison.py"]
                and obj["input_sha256"]["table"] == hashes["SPARC_table.txt"]
                and all(stored[k] == v for k, v in settings.items()
                        if k not in ("integration", "boundary_log_tolerance", "boundary_max_refinements"))
                and obj["plot"]["status"] == "created"
                and obj["plot"]["path"] == str(path.with_suffix(".png"))):
            return False
        if settings["integration"] == "boundary":
            numeric = obj["probability_numerics"]
            return (p["boundary_entry_sha256"] == hashes["inverse_rar_boundary.py"]
                    and p["probability_recursion_sha256"] == hashes["probability_recursion.py"]
                    and p["utils_sha256"] == hashes["utils.py"]
                    and numeric["method"] == "boundary_log_cumulative_v1"
                    and numeric["log_tolerance"] == settings["boundary_log_tolerance"]
                    and numeric["max_refinements"] == settings["boundary_max_refinements"])
        return "probability_numerics" not in obj and p["utils_sha256"] == hashes["utils.py"]
    except (OSError, ValueError, KeyError, TypeError):
        return False


def save_report(path, report):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shard", type=int, help="one-based shard number")
    parser.add_argument("--shards", type=int, default=10)
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--galaxies", nargs="+", help="restrict to named comparison galaxies")
    parser.add_argument("--integration", choices=("boundary", "legacy"), default="boundary")
    parser.add_argument("--prob-grid-size", type=int, default=400,
                        help="starting resolution for boundary, fixed resolution for legacy")
    parser.add_argument("--ml-samples", type=int, default=9)
    parser.add_argument("--zpadding", type=float, default=8.)
    parser.add_argument("--boundary-log-tolerance", type=float, default=2e-5)
    parser.add_argument("--boundary-max-refinements", type=int, default=8)
    parser.add_argument("--qp-backend", choices=("slsqp", "osqp"), default="slsqp")
    parser.add_argument("--osqp-eps", type=float, default=1e-9)
    parser.add_argument("--osqp-max-iter", type=int, default=100000)
    parser.add_argument("--osqp-time-limit", type=float, default=10.0)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--audit", action="store_true", help="read-only selected-output audit")
    args = parser.parse_args()
    if not 1 <= args.shards <= 175 or (args.shard is not None and not 1 <= args.shard <= args.shards):
        parser.error("shard/shards must define a valid allocation")
    if args.prob_grid_size < 2 or args.ml_samples < 1 or args.boundary_max_refinements < 1:
        parser.error("invalid probability grid, quadrature or refinements")
    if not (0 < args.zpadding < float("inf") and 0 < args.boundary_log_tolerance < float("inf")):
        parser.error("finite positive padding and tolerance required")
    if not (args.dry_run or args.audit) and args.shard is None:
        parser.error("--shard is required for computation")
    if (args.osqp_eps <= 0 or args.osqp_eps >= float("inf") or args.osqp_max_iter < 1
            or not 0 < args.osqp_time_limit < float("inf")):
        parser.error("positive finite OSQP controls required")
    if not (args.dry_run or args.audit) and args.qp_backend == "osqp":
        import importlib.util
        if importlib.util.find_spec("osqp") is None:
            parser.error("OSQP unavailable in this interpreter; use the isolated .venv-osqp environment")
    if args.output_root is None:
        tag = (f"{args.integration}_p{args.prob_grid_size}_q{args.ml_samples}"
               f"_z{args.zpadding:g}_tol{args.boundary_log_tolerance:g}_ref{args.boundary_max_refinements}"
               f"_qp{args.qp_backend}")
        args.output_root = DEFAULT_ROOT / tag
    bins, costs = assignments(args.shards)
    if args.galaxies:
        requested = set(args.galaxies)
        known = {x["galaxy"] for group in bins for x in group}
        if len(requested) != len(args.galaxies) or requested - known:
            parser.error("duplicate or unknown comparison galaxies")
        bins = [[x for x in group if x["galaxy"] in requested] for group in bins]
        costs = [sum(x["weight"] for x in group) for group in bins]
    selected = range(args.shards) if args.shard is None else [args.shard - 1]
    settings, hashes = run_settings(args), source_hashes()
    if args.dry_run:
        print(json.dumps({"output_root": str(args.output_root), "settings": settings,
                          "source_sha256": hashes}, indent=2))
        for index in selected:
            print(f"shard {index + 1:02d}: {len(bins[index])} galaxies, weight {costs[index]}")
            for item in bins[index]:
                path = args.output_root / "galaxies" / (item["galaxy"] + "_inverse.json")
                print(json.dumps(command_for(item["galaxy"], path, args)))
        return 0
    if args.audit:
        found = {"matching": [], "missing_or_mismatched": [], "failed_diagnostics": []}
        for index in selected:
            for item in bins[index]:
                galaxy = item["galaxy"]
                path = args.output_root / "galaxies" / (galaxy + "_inverse.json")
                ok = existing_ok(path, galaxy, settings, hashes)
                found["matching" if ok else "missing_or_mismatched"].append(galaxy)
                if ok and json.loads(path.read_text())["status"] not in ("feasible_weak_order", "infeasible"):
                    found["failed_diagnostics"].append(galaxy)
        print(json.dumps({k: {"count": len(v), "galaxies": v} for k, v in found.items()}, indent=2))
        return int(bool(found["missing_or_mismatched"] or found["failed_diagnostics"]))
    index = args.shard - 1
    report_path = args.output_root / f"shard_{args.shard:02d}.json"
    names = [x["galaxy"] for x in bins[index]]
    if report_path.exists():
        old = json.loads(report_path.read_text())
        if old.get("settings") != settings or old.get("source_sha256") != hashes or old.get("galaxies") != names:
            raise ValueError("existing shard has different settings/source/selection; use a new output root")
    report = {"created_utc": datetime.now(timezone.utc).isoformat(), "shard": args.shard,
              "shards": args.shards, "galaxies": names, "estimated_weight": costs[index],
              "source_sha256": hashes, "settings": settings, "results": {}}
    failures = 0
    for item in bins[index]:
        if source_hashes() != hashes:
            raise RuntimeError("source changed during shard; stopping to prevent mixed algorithms")
        galaxy = item["galaxy"]
        path = args.output_root / "galaxies" / (galaxy + "_inverse.json")
        if existing_ok(path, galaxy, settings, hashes):
            status = json.loads(path.read_text())["status"]
            failures += int(status not in ("feasible_weak_order", "infeasible"))
            report["results"][galaxy] = {"status": status, "existing": True, "path": str(path)}
        elif path.exists() or path.with_suffix(".png").exists():
            failures += 1
            report["results"][galaxy] = {"status": "collision_needs_review", "path": str(path)}
        else:
            path.parent.mkdir(parents=True, exist_ok=True)
            started = time.monotonic()
            completed = subprocess.run(command_for(galaxy, path, args), cwd=HERE, text=True,
                                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            elapsed = round(time.monotonic() - started, 2)
            if completed.returncode == 0 and existing_ok(path, galaxy, settings, hashes):
                status = json.loads(path.read_text())["status"]
                failures += int(status not in ("feasible_weak_order", "infeasible"))
                report["results"][galaxy] = {"status": status, "seconds": elapsed, "path": str(path)}
            else:
                failures += 1
                report["results"][galaxy] = {"status": "error", "seconds": elapsed,
                    "returncode": completed.returncode, "stdout": completed.stdout[-4000:],
                    "stderr": completed.stderr[-8000:], "path": str(path)}
        print(f"{galaxy}: {report['results'][galaxy]['status']}", flush=True)
        save_report(report_path, report)
    print(f"shard {args.shard:02d}: {len(names)} attempted, {failures} errors", flush=True)
    return int(bool(failures))


if __name__ == "__main__":
    sys.exit(main())
