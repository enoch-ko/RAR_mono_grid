#!/usr/bin/env python3
"""Bounded, reproducible OSQP/SLSQP fit comparison; no probability calculation.
Run one galaxy at a time in .venv-osqp; full-sample runs require queue approval.
"""
import argparse
import importlib.metadata
import json
import signal
import sys
import time
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import inverse_rar_diagnostic as inv
import utils


def expired(signum, frame):
    raise TimeoutError("bounded login-node time budget exceeded")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--galaxy", required=True)
    parser.add_argument("--seconds", type=int, default=40)
    parser.add_argument("--osqp-eps", type=float, default=1e-9)
    parser.add_argument("--osqp-max-iter", type=int, default=100000)
    parser.add_argument("--osqp-time-limit", type=float, default=10.)
    parser.add_argument("--output-root", type=Path, default=Path(utils.DEFAULT_PROJECT_ROOT) /
                        "outputs/binned_ML_opt/osqp_comparison_20261006")
    args = parser.parse_args()
    if args.seconds < 1:
        parser.error("positive time budget required")
    if Path(sys.prefix).name != ".venv-osqp":
        parser.error("use .venv-osqp/bin/python to retain the verified dependency environment")
    path = args.output_root / (args.galaxy + "_comparison.json")
    basis_path = args.output_root / (args.galaxy + "_response.npz")
    if path.exists() or basis_path.exists():
        raise FileExistsError(path)
    args.output_root.mkdir(parents=True, exist_ok=True)
    result = {"galaxy": args.galaxy, "diagnostic_only": True, "fit_only": True,
              "created_unix": time.time(), "budget_seconds": args.seconds,
              "versions": {k: importlib.metadata.version(k) for k in
                           ("numpy", "scipy", "osqp", "vcdisk")},
              "python": sys.version, "executable": sys.executable,
              "source_sha256": {Path(p).name: inv.sha256(p) for p in
                                (__file__, inv.__file__, inv.mono.__file__, utils.__file__)},
              "settings": {"min_ml": .01, "max_ml": 10., "retry_max_ml": 20.,
                           "length_fraction": .3, "knots_from_dens": True,
                           "gas_ml": 1., "bulge_ml": .7, "gbar_floor_relative": 1e-8,
                           "osqp_eps": args.osqp_eps, "osqp_max_iter": args.osqp_max_iter,
                           "osqp_time_limit": args.osqp_time_limit},
              "deferred": ["post-fit probabilities", "full sample", "compute-node environment smoke"]}
    start = time.monotonic()
    signal.signal(signal.SIGALRM, expired)
    signal.alarm(args.seconds)
    try:
        result["cache_gate"] = inv.cache_gate(args.galaxy, Path(utils.DEFAULT_PROJECT_ROOT) / "outputs")
        rad, sb, r, v, e, gas, bulge, zd, excluded = inv.galaxy_data(
            args.galaxy, inv.mono.DEFAULT_DATA_ROOT, utils.DEFAULT_DENS_ROOT)
        result["input_sha256"] = {
            "rotmod": inv.sha256(Path(inv.mono.DEFAULT_DATA_ROOT) / (args.galaxy + "_rotmod.dat")),
            "dens": inv.sha256(Path(utils.DEFAULT_DENS_ROOT) / (args.galaxy + ".dens")),
            "table": inv.sha256(utils.DEFAULT_TABLE_PATH)}
        result.update(knots_kpc=rad, retained_radii_kpc=r, excluded_rows=excluded)
        H, D, seen = inv.build_response(rad, sb, r, zd, rad)
        result["basis_seconds"] = time.monotonic() - start
        result["prior_direct_cache_qc"] = inv.direct_check(rad, sb, r, zd, H, D, np.full(len(rad), .5))
        c = utils.vel2acc(gas * np.abs(gas) + .7 * bulge ** 2, r)
        observed = utils.vel2acc(v ** 2, r)
        # A single immutable response for both solvers also makes follow-up
        # numerical investigations inexpensive, without new Casertano solves.
        np.savez_compressed(basis_path, H=H, D=D, c=c, observed=observed, knots=rad,
                            r=r, brightness=sb, zd=zd)
        result["response"] = {"path": str(basis_path), "sha256": inv.sha256(basis_path)}
        result["backends"] = {}
        for backend in ("slsqp", "osqp"):
            options = SimpleNamespace(**dict(result["settings"], qp_backend=backend))
            before = time.monotonic()
            fit, attempts, upper = inv.solve_with_retry(D, c, observed, rad, options)
            record = {"fit": fit, "attempts": attempts, "upper_ml": upper,
                      "seconds": time.monotonic() - before}
            result["backends"][backend] = record
            if fit["status"] == "feasible_weak_order":
                record["direct_cache_qc"] = inv.direct_check(rad, sb, r, zd, H, D, fit["xi"])
                record["profile_min_max"] = [float((H @ fit["xi"]).min()), float((H @ fit["xi"]).max())]
                record["strict_positive_gbar"] = bool(np.all(fit["g_bar"] > 0))
            print(args.galaxy, backend, fit["status"],
                  fit.get("objective"), "seconds", round(record["seconds"], 3), flush=True)
        old_path = (Path(utils.DEFAULT_PROJECT_ROOT) / "outputs/rerun_boundary_20261006/"
                    "binned_ML_opt/galaxies" / (args.galaxy + "_inverse.json"))
        if old_path.exists():
            old = json.loads(old_path.read_text())
            result["previous_artifact"] = {"path": str(old_path), "sha256": inv.sha256(old_path),
                "status": old["status"], "objective": old["optimization"].get("objective")}
        fits = [result["backends"][b]["fit"] for b in ("slsqp", "osqp")]
        if all(f["status"] == "feasible_weak_order" for f in fits):
            result["comparison"] = {
                "objective_osqp_minus_slsqp": fits[1]["objective"] - fits[0]["objective"],
                "max_abs_xi_delta": float(np.max(np.abs(fits[1]["xi"] - fits[0]["xi"]))),
                "max_abs_gbar_delta_m_s2": float(np.max(np.abs(fits[1]["g_bar"] - fits[0]["g_bar"])))}
        result["status"] = "completed_fit_comparison"
    except Exception as exc:
        result.update(status="error", error=type(exc).__name__ + ": " + str(exc))
    finally:
        signal.alarm(0)
    result["seconds"] = time.monotonic() - start
    path.write_text(json.dumps(inv.json_safe(result), indent=2, allow_nan=False) + "\n")
    print(result["status"], path, "total_seconds", round(result["seconds"], 3), flush=True)
    return int(result["status"] == "error")


if __name__ == "__main__":
    sys.exit(main())
