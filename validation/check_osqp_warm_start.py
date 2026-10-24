#!/mnt/users/koe/RAR_mono_grid/.venv-osqp/bin/python
"""Diagnostic-only saved-response cold/warm OSQP comparison; no probabilities.
The saved SLSQP candidate is a seed only. Existing acceptance checks are retained.
"""
import os
for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[name] = "1"
import argparse
import importlib.metadata
import json
import signal
import sys
import time
from pathlib import Path
from unittest.mock import patch
import numpy as np
from scipy.linalg import solve_triangular
import inverse_rar_diagnostic as inv
import utils


def solve_seeded(D, c, observed, knots, settings, seed, budget, limit, eps=1e-9):
    """Substitute the initial iterate only; delegate all audits to solve_qp."""
    original = inv.osqp_solve
    initial = {}
    def seeded(C, rhs, z0, eps, max_iter, time_limit):
        z = z0 if seed is None else np.asarray(seed, float)
        if z.shape != z0.shape or not np.all(np.isfinite(z)):
            raise ValueError("invalid saved warm-start candidate")
        initial.update(kind="lp" if seed is None else "saved_slsqp_candidate",
                       objective=float(z @ z), max_scaled_slack_violation=float(max(0., -np.min(C @ z - rhs))))
        return original(C, rhs, z, eps, max_iter, time_limit)
    with patch.object(inv, "osqp_solve", seeded):
        fit = inv.solve_qp(D, c, observed, knots, settings["min_ml"], settings["max_ml"],
                           settings["length_fraction"], gbar_floor_relative=settings["gbar_floor_relative"],
                           qp_backend="osqp", osqp_eps=eps, osqp_max_iter=budget, osqp_time_limit=limit)
    return fit, initial


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--galaxy", required=True)
    p.add_argument("--budgets", type=int, nargs="+", default=[1000000])
    p.add_argument("--osqp-eps", type=float, default=1e-9)
    p.add_argument("--time-limit", type=float, default=10.)
    p.add_argument("--seconds", type=int, default=40)
    p.add_argument("--input-root", type=Path, default=Path(utils.DEFAULT_PROJECT_ROOT) / "outputs/binned_ML_opt/osqp_refined_20261006")
    p.add_argument("--output-root", type=Path, required=True)
    args = p.parse_args()
    if Path(sys.prefix).name != ".venv-osqp" or min(args.budgets) < 1 or not np.isfinite(args.osqp_eps) or args.osqp_eps <= 0 or not np.isfinite(args.time_limit) or args.time_limit <= 0 or args.seconds < 1:
        p.error("dedicated .venv-osqp interpreter and positive budgets required")
    path = args.output_root / (args.galaxy + "_warm_start.json")
    if path.exists():
        raise FileExistsError(path)
    args.output_root.mkdir(parents=True, exist_ok=True)
    comparison = args.input_root / (args.galaxy + "_comparison.json")
    basis = args.input_root / (args.galaxy + "_response.npz")
    saved = json.loads(comparison.read_text())
    if saved["galaxy"] != args.galaxy or inv.sha256(basis) != saved["response"]["sha256"]:
        raise ValueError("wrong galaxy or stale response")
    for source in (inv.__file__, inv.mono.__file__, utils.__file__):
        if inv.sha256(source) != saved["source_sha256"][Path(source).name]:
            raise ValueError("saved seed/model source mismatch: " + source)
    current_inputs = {"rotmod": inv.sha256(Path(inv.mono.DEFAULT_DATA_ROOT) / (args.galaxy + "_rotmod.dat")),
                      "dens": inv.sha256(Path(utils.DEFAULT_DENS_ROOT) / (args.galaxy + ".dens")),
                      "table": inv.sha256(utils.DEFAULT_TABLE_PATH)}
    if current_inputs != saved["input_sha256"]:
        raise ValueError("saved input mismatch")
    result = {"galaxy": args.galaxy, "diagnostic_only": True, "fit_only": True,
              "created_unix": time.time(), "arguments": vars(args), "input_sha256": current_inputs,
              "response_sha256": inv.sha256(basis), "seed_json_sha256": inv.sha256(comparison),
              "source_sha256": {Path(f).name: inv.sha256(f) for f in (__file__, inv.__file__, inv.mono.__file__, utils.__file__)},
              "versions": {k: importlib.metadata.version(k) for k in ("numpy", "scipy", "osqp", "vcdisk")},
              "python": sys.version, "executable": sys.executable, "runs": []}
    result["arguments"] = {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()}
    start = time.monotonic()
    def expired(signum, frame):
        raise TimeoutError("overall diagnostic budget exceeded")
    signal.signal(signal.SIGALRM, expired)
    signal.alarm(args.seconds)
    def checkpoint():
        result["seconds"] = time.monotonic() - start
        path.write_text(json.dumps(inv.json_safe(result), indent=2, allow_nan=False) + "\n")
    try:
        result["cache_gate"] = inv.cache_gate(args.galaxy, Path(utils.DEFAULT_PROJECT_ROOT) / "outputs")
        data = np.load(basis, allow_pickle=False)
        D, c, observed, knots = [data[k] for k in ("D", "c", "observed", "knots")]
        candidate = saved["backends"]["slsqp"]["fit"]
        diagnostics = candidate["solver_diagnostics"]
        seed = diagnostics.get("candidate_z")
        if seed is None and "xi" in candidate:
            seed = solve_triangular(np.linalg.cholesky(inv.kernel(knots, saved["settings"]["length_fraction"])), np.asarray(candidate["xi"]) - .5, lower=True)
        if seed is None:
            raise ValueError("no actual SLSQP candidate available")
        seed = np.asarray(seed, float)
        result["seed_status"] = candidate["status"]
        result["seed_optimality"] = candidate.get("optimality", diagnostics.get("optimality"))
        result["prior_direct_cache_qc"] = inv.direct_check(knots, data["brightness"], data["r"], float(data["zd"]), data["H"], D, np.full(len(knots), .5))
        for budget in args.budgets:
            for mode in ("cold_lp", "warm_slsqp"):
                before = time.monotonic()
                fit, initial = solve_seeded(D, c, observed, knots, saved["settings"], None if mode == "cold_lp" else seed, budget, args.time_limit, args.osqp_eps)
                record = {"mode": mode, "iteration_budget": budget, "time_limit": args.time_limit,
                          "initial": initial, "fit": fit, "seconds": time.monotonic() - before}
                result["runs"].append(record)
                if fit["status"] == "feasible_weak_order":
                    record["direct_cache_qc"] = inv.direct_check(knots, data["brightness"], data["r"], float(data["zd"]), data["H"], D, fit["xi"])
                checkpoint()
                sd = fit["solver_diagnostics"]
                print(args.galaxy, mode, budget, fit["status"], sd.get("solver", {}).get("iterations"), sd.get("objective"), round(record["seconds"], 3), flush=True)
        result["status"] = "completed_fit_comparison"
    except Exception as exc:
        result.update(status="error", error=type(exc).__name__ + ": " + str(exc))
    finally:
        signal.alarm(0)
        checkpoint()
    print(result["status"], path, flush=True)
    return int(result["status"] == "error")


if __name__ == "__main__":
    sys.exit(main())