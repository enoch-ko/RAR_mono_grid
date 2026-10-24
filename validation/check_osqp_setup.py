#!/usr/bin/env python3
"""Bounded OSQP environment checks; no galaxy fits or production solver changes.

Run with .venv-osqp/bin/python -B check_osqp_setup.py --output <new JSON>.
An installation smoke test is not validation of the inverse-RAR optimizer.
"""
import argparse
import hashlib
import importlib
import importlib.metadata as metadata
import json
import platform
import signal
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import osqp
from scipy import sparse

HERE = Path(__file__).resolve().parent
SETTINGS = dict(verbose=False, eps_abs=1e-9, eps_rel=1e-9,
                max_iter=10000, polishing=True, time_limit=5.0)
SOURCES = ("mono_SPARC_RAR.py", "utils.py", "probability_recursion.py",
           "inverse_rar_diagnostic.py", "inverse_rar_boundary.py",
           "run_inverse_rar_batch.py", "run_binned_ML_opt.py",
           "plot_inverse_rar_comparison.py")


def hashes():
    return {name: hashlib.sha256((HERE / name).read_bytes()).hexdigest()
            for name in SOURCES}


def deadline(signum, frame):
    raise TimeoutError("OSQP setup checks exceeded 30 seconds")


def check_case(name, P, q, A, lower, upper, expected=None):
    P, A = np.asarray(P, float), np.asarray(A, float)
    q, lower, upper = [np.asarray(x, float) for x in (q, lower, upper)]
    solver = osqp.OSQP()
    solver.setup(P=sparse.csc_matrix(P), q=q, A=sparse.csc_matrix(A),
                 l=lower, u=upper, **SETTINGS)
    result = solver.solve(raise_error=False)
    info = result.info
    record = {"case": name, "status": info.status,
              "status_val": int(info.status_val), "iterations": int(info.iter),
              "solver_primal_residual": float(info.prim_res),
              "solver_dual_residual": float(info.dual_res)}
    if expected is None:
        if info.status != "primal infeasible":
            raise AssertionError(f"{name}: expected primal infeasibility")
        record["passed"] = True
        return record
    if info.status != "solved" or not np.all(np.isfinite(result.x)):
        raise AssertionError(f"{name}: expected a finite solved result")
    x, y = np.asarray(result.x), np.asarray(result.y).copy()
    np.testing.assert_allclose(x, expected, rtol=0, atol=1e-7)
    ax = A @ x
    violation = max(0.0, float(np.max(lower - ax)), float(np.max(ax - upper)))
    # Infinite bounds cannot carry a multiplier. Enforce valid dual signs
    # before deriving an independent lower bound, rather than trusting status.
    y[(y > 0) & ~np.isfinite(upper)] = 0.0
    y[(y < 0) & ~np.isfinite(lower)] = 0.0
    t = q + A.T @ y
    stationarity = float(np.max(np.abs(P @ x + t)))
    fu, fl = np.isfinite(upper), np.isfinite(lower)
    bounds_term = (np.maximum(y[fu], 0) @ upper[fu]
                   + np.minimum(y[fl], 0) @ lower[fl])
    objective = float(0.5 * x @ P @ x + q @ x)
    dual_lower = float(-0.5 * t @ np.linalg.solve(P, t) - bounds_term)
    gap = objective - dual_lower
    if violation > 1e-7 or stationarity > 1e-7 or abs(gap) > 1e-7 * max(1, abs(objective)):
        raise AssertionError(f"{name}: independent feasibility/optimality check failed")
    record.update(passed=True, x=x.tolist(), objective=objective,
                  max_constraint_violation=violation, stationarity_norm=stationarity,
                  dual_lower_bound=dual_lower, primal_dual_gap=gap)
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    if Path(sys.prefix).resolve() != HERE / ".venv-osqp":
        raise RuntimeError("Use the dedicated .venv-osqp/bin/python interpreter")
    signal.signal(signal.SIGALRM, deadline)
    signal.alarm(30)
    before = hashes()
    imports = {}
    for name in ("numpy", "scipy", "matplotlib", "pandas", "vcdisk",
                 "osqp", "inverse_rar_diagnostic", "run_inverse_rar_batch",
                 "probability_recursion", "plot_inverse_rar_comparison"):
        module = importlib.import_module(name)
        imports[name] = str(Path(module.__file__).resolve())
    cases = [
        check_case("zero_deviation_reference", 2 * np.eye(2), [0, 0],
                   np.eye(2), [-1, -1], [1, 1], [0, 0]),
        check_case("nonzero_projection", 2 * np.eye(2), [0, 0],
                   [[1, 1], [1, 0], [0, 1]], [1, 0, 0],
                   [np.inf, np.inf, np.inf], [0.5, 0.5]),
        check_case("official_two_variable_example", [[4, 1], [1, 2]], [1, 1],
                   [[1, 1], [1, 0], [0, 1]], [1, 0, 0], [1, 0.7, 0.7], [0.3, 0.7]),
        check_case("infeasible_constraints", [[2]], [0], [[1], [1]],
                   [1, -np.inf], [np.inf, 0]),
    ]
    if hashes() != before:
        raise RuntimeError("Analysis source changed during setup check")
    versions = {name: metadata.version(name) for name in (
        "osqp", "joblib", "numpy", "scipy", "matplotlib", "pandas", "vcdisk",
        "jinja2", "MarkupSafe", "cloudpickle", "setuptools")}
    freeze = subprocess.run([sys.executable, "-m", "pip", "freeze"],
                            check=True, capture_output=True, text=True, timeout=10)
    report = dict(status="setup_checks_passed", diagnostic_only=True,
                  scope="imports and analytic QPs only; no galaxy fits",
                  created_utc=datetime.now(timezone.utc).isoformat(),
                  python=sys.version, executable=sys.executable, base_prefix=sys.base_prefix,
                  system_site_packages=True, platform=platform.platform(),
                  versions=versions, import_paths=imports, settings=SETTINGS,
                  cases=cases, analysis_source_sha256=before,
                  check_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  pip_freeze=freeze.stdout.splitlines(),
                  existing_plot_directory=str(HERE / "outputs/rerun_boundary_20261006/binned_ML_opt/galaxies"),
                  compute_node_validation="not performed; separately approved queue check required")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write("\n")
    signal.alarm(0)
    print(f"PASS: {len(cases)} analytic QP checks and {len(imports)} imports")
    print(args.output)


if __name__ == "__main__":
    main()
