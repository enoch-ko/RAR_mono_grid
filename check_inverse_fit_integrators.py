#!/usr/bin/env python3
"""Bounded checks 1 and 4 only: corrected fits and same-profile probabilities.
No nuisance marginalisation or tail/quadrature convergence claim.
"""
import argparse
import json
import signal
import time
import warnings
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import inverse_rar_diagnostic as inv
import probability_recursion as boundary
import utils


def safe(value):
    if isinstance(value, dict):
        return {k: safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [safe(v) for v in value]
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.bool_):
        return bool(value)
    return value


def expired(signum, frame):
    raise TimeoutError("bounded login-node time budget exceeded")


def hashes():
    return {Path(p).name: inv.sha256(p) for p in
            (__file__, inv.__file__, inv.mono.__file__, utils.__file__, boundary.__file__)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--galaxy", required=True)
    parser.add_argument("--compare-only", action="store_true")
    parser.add_argument("--seconds", type=int, default=40)
    parser.add_argument("--output-root", type=Path, default=Path(utils.DEFAULT_PROJECT_ROOT) /
                        "outputs/inverse_rar_diagnostic/comparison_cholesky/checks_1_4_20261006")
    args = parser.parse_args()
    if args.seconds < 1:
        parser.error("positive time budget required")
    args.output_root.mkdir(parents=True, exist_ok=True)
    path = args.output_root / (args.galaxy + ("_probabilities.json" if args.compare_only else "_fit.json"))
    if path.exists():
        raise FileExistsError(path)
    result = {"galaxy": args.galaxy, "source_sha256": hashes(), "diagnostic_only": True,
              "deferred": ["nuisance marginalisation", "tail and quadrature convergence"],
              "created_unix": time.time(), "budget_seconds": args.seconds}
    start = time.monotonic()
    signal.signal(signal.SIGALRM, expired)
    signal.alarm(args.seconds)
    try:
        if not args.compare_only:
            result["cache_gate"] = inv.cache_gate(args.galaxy, Path(utils.DEFAULT_PROJECT_ROOT) / "outputs")
            rad, sb, r, v, e, gas, bulge, zd, excluded = inv.galaxy_data(
                args.galaxy, inv.mono.DEFAULT_DATA_ROOT, utils.DEFAULT_DENS_ROOT)
            result["input_sha256"] = {
                "rotmod": inv.sha256(Path(inv.mono.DEFAULT_DATA_ROOT) / (args.galaxy + "_rotmod.dat")),
                "dens": inv.sha256(Path(utils.DEFAULT_DENS_ROOT) / (args.galaxy + ".dens")),
                "table": inv.sha256(utils.DEFAULT_TABLE_PATH)}
            H, D, seen = inv.build_response(rad, sb, r, zd, rad)
            result["prior_direct_cache_qc"] = inv.direct_check(rad, sb, r, zd, H, D, np.full(len(rad), .5))
            c = utils.vel2acc(gas * np.abs(gas) + .7 * bulge ** 2, r)
            observed = utils.vel2acc(v ** 2, r)
            fit, attempts, upper = inv.solve_with_retry(D, c, observed, rad,
                SimpleNamespace(min_ml=.01, max_ml=10., retry_max_ml=20., length_fraction=.3))
            result.update(status=fit["status"], fit=fit, attempts=attempts,
                          upper_ml=upper, fixed_gas_ml=1., fixed_bulge_ml=.7,
                          knots_kpc=rad, r=r, vobs=v, errv=e, excluded_rows=excluded)
            if fit["status"] == "feasible_weak_order":
                result["optimized_direct_cache_qc"] = inv.direct_check(rad, sb, r, zd, H, D, fit["xi"])
                bad = np.flatnonzero(~np.isfinite(fit["g_bar"]) | (fit["g_bar"] <= 0))
                result["nonpositive_gbar_indices"] = bad
                result["nonpositive_gbar_radii_kpc"] = r[bad]
                result["physical_positive"] = len(bad) == 0
                result["profile_min_max"] = [np.min(H @ fit["xi"]), np.max(H @ fit["xi"])]
        else:
            fit_path = args.output_root / (args.galaxy + "_fit.json")
            saved = json.loads(fit_path.read_text())
            if saved["source_sha256"] != hashes():
                raise ValueError("fit provenance does not match current code")
            result["fit_sha256"] = inv.sha256(fit_path)
            if saved.get("status") != "feasible_weak_order" or not saved.get("physical_positive"):
                result["status"] = "skipped_invalid_fit"
            else:
                arrays = [np.asarray(saved["fit"]["g_bar"]), np.asarray(saved["vobs"]),
                          np.asarray(saved["errv"]), np.asarray(saved["r"])]
                result.update(status="comparison", ngrid=400, zpadding=8., tie_atol=0.,
                              max_tie_permutations=720, same_profile=True,
                              probability_kind="conditional; gas=1, bulge=0.7")
                for name, fun in [("legacy", utils.log_monotonic_probability_for_gbar),
                                  ("boundary", boundary.boundary_log_probability_for_gbar)]:
                    call_start = time.monotonic()
                    try:
                        with warnings.catch_warnings(record=True) as seen:
                            warnings.simplefilter("always")
                            lp = float(fun(*arrays, ngrid=400, zpadding=8., tie_atol=0., max_tie_permutations=720))
                        result[name] = {"logp": lp, "finite": np.isfinite(lp),
                                        "warnings": [str(w.message) for w in seen]}
                    except Exception as exc:
                        result[name] = {"finite": False, "error": type(exc).__name__ + ": " + str(exc)}
                    result[name]["seconds"] = time.monotonic() - call_start
                if all(result[n]["finite"] for n in ("legacy", "boundary")):
                    result["delta_logp_boundary_minus_legacy"] = result["boundary"]["logp"] - result["legacy"]["logp"]
    except Exception as exc:
        result.update(status="error", error=type(exc).__name__ + ": " + str(exc))
    finally:
        signal.alarm(0)
    result["seconds"] = time.monotonic() - start
    path.write_text(json.dumps(safe(result), indent=2, allow_nan=False) + "\n")
    print(json.dumps(safe({k: v for k, v in result.items() if k in
          ("galaxy", "status", "error", "seconds", "physical_positive", "nonpositive_gbar_radii_kpc",
           "profile_min_max", "legacy", "boundary", "delta_logp_boundary_minus_legacy")})), flush=True)


if __name__ == "__main__":
    main()
