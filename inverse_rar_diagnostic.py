#!/usr/bin/env python3
"""First-order inverse RAR diagnostic; not a validated uncertainty-aware optimum.

Find the bounded piecewise-linear disc M/L closest to 0.5 in a normalized
squared-exponential-kernel metric, subject to weak central-value RAR ordering.
The common covariance amplitude cancels from this hard-constrained optimum.
Gas/bulge M/L are fixed in the QP and marginalised only after optimization.
"""
import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, linprog, minimize

import mono_SPARC_RAR as mono
from utils import (DEFAULT_DENS_ROOT, DEFAULT_PROJECT_ROOT, DEFAULT_TABLE_PATH, get_SBdisk,
                   load_sparc_table, lognormal_quadrature, vel2acc)


VALIDATION_GLOB = "real_disk_cache_all175_p101_*_20261005.json"


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def profile_knots(rad, args):
    """Choose every .dens radius, or the requested uniform coarse grid."""
    if getattr(args, "knots_from_dens", False):
        return np.asarray(rad, float).copy()
    return np.linspace(rad[0], rad[-1], args.knots)


def hat_matrix(radii, knots):
    """Nonnegative linear hats, summing to one throughout the solver domain."""
    radii, knots = np.asarray(radii, float), np.asarray(knots, float)
    if (knots.ndim != 1 or len(knots) < 2 or np.any(np.diff(knots) <= 0)
            or np.any(radii < knots[0]) or np.any(radii > knots[-1])):
        raise ValueError("invalid knots or radii outside knot domain")
    return np.column_stack([np.interp(radii, knots, np.eye(len(knots))[j])
                            for j in range(len(knots))])


def kernel(knots, length_fraction, nugget=1e-8):
    """Unit-diagonal squared-exponential correlation plus numerical nugget."""
    if not 0 < length_fraction <= 2 or not 0 < nugget < 1e-2:
        raise ValueError("invalid kernel hyperparameters")
    x = (knots - knots[0]) / (knots[-1] - knots[0])
    K = np.exp(-0.5 * ((x[:, None] - x[None, :]) / length_fraction) ** 2)
    K += nugget * np.eye(len(knots))
    return K


def ordering_system(D, c, g_obs):
    """Adjacent central g_obs order constraints A xi >= b, dimensionless."""
    order = np.argsort(g_obs, kind="mergesort")
    scale = max(float(np.max(np.abs(c + D @ np.full(D.shape[1], 0.5)))), 1e-30)
    A = np.diff(D[order], axis=0) / scale
    b = -np.diff(c[order]) / scale
    return A, b, order, scale


def solve_qp(D, c, g_obs, knots, lower, upper, length_fraction,
             feasibility_tol=1e-8, active_tol=1e-7):
    """Solve a convex QP, with an LP feasibility check and audited residuals."""
    if not (0 < lower < 0.5 < upper):
        raise ValueError("bounds must contain the reference M/L=0.5")
    A, b, order, scale = ordering_system(D, c, g_obs)
    n = len(knots)
    feasible = linprog(np.zeros(n), A_ub=-A, b_ub=-b,
                       bounds=[(lower, upper)] * n, method="highs")
    if feasible.status == 2:
        return {"status": "infeasible", "lp_message": feasible.message}
    if not feasible.success:
        return {"status": "feasibility_unknown", "lp_message": feasible.message}
    invK = np.linalg.inv(kernel(knots, length_fraction))
    metric_scale = float(np.max(np.abs(invK)))
    scaled_invK = invK / metric_scale
    m = np.full(n, 0.5)
    def objective(x):
        delta = x - m
        return float(delta @ scaled_invK @ delta)
    def gradient(x):
        return 2 * scaled_invK @ (x - m)
    result = minimize(objective, feasible.x, jac=gradient, method="SLSQP",
                      bounds=Bounds(lower, upper),
                      constraints=[LinearConstraint(A, b, np.inf)],
                      options={"ftol": 1e-11, "maxiter": 3000})
    x = np.asarray(result.x)
    residual = A @ x - b
    max_violation = float(max(0.0, -np.min(residual)))
    if not result.success or max_violation > feasibility_tol:
        return {"status": "qp_failure", "message": result.message,
                "max_scaled_violation": max_violation}
    g_bar = c + D @ x
    gaps = np.diff(g_bar[order])
    return {"status": "feasible_weak_order", "xi": x, "g_bar": g_bar,
            "order": order, "gaps": gaps, "objective": objective(x) * metric_scale,
            "max_scaled_violation": max_violation,
            "active_order_pairs": [list(map(int, order[j:j + 2]))
                                   for j in np.flatnonzero(residual <= active_tol)],
            "active_lower_knots": np.flatnonzero(x - lower <= active_tol).tolist(),
            "active_upper_knots": np.flatnonzero(upper - x <= active_tol).tolist(),
            "scale_m_s2": scale, "qp_message": result.message}


def cache_gate(galaxy, root):
    """Require the completed all-175 diagnostic for this exact driver version."""
    files = sorted(Path(root).glob(VALIDATION_GLOB))
    if len(files) != 5:
        raise ValueError("five completed all-175 cache-validation shards required")
    code_hash = sha256(mono.__file__)
    records, qc = [], []
    for path in files:
        data = json.loads(path.read_text())
        if data.get("status") != "diagnostic_only" or data.get("code_sha256") != code_hash:
            raise ValueError(f"stale or incomplete cache validation: {path}")
        records += [x for x in data["records"] if x["galaxy"] == galaxy]
        qc += [x for x in data["galaxies"] if x["galaxy"] == galaxy]
    if (len(qc) != 1 or not records or qc[0]["basis_warning_count"]
            or any(not x["force_pass"] or not x["rank_agrees"]
                   or not x["prob_valid_agrees"] for x in records)):
        raise ValueError(f"cache validation failed or missing for {galaxy}")
    return {"files": [str(p) for p in files], "records": len(records),
            "zero_start": qc[0]["photometric_zero_start"],
            "retained_kinematic_points": qc[0]["kinematic_retained"]}


def galaxy_data(galaxy, data_root, dens_root):
    table = load_sparc_table()
    row = table.loc[table.Galaxy == galaxy]
    if len(row) != 1:
        raise ValueError(f"expected one catalogue row for {galaxy}")
    raw = np.loadtxt(Path(data_root) / f"{galaxy}_rotmod.dat", ndmin=2)
    if raw.shape[1] < 6:
        raise ValueError("rotmod needs at least six columns")
    rad, brightness = get_SBdisk(galaxy, dens_root=dens_root, table=table)
    if (len(rad) < 3 or not np.all(np.isfinite(rad)) or
            np.any(np.diff(rad) <= 0) or np.any(rad < 0) or
            not np.all(np.isfinite(brightness))):
        raise ValueError("invalid full photometric grid")
    mask = np.isfinite(raw[:, 0]) & (raw[:, 0] >= rad[0]) & (raw[:, 0] <= rad[-1])
    r, vobs, errv, vgas, vbul = [raw[mask, j] for j in (0, 1, 2, 4, 5)]
    if (len(r) < 2 or np.any(r <= 0) or np.any(errv <= 0)
            or not np.all(np.isfinite([r, vobs, errv, vgas, vbul]))):
        raise ValueError("invalid retained kinematic data")
    zd = 0.196 * float(row.iloc[0]["Rdisk"]) ** 0.633
    return rad, brightness, r, vobs, errv, vgas, vbul, zd, len(raw) - len(r)


def build_response(rad, brightness, r, zd, knots):
    H = hat_matrix(rad, knots)
    columns, warnings = [], []
    for j in range(len(knots)):
        force, seen, error = mono.solve_disk_basis(rad, brightness * H[:, j], zd, r)
        if error or force is None or not np.all(np.isfinite(force)) or seen:
            raise ValueError(f"hat basis {j} failed: {error}; warnings={seen}")
        columns.append(vel2acc(force, r))
        warnings.extend(seen)
    return H, np.column_stack(columns), warnings


def direct_check(rad, brightness, r, zd, H, D, xi, tolerance=1e-8):
    force, seen, error = mono.solve_disk_basis(rad, brightness * (H @ xi), zd, r)
    if error or force is None or seen:
        raise ValueError(f"direct force solve failed: {error}; warnings={seen}")
    direct = vel2acc(force, r)
    cached = D @ xi
    scale = max(float(np.max(np.abs(direct))), 1e-30)
    delta = float(np.max(np.abs(direct - cached)))
    if not np.all(np.isfinite(direct)) or delta > tolerance * scale:
        raise ValueError(f"direct/cache force mismatch: {delta / scale}")
    return {"max_abs_accel_delta_m_s2": delta, "relative_to_peak": delta / scale,
            "tolerance": tolerance}


def solve_with_retry(D, c, observed, knots, args):
    """Retry only a proven LP-infeasible fit at a wider diagnostic bound."""
    upper = float(args.max_ml)
    fit = solve_qp(D, c, observed, knots, args.min_ml, upper,
                   args.length_fraction)
    attempts = [{"upper_ml": upper, "status": fit["status"],
                 "message": fit.get("lp_message", fit.get("message", ""))}]
    retry_upper = float(getattr(args, "retry_max_ml", 20.0))
    if fit["status"] == "infeasible" and retry_upper > upper:
        upper = retry_upper
        fit = solve_qp(D, c, observed, knots, args.min_ml, upper,
                       args.length_fraction)
        attempts.append({"upper_ml": upper, "status": fit["status"],
                         "message": fit.get("lp_message", fit.get("message", ""))})
    return fit, attempts, upper


def analyze(args):
    gate = cache_gate(args.galaxy, args.validation_root)
    rad, brightness, r, vobs, errv, vgas, vbul, zd, excluded = galaxy_data(
        args.galaxy, args.data_root, args.dens_root)
    knots = profile_knots(rad, args)
    H, D, _ = build_response(rad, brightness, r, zd, knots)
    prior = np.full(len(knots), 0.5)
    prior_qc = direct_check(rad, brightness, r, zd, H, D, prior)
    c = vel2acc(vgas * np.abs(vgas) * args.gas_ml +
                vbul ** 2 * args.bulge_ml, r)
    observed = vel2acc(vobs ** 2, r)
    fit, attempts, effective_upper = solve_with_retry(D, c, observed, knots, args)
    output = {"status": fit["status"], "galaxy": args.galaxy,
              "diagnostic_only": True, "reference_ml": 0.5,
              "gas_ml_qp": args.gas_ml, "bulge_ml_qp": args.bulge_ml,
              "knots_kpc": knots.tolist(),
              "knot_source": "dens_radii" if getattr(args, "knots_from_dens", False) else "uniform",
              "photometric_radii_kpc": rad.tolist(),
              "retained_radii_kpc": r.tolist(), "excluded_kinematic_rows": excluded,
              "kernel": {"kind": "squared_exponential_correlation",
                         "length_fraction": args.length_fraction, "nugget": 1e-8,
                         "amplitude": "cancels from hard-constrained optimum"},
              "requested_bounds_ml": [args.min_ml, args.max_ml],
              "bounds_ml": [args.min_ml, effective_upper],
              "fit_attempts": attempts, "cache_gate": gate,
              "numerical_tolerances": {"scaled_order_feasibility": 1e-8,
                                       "scaled_active_constraint": 1e-7,
                                       "direct_force_relative_to_peak": 1e-8},
              "input_sha256": {"rotmod": sha256(Path(args.data_root) / f"{args.galaxy}_rotmod.dat"),
                               "dens": sha256(Path(args.dens_root) / f"{args.galaxy}.dens"),
                               "table": sha256(DEFAULT_TABLE_PATH)},
              "prior_direct_cache_qc": prior_qc,
              "provenance": {"created_utc": datetime.now(timezone.utc).isoformat(),
                             "script_sha256": sha256(__file__),
                             "driver_sha256": sha256(mono.__file__),
                             "python": sys.version,
                             "args": {k: str(v) for k, v in vars(args).items()}}}
    if fit["status"] != "feasible_weak_order":
        output["optimization"] = fit
        return output
    xi = fit.pop("xi")
    g_bar = fit.pop("g_bar")
    order = fit.pop("order")
    gaps = fit.pop("gaps")
    output["optimized_direct_cache_qc"] = direct_check(
        rad, brightness, r, zd, H, D, xi)
    profile = H @ xi
    if np.min(profile) < args.min_ml - 1e-8 or np.max(profile) > effective_upper + 1e-8:
        raise ValueError("physical bounds violated on solver grid")
    probability_args = SimpleNamespace(
        gas_ml=args.gas_ml, bulge_ml=args.bulge_ml,
        marginalise_gas_bulge=True, prob_grid_size=args.prob_grid_size,
        zpadding=8.0, tie_atol=0.0, max_tie_permutations=720)
    gas_nodes, gas_weights = lognormal_quadrature(args.gas_ml, args.gas_ml_dex,
                                                   args.ml_samples)
    bulge_nodes, bulge_weights = lognormal_quadrature(args.bulge_ml,
                                                      args.bulge_ml_dex, args.ml_samples)
    conditional, marginal = mono.evaluate_disk_term(
        D @ xi / (1e6 / (r * 3.086e19)), vgas, vbul, r, vobs, errv,
        probability_args, gas_nodes, gas_weights, bulge_nodes, bulge_weights)
    output["optimization"] = fit
    output["optimization"].update(
        xi=xi.tolist(), profile_on_photometric_grid=profile.tolist(),
        max_abs_deviation_ml=float(np.max(np.abs(xi - 0.5))),
        mahalanobis_distance=float(np.sqrt(fit["objective"])),
        gbar_m_s2=g_bar.tolist(), observed_order=order.tolist(),
        adjacent_gbar_gaps_m_s2=gaps.tolist(),
        boundary_ties_present=bool(fit["active_order_pairs"]))
    output["post_optimization_probability"] = {
        "conditional_logp": conditional[0], "conditional_valid": conditional[1],
        "conditional_reason": conditional[2], "marginal_logp": marginal[0],
        "marginal_valid": marginal[1], "marginal_reason": marginal[2],
        "gas_ml_dex": args.gas_ml_dex, "bulge_ml_dex": args.bulge_ml_dex,
        "ml_samples_per_term": args.ml_samples, "prob_grid_size": args.prob_grid_size,
        "no_reoptimization": True}
    if not conditional[1] or not marginal[1]:
        output["status"] = "probability_invalid"
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--galaxy", required=True)
    parser.add_argument("--output", type=Path,
                        help="JSON path; defaults to outputs/inverse_rar_diagnostic/<galaxy>_inverse.json")
    parser.add_argument("--data-root", type=Path,
                        default=mono.DEFAULT_DATA_ROOT)
    parser.add_argument("--dens-root", type=Path, default=Path(DEFAULT_DENS_ROOT))
    parser.add_argument("--validation-root", type=Path,
                        default=Path(DEFAULT_PROJECT_ROOT) / "outputs")
    parser.add_argument("--knots", type=int, default=8)
    knot_group = parser.add_mutually_exclusive_group()
    knot_group.add_argument("--knots-from-dens", dest="knots_from_dens",
                            action="store_true", default=True,
                            help="Use every .dens radius as a knot (default)")
    knot_group.add_argument("--uniform-knots", dest="knots_from_dens",
                            action="store_false", help="Use --knots uniform knots")
    parser.add_argument("--length-fraction", type=float, default=0.3)
    parser.add_argument("--min-ml", type=float, default=0.01)
    parser.add_argument("--max-ml", type=float, default=10.0,
                        help="Exploratory first-attempt upper bound")
    parser.add_argument("--retry-max-ml", type=float, default=20.0,
                        help="Retry LP-infeasible fits at this upper bound")
    parser.add_argument("--gas-ml", type=float, default=1.0)
    parser.add_argument("--bulge-ml", type=float, default=0.7)
    parser.add_argument("--gas-ml-dex", type=float, default=0.04)
    parser.add_argument("--bulge-ml-dex", type=float, default=0.1)
    parser.add_argument("--ml-samples", type=int, default=9)
    parser.add_argument("--prob-grid-size", type=int, default=400)
    parser.add_argument("--replace-existing", action="store_true",
                        help="Intentionally regenerate this diagnostic JSON/PNG pair")
    args = parser.parse_args()
    if args.knots < 2 or args.ml_samples < 1 or args.prob_grid_size < 2:
        parser.error("knots, ml samples and probability grid must be valid")
    if args.output is None:
        args.output = (Path(DEFAULT_PROJECT_ROOT) / "outputs" /
                       "inverse_rar_diagnostic" / f"{args.galaxy}_inverse.json")
    if not (0 < args.min_ml < 0.5 < args.max_ml) or args.retry_max_ml <= 0.5:
        parser.error("M/L bounds must contain the reference 0.5")
    plot_path = args.output.with_suffix(".png")
    if (args.output.exists() or plot_path.exists()) and not args.replace_existing:
        raise FileExistsError(f"result or plot already exists: {args.output}, {plot_path}")
    result = analyze(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    import plot_inverse_rar_comparison as plotter
    result["provenance"]["plot_script_sha256"] = sha256(plotter.__file__)
    try:
        result["plot"] = plotter.render_plot(
            result, plot_path, args.data_root, args.dens_root,
            sys.modules[__name__], args.replace_existing)
    except Exception as exc:
        result["plot"] = {"status": "failed", "path": str(plot_path),
                          "reason": f"{type(exc).__name__}: {exc}"}
        args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
        raise
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(result["status"], args.output, plot_path, flush=True)


if __name__ == "__main__":
    main()


