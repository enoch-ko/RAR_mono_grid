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
import time
from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.linalg import solve_triangular
from scipy.optimize import LinearConstraint, linprog, minimize, nnls

import mono_SPARC_RAR as mono
from sparc_rotmod_columns import (ROTMOD_PARSER_ID, read_rotmod, read_column, parser_provenance)
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


def json_safe(value):
    """Failure candidates are diagnostics, with non-finite values encoded as null."""
    if isinstance(value, dict):
        return {k: json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [json_safe(v) for v in value]
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.bool_):
        return bool(value)
    return value


def qp_options(args):
    """Keep the historical backend unless OSQP is explicitly selected."""
    return {"qp_backend": getattr(args, "qp_backend", "slsqp"),
            "osqp_eps": getattr(args, "osqp_eps", 1e-9),
            "osqp_max_iter": getattr(args, "osqp_max_iter", 100000),
            "osqp_time_limit": getattr(args, "osqp_time_limit", 10.0)}


def osqp_solve(C, rhs, z0, eps, max_iter, time_limit):
    """min z.T z, C z >= rhs; no change to physical model or constraints."""
    import osqp
    from scipy import sparse
    settings = {"eps_abs": eps, "eps_rel": eps, "max_iter": max_iter,
                "time_limit": time_limit, "polishing": True,
                "adaptive_rho_interval": 50, "verbose": False,
                "scaling": 0, "sigma": 1e-9,
                "eps_prim_inf": 1e-12, "eps_dual_inf": 1e-12}
    solver = osqp.OSQP()
    solver.setup(P=2.0 * sparse.eye(C.shape[1], format="csc"),
                 q=np.zeros(C.shape[1]), A=sparse.csc_matrix(C),
                 l=rhs, u=np.full(len(rhs), np.inf), **settings)
    solver.warm_start(x=z0)
    result = solver.solve(raise_error=False)
    info = result.info
    diagnostics = {"version": osqp.__version__, "settings": settings,
                   "status": info.status, "status_code": int(info.status_val),
                   "iterations": int(info.iter), "polish_status": int(info.status_polish),
                   "primal_residual": float(info.prim_res),
                   "dual_residual": float(info.dual_res),
                   "run_time_seconds": float(info.run_time),
                   "rho_estimate": float(info.rho_estimate)}
    if info.status_val in (3, 4, 5, 6):
        # OSQP returns sentinel x/y values for infeasibility certificates,
        # not an iterate. LP success means this is still a NUMERICAL failure.
        diagnostics["candidate_available"] = False
        diagnostics["primal_infeasibility_certificate"] = json_safe(result.prim_inf_cert)
        diagnostics["dual_infeasibility_certificate"] = json_safe(result.dual_inf_cert)
        return None, None, False, diagnostics
    # For a lower-bound-only row OSQP y is negative: lambda=-y >=0.
    return result.x, result.y, info.status_val == 1, diagnostics


def solve_qp(D, c, g_obs, knots, lower, upper, length_fraction,
             feasibility_tol=1e-8, active_tol=1e-7, optimality_tol=1e-6,
             gbar_floor_relative=1e-8, qp_backend="slsqp", osqp_eps=1e-9,
             osqp_max_iter=100000, osqp_time_limit=10.0):
    """Same whitened QP, selectable solver, independent physical/dual checks."""
    if qp_backend not in ("slsqp", "osqp"):
        raise ValueError("unknown QP backend")
    if not (0 < lower < 0.5 < upper):
        raise ValueError("bounds must contain the reference M/L=0.5")
    if optimality_tol <= 0 or not np.isfinite(optimality_tol):
        raise ValueError("positive finite optimality tolerance required")
    if (not np.isfinite(osqp_eps) or osqp_eps <= 0 or osqp_max_iter < 1
            or not np.isfinite(osqp_time_limit) or osqp_time_limit <= 0):
        raise ValueError("positive finite OSQP controls required")
    A, b, order, scale = ordering_system(D, c, g_obs)
    if not np.isfinite(gbar_floor_relative) or gbar_floor_relative <= 0:
        raise ValueError("positive finite relative gbar floor required")
    # Strictly positive signed net gbar: no absolute values or force clipping.
    order_rows = len(b)
    A = np.vstack([A, D / scale])
    b = np.concatenate([b, gbar_floor_relative - c / scale])
    n = len(knots)
    m = np.full(n, 0.5)
    K = kernel(knots, length_fraction)
    L = np.linalg.cholesky(K)
    diagnostics = {"backend_requested": qp_backend, "knots": n,
                   "order_rows": order_rows, "positivity_rows": len(c),
                   "kernel_condition_number": float(np.linalg.cond(K)),
                   "feasibility_tolerance": feasibility_tol,
                   "active_tolerance": active_tol, "optimality_tolerance": optimality_tol}

    def failed(stage, message, status="qp_failure", **extra):
        # No top-level xi on a failed result: plotter must not treat it as a fit.
        diagnostics["failure_stage"] = stage
        return json_safe(dict(status=status, message=str(message),
                              qp_backend=qp_backend, solver_diagnostics=diagnostics, **extra))

    reference_residual = A @ m - b
    reference_feasible = bool(np.all(reference_residual >= 0))
    if reference_feasible:
        x, z = m.copy(), np.zeros(n)
        message = "Reference M/L=0.5 is feasible; exact zero-penalty optimum"
        optimality = {"method": "exact_reference", "gap": 0.0,
                      "tolerance": optimality_tol, "passed": True}
        diagnostics.update(backend_executed="exact_reference", objective=0.0)
        iterations = 0
    else:
        feasible = linprog(np.zeros(n), A_ub=-A, b_ub=-b,
                           bounds=[(lower, upper)] * n, method="highs")
        diagnostics["lp"] = {"status_code": int(feasible.status),
                             "success": bool(feasible.success), "message": feasible.message}
        if feasible.status == 2:
            return failed("lp", feasible.message, "infeasible", lp_message=feasible.message)
        if not feasible.success:
            return failed("lp", feasible.message, "feasibility_unknown", lp_message=feasible.message)
        C = np.vstack([A @ L, L, -L])
        rhs = np.concatenate([b - A @ m, np.full(n, lower) - m,
                              m - np.full(n, upper)])
        row_scale = np.maximum(np.maximum(np.linalg.norm(C, axis=1), np.abs(rhs)), 1e-12)
        C, rhs = C / row_scale[:, None], rhs / row_scale
        z0 = solve_triangular(L, feasible.x - m, lower=True)
        diagnostics.update(backend_executed=qp_backend,
                           row_scale_range=[float(row_scale.min()), float(row_scale.max())],
                           initial_objective=float(z0 @ z0))
        started = time.monotonic()
        dual = None
        try:
            if qp_backend == "osqp":
                z, dual, success, details = osqp_solve(
                    C, rhs, z0, osqp_eps, osqp_max_iter, osqp_time_limit)
                diagnostics["solver"] = details
                message, iterations = details["status"], details["iterations"]
            else:
                result = minimize(lambda z: float(z @ z), z0,
                                  jac=lambda z: 2 * z, method="SLSQP",
                                  constraints=[LinearConstraint(C, rhs, np.inf)],
                                  options={"ftol": 1e-11, "maxiter": 3000})
                z, success = result.x, bool(result.success)
                message, iterations = str(result.message), int(result.nit)
                diagnostics["solver"] = {"status": message,
                    "status_code": getattr(result, "status", None), "iterations": iterations,
                    "settings": {"ftol": 1e-11, "maxiter": 3000}}
        except (ImportError, ValueError, RuntimeError) as exc:
            diagnostics["elapsed_seconds"] = time.monotonic() - started
            return failed("solver_exception", str(exc))
        diagnostics["elapsed_seconds"] = time.monotonic() - started
        if z is None:
            return failed("solver", "solver returned no candidate: " + message)
        z = np.asarray(z, float)
        with np.errstate(invalid="ignore", over="ignore"):
            x = m + L @ z
        diagnostics.update(candidate_z=z, candidate_xi=x,
                           candidate_finite=bool(np.all(np.isfinite(x))))
        if not np.all(np.isfinite(x)):
            return failed("candidate_nonfinite", message)
        residual = A @ x - b
        slack = C @ z - rhs
        diagnostics.update(objective=float(z @ z),
            max_order_violation=float(max(0.0, -np.min(residual[:order_rows]))),
            max_positivity_violation=float(max(0.0, -np.min(residual[order_rows:]))),
            max_bound_violation=float(max(0.0, lower - x.min(), x.max() - upper)),
            max_whitened_violation=float(max(0.0, -slack.min())),
            candidate_minimum_gbar_m_s2=float(np.min(c + D @ x)))
        # Independently recompute a valid dual lower bound. OSQP multipliers
        # are checked in OUR matrices, not accepted solely from solver status.
        multipliers = np.zeros(len(rhs))
        if dual is not None and np.all(np.isfinite(dual)):
            multipliers = np.maximum(-np.asarray(dual), 0.0)
            multiplier_source = "osqp_lower_bound_dual"
        else:
            active = slack <= active_tol
            multiplier_source = "active_set_nnls"
            try:
                if np.any(active):
                    multipliers[active], _ = nnls(C[active].T, 2 * z,
                                                  maxiter=30 * len(rhs))
            except RuntimeError as exc:
                return failed("dual_audit", str(exc))
        stationarity = 2 * z - C.T @ multipliers
        gap = float(stationarity @ stationarity / 4 + multipliers @ slack)
        tolerance = optimality_tol * max(1.0, float(z @ z))
        optimality = {"method": "nonnegative_multiplier_duality_gap", "gap": gap,
                      "multiplier_source": multiplier_source,
                      "dual_lower_bound": float(z @ z - gap),
                      "stationarity_norm": float(np.linalg.norm(stationarity)),
                      "tolerance": tolerance,
                      "passed": bool(np.isfinite(gap) and abs(gap) <= tolerance)}
        diagnostics["optimality"] = optimality
        if not success:
            return failed("solver", message, optimality=optimality)
        if not optimality["passed"]:
            return failed("optimality", "optimality check failed", optimality=optimality)
    residual = A @ x - b
    max_violation = float(max(0.0, -np.min(residual)))
    bound_violation = float(max(0.0, lower - np.min(x), np.max(x) - upper))
    if max_violation > feasibility_tol or bound_violation > feasibility_tol:
        return failed("physical_feasibility", "physical/order feasibility failed",
                      max_scaled_violation=max_violation, max_bound_violation=bound_violation)
    whitened = solve_triangular(L, x - m, lower=True)
    g_bar = c + D @ x
    if not np.all(np.isfinite(g_bar)) or np.any(g_bar <= 0):
        return failed("positive_gbar", "strict positive gbar check failed")
    # Successful xi is already saved below; avoid duplicating profile arrays.
    diagnostics.pop("candidate_xi", None)
    diagnostics.pop("candidate_z", None)
    return {"status": "feasible_weak_order", "xi": x, "g_bar": g_bar,
            "order": order, "gaps": np.diff(g_bar[order]),
            "objective": float(whitened @ whitened),
            "optimizer_method": "cholesky_whitened_v1", "qp_backend": qp_backend,
            "solver_diagnostics": json_safe(diagnostics),
            "constraint_version": "positive_net_gbar_v1",
            "gbar_floor_relative": gbar_floor_relative,
            "gbar_floor_m_s2": gbar_floor_relative * scale,
            "minimum_gbar_m_s2": float(np.min(g_bar)),
            "active_positivity_radii": np.flatnonzero(residual[order_rows:] <= active_tol).tolist(),
            "reference_feasible": reference_feasible,
            "iterations": iterations, "optimality": optimality,
            "max_scaled_violation": max_violation, "max_bound_violation": bound_violation,
            "active_order_pairs": [list(map(int, order[j:j + 2]))
                                   for j in np.flatnonzero(residual[:order_rows] <= active_tol)],
            "active_lower_knots": np.flatnonzero(x - lower <= active_tol).tolist(),
            "active_upper_knots": np.flatnonzero(upper - x <= active_tol).tolist(),
            "scale_m_s2": scale, "qp_message": message}


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
    raw = read_rotmod(Path(data_root) / f"{galaxy}_rotmod.dat")
    rad, brightness = get_SBdisk(galaxy, dens_root=dens_root, table=table)
    if (len(rad) < 3 or not np.all(np.isfinite(rad)) or
            np.any(np.diff(rad) <= 0) or np.any(rad < 0) or
            not np.all(np.isfinite(brightness))):
        raise ValueError("invalid full photometric grid")
    mask = np.isfinite(raw[:, 0]) & (raw[:, 0] >= rad[0]) & (raw[:, 0] <= rad[-1])
    # Select by the verified header: gas is column3, NOT the disc template in4.
    r, vobs, errv, vgas, vbul = [read_column(raw, name, mask)
        for name in ("Rad", "Vobs", "errV", "Vgas", "Vbul")]
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
                   args.length_fraction,
                   gbar_floor_relative=getattr(args, "gbar_floor_relative", 1e-8),
                   **qp_options(args))
    attempts = [{"upper_ml": upper, "status": fit["status"],
                 "message": fit.get("lp_message", fit.get("message", ""))}]
    retry_upper = float(getattr(args, "retry_max_ml", 20.0))
    if fit["status"] == "infeasible" and retry_upper > upper:
        upper = retry_upper
        fit = solve_qp(D, c, observed, knots, args.min_ml, upper,
                       args.length_fraction,
                   gbar_floor_relative=getattr(args, "gbar_floor_relative", 1e-8),
                   **qp_options(args))
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
              "positivity_constraint": {"kind": "signed_net_gbar",
                                        "relative_numerical_floor": getattr(args, "gbar_floor_relative", 1e-8),
                                        "strict_positive_postcheck": True},
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
                             "script_file": Path(__file__).name,
                             "input_parser": parser_provenance(),
                             "driver_sha256": sha256(mono.__file__),
                             "utils_sha256": sha256(Path(DEFAULT_PROJECT_ROOT) / "utils.py"),
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
        zpadding=getattr(args, "zpadding", 8.0), tie_atol=0.0, max_tie_permutations=720)
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
        "zpadding": getattr(args, "zpadding", 8.0), "no_reoptimization": True}
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
    parser.add_argument("--gbar-floor-relative", type=float, default=1e-8,
                        help="Numerical positivity margin relative to reference peak acceleration")
    parser.add_argument("--qp-backend", choices=("slsqp", "osqp"), default="slsqp")
    parser.add_argument("--osqp-eps", type=float, default=1e-9)
    parser.add_argument("--osqp-max-iter", type=int, default=100000)
    parser.add_argument("--osqp-time-limit", type=float, default=10.0)
    parser.add_argument("--gas-ml", type=float, default=1.0)
    parser.add_argument("--bulge-ml", type=float, default=0.7)
    parser.add_argument("--gas-ml-dex", type=float, default=0.04)
    parser.add_argument("--bulge-ml-dex", type=float, default=0.1)
    parser.add_argument("--ml-samples", type=int, default=9)
    parser.add_argument("--prob-grid-size", type=int, default=400)
    parser.add_argument("--zpadding", type=float, default=8.0,
                        help="Gaussian tail padding; vary independently for convergence")
    parser.add_argument("--replace-existing", action="store_true",
                        help="Intentionally regenerate this diagnostic JSON/PNG pair")
    args = parser.parse_args()
    if not np.isfinite(args.zpadding) or args.zpadding <= 0:
        parser.error("positive finite zpadding required")
    if args.knots < 2 or args.ml_samples < 1 or args.prob_grid_size < 2:
        parser.error("knots, ml samples and probability grid must be valid")
    if args.output is None:
        args.output = (Path(DEFAULT_PROJECT_ROOT) / "outputs" /
                       "rotmod_columns_corrected_20261015" / "qp" / f"{args.galaxy}_inverse.json")
    if not (0 < args.min_ml < 0.5 < args.max_ml) or args.retry_max_ml <= 0.5:
        parser.error("M/L bounds must contain the reference 0.5")
    plot_path = args.output.with_suffix(".png")
    if (args.output.exists() or plot_path.exists()) and not args.replace_existing:
        raise FileExistsError(f"result or plot already exists: {args.output}, {plot_path}")
    result = analyze(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    import plot_inverse_rar_comparison_corrected as plotter
    result["provenance"]["plot_script_sha256"] = sha256(plotter.__file__)
    try:
        result["plot"] = plotter.render_plot(
            result, plot_path, args.data_root, args.dens_root,
            sys.modules[__name__], args.replace_existing)
    except Exception as exc:
        result["plot"] = {"status": "failed", "path": str(plot_path),
                          "reason": f"{type(exc).__name__}: {exc}"}
        args.output.write_text(json.dumps(json_safe(result), indent=2, allow_nan=False) + "\n")
        raise
    args.output.write_text(json.dumps(json_safe(result), indent=2, allow_nan=False) + "\n")
    print(result["status"], args.output, plot_path, flush=True)


if __name__ == "__main__":
    main()
