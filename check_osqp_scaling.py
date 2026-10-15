#!/usr/bin/env python3
"""Bounded OSQP objective-scaling probes on a saved immutable disc response.
Multiplying the whole objective by a positive scalar leaves the optimum unchanged.
This is numerical investigation only; never accepts an unaudited solver status.
"""
import argparse
import json
import signal
from pathlib import Path
import numpy as np
import osqp
from scipy import sparse
import inverse_rar_diagnostic as inv
from check_osqp_fits import expired


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparison", required=True, type=Path)
    parser.add_argument("--seconds", type=int, default=35)
    args = parser.parse_args()
    old = json.loads(args.comparison.read_text())
    response = Path(old["response"]["path"])
    if inv.sha256(response) != old["response"]["sha256"]:
        raise ValueError("saved response hash mismatch")
    output = args.comparison.with_name(old["galaxy"] + "_scaling.json")
    if output.exists():
        raise FileExistsError(output)
    result = {"galaxy": old["galaxy"], "diagnostic_only": True,
              "script_sha256": inv.sha256(__file__), "osqp_version": osqp.__version__,
              "input_comparison_sha256": inv.sha256(args.comparison),
              "response_sha256": inv.sha256(response), "cases": []}
    signal.signal(signal.SIGALRM, expired)
    signal.alarm(args.seconds)
    try:
        data = np.load(response)
        D, c, observed, knots = [data[k] for k in ("D", "c", "observed", "knots")]
        A, b, order, scale = inv.ordering_system(D, c, observed)
        A = np.vstack([A, D / scale])
        b = np.r_[b, 1e-8 - c / scale]
        n, m = len(knots), np.full(len(knots), .5)
        L = np.linalg.cholesky(inv.kernel(knots, .3))
        C = np.vstack([A @ L, L, -L])
        rhs = np.r_[b - A @ m, .01 - m, m - 10.]
        row_scale = np.maximum(np.maximum(np.linalg.norm(C, axis=1), np.abs(rhs)), 1e-12)
        C, rhs = C / row_scale[:, None], rhs / row_scale
        feasible = inv.linprog(np.zeros(n), A_ub=-A, b_ub=-b,
                              bounds=[(.01, 10.)] * n, method="highs")
        if not feasible.success:
            raise ValueError("LP did not verify feasibility")
        z0 = inv.solve_triangular(L, feasible.x - m, lower=True)
        for factor in (1., 1e-4, 1e-6, 1e-8):
            settings = dict(eps_abs=1e-9, eps_rel=1e-9, eps_prim_inf=1e-12,
                            eps_dual_inf=1e-12, max_iter=100000, time_limit=5.,
                            polishing=True, scaling=0, sigma=1e-9,
                            adaptive_rho_interval=50, verbose=False)
            solver = osqp.OSQP()
            solver.setup(P=2. * factor * sparse.eye(n, format="csc"), q=np.zeros(n),
                         A=sparse.csc_matrix(C), l=rhs, u=np.full(len(rhs), np.inf), **settings)
            solver.warm_start(x=z0)
            solved = solver.solve(raise_error=False)
            rec = {"objective_factor": factor, "settings": settings,
                   "status": solved.info.status, "status_code": solved.info.status_val,
                   "iterations": solved.info.iter, "seconds": solved.info.run_time,
                   "primal_residual": solved.info.prim_res, "dual_residual": solved.info.dual_res,
                   "rho_estimate": solved.info.rho_estimate, "polish_status": solved.info.status_polish}
            if solved.info.status_val not in (3, 4, 5, 6) and solved.x is not None:
                z = solved.x
                x = m + L @ z
                multipliers = np.maximum(-solved.y / factor, 0.)
                slack = C @ z - rhs
                stationarity = 2. * z - C.T @ multipliers
                gap = float(stationarity @ stationarity / 4. + multipliers @ slack)
                J = float(z @ z)
                violation = float(max(0., -np.min(A @ x - b)))
                bound_violation = float(max(0., .01 - x.min(), x.max() - 10.))
                rec.update(objective=J, gap=gap, gap_tolerance=1e-6 * max(1., J),
                           violation=violation, bound_violation=bound_violation,
                           candidate_xi=x, candidate_z=z,
                           audited=bool(solved.info.status_val == 1 and violation <= 1e-8
                                        and bound_violation <= 1e-8 and abs(gap) <= 1e-6 * max(1., J)
                                        and np.all(c + D @ x > 0)))
            result["cases"].append(rec)
            print(old["galaxy"], factor, rec["status"], "J", rec.get("objective"),
                  "gap", rec.get("gap"), "violation", rec.get("violation"),
                  "audited", rec.get("audited", False), flush=True)
    except Exception as exc:
        result["error"] = type(exc).__name__ + ": " + str(exc)
    finally:
        signal.alarm(0)
    output.write_text(json.dumps(inv.json_safe(result), indent=2, allow_nan=False) + "\n")
    print(output, flush=True)


if __name__ == "__main__":
    main()
