#!/usr/bin/env python3
"""Light, dependency-free unittest checks for the first-order inverse diagnostic."""
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

import inverse_rar_diagnostic as inverse
import plot_inverse_rar_comparison as plotter


class InverseDiagnosticTests(unittest.TestCase):
    def test_hats_partition_and_physical_bounds(self):
        knots = np.array([0., 1., 2., 3.])
        radii = np.linspace(0., 3., 13)
        H = inverse.hat_matrix(radii, knots)
        np.testing.assert_allclose(H.sum(axis=1), 1.)
        self.assertGreaterEqual(H.min(), 0.)
        xi = np.array([0.1, 0.4, 1.2, 0.8])
        self.assertGreaterEqual((H @ xi).min(), xi.min())
        self.assertLessEqual((H @ xi).max(), xi.max())

    def test_exact_dens_radii_as_knots(self):
        radii = np.array([0., 0.17, 0.43, 1.2, 3.0])
        knots = inverse.profile_knots(
            radii, SimpleNamespace(knots=4, knots_from_dens=True))
        np.testing.assert_array_equal(knots, radii)
        np.testing.assert_allclose(inverse.hat_matrix(radii, knots),
                                   np.eye(len(radii)), atol=1e-15)

    def test_reference_already_feasible(self):
        D = np.eye(2)
        fit = inverse.solve_qp(D, np.array([0., 0.2]),
                               np.array([1., 2.]), np.array([0., 1.]),
                               0.01, 2., 0.3)
        self.assertEqual(fit["status"], "feasible_weak_order")
        np.testing.assert_allclose(fit["xi"], 0.5, atol=1e-6)
        self.assertLess(fit["max_scaled_violation"], 1e-8)

    def test_boundary_tie_is_not_strict(self):
        D = np.eye(2)
        fit = inverse.solve_qp(D, np.array([0.4, 0.]),
                               np.array([1., 2.]), np.array([0., 1.]),
                               0.01, 2., 0.3)
        self.assertEqual(fit["status"], "feasible_weak_order")
        self.assertTrue(fit["active_order_pairs"])
        self.assertLess(abs(fit["gaps"][0]), 1e-7)

    def test_dense_reference_is_exact_without_optimizer(self):
        knots = np.linspace(0., 1., 73)
        D = np.outer(np.arange(1., 5.), np.ones(73) / 73)
        with patch.object(inverse, "minimize", side_effect=AssertionError("no QP needed")):
            fit = inverse.solve_qp(D, np.zeros(4), np.arange(1., 5.),
                                   knots, .01, 10., .3)
        np.testing.assert_array_equal(fit["xi"], np.full(73, .5))
        self.assertEqual(fit["objective"], 0.)
        self.assertEqual(fit["optimality"]["method"], "exact_reference")

    def test_whitened_solution_matches_analytic_projection(self):
        knots = np.linspace(0., 1., 19)
        D = np.zeros((2, 19))
        D[0, 0], D[1, -1] = 1., 1.
        c = np.array([.4, 0.])
        fit = inverse.solve_qp(D, c, np.array([1., 2.]), knots, .01, 10., .3)
        self.assertEqual(fit["status"], "feasible_weak_order")
        K = inverse.kernel(knots, .3)
        a = D[1] - D[0]
        expected = np.full(19, .5) + K @ a * (.4 / (a @ K @ a))
        np.testing.assert_allclose(fit["xi"], expected, atol=2e-6, rtol=0)
        self.assertAlmostEqual(fit["objective"], .4**2 / (a @ K @ a), places=7)
        self.assertTrue(fit["optimality"]["passed"])
        self.assertFalse(fit["reference_feasible"])

    def test_feasible_but_nonoptimal_solver_result_rejected(self):
        knots = np.array([0., 1.])
        L = np.linalg.cholesky(inverse.kernel(knots, .3))
        bad_x = np.array([.01, 1.])
        bad_z = inverse.solve_triangular(L, bad_x - .5, lower=True)
        fake = SimpleNamespace(success=True, x=bad_z, nit=1, message="fake success")
        with patch.object(inverse, "minimize", return_value=fake):
            fit = inverse.solve_qp(np.eye(2), np.array([.4, 0.]),
                                   np.array([1., 2.]), knots, .01, 2., .3)
        self.assertEqual(fit["status"], "qp_failure")
        self.assertFalse(fit["optimality"]["passed"])

    def test_whitened_bounds_are_enforced(self):
        fit = inverse.solve_qp(np.eye(2), np.array([.7, 0.]),
                               np.array([1., 2.]), np.array([0., 1.]),
                               .01, .8, .3)
        self.assertEqual(fit["status"], "feasible_weak_order")
        self.assertGreaterEqual(min(fit["xi"]), .01 - 1e-8)
        self.assertLessEqual(max(fit["xi"]), .8 + 1e-8)
        self.assertTrue(fit["active_upper_knots"])
        self.assertTrue(fit["optimality"]["passed"])

    def test_infeasible(self):
        fit = inverse.solve_qp(np.zeros((2, 2)), np.array([1., 0.]),
                               np.array([1., 2.]), np.array([0., 1.]),
                               0.01, 2., 0.3)
        self.assertEqual(fit["status"], "infeasible")

    def test_infeasible_retries_at_twenty(self):
        args = SimpleNamespace(min_ml=0.01, max_ml=10., retry_max_ml=20.,
                               length_fraction=0.3)
        with patch.object(inverse, "solve_qp", side_effect=[
                {"status": "infeasible", "lp_message": "first"},
                {"status": "infeasible", "lp_message": "second"}]) as solve:
            fit, attempts, upper = inverse.solve_with_retry(
                np.zeros((2, 2)), np.zeros(2), np.array([1., 2.]),
                np.array([0., 1.]), args)
        self.assertEqual(solve.call_count, 2)
        self.assertEqual([call.args[5] for call in solve.call_args_list], [10., 20.])
        self.assertEqual([a["upper_ml"] for a in attempts], [10., 20.])
        self.assertEqual(upper, 20.)
        self.assertEqual(fit["status"], "infeasible")

    def test_qp_failure_does_not_trigger_bound_retry(self):
        args = SimpleNamespace(min_ml=0.01, max_ml=10., retry_max_ml=20.,
                               length_fraction=0.3)
        with patch.object(inverse, "solve_qp", return_value={
                "status": "qp_failure", "message": "optimizer error"}) as solve:
            fit, attempts, upper = inverse.solve_with_retry(
                np.zeros((2, 2)), np.zeros(2), np.array([1., 2.]),
                np.array([0., 1.]), args)
        solve.assert_called_once()
        self.assertEqual(upper, 10.)
        self.assertEqual(len(attempts), 1)
        self.assertEqual(fit["status"], "qp_failure")

    def test_infeasible_plot_is_sparc_only(self):
        rad = np.array([0., 1., 2.])
        r = np.array([0.5, 1.5])
        galaxy = (rad, np.ones(3), r, np.array([30., 40.]),
                  np.ones(2), np.array([10., 12.]), np.zeros(2), 0.2, 0)
        result = {"galaxy": "synthetic", "status": "infeasible",
                  "reference_ml": 0.5, "gas_ml_qp": 1., "bulge_ml_qp": 0.7,
                  "bounds_ml": [0.01, 20.], "retained_radii_kpc": r.tolist(),
                  "optimization": {"status": "infeasible"}}
        with tempfile.TemporaryDirectory() as root:
            root = Path(root)
            raw = np.column_stack([r, galaxy[3], galaxy[4],
                                   np.array([15., 20.]), galaxy[5], galaxy[6]])
            np.savetxt(root / "synthetic_rotmod.dat", raw)
            with (patch.object(inverse, "galaxy_data", return_value=galaxy),
                  patch.object(inverse.mono, "solve_disk_basis",
                               side_effect=AssertionError("SPARC-only plot called Casertano"))):
                qc = plotter.render_plot(result, root / "synthetic.png",
                                         root, root, inverse)
            self.assertEqual(qc["kind"], "sparc_reference_only")
            self.assertTrue((root / "synthetic.png").is_file())

    def test_current_cache_gate(self):
        root = Path(inverse.DEFAULT_PROJECT_ROOT) / "outputs"
        gate = inverse.cache_gate("UGC07399", root)
        self.assertGreater(gate["records"], 0)

    def test_zero_start_direct_vs_cached(self):
        rad = np.linspace(0., 4., 17)
        brightness = 1e8 * np.exp(-rad)
        r = rad[2::2]
        knots = np.linspace(0., 4., 3)
        H, D, warnings = inverse.build_response(rad, brightness, r, 0.2, knots)
        self.assertFalse(warnings)
        qc = inverse.direct_check(rad, brightness, r, 0.2, H, D,
                                  np.array([0.3, 0.7, 0.5]))
        self.assertLess(qc["relative_to_peak"], 1e-8)

    def test_feasible_path_evaluates_probability_after_optimization(self):
        rad = np.array([0., 1., 2.])
        r = np.array([0.5, 1.5])
        H = inverse.hat_matrix(rad, np.array([0., 2.]))
        D = np.array([[1e-10, 0.], [0., 2e-10]])
        args = SimpleNamespace(
            galaxy="synthetic", validation_root=Path("."), data_root=Path("."),
            dens_root=Path("."), knots=2, min_ml=0.01, max_ml=2.,
            length_fraction=0.3, gas_ml=1., bulge_ml=0.7,
            gas_ml_dex=0.04, bulge_ml_dex=0.1, ml_samples=1,
            prob_grid_size=101)
        galaxy = (rad, np.ones(3), r, np.array([30., 80.]),
                  np.ones(2), np.zeros(2), np.zeros(2), 0.2, 0)
        with (patch.object(inverse, "cache_gate", return_value={"records": 1}),
              patch.object(inverse, "galaxy_data", return_value=galaxy),
              patch.object(inverse, "build_response", return_value=(H, D, [])),
              patch.object(inverse, "direct_check",
                           return_value={"relative_to_peak": 0.}) as qc,
              patch.object(inverse, "sha256", return_value="hash"),
              patch.object(inverse, "lognormal_quadrature",
                           return_value=(np.ones(1), np.ones(1))),
              patch.object(inverse.mono, "evaluate_disk_term",
                           return_value=((-0.2, True, None, []),
                                         (-0.3, True, None, []))) as probability):
            result = inverse.analyze(args)
        self.assertEqual(result["status"], "feasible_weak_order")
        np.testing.assert_allclose(result["optimization"]["xi"], 0.5, atol=1e-6)
        self.assertEqual(result["post_optimization_probability"]["marginal_logp"], -0.3)
        self.assertTrue(result["post_optimization_probability"]["no_reoptimization"])
        self.assertEqual(qc.call_count, 2)
        probability.assert_called_once()
        self.assertEqual(set(result["input_sha256"]), {"rotmod", "dens", "table"})



class OSQPDiagnosticTests(unittest.TestCase):
    def case(self, backend="osqp", **kwargs):
        return inverse.solve_qp(np.eye(2), np.array([.4, 0.]),
            np.array([1., 2.]), np.array([0., 1.]), .01, 2., .3,
            qp_backend=backend, **kwargs)

    def test_osqp_matches_analytic_projection(self):
        import importlib.util
        if importlib.util.find_spec("osqp") is None:
            self.skipTest("optional OSQP environment not installed")
        fit = self.case()
        self.assertEqual(fit["status"], "feasible_weak_order")
        np.testing.assert_allclose(fit["xi"], [.3, .7], atol=1e-8)
        self.assertTrue(fit["optimality"]["passed"])
        self.assertLessEqual(fit["max_scaled_violation"], 1e-8)
        self.assertTrue(fit["solver_diagnostics"]["solver"]["version"])
        np.testing.assert_allclose(fit["xi"], self.case("slsqp")["xi"], atol=1e-8)

    def test_osqp_retains_signed_component_and_positive_net_force(self):
        import importlib.util
        if importlib.util.find_spec("osqp") is None:
            self.skipTest("optional OSQP environment not installed")
        fit = inverse.solve_qp(np.diag([1., -1.]), np.array([0., .7]),
            np.array([1., 2.]), np.array([0., 1.]), .01, 2., .3, qp_backend="osqp")
        self.assertEqual(fit["status"], "feasible_weak_order")
        np.testing.assert_allclose(fit["xi"], [.35, .35], atol=1e-8)
        self.assertTrue(np.all(fit["g_bar"] > 0))
        self.assertTrue(fit["optimality"]["passed"])

    def test_osqp_reference_stays_exact_without_solver(self):
        with patch.object(inverse, "osqp_solve", side_effect=AssertionError("no solver needed")):
            fit = inverse.solve_qp(np.eye(2), np.array([0., .2]),
                np.array([1., 2.]), np.array([0., 1.]), .01, 2., .3, qp_backend="osqp")
        np.testing.assert_array_equal(fit["xi"], [.5, .5])
        self.assertEqual(fit["objective"], 0.)
        self.assertEqual(fit["solver_diagnostics"]["backend_executed"], "exact_reference")

    def test_invalid_backend_and_controls(self):
        for kw in ({"qp_backend": "unknown"}, {"osqp_eps": np.nan},
                   {"osqp_max_iter": 0}, {"osqp_time_limit": 0.}):
            with self.assertRaises(ValueError):
                inverse.solve_qp(np.eye(2), np.array([.4, 0.]), np.array([1., 2.]),
                    np.array([0., 1.]), .01, 2., .3, **kw)

    def test_failed_candidate_is_not_an_accepted_profile(self):
        fake = SimpleNamespace(success=False, x=np.array([0., 0.]), nit=2,
                               message="line search failed", status=8)
        with patch.object(inverse, "minimize", return_value=fake):
            fit = self.case("slsqp")
        self.assertEqual(fit["status"], "qp_failure")
        self.assertNotIn("xi", fit)
        detail = fit["solver_diagnostics"]
        self.assertEqual(detail["failure_stage"], "solver")
        self.assertEqual(detail["candidate_xi"], [.5, .5])
        self.assertGreater(detail["max_order_violation"], 0.)
        import json
        json.dumps(fit, allow_nan=False)

    def test_nonfinite_candidate_has_valid_failure_json(self):
        fake = SimpleNamespace(success=False, x=np.array([np.nan, np.inf]), nit=1,
                               message="bad candidate", status=8)
        with patch.object(inverse, "minimize", return_value=fake):
            fit = self.case("slsqp")
        self.assertEqual(fit["solver_diagnostics"]["failure_stage"], "candidate_nonfinite")
        self.assertNotIn("xi", fit)
        import json
        json.dumps(fit, allow_nan=False)

    def test_osqp_success_status_does_not_bypass_optimality(self):
        L = np.linalg.cholesky(inverse.kernel(np.array([0., 1.]), .3))
        z = inverse.solve_triangular(L, np.array([.01, 1.]) - .5, lower=True)
        details = {"status": "solved", "iterations": 1}
        with patch.object(inverse, "osqp_solve", return_value=(z, None, True, details)):
            fit = self.case()
        self.assertEqual(fit["status"], "qp_failure")
        self.assertEqual(fit["solver_diagnostics"]["failure_stage"], "optimality")
        self.assertFalse(fit["optimality"]["passed"])

    def test_missing_osqp_is_dependency_failure_not_infeasibility(self):
        with patch.object(inverse, "osqp_solve", side_effect=ImportError("OSQP unavailable")):
            fit = self.case()
        self.assertEqual(fit["status"], "qp_failure")
        self.assertEqual(fit["solver_diagnostics"]["failure_stage"], "solver_exception")
        self.assertTrue(fit["solver_diagnostics"]["lp"]["success"])

    def test_batch_forwarding_and_backend_provenance(self):
        import run_inverse_rar_batch as batch
        args = InverseBatchTests().args()
        args.qp_backend, args.osqp_eps = "osqp", 1e-9
        args.osqp_max_iter, args.osqp_time_limit = 100000, 10.
        command = batch.command_for("D512-2", Path("result.json"), args)
        self.assertEqual(command[command.index("--qp-backend") + 1], "osqp")
        self.assertEqual(command[command.index("--osqp-time-limit") + 1], "10.0")
        self.assertEqual(batch.run_settings(args)["qp_backend"], "osqp")


class InverseBatchTests(unittest.TestCase):
    def args(self, integration="boundary"):
        return SimpleNamespace(integration=integration, prob_grid_size=400,
                               ml_samples=9, zpadding=8., boundary_log_tolerance=2e-5,
                               boundary_max_refinements=8)

    def test_boundary_and_legacy_commands_share_fit_settings(self):
        import run_inverse_rar_batch as batch
        new = batch.command_for("D512-2", Path("result.json"), self.args())
        old = batch.command_for("D512-2", Path("result.json"), self.args("legacy"))
        self.assertEqual(Path(new[1]).name, "inverse_rar_boundary.py")
        self.assertEqual(Path(old[1]).name, "inverse_rar_diagnostic.py")
        self.assertEqual(new[2:new.index("--boundary-log-tolerance")], old[2:])
        self.assertIn("--zpadding", new)

    def test_new_artifacts_require_probability_hash_and_settings(self):
        import json
        import run_inverse_rar_batch as batch
        settings = batch.run_settings(self.args())
        hashes = batch.source_hashes()
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "synthetic.json"
            p.with_suffix(".png").touch()
            obj = {"galaxy": "synthetic", "provenance": {
                "args": {k: v for k, v in settings.items() if k not in
                         ("integration", "boundary_log_tolerance", "boundary_max_refinements")},
                "script_sha256": hashes["inverse_rar_diagnostic.py"],
                "driver_sha256": hashes["mono_SPARC_RAR.py"],
                "plot_script_sha256": hashes["plot_inverse_rar_comparison.py"],
                "utils_sha256": hashes["utils.py"],
                "boundary_entry_sha256": hashes["inverse_rar_boundary.py"],
                "probability_recursion_sha256": hashes["probability_recursion.py"]},
                "input_sha256": {"table": hashes["SPARC_table.txt"]},
                "probability_numerics": {"method": "boundary_log_cumulative_v1",
                    "log_tolerance": 2e-5, "max_refinements": 8},
                "plot": {"status": "created", "path": str(p.with_suffix(".png"))}}
            p.write_text(json.dumps(obj))
            self.assertTrue(batch.existing_ok(p, "synthetic", settings, hashes))
            obj["provenance"]["probability_recursion_sha256"] = "stale"
            p.write_text(json.dumps(obj))
            self.assertFalse(batch.existing_ok(p, "synthetic", settings, hashes))
            obj["provenance"]["probability_recursion_sha256"] = hashes["probability_recursion.py"]
            obj["probability_numerics"]["log_tolerance"] = 1e-4
            p.write_text(json.dumps(obj))
            self.assertFalse(batch.existing_ok(p, "synthetic", settings, hashes))

    def test_comparison_dry_run_is_read_only(self):
        import io
        import sys
        from contextlib import redirect_stdout
        import run_inverse_rar_batch as batch
        bins = [[{"galaxy": "D512-2", "weight": 1, "knots": 19}]]
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "not_created"
            argv = ["run_inverse_rar_batch.py", "--dry-run", "--shards", "1",
                    "--galaxies", "D512-2", "--output-root", str(output)]
            with (patch.object(sys, "argv", argv),
                  patch.object(batch, "assignments", return_value=(bins, [1])),
                  patch.object(batch.subprocess, "run", side_effect=AssertionError("dry run launched work")),
                  redirect_stdout(io.StringIO())):
                self.assertEqual(batch.main(), 0)
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
