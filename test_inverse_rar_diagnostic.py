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


if __name__ == "__main__":
    unittest.main()


