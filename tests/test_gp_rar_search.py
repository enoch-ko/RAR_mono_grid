"""Focused invariants for the exploratory GP proposal search."""
import unittest
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

import gp_rar_search as gp


class GPProposalTests(unittest.TestCase):
    def test_proposal_lengthscale_is_in_rd_not_reference_span(self):
        knots = np.array([1.0, 3.0, 7.0])
        L, ell = gp.proposal_cholesky(knots, 2.0, 0.5)
        self.assertEqual(ell, 1.0)
        np.testing.assert_allclose(L @ L.T,
            np.exp(-0.5 * (knots[:, None] - knots[None, :])**2)
            + 1e-8 * np.eye(3), rtol=1e-12, atol=1e-12)

    def test_signed_nuisance_and_knot_bounds_limit_amplitude(self):
        mean = np.array([0.5, 0.5])
        direction = np.array([1.0, -1.0])
        D = np.eye(2)
        nuisance = np.array([[-0.6, 0.4]])
        lo, hi = gp.physical_interval(mean, direction, 0.01, 10.0,
                                      D, nuisance, 0.01)
        self.assertAlmostEqual(lo, 0.11)
        self.assertAlmostEqual(hi, 0.49)
        amplitudes = gp.amplitude_grid((lo, hi), 6)
        self.assertAlmostEqual(amplitudes[0], lo)
        self.assertAlmostEqual(amplitudes[-1], hi)
        for amplitude in amplitudes:
            xi = mean + amplitude * direction
            self.assertGreaterEqual(float(np.min(nuisance + D @ xi)), 0.01 - 1e-12)
            self.assertGreaterEqual(float(np.min(xi)), 0.01 - 1e-12)

    def test_frontier_uses_same_fixed_J_across_proposals(self):
        rows = [dict(valid=True, J_ref=0.0, marginal_logp=-2.0),
                dict(valid=True, J_ref=1.0, marginal_logp=-1.0),
                dict(valid=True, J_ref=0.5, marginal_logp=-3.0),
                dict(valid=True, J_ref=2.0, marginal_logp=-1.0),
                dict(valid=False, J_ref=0.1)]
        self.assertEqual(gp.frontier_indices(rows), [0, 1])

    def test_probability_cache_depends_on_ordering(self):
        options = SimpleNamespace(boundary_log_tolerance=2e-5,
                                  boundary_max_refinements=8)
        evaluator = gp.ProbabilityEvaluator(
            np.eye(2), np.ones(2), np.ones(2), np.ones(2),
            np.zeros(2), np.zeros(2), options,
            np.array([1.0]), np.array([1.0]),
            np.array([0.7]), np.array([1.0]),
            np.ones((1, 2)), 1e-8)
        with patch.object(gp.boundary, 'boundary_log_probability_for_gbar',
                          return_value=-0.7) as solver:
            x = dict(ngrid=101, zpadding=8.0, tie_atol=0.0,
                     max_tie_permutations=720)
            self.assertEqual(evaluator._ordered_logp(np.array([1., 2.]), **x), -0.7)
            self.assertEqual(evaluator._ordered_logp(np.array([3., 4.]), **x), -0.7)
            self.assertEqual(evaluator._ordered_logp(np.array([4., 3.]), **x), -0.7)
            self.assertEqual(solver.call_count, 2)
            self.assertEqual(evaluator.cache_hits, 1)

    def test_nonpositive_signed_gbar_is_rejected_before_probability(self):
        options = SimpleNamespace()
        evaluator = gp.ProbabilityEvaluator(
            np.array([[1.0]]), np.array([1.0]), np.array([2.0]),
            np.array([0.1]), np.array([0.0]), np.array([0.0]), options,
            np.array([1.0]), np.array([1.0]),
            np.array([0.7]), np.array([1.0]),
            np.array([[-1.0]]), 1e-8)
        with patch.object(gp.mono, 'evaluate_disk_term') as probability:
            result = evaluator.evaluate(np.array([0.5]))
            self.assertFalse(result['valid'])
            probability.assert_not_called()


    def test_real_boundary_probability_interface_for_two_points(self):
        options = SimpleNamespace(gas_ml=1.0, bulge_ml=0.7,
                                  marginalise_gas_bulge=True,
                                  prob_grid_size=51, zpadding=8.0,
                                  tie_atol=0.0, max_tie_permutations=720,
                                  boundary_log_tolerance=2e-5,
                                  boundary_max_refinements=8)
        evaluator = gp.ProbabilityEvaluator(
            np.eye(2) * 1e-11, np.array([1.0, 2.0]),
            np.array([20.0, 25.0]), np.array([2.0, 2.0]),
            np.zeros(2), np.zeros(2), options,
            np.array([1.0]), np.array([1.0]),
            np.array([0.7]), np.array([1.0]),
            np.zeros((1, 2)), 1e-14)
        result = evaluator.evaluate(np.array([0.5, 0.5]))
        self.assertTrue(result["valid"], result)
        self.assertLessEqual(result["marginal_logp"], 0.0)
        self.assertAlmostEqual(result["marginal_logp"],
                               result["conditional_logp"], places=12)

    def test_warning_is_not_erased_by_order_cache(self):
        import warnings
        options = SimpleNamespace(boundary_log_tolerance=2e-5,
                                  boundary_max_refinements=8)
        evaluator = gp.ProbabilityEvaluator(
            np.eye(2), np.ones(2), np.ones(2), np.ones(2),
            np.zeros(2), np.zeros(2), options,
            np.array([1.0]), np.array([1.0]),
            np.array([0.7]), np.array([1.0]),
            np.ones((1, 2)), 1e-8)
        def warn_once(*values, **options):
            warnings.warn("tail truncation", RuntimeWarning)
            return -0.7
        with patch.object(gp.boundary, 'boundary_log_probability_for_gbar',
                          side_effect=warn_once) as solver:
            for _ in range(2):
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter('always')
                    evaluator._ordered_logp(np.array([1., 2.]),
                                            ngrid=101, zpadding=8.0,
                                            tie_atol=0.0, max_tie_permutations=720)
                self.assertEqual(len(caught), 1)
                self.assertIn('tail truncation', str(caught[0].message))
            self.assertEqual(solver.call_count, 2)
            self.assertEqual(len(evaluator.memo), 0)

    def test_imported_qp_requires_saved_acceptance_audits(self):
        import json
        import tempfile
        from pathlib import Path
        inputs = {"rotmod": "a", "dens": "b", "table": "c"}
        data = {"galaxy": "UGC07399", "input_sha256": inputs,
                "knots_kpc": [1.0, 2.0], "bounds_ml": [0.01, 10.0],
                "kernel": {"length_fraction": 0.3, "nugget": 1e-8},
                "gas_ml_qp": 1.0, "bulge_ml_qp": 0.7,
                "optimized_direct_cache_qc": {"relative_to_peak": 0.0},
                "provenance": {"driver_sha256": gp.inverse.sha256(gp.mono.__file__),
                               "utils_sha256": gp.inverse.sha256(
                                   Path(gp.DEFAULT_PROJECT_ROOT) / "utils.py")},
                "optimization": {"status": "feasible_weak_order", "xi": [0.5, 0.5],
                                 "constraint_version": "positive_net_gbar_v1",
                                 "gbar_floor_relative": 1e-8,
                                 "optimality": {"passed": False},
                                 "max_scaled_violation": 0.0,
                                 "max_bound_violation": 0.0}}
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "old.json"
            path.write_text(json.dumps(data))
            args = SimpleNamespace(length_fraction=0.3, gas_ml=1.0, bulge_ml=0.7,
                                   gbar_floor_relative=1e-8)
            with self.assertRaisesRegex(ValueError, "audit failed"):
                gp.reference_profile(path, "UGC07399", np.array([1.0, 2.0]),
                                     0.01, 10.0, args, inputs)
            data["optimization"]["optimality"]["passed"] = True
            path.write_text(json.dumps(data))
            xi, info = gp.reference_profile(path, "UGC07399",
                np.array([1.0, 2.0]), 0.01, 10.0, args, inputs)
            np.testing.assert_array_equal(xi, np.array([0.5, 0.5]))
            self.assertEqual(info["status"], "saved_audits_passed_current_model_checked")
            data["optimization"]["max_scaled_violation"] = None
            path.write_text(json.dumps(data))
            with self.assertRaisesRegex(ValueError, "audit failed"):
                gp.reference_profile(path, "UGC07399", np.array([1.0, 2.0]),
                                     0.01, 10.0, args, inputs)
            data["optimization"]["max_scaled_violation"] = 0.0
            data["optimization"]["gbar_floor_relative"] = 1e-9
            path.write_text(json.dumps(data))
            with self.assertRaisesRegex(ValueError, "model/settings differ"):
                gp.reference_profile(path, "UGC07399", np.array([1.0, 2.0]),
                                     0.01, 10.0, args, inputs)

if __name__ == '__main__':
    unittest.main()
