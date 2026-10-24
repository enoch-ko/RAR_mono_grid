"""Warm-start changes initial iterate only; acceptance cannot be bypassed."""
import unittest
from unittest.mock import patch
import numpy as np
import check_osqp_warm_start as warm
import inverse_rar_diagnostic as inv


class WarmStartTests(unittest.TestCase):
    def run_fit(self, seed, eps=1e-9):
        return warm.solve_seeded(np.eye(2), np.array([1., 0.]), np.array([1., 2.]),
                                 np.array([0., 1.]), {"min_ml": .01, "max_ml": 10.,
                                 "length_fraction": .3, "gbar_floor_relative": 1e-8}, seed, 1234567, 10., eps)

    def test_seed_forwarded_not_accepted(self):
        seed = np.array([3., 4.])
        seen = []
        def fake(C, rhs, z0, eps, budget, limit):
            seen.append(z0.copy())
            self.assertEqual(budget, 1234567)
            return z0, np.zeros(len(rhs)), True, {"status": "solved", "iterations": 1}
        with patch.object(inv, "osqp_solve", fake):
            fit, initial = self.run_fit(seed)
        np.testing.assert_array_equal(seen[0], seed)
        self.assertEqual(initial["objective"], 25.)
        self.assertEqual(fit["status"], "qp_failure")
        self.assertNotIn("xi", fit)

    def test_accuracy_forwarded(self):
        seen = []
        def fake(C, rhs, z0, eps, budget, limit):
            seen.append(eps)
            return z0, None, False, {"status": "maximum iterations reached", "iterations": 1}
        with patch.object(inv, "osqp_solve", fake):
            self.run_fit(None, eps=1e-11)
        self.assertEqual(seen, [1e-11])

    def test_invalid_seed_rejected(self):
        for seed in ([1.], [np.nan, 1.]):
            with self.subTest(seed=seed):
                fit, initial = self.run_fit(seed)
                self.assertEqual(fit["solver_diagnostics"]["failure_stage"], "solver_exception")
                self.assertNotIn("xi", fit)

    def test_cold_start_retained(self):
        def fake(C, rhs, z0, eps, budget, limit):
            return z0, None, False, {"status": "maximum iterations reached", "iterations": 1}
        with patch.object(inv, "osqp_solve", fake):
            fit, initial = self.run_fit(None)
        self.assertEqual(initial["kind"], "lp")
        self.assertEqual(fit["status"], "qp_failure")


if __name__ == "__main__":
    unittest.main()