"""Convex signed-net-acceleration positivity regressions."""
import unittest
import numpy as np
from inverse_rar_diagnostic import solve_qp

class PositivityTests(unittest.TestCase):
    def fit(self, D, c):
        return solve_qp(np.array(D, float), np.array(c, float),
                        np.arange(len(c)), np.array([0., 1.]), .01, 10., .3)

    def test_reference_retained(self):
        fit = self.fit([[1, 0], [0, 2]], [0, 0])
        np.testing.assert_array_equal(fit["xi"], [.5, .5])
        self.assertEqual(fit["objective"], 0.)
        self.assertGreater(min(fit["g_bar"]), 0.)

    def test_negative_ordered_reference_requires_positive_fit(self):
        fit = self.fit([[1, 0], [1, 0]], [-1, -.9])
        self.assertEqual(fit["status"], "feasible_weak_order")
        self.assertGreater(min(fit["g_bar"]), 0.)
        self.assertGreater(fit["objective"], 0.)
        self.assertTrue(fit["optimality"]["passed"])
        self.assertEqual(fit["constraint_version"], "positive_net_gbar_v1")

    def test_positive_infeasibility_is_not_accepted(self):
        fit = self.fit([[-1, 0], [-1, 0]], [0, 0])
        self.assertEqual(fit["status"], "infeasible")

    def test_nonpositive_floor_rejected(self):
        with self.assertRaises(ValueError):
            solve_qp(np.eye(2), np.ones(2), np.arange(2),
                     np.arange(2), .01, 10., .3, gbar_floor_relative=0.)

if __name__ == "__main__":
    unittest.main()
