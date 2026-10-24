"""Synthetic-only tests for the bounded symmetric-velocity prototype."""
import math
import unittest

import numpy as np

from bounded_uniform_bayes_factor_20261021 import (
    ProbabilityResult,
    _finite_support_probability,
    estimate_bayes_factor,
    estimate_broad_limit_bayes_factor,
)


class BoundedUniformBayesFactorTests(unittest.TestCase):
    def test_two_bin_unequal_radius_prior_matches_analytic_value(self):
        # X1~U(0,Vmax), X2~U(0,Vmax/2); P(X1>X2)=3/4.
        result = _finite_support_probability(
            np.array([1.0, 4.0]), np.array([2.0, 1.0]), np.zeros(2),
            np.ones(2), 10.0, 8001, 0.0, posterior=False,
        )
        self.assertAlmostEqual(result.probability, 0.75, delta=2e-4)

    def test_prior_quadrature_refinement(self):
        args = (np.array([1.0, 4.0]), np.array([2.0, 1.0]), np.zeros(2), np.ones(2))
        coarse = _finite_support_probability(*args, 10.0, 1001, 0.0, posterior=False)
        fine = _finite_support_probability(*args, 10.0, 16001, 0.0, posterior=False)
        self.assertLess(abs(fine.probability - 0.75), abs(coarse.probability - 0.75))
        self.assertAlmostEqual(fine.probability, 0.75, delta=2e-4)

    def test_prior_scale_invariance(self):
        args = (np.array([1.0, 4.0]), np.array([2.0, 1.0]), np.zeros(2), np.ones(2))
        a = _finite_support_probability(*args, 2.0, 2001, 0.0, posterior=False)
        b = _finite_support_probability(*args, 200.0, 2001, 0.0, posterior=False)
        self.assertAlmostEqual(a.log_probability, b.log_probability, delta=2e-10)

    def test_identical_posterior_bins_give_half_ordering_probability(self):
        result = _finite_support_probability(
            np.ones(2), np.array([2.0, 1.0]), np.array([30.0, 30.0]),
            np.array([7.0, 7.0]), 100.0, 8001, 0.0, posterior=True,
        )
        self.assertAlmostEqual(result.probability, 0.5, delta=2e-4)

    def test_symmetric_posterior_is_invariant_to_velocity_sign_flip(self):
        args = (np.array([1.0, 2.0, 3.0]), np.array([3.0, 2.0, 1.0]),
                np.array([18.0, 22.0, 9.0]), np.array([7.0, 5.0, 4.0]))
        a = _finite_support_probability(*args, 40.0, 4001, 0.0, posterior=True)
        flipped = (args[0], args[1], -args[2], args[3])
        b = _finite_support_probability(*flipped, 40.0, 4001, 0.0, posterior=True)
        self.assertAlmostEqual(a.log_probability, b.log_probability, delta=2e-10)

    def test_tied_gbar_uses_average_over_orderings(self):
        result = _finite_support_probability(
            np.ones(2), np.ones(2), np.array([20.0, 20.0]), np.array([3.0, 3.0]),
            40.0, 4001, 0.0, posterior=True,
        )
        self.assertEqual(result.orderings_averaged, 2)
        self.assertAlmostEqual(result.probability, 0.5, delta=2e-4)

    def test_mc_crosscheck_is_reproducible_and_reports_uncertainty(self):
        args = dict(radii=[1.0, 4.0], g_bar=[2.0, 1.0], v_obs=[3.0, 2.0],
                    sigma_v=[2.0, 2.0], vmax=10.0, seed=1729, ngrid=2001, mc_draws=30000)
        a = estimate_bayes_factor(**args)
        b = estimate_bayes_factor(**args)
        self.assertEqual(a.prior_monte_carlo.successes, b.prior_monte_carlo.successes)
        self.assertAlmostEqual(a.prior_monte_carlo.probability, 0.75, delta=0.015)
        self.assertGreater(a.prior_monte_carlo.standard_error, 0.0)
        self.assertAlmostEqual(a.log_bayes_factor,
                               a.posterior.log_probability - a.prior.log_probability)

    def test_zero_success_upper_bound_formula(self):
        estimate = estimate_bayes_factor(
            [1.0, 4.0, 16.0], [3.0, 2.0, 1.0], [100.0, 100.0, 100.0],
            [0.01, 0.01, 0.01], vmax=1.0, seed=7, ngrid=1001, mc_draws=5,
        )
        mc = estimate.prior_monte_carlo
        if mc.successes == 0:
            self.assertAlmostEqual(mc.zero_success_upper_95, 1.0 - 0.05 ** (1.0 / mc.draws))

    def test_iid_uniform_prior_ordering_matches_factorial(self):
        n = 5
        result = _finite_support_probability(
            np.ones(n), np.arange(n, 0, -1, dtype=float), np.zeros(n), np.ones(n),
            10.0, 4001, 0.0, posterior=False,
        )
        self.assertAlmostEqual(result.probability, 1.0 / math.factorial(n), delta=2e-5)

    def test_log_probability_is_retained_when_float_probability_underflows(self):
        result = ProbabilityResult(-800.0, 129, 1)
        self.assertEqual(result.log_probability, -800.0)
        self.assertEqual(result.probability, 0.0)

    def test_broad_limit_equal_gaussian_bins_give_half(self):
        result = estimate_broad_limit_bayes_factor(
            [1.0, 1.0], [2.0, 1.0], [30.0, 30.0], [7.0, 7.0],
            seed=19, ngrid=8001, mc_draws=20000, tail_tolerance=1e-10,
        )
        self.assertAlmostEqual(result["prior"].probability, 0.5, delta=2e-4)
        self.assertAlmostEqual(result["posterior_probability_lower"], 0.5, delta=2e-4)
        self.assertLessEqual(result["posterior_probability_lower"], result["posterior_probability_upper"])
        self.assertLessEqual(result["posterior_probability_upper"], 0.5 + 2e-4)

    def test_broad_limit_unbounded_gaussian_matches_large_finite_bound(self):
        radii = [1.0, 2.0]
        g_bar = [2.0, 1.0]
        v_obs = [25.0, 15.0]
        sigma_v = [2.0, 3.0]
        broad = estimate_broad_limit_bayes_factor(
            radii, g_bar, v_obs, sigma_v, seed=23, ngrid=12001,
            mc_draws=10000, tail_tolerance=1e-10,
        )
        finite = _finite_support_probability(
            np.asarray(radii), np.asarray(g_bar), np.asarray(v_obs), np.asarray(sigma_v),
            100.0, 12001, 0.0, posterior=True,
        )
        self.assertAlmostEqual(broad["posterior_probability_lower"], finite.probability, delta=3e-4)

    def test_explicit_bound_and_input_validation(self):
        with self.assertRaises(TypeError):
            estimate_bayes_factor([1, 2], [2, 1], [1, 1], [1, 1], seed=1)
        with self.assertRaisesRegex(ValueError, "Vmax"):
            estimate_bayes_factor([1, 2], [2, 1], [1, 1], [1, 1], vmax=0, seed=1)
        with self.assertRaisesRegex(ValueError, "strictly positive"):
            estimate_bayes_factor([0, 2], [2, 1], [1, 1], [1, 1], vmax=10, seed=1)


if __name__ == "__main__":
    unittest.main()
