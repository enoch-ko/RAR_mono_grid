"""Lightweight exact/reference tests: python3 -m unittest -v test_probability_recursion."""
import math
import unittest
import numpy as np
from scipy.integrate import quad
from scipy.special import ndtr
from probability_recursion import (boundary_log_probability, ProbabilityConvergenceError,
                                   boundary_log_probability_for_gbar)


class BoundaryRecursionTests(unittest.TestCase):
    def test_single(self):
        self.assertEqual(boundary_log_probability([1], [1], [1]), 0.0)

    def test_exchangeable(self):
        for n in (2, 3, 5, 8):
            with self.subTest(n=n):
                p = boundary_log_probability([1]*n, [1]*n, [1]*n)
                self.assertAlmostEqual(p, -math.lgamma(n+1), delta=2e-5)

    def test_two_zero_mean_exact(self):
        for ratio in (.001, .2, 1., 5., 1000.):
            with self.subTest(ratio=ratio):
                expected = np.log(2/np.pi*np.arctan(ratio))
                actual = boundary_log_probability([0, 0], [ratio, 1], [1, 1])
                self.assertAlmostEqual(actual, expected, delta=2e-5)

    def test_two_nonzero_means_reference(self):
        mu = [1.2, 2.3]
        scale = [2., .7]
        f = lambda x: (np.exp(-.5*(x/scale[0]-mu[0])**2)
                        +np.exp(-.5*(-x/scale[0]-mu[0])**2))/(scale[0]*np.sqrt(2*np.pi))
        cdf = lambda x: ndtr(x/scale[1]-mu[1])-ndtr(-x/scale[1]-mu[1])
        expected = np.log(quad(lambda x: f(x)*cdf(x), 0, 30, epsabs=1e-13)[0])
        actual = boundary_log_probability(np.array(mu)*scale, scale, [1, 1])
        self.assertAlmostEqual(actual, expected, delta=2e-5)

    def test_narrow_interval_reference(self):
        expected = quad(lambda z: .5*np.exp(-.5*(z-1)**2)/np.sqrt(2*np.pi)
                        *(ndtr(.001*abs(z)-1)-ndtr(-.001*abs(z)-1))**2,
                        -10, 12, epsabs=1e-16)[0]
        for resolution in (400, 401, 800):
            actual = boundary_log_probability([.001, 1, 1], [.001, 1, 1],
                                               [1, 1, 1], ngrid=resolution)
            self.assertAlmostEqual(actual, np.log(expected), delta=2e-5)

    def test_below_old_floor(self):
        actual = boundary_log_probability([1, 80], [1, 1], [1, 1])
        self.assertTrue(np.isfinite(actual))
        self.assertLess(actual, -1000)

    def test_sign_and_scale_invariance(self):
        v = np.array([1., 2., 3.])
        e = np.array([.2, 1., 2.])
        r = np.array([1., 2., 3.])
        expected = boundary_log_probability(v, e, r)
        self.assertAlmostEqual(expected, boundary_log_probability(-v, e, r), delta=1e-10)
        self.assertAlmostEqual(expected, boundary_log_probability(v*1e4, e*1e4, r*7), delta=1e-9)

    def test_padding_convergence(self):
        a = boundary_log_probability([1, 2, 3], [1, 1, 1], [1, 2, 3], zpadding=8)
        b = boundary_log_probability([1, 2, 3], [1, 1, 1], [1, 2, 3], zpadding=10)
        self.assertAlmostEqual(a, b, delta=2e-5)

    def test_refinement_recorded(self):
        d = boundary_log_probability([1, 1], [1, 1], [1, 1], return_diagnostics=True)
        self.assertGreaterEqual(len(d['levels']), 2)
        self.assertLessEqual(d['levels'][-1]['log_delta'], d['log_tolerance'])

    def test_nonconvergence_not_zero(self):
        with self.assertRaises(ProbabilityConvergenceError):
            boundary_log_probability([1]*5, [1]*5, [1]*5, ngrid=32,
                                     max_refinements=1, log_tolerance=1e-14)

    def test_baryonic_sort(self):
        g = np.array([1., 3., 2.])
        v, e, r = np.array([1., 2., 3.]), np.ones(3), np.array([1., 2., 3.])
        o = np.argsort(-g)
        expected = boundary_log_probability(v[o], e[o], r[o])
        actual = boundary_log_probability_for_gbar(g, v, e, r)
        self.assertAlmostEqual(actual, expected, delta=1e-10)

    def test_tied_order_average(self):
        actual = boundary_log_probability_for_gbar([1., 1.], [1., 3.], [1., 1.], [1., 1.])
        self.assertAlmostEqual(actual, -np.log(2), delta=2e-5)

    def test_tie_cap_checked_before_enumeration(self):
        from unittest.mock import patch
        with patch('probability_recursion.boundary_log_probability', return_value=-1.) as solver:
            with self.assertWarnsRegex(RuntimeWarning, '5040'):
                boundary_log_probability_for_gbar([1]*7, [1]*7, [1]*7, [1]*7)
            self.assertEqual(solver.call_count, 1)

    def test_bad_inputs(self):
        for args in (([], [], []), ([1], [0], [1]), ([1], [1], [0]),
                     ([np.nan], [1], [1]), ([1, 2], [1], [1])):
            with self.subTest(args=args), self.assertRaises(ValueError):
                boundary_log_probability(*args)


if __name__ == '__main__':
    unittest.main()
