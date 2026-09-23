"""Lightweight unit tests for the prior-marginalisation machinery.

These tests do not read SPARC files or launch a galaxy analysis. They isolate the
numerical quadrature and output-contract logic so they can run on a login node.
"""

import unittest
from types import SimpleNamespace

import numpy as np

from mono_SPARC_RAR import (
    fully_marginalise_grid,
    grid_prior_weights,
    trapezoid_weights,
)


class MarginalisationTests(unittest.TestCase):
    """Check normalisation, prior shape, and deterministic integration cases."""

    @staticmethod
    def args(profile="linear"):
        return SimpleNamespace(
            marginalise_disk=True,
            profile=profile,
            disk_ml=0.5,
            disk_ml_dex=0.1,
            slope_min=-5.0,
            slope_max=5.0,
            slope_prior_sigma=0.5,
            grid_size=101,
        )

    def test_trapezoid_weights_cover_interval(self):
        coordinates = np.linspace(-2.0, 3.0, 11)
        weights = trapezoid_weights(coordinates)
        self.assertAlmostEqual(float(np.sum(weights)), 5.0)
        self.assertTrue(np.all(weights > 0.0))

    def test_declared_priors_are_normalised_and_symmetric(self):
        args = self.args()
        log10_ml = np.linspace(np.log10(0.5) - 1.0, np.log10(0.5) + 1.0, 101)
        slopes = np.linspace(-5.0, 5.0, 101)
        disk, slope, disk_mass, slope_mass = grid_prior_weights(log10_ml, slopes, args)

        self.assertAlmostEqual(float(np.sum(disk)), 1.0)
        self.assertAlmostEqual(float(np.sum(slope["broad_uniform"])), 1.0)
        self.assertAlmostEqual(float(np.sum(slope["sparc_gaussian"])), 1.0)
        np.testing.assert_allclose(slope["sparc_gaussian"], slope["sparc_gaussian"][::-1])
        self.assertGreater(disk_mass, 0.999)
        self.assertGreater(slope_mass["sparc_gaussian"], 0.999)
        broad_mean_abs = np.sum(np.abs(slopes) * slope["broad_uniform"])
        gaussian_mean_abs = np.sum(np.abs(slopes) * slope["sparc_gaussian"])
        self.assertLess(gaussian_mean_abs, broad_mean_abs)

    def test_constant_profile_preserves_constant_probability(self):
        args = self.args(profile="constant")
        log10_ml = np.linspace(np.log10(0.5) - 1.0, np.log10(0.5) + 1.0, 101)
        slopes = np.linspace(-5.0, 5.0, 101)
        logp = np.full((1, 101), np.log(0.25))
        valid = np.ones_like(logp, dtype=bool)
        negative = np.zeros_like(logp, dtype=bool)

        rows = fully_marginalise_grid("synthetic", args, log10_ml, slopes, logp, valid, negative)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["probability_kind"], "fully_marginalised_sparc_disk_prior")
        self.assertAlmostEqual(rows[0]["log_probability"], np.log(0.25))

    def test_radial_profiles_emit_both_prior_cases(self):
        args = self.args(profile="linear")
        log10_ml = np.linspace(np.log10(0.5) - 1.0, np.log10(0.5) + 1.0, 101)
        slopes = np.linspace(-5.0, 5.0, 101)
        logp = np.full((101, 101), np.log(0.4))
        valid = np.ones_like(logp, dtype=bool)
        negative = np.zeros_like(logp, dtype=bool)

        rows = fully_marginalise_grid("synthetic", args, log10_ml, slopes, logp, valid, negative)
        self.assertEqual(
            {row["probability_kind"] for row in rows},
            {
                "fully_marginalised_broad_slope_prior",
                "fully_marginalised_sparc_slope_prior",
            },
        )
        for row in rows:
            self.assertTrue(row["valid"])
            self.assertAlmostEqual(row["log_probability"], np.log(0.4))


if __name__ == "__main__":
    unittest.main()
