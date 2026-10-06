"""Lightweight standard-library tests for marginalisation and output contracts."""

import unittest
from types import SimpleNamespace

import numpy as np

from mono_SPARC_RAR import (
    append_row,
    fully_marginalise_grid,
    grid_prior_weights,
    trapezoid_weights,
)


class MarginalisationTests(unittest.TestCase):
    @staticmethod
    def args(profile="linear"):
        return SimpleNamespace(marginalise_disk=True, profile=profile, disk_ml=0.5, disk_ml_dex=0.1, slope_min=-5.0, slope_max=5.0, grid_size=101)

    def test_trapezoid_weights_cover_interval(self):
        weights = trapezoid_weights(np.linspace(-2.0, 3.0, 11))
        self.assertAlmostEqual(float(np.sum(weights)), 5.0)
        self.assertTrue(np.all(weights > 0.0))

    def test_declared_priors_are_normalised(self):
        args = self.args()
        x = np.linspace(np.log10(0.5) - 1.0, np.log10(0.5) + 1.0, 101)
        slopes = np.linspace(-5.0, 5.0, 101)
        disk, slope, covered_mass = grid_prior_weights(x, slopes, args)
        self.assertAlmostEqual(float(np.sum(disk)), 1.0)
        self.assertAlmostEqual(float(np.sum(slope)), 1.0)
        self.assertGreater(covered_mass, 0.999)

    def test_constant_profile_preserves_probability(self):
        args = self.args("constant")
        x = np.linspace(np.log10(0.5) - 1.0, np.log10(0.5) + 1.0, 101)
        slopes = np.linspace(-5.0, 5.0, 101)
        logp = np.full((1, 101), np.log(0.25))
        valid = np.ones_like(logp, dtype=bool)
        result = fully_marginalise_grid("synthetic", args, x, slopes, logp, valid, np.zeros_like(valid))
        self.assertTrue(result["valid"])
        self.assertEqual(result["probability_kind"], "fully_marginalised_sparc_prior")
        self.assertAlmostEqual(result["log_probability"], np.log(0.25))

    def test_linear_profile_uses_uniform_slope_prior(self):
        args = self.args("linear")
        x = np.linspace(np.log10(0.5) - 1.0, np.log10(0.5) + 1.0, 101)
        slopes = np.linspace(-5.0, 5.0, 101)
        logp = np.full((101, 101), np.log(0.4))
        valid = np.ones_like(logp, dtype=bool)
        result = fully_marginalise_grid("synthetic", args, x, slopes, logp, valid, np.zeros_like(valid))
        self.assertTrue(result["valid"])
        self.assertEqual(result["slope_prior"], "uniform(-5.0,5.0)")
        self.assertAlmostEqual(result["log_probability"], np.log(0.4))
        self.assertAlmostEqual(result["retained_profile_prior_mass"], 1.0)

    def test_grid_writer_omits_repeated_run_settings(self):
        rows = []
        append_row(rows, "synthetic", "linear", 0, 0, np.log10(0.5), 0.5, 0.0, np.log(0.5), True, "", [], self.args(), "conditional")
        self.assertEqual(set(rows[0]), {"galaxy", "profile", "probability_kind", "grid_i", "grid_j", "log10_upsilon_disk", "upsilon_disk_mean", "disk_slope", "log_probability", "valid", "invalid_reason", "warnings"})


if __name__ == "__main__":
    unittest.main()
