#!/usr/bin/env python3
"""Lightweight regression tests for sparc_scale_proxies.py (no pytest required)."""
from pathlib import Path
import tempfile
import unittest

import numpy as np

from sparc_scale_proxies import (
    A0_DEFAULT_M_S2,
    G_KPC_KMS2_PER_MSUN,
    KPC_METRES,
    OUTPUT_FIELDS,
    disk_mass_for_profile,
    mass_distribution,
    prepare_light_profile,
)
from utils import compute_A_B, lognormal_quadrature


class ScaleProxyTests(unittest.TestCase):
    def setUp(self):
        self.radius = np.linspace(0.0, 6.0, 121)
        self.disk_sb = 2.0e6 * np.exp(-self.radius / 1.7)
        self.bulge_sb = 1.0e6 * np.exp(-self.radius / 0.45)
        self.disk_raw_lum, _ = compute_A_B(self.radius, self.disk_sb)
        self.bulge_raw_lum, _ = compute_A_B(self.radius, self.bulge_sb)
        self.catalog_lum_1e9 = 8.0
        scale = self.catalog_lum_1e9 * 1.0e9 / (self.disk_raw_lum + self.bulge_raw_lum)
        self.light = {
            "radius": self.radius,
            "raw_disk_sb": self.disk_sb,
            "scale": scale,
            "disk_luminosity_1e9": self.disk_raw_lum * scale / 1.0e9,
            "bulge_luminosity_1e9": self.bulge_raw_lum * scale / 1.0e9,
        }

    def test_all_radial_profiles_conserve_luminosity_weighted_disk_mass(self):
        mean_upsilon = 0.5
        expected = mean_upsilon * self.light["disk_luminosity_1e9"]
        cases = [("constant", 0.0), ("linear", 0.025), ("linear", -0.025),
                 ("log_linear", 0.08), ("log_linear", -0.08)]
        for profile, slope in cases:
            with self.subTest(profile=profile, slope=slope):
                mass, radial_upsilon = disk_mass_for_profile(profile, mean_upsilon, slope, self.light)
                self.assertTrue(np.all(radial_upsilon > 0.0))
                self.assertAlmostEqual(mass, expected, delta=1e-11 * expected)

    def test_bulge_light_is_subtracted_from_catalogue_total(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "synthetic.dens"
            rows = ["# Rad[kpc]\tSBdisk[Lsun/pc^2]\tSBbul[Lsun/pc^2]"]
            rows.extend(f"{r:.8g}\t{d:.8g}\t{b:.8g}" for r, d, b in zip(self.radius[::30], self.disk_sb[::30] / 1e6, self.bulge_sb[::30] / 1e6))
            path.write_text("\n".join(rows) + "\n")
            light = prepare_light_profile(path, self.catalog_lum_1e9)
        self.assertGreater(light["bulge_luminosity_1e9"], 0.0)
        self.assertAlmostEqual(light["disk_luminosity_1e9"] + light["bulge_luminosity_1e9"], self.catalog_lum_1e9, places=12)
        self.assertAlmostEqual(light["disk_luminosity_1e9"], self.catalog_lum_1e9 - light["bulge_luminosity_1e9"], places=12)

    def test_mass_budget_and_proxy_units(self):
        gas_nodes, gas_weights = lognormal_quadrature(1.0, 0.04, 9)
        bulge_nodes, bulge_weights = lognormal_quadrature(0.7, 0.1, 9)
        m_disk = 0.5 * self.light["disk_luminosity_1e9"]
        mhi = 2.0
        lbulge = self.light["bulge_luminosity_1e9"]
        a_proxy = 2.4 / 0.85
        result = mass_distribution(m_disk, mhi, lbulge, gas_nodes, gas_weights,
                                   bulge_nodes, bulge_weights, 1.0, 0.7, 1.33,
                                   A0_DEFAULT_M_S2, a_proxy)
        expected_mbar = m_disk + 1.33 * mhi + 0.7 * lbulge
        a0_kms2_kpc = A0_DEFAULT_M_S2 * KPC_METRES / 1.0e6
        expected_rm = np.sqrt(G_KPC_KMS2_PER_MSUN * expected_mbar * 1.0e9 / a0_kms2_kpc)
        self.assertAlmostEqual(result["M_bar_median_1e9_Msun"], expected_mbar, places=12)
        self.assertAlmostEqual(result["r_M_median_kpc"], expected_rm, places=12)
        self.assertAlmostEqual(result["eta_median"], a_proxy / expected_rm, places=12)


    def test_output_schema_is_compact_and_unique(self):
        self.assertEqual(len(OUTPUT_FIELDS), len(set(OUTPUT_FIELDS)))
        self.assertTrue({
            "galaxy", "grid_i", "grid_j", "upsilon_disk_mean", "disk_slope",
            "source_grid_valid", "mass_status", "M_disk_1e9_Msun",
            "M_bar_median_1e9_Msun", "r_M_median_kpc", "eta_median",
        }.issubset(OUTPUT_FIELDS))
        repeated_or_derivable = {
            "profile", "source_probability_kind", "log10_upsilon_disk", "mass_valid",
            "L_total_catalog_1e9_Lsun", "L_disk_1e9_Lsun", "L_bulge_1e9_Lsun",
            "bulge_light_fraction", "a_proxy_kpc", "Rdisk_kpc", "MHI_raw_1e9_Msun",
            "gas_helium_factor", "M_gas_median_1e9_Msun", "M_bulge_median_1e9_Msun",
        }
        self.assertTrue(repeated_or_derivable.isdisjoint(OUTPUT_FIELDS))


if __name__ == "__main__":
    unittest.main(verbosity=2)
