#!/usr/bin/env python3
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent


class TestMLProfiles(unittest.TestCase):
    def test_known_profiles_and_exclusions(self):
        with tempfile.TemporaryDirectory() as temp:
            temp = Path(temp)
            source = temp / 'input.csv'
            pd.DataFrame({
                'Galaxy': ['CONST'] * 3 + ['LINEAR'] * 3 + ['LOG'] * 3 + ['BAD'] * 2,
                'Sample': ['XMM'] * 3 + ['XMM'] * 3 + ['COSMOS_kin'] * 3 + ['LADUMA'] * 2,
                'Radius_kpc': [1., 2., 3., 1., 2., 3., 1., 2., 3., 1., 2.],
                'ml_K': [2., 2., 2., 1., 2., 3., 1., 2., 4., 0., 2.],
            }).to_csv(source, index=False)
            output = temp / 'out'
            subprocess.run(
                [sys.executable, str(ROOT / 'fit_ml_profiles.py'), '--input', str(source), '--output-dir', str(output)],
                check=True,
            )
            fits = pd.read_csv(output / 'ml_profile_fits.csv').set_index('galaxy')
            self.assertEqual(len(fits), 4)
            self.assertAlmostEqual(fits.loc['CONST', 'constant_c'], 2.)
            self.assertAlmostEqual(fits.loc['CONST', 'linear_slope'], 0.)
            self.assertAlmostEqual(fits.loc['LINEAR', 'linear_c_at_Rp'], 2.)
            self.assertAlmostEqual(fits.loc['LINEAR', 'linear_slope'], 1.)
            self.assertTrue(fits.loc['LOG', 'log_linear_success'])
            self.assertEqual(fits.loc['BAD', 'n_excluded'], 1)
            self.assertFalse(fits.loc['BAD', 'linear_success'])
            self.assertEqual(len(pd.read_csv(output / 'ml_profile_exclusions.csv')), 1)
            metadata = json.loads((output / 'ml_profile_metadata.json').read_text())
            self.assertEqual(metadata['n_unique_sample_galaxy_groups'], 4)


if __name__ == '__main__':
    unittest.main()
