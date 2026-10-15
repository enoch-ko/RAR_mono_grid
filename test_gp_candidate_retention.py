"""Run with python3 -B -m unittest test_gp_candidate_retention -v; no pytest."""
import ast
import json
import math
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import numpy as np
import pandas as pd
from gp_candidate_retention import CandidateRetention, DEFAULT_THRESHOLDS, limited_json, RetentionStorageLimitError
import gp_rar_search as old
import gp_rar_search_retained as new
import plot_ml_profile_clouds as clouds
import plot_gp_rar_search as plots
import summarize_gp_thresholds as thresholds
import run_gp_rar_retained_batch as batch


def row(j, p, source='GP_direction', valid=True, **extra):
    return dict(J_ref=j, conditional_logp=math.log(p), marginal_logp=math.log(p),
                source=source, valid=valid, **extra)


class RetentionTests(unittest.TestCase):
    def test_random_offline_minima_and_frontier(self):
        rng = np.random.default_rng(42)
        full = [row(float(j), float(p)) for j, p in zip(rng.random(2000)*100,
                                                       rng.uniform(1e-8, .9, 2000))]
        full[0]['source'] = 'SPARC_constant_0.5'
        keeper = CandidateRetention()
        for i, r in enumerate(full):
            keeper.add(r, [i, i+1])
        records, profiles, frontier, summary = keeper.export()
        for name, probability in DEFAULT_THRESHOLDS.items():
            eligible = [i for i, r in enumerate(full) if r['marginal_logp'] >= math.log(probability)]
            best = min(eligible, key=lambda i: (full[i]['J_ref'], i))
            saved = records[summary['threshold_winners'][name]['record_index']]
            self.assertEqual(saved['original_candidate_index'], best)
        sorted_ids = sorted(range(len(full)), key=lambda i: (full[i]['J_ref'], -full[i]['marginal_logp'], i))
        expected, best_lp = [], -math.inf
        for i in sorted_ids:
            if full[i]['marginal_logp'] > best_lp:
                expected.append(i)
                best_lp = full[i]['marginal_logp']
        self.assertEqual([records[i]['original_candidate_index'] for i in frontier], expected)
        for saved, xi in zip(records, profiles):
            self.assertEqual(xi.tolist(), [saved['original_candidate_index'], saved['original_candidate_index']+1])
        self.assertEqual(summary['evaluated'], 2000)

    def test_threshold_tie_pin_survives_dominance(self):
        keeper = CandidateRetention({'cut': .05})
        keeper.add(row(0., .01, 'SPARC_constant_0.5'), [.5])
        keeper.add(row(1., .06), [.6])
        keeper.add(row(1., .08), [.7])
        records, _, frontier, summary = keeper.export()
        self.assertEqual(records[summary['threshold_winners']['cut']['record_index']]['original_candidate_index'], 1)
        self.assertEqual([records[i]['original_candidate_index'] for i in frontier], [0, 2])

    def test_near_cut_no_epsilon_loss(self):
        keeper = CandidateRetention({'cut': .05})
        keeper.add(row(1., .05-1e-15), [.5])
        keeper.add(row(2., .05+1e-15), [.6])
        records, _, frontier, summary = keeper.export()
        self.assertEqual(len(frontier), 2)
        self.assertEqual(records[summary['threshold_winners']['cut']['record_index']]['J_ref'], 2.)

    def test_invalid_baseline_counts_and_unreachable(self):
        keeper = CandidateRetention({'cut': .9})
        for i in range(100):
            keeper.add(row(0., .1, 'SPARC_constant_0.5' if i == 0 else 'GP_direction',
                           False, reason='signed_gbar_below_floor'), [.5])
        records, _, frontier, summary = keeper.export()
        self.assertEqual(len(records), 1)
        self.assertEqual(frontier, [])
        self.assertEqual(summary['invalid'], 100)
        self.assertEqual(summary['failure_counts'], {'signed_gbar_below_floor': 100})
        self.assertEqual(len(summary['failure_examples']), 8)
        self.assertIsNone(summary['threshold_winners']['cut']['record_index'])

    def test_dominated_collection_is_bounded(self):
        keeper = CandidateRetention()
        keeper.add(row(0., .5, 'SPARC_constant_0.5'), [.5])
        for i in range(10000):
            keeper.add(row(float(i+1), .1), [.6])
        records, _, _, summary = keeper.export()
        self.assertEqual(len(records), 1)
        self.assertEqual(summary['discarded'], 10000)

    def test_extra_probability_and_invalid_scores(self):
        keeper = CandidateRetention({'4.5percent': .045})
        keeper.add(row(1., .045), [.5])
        self.assertEqual(keeper.export()[3]['threshold_winners']['4.5percent']['status'], 'reached')
        for probability in [0., -1., float('nan'), 1.01]:
            with self.assertRaises(ValueError):
                CandidateRetention({'bad': probability})
        with self.assertRaises(ValueError):
            CandidateRetention().add(row(float('nan'), .2), [.5])
        with self.assertRaises(ValueError):
            CandidateRetention().add(row(1., .2), [float('nan')])

    def test_numerical_functions_identical_to_frozen_parent(self):
        def definitions(path):
            return {node.name: ast.dump(node, include_attributes=False)
                    for node in ast.parse(Path(path).read_text()).body
                    if isinstance(node, (ast.FunctionDef, ast.ClassDef))}
        parent, retained = definitions(old.__file__), definitions(new.__file__)
        for name in ['nested_direction_rng', 'proposal_cholesky', 'physical_interval',
                     'amplitude_grid', 'ProbabilityEvaluator', 'reference_profile',
                     'warm_reference_profile']:
            self.assertEqual(parent[name], retained[name], name)

    def test_mocked_driver_matches_full_storage(self):
        """Synthetic dependency mocks: no galaxy files, force integrals or probability runs."""
        radius = np.array([1., 2.])
        frame = pd.DataFrame({'Galaxy': ['D512-2'], 'Rdisk': [1.]})
        args = SimpleNamespace(sampling_scheme='nested-v1', replicate=0, seed=20261013,
             directions=8, amplitude_steps=6, ell_over_rd=[1.], length_fraction=.3,
             min_ml=.01, max_ml=10., gbar_floor_relative=1e-8, gas_ml=1., gas_ml_dex=.04,
             bulge_ml=.7, bulge_ml_dex=.1, ml_samples=1, prob_grid_size=32, zpadding=8.,
             boundary_log_tolerance=2e-5, boundary_max_refinements=8,
             data_root=Path('/synthetic'), dens_root=Path('/synthetic'),
             validation_root=Path('/synthetic'), retain_probability=[])
        def evaluate(self, xi):
            return dict(valid=True, conditional_logp=-float(np.sum((xi-.8)**2)),
                        marginal_logp=-float(np.sum((xi-.8)**2)), minimum_nuisance_gbar_m_s2=.1)
        with tempfile.TemporaryDirectory(prefix='rar-retention-test-') as tmp:
            results = []
            for module, mode in [(old, 'compact'), (new, 'retained-only')]:
                args.output_root = Path(tmp) / mode
                args.checkpoint_mode = mode
                with patch.object(module.inverse, 'sha256', return_value='synthetic-hash'), \
                     patch.object(module.inverse, 'cache_gate', return_value={'passed': True}), \
                     patch.object(module.inverse, 'galaxy_data', return_value=(radius, np.ones(2), radius,
                          np.ones(2), np.ones(2), np.zeros(2), np.zeros(2), .1, [])), \
                     patch.object(module.inverse, 'build_response', return_value=(np.eye(2), np.eye(2), [])), \
                     patch.object(module.inverse, 'direct_check', return_value={'relative_to_peak': 0.}), \
                     patch.object(module.inverse, 'kernel', return_value=np.eye(2)), \
                     patch.object(module, 'load_sparc_table', return_value=frame), \
                     patch.object(module.ProbabilityEvaluator, 'evaluate', evaluate):
                    results.append(module.galaxy_run('D512-2', args))
            full, saved = results
            self.assertEqual(saved['checkpoint_storage']['summary']['evaluated'], len(full['records']))
            for r in saved['records']:
                i = r['original_candidate_index']
                self.assertEqual({k:v for k,v in r.items() if k != 'original_candidate_index'}, full['records'][i])
            selected = saved['retained_profiles']
            original_array = np.load(full['checkpoint_storage']['profiles_file'])
            for item in selected:
                i = saved['records'][item['record_index']]['original_candidate_index']
                np.testing.assert_array_equal(item['xi'], original_array[i])
            paths = list((Path(tmp) / 'retained-only').rglob('*'))
            self.assertFalse(any(p.suffix == '.npy' or p.name.startswith('candidate_') for p in paths))
            progress = json.loads((Path(tmp)/'retained-only'/'D512-2_gp_checkpoints'/'progress.json').read_text())
            self.assertEqual(progress['records_so_far'], len(full['records']))
            self.assertEqual(len(progress['retained_snapshot']['records']), len(saved['records']))
            json.dumps(saved, allow_nan=False)
            path = Path(tmp)/'retained-only'/'D512-2_gp_search.json'
            path.write_text(json.dumps(saved, allow_nan=False))
            report, _ = thresholds.audit(path)
            self.assertEqual(report['records'], len(full['records']))
            self.assertEqual(report['retained_records'], len(saved['records']))
            for name, p in thresholds.THRESHOLDS.items():
                expected = thresholds.minimum(full['records'], p)
                actual = report['minimum_found_J'][name]
                self.assertEqual(None if expected is None else expected['J_ref'],
                                 None if actual is None else actual['J_ref'])
            panels, _ = clouds.load_gp([Path(tmp)/'retained-only'], [2., 5.], 'one-sided')
            self.assertTrue(all(len(values) == 1 for values in panels.values()))
            self.assertTrue(plots.targets(path, saved, Path(tmp)/'plots'))
            with patch.object(plots.inverse, 'sha256', return_value='synthetic-hash'):
                plots.check_provenance(saved, Path('/synthetic'), Path('/synthetic'))
                bad = dict(saved, provenance=dict(saved['provenance'], script_file='../outside.py'))
                with self.assertRaises(ValueError):
                    plots.check_provenance(bad, Path('/synthetic'), Path('/synthetic'))

    def test_byte_guard_never_truncates(self):
        value = {"xi": [0.5] * 1000}
        with self.assertRaises(RetentionStorageLimitError):
            limited_json(value, 1024)
        self.assertEqual(json.loads(limited_json(value, 10000)), value)

    def test_future_launcher_cannot_select_full_storage(self):
        command = batch.payload('D512-2', Path('/synthetic/output'))
        self.assertTrue(command[0].endswith('/gp_rar_search_retained.py'))
        self.assertEqual(command[command.index('--checkpoint-mode')+1], 'retained-only')
        self.assertTrue(str(batch.DEFAULT_OUTPUT).endswith('gp_rar_full_retained_20261014_d5000_a16'))



if __name__ == '__main__':
    unittest.main()
