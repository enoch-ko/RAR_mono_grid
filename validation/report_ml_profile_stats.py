#!/usr/bin/env python3
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd

SPECS = {
    'constant': {'normalization': 'constant_c', 'success': 'constant_success', 'intercept': None, 'slope': None},
    'linear': {'normalization': 'linear_c_at_Rp', 'success': 'linear_success', 'intercept': 'linear_intercept_c0', 'slope': 'linear_slope'},
    'log_linear': {'normalization': 'log_linear_c_at_Rp', 'success': 'log_linear_success', 'intercept': 'log_linear_intercept_c0', 'slope': 'log_linear_slope'},
}
GROUPS = ['COSMOS', 'LADUMA', 'XMM', 'ALL']

def bootstrap(values, n_boot, rng):
    if len(values) == 0:
        return {}
    draws = rng.choice(values, size=(n_boot, len(values)), replace=True)
    return {
        'mean_ci_low': float(np.percentile(draws.mean(axis=1), 2.5)),
        'mean_ci_high': float(np.percentile(draws.mean(axis=1), 97.5)),
        'median_ci_low': float(np.percentile(np.median(draws, axis=1), 2.5)),
        'median_ci_high': float(np.percentile(np.median(draws, axis=1), 97.5)),
    }

def metric_row(values, group, profile, metric, n_boot, rng):
    values = pd.to_numeric(values, errors='coerce').dropna().to_numpy(float)
    row = {'summary_group': group, 'profile': profile, 'metric': metric, 'n': len(values)}
    if len(values) == 0:
        return row
    row.update({
        'mean': float(values.mean()),
        'sd': float(values.std(ddof=1)) if len(values) > 1 else np.nan,
        'median': float(np.median(values)),
        'iqr': float(np.percentile(values, 75) - np.percentile(values, 25)),
        'p02_5': float(np.percentile(values, 2.5)),
        'p16': float(np.percentile(values, 16)),
        'p84': float(np.percentile(values, 84)),
        'p97_5': float(np.percentile(values, 97.5)),
        'min': float(values.min()),
        'max': float(values.max()),
    })
    row.update(bootstrap(values, n_boot, rng))
    if metric == 'slope':
        row.update({
            'fraction_negative': float(np.mean(values < 0)),
            'fraction_zero': float(np.mean(values == 0)),
            'fraction_positive': float(np.mean(values > 0)),
        })
    return row

def run(input_csv, output_dir, n_boot, seed):
    fits = pd.read_csv(input_csv)
    output_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    rows = []
    quality = []
    for group in GROUPS:
        subset = fits if group == 'ALL' else fits[fits['sample_group'] == group]
        quality.append({
            'summary_group': group,
            'n_groups': len(subset),
            'n_constant_success': int(subset['constant_success'].fillna(False).sum()),
            'n_linear_success': int(subset['linear_success'].fillna(False).sum()),
            'n_log_linear_success': int(subset['log_linear_success'].fillna(False).sum()),
            'n_excluded_rows': int(subset['n_excluded'].fillna(0).sum()),
            'n_duplicate_radius_groups': int(subset['duplicate_radius_flag'].fillna(False).sum()),
            'n_profile_count_mismatch_groups': int(subset['profile_count_mismatch_flag'].fillna(False).sum()),
        })
        for profile, spec in SPECS.items():
            success = subset[spec['success']].fillna(False).astype(bool)
            rows.append(metric_row(subset.loc[success, spec['normalization']], group, profile, 'normalization_at_Rp', n_boot, rng))
            if spec['intercept'] is not None:
                rows.append(metric_row(subset.loc[success, spec['intercept']], group, profile, 'intercept_R0', n_boot, rng))
            if spec['slope'] is not None:
                rows.append(metric_row(subset.loc[success, spec['slope']], group, profile, 'slope', n_boot, rng))
    pd.DataFrame(rows).to_csv(output_dir / 'ml_profile_statistics.csv', index=False)
    pd.DataFrame(quality).to_csv(output_dir / 'ml_profile_statistics_quality.csv', index=False)
    metadata = {
        'input_csv': str(input_csv.resolve()),
        'input_sha256': hashlib.sha256(input_csv.read_bytes()).hexdigest(),
        'n_bootstrap': n_boot,
        'seed': seed,
        'weighting': 'equal galaxy/sample-group weighting; one fitted row per group',
        'groups': GROUPS,
        'cosmos_definition': 'COSMOS_DR1 and COSMOS_kin combined in sample_group=COSMOS',
        'bootstrap': 'resample fitted groups with replacement; percentile 95% intervals',
        'normalization_definition': 'constant_c, or fitted c at mean radius Rp for linear/log-linear profiles',
        'intercept_definition': 'linear_intercept_c0 and log_linear_intercept_c0 are the fitted R=0 intercepts in the a+bR convention',
        'slope_units': {'linear': 'native M/L per kpc', 'log_linear': 'inverse kpc'},
    }
    (output_dir / 'ml_profile_statistics_metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, default=Path('outputs/ml_profile_fits/ml_profile_fits.csv'))
    parser.add_argument('--output-dir', type=Path, default=Path('outputs/ml_profile_stats'))
    parser.add_argument('--n-bootstrap', type=int, default=10000)
    parser.add_argument('--seed', type=int, default=20260922)
    args = parser.parse_args()
    run(args.input, args.output_dir, args.n_bootstrap, args.seed)
