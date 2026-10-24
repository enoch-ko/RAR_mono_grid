#!/usr/bin/env python3
"""Denominator-only broad-limit pilot with an explicit, separate order-receipt preparation mode."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from bounded_uniform_bayes_factor_20261021 import _finite_support_probability, _prior_mc

DEFAULT_GRIDS = (8001, 16001, 32001)
SEED = 1729
DRAWS = 100000


def _sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def _canonical_conditional_row(csv_path, galaxy):
    with Path(csv_path).open(newline='') as f:
        rows = [r for r in csv.DictReader(f)
                if r.get('galaxy') == galaxy and r.get('profile') == 'constant'
                and r.get('probability_kind') == 'conditional'
                and math.isclose(float(r['upsilon_disk_mean']), .5, abs_tol=1e-10)
                and (r.get('disk_slope', '').lower() == 'nan' or r.get('disk_slope', '') == '')]
    if len(rows) != 1 or rows[0].get('valid', '').lower() != 'true':
        raise ValueError(f'expected exactly one valid constant conditional M/L=0.5 row for {galaxy}; found {len(rows)}')
    value = float(rows[0]['log_probability'])
    if not math.isfinite(value):
        raise ValueError('canonical conditional numerator is non-finite')
    return rows[0], value


def _prepare_order_receipt(galaxy, canonical_csv, metadata_path, receipt_path,
                            data_root=None, dens_root=None, table_path=None):
    """Prepare a fixed-theta order receipt. Called only by explicit CLI mode; computes no Q/P."""
    from mono_SPARC_RAR import (DEFAULT_DATA_ROOT, DEFAULT_DENS_ROOT, DEFAULT_TABLE_PATH,
                                solve_disk_basis)
    from utils import get_SBdisk, load_sparc_table, vel2acc

    out = Path(receipt_path)
    if out.exists():
        raise FileExistsError(f'refusing to overwrite receipt: {out}')
    metadata_path = Path(metadata_path)
    canonical_csv = Path(canonical_csv)
    meta = json.loads(metadata_path.read_text())
    row, logp = _canonical_conditional_row(canonical_csv, galaxy)
    if meta.get('conditional_probability') is not True:
        raise ValueError('run metadata does not confirm conditional probability output')
    if not math.isclose(float(meta.get('disk_ml_median', -1)), .5):
        raise ValueError('run metadata disk M/L does not match the pinned baseline')
    if not math.isclose(float(meta.get('gas_ml_median', -1)), 1.0):
        raise ValueError('run metadata gas M/L does not match the pinned central value')

    data_root = Path(data_root) if data_root else Path(DEFAULT_DATA_ROOT)
    dens_root = Path(dens_root) if dens_root else Path(DEFAULT_DENS_ROOT)
    table_path = Path(table_path) if table_path else Path(DEFAULT_TABLE_PATH)
    table = load_sparc_table(table_path)
    matches = table.loc[table['Galaxy'].astype(str) == galaxy]
    if len(matches) != 1:
        raise ValueError(f'expected one SPARC master row for {galaxy}; found {len(matches)}')
    master = matches.iloc[0]
    if float(master.get('Lbul', 0.0)) != 0.0:
        raise ValueError('this pilot receipt is restricted to the verified bulgeless baseline')

    rotmod = data_root / f'{galaxy}_rotmod.dat'
    raw = np.atleast_2d(np.loadtxt(rotmod))
    if raw.shape[1] != 8:
        raise ValueError(f'{rotmod} has {raw.shape[1]} columns; expected the canonical 8-column rotmod layout')
    cols = {name: raw[:, i] for i, name in enumerate(
        ('Rad', 'Vobs', 'errV', 'Vgas', 'Vdisk', 'Vbul', 'SBdisk', 'SBbul'))}
    r_input = cols['Rad']
    rad, sb_disk = get_SBdisk(galaxy, dens_root=dens_root, table=table)
    if len(rad) < 2 or len(rad) != len(sb_disk) or not np.all(np.isfinite(rad)):
        raise ValueError('invalid full .dens disk profile returned by the canonical loader')
    overlap = (
        np.isfinite(r_input) & (r_input >= float(rad[0])) & (r_input <= float(rad[-1]))
        & np.isfinite(cols['Vobs']) & np.isfinite(cols['errV'])
        & np.isfinite(cols['Vgas']) & np.isfinite(cols['Vbul'])
        & (r_input > 0) & (cols['errV'] > 0)
    )
    dens_file = dens_root / f'{galaxy}.dens'
    if not dens_file.is_file():
        raise FileNotFoundError(f'expected canonical full .dens source: {dens_file}')
    if np.count_nonzero(overlap) < 2:
        raise ValueError('fewer than two kinematic radii overlap the full .dens radial range')
    r = r_input[overlap]
    vobs, errv = cols['Vobs'][overlap], cols['errV'][overlap]
    vgas, vbul = cols['Vgas'][overlap], cols['Vbul'][overlap]
    arrays = (r, vobs, errv, vgas, vbul)
    if not all(np.all(np.isfinite(a)) for a in arrays) or np.any(r <= 0) or np.any(errv <= 0):
        raise ValueError('scorer-retained overlapping inputs must be finite with positive radii/errors')
    rdisk = float(master['Rdisk'])
    zd = 0.196 * rdisk ** 0.633
    # Exactly one constant-disk basis solve. The helper receives the full .dens profile
    # and applies the established zero-radius solver treatment internally.
    basis, warnings, error = solve_disk_basis(rad, sb_disk, zd, r)
    if error:
        raise RuntimeError('single baseline force/order solve failed: ' + str(error))
    disk_term = .5 * np.asarray(basis, dtype=float)
    vbar2 = disk_term + vgas * np.abs(vgas) + vbul ** 2 * float(meta.get('bulge_ml_median', .7))
    if vbar2.shape != r.shape or not np.all(np.isfinite(vbar2)):
        raise ValueError('non-finite or misaligned fixed-theta baryonic velocity-squared field')
    gbar = np.asarray(vel2acc(vbar2, r), dtype=float)
    if not np.all(np.isfinite(gbar)) or np.any(gbar <= 0):
        raise ValueError('scorer-retained gbar must be finite and positive')
    order = np.argsort(-gbar, kind='stable')
    sorted_values = gbar[order]
    groups = []
    for i, value in enumerate(sorted_values):
        if not groups or value != sorted_values[groups[-1][0]]:
            groups.append([i])
        else:
            groups[-1].append(i)
    receipt = {
        'galaxy': galaxy, 'probability_kind': 'conditional',
        'saved_logp_conditional': logp,
        'numerator_source': {'csv': str(canonical_csv), 'row': row,
                             'csv_sha256': _sha256(canonical_csv),
                             'run_metadata': str(metadata_path),
                             'run_metadata_sha256': _sha256(metadata_path)},
        'conditioning': {'profile': 'constant', 'disk_ml': .5, 'gas_ml_factor': 1.0,
                         'bulge': 'absent'},
        'radii_kpc': r.tolist(), 'gbar_km2_s2_kpc': gbar.tolist(),
        'vobs_km_s': vobs.tolist(), 'sigma_v_km_s': errv.tolist(),
        'descending_gbar_order_zero_based': order.tolist(),
        'exact_tie_groups_in_sorted_order': [g for g in groups if len(g) > 1],
        'tie_atol': 0.0, 'scorer_retained_mask': overlap.tolist(),
        'scorer_mask_definition': 'finite rotmod radius inside full .dens [Rmin,Rmax]; all retained Vobs, errV, Vgas, Vbul and radius finite; errV>0; r>0',
        'force_solver_call_count': 1,
        'force_source_sha256': _sha256(Path(__file__).with_name('mono_SPARC_RAR.py')),
        'utils_source_sha256': _sha256(Path(__file__).with_name('utils.py')),
        'rotmod_source': str(rotmod), 'rotmod_sha256': _sha256(rotmod),
        'dens_source': str(dens_file), 'dens_sha256': _sha256(dens_file),
        'dens_root': str(dens_root), 'table_source': str(table_path),
        'table_sha256': _sha256(table_path), 'solver_warnings': list(warnings),
        'run_settings': {k: meta.get(k) for k in ('source_sha256', 'data_root', 'dens_root',
                          'disk_ml_median', 'gas_ml_median', 'bulge_ml_median',
                          'conditional_probability', 'gas_bulge_marginalisation')},
        'note': 'Receipt preparation computes only the single constant disk force response and gbar ordering; it does not evaluate or recompute Q/P.'}
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({'prepared_receipt': str(out), 'galaxy': galaxy,
                      'retained_radii': len(r), 'force_solver_calls': 1,
                      'conditional_logP_reused': logp, 'Q_or_P_recomputed': False}, indent=2))


def _load_receipt(path):
    data = json.loads(Path(path).read_text())
    required = ('galaxy', 'radii_kpc', 'gbar_km2_s2_kpc', 'vobs_km_s', 'sigma_v_km_s',
                'saved_logp_conditional', 'numerator_source', 'conditioning', 'probability_kind')
    missing = [key for key in required if key not in data]
    if missing:
        raise ValueError('receipt missing required fields: ' + ', '.join(missing))
    if data['probability_kind'] != 'conditional':
        raise ValueError('only an explicitly identified fixed-theta conditional numerator is accepted')
    c = data['conditioning']
    if c.get('profile') != 'constant' or not math.isclose(float(c.get('disk_ml', -1)), .5):
        raise ValueError('receipt is not pinned to the constant disk_ML=0.5 baseline')
    if not math.isclose(float(c.get('gas_ml_factor', -1)), 1.0) or c.get('bulge') not in ('absent', 'none'):
        raise ValueError('receipt conditioning must be gas factor 1 and bulgeless')
    r = np.asarray(data['radii_kpc'], dtype=float)
    gb = np.asarray(data['gbar_km2_s2_kpc'], dtype=float)
    vo = np.asarray(data['vobs_km_s'], dtype=float)
    ev = np.asarray(data['sigma_v_km_s'], dtype=float)
    n = len(r)
    if n < 2 or any(x.shape != (n,) for x in (gb, vo, ev)):
        raise ValueError('radii, gbar, Vobs and sigma arrays must be aligned 1-D arrays')
    if not all(np.all(np.isfinite(x)) for x in (r, gb, vo, ev)) or np.any(r <= 0) or np.any(ev <= 0):
        raise ValueError('inputs must be finite, with positive radii and velocity errors')
    order = np.argsort(-gb, kind='stable')
    if data.get('descending_gbar_order_zero_based') != order.tolist():
        raise ValueError('saved ordering does not exactly match supplied gbar values')
    lp = float(data['saved_logp_conditional'])
    if not math.isfinite(lp):
        raise ValueError('saved conditional logP must be finite')
    return data, r, gb, vo, ev, lp


def denominator_only(receipt, grid_sizes=DEFAULT_GRIDS, draws=DRAWS, seed=SEED):
    data, r, gb, vo, ev, saved_logp = _load_receipt(receipt)
    recursions = []
    for ngrid in grid_sizes:
        prior = _finite_support_probability(r, gb, vo, ev, 1.0, ngrid, 0.0, posterior=False)
        recursions.append({'ngrid': ngrid, 'prior_log_probability': prior.log_probability,
                           'prior_probability': math.exp(prior.log_probability)})
    prior_logp = recursions[-1]['prior_log_probability']
    mc = _prior_mc(r, gb, 1.0, draws, 0.0, seed, prior_logp)
    return {'scope': 'denominator-only conditional broad-limit ratio; no force, Q, or numerator recomputation',
            'galaxy': data['galaxy'], 'numerator_source': data['numerator_source'],
            'saved_logp_conditional': saved_logp, 'prior_velocity_support_unit_interval': 'Uniform(-1,1)',
            'prior_log_probability': prior_logp, 'log_bayes_factor_using_saved_numerator': saved_logp-prior_logp,
            'bayes_factor_using_saved_numerator': math.exp(saved_logp-prior_logp),
            'recursion_grid_refinement': recursions, 'prior_mc_crosscheck': mc.__dict__,
            'seed': seed, 'mc_draws': draws, 'tie_atol': 0.0,
            'conditioning': data['conditioning'],
            'notes': ['Common velocity scaling leaves the ordering denominator unchanged.',
                      'Gaussian Q is not recalculated; numerator compatibility must be established before interpreting the ratio.',
                      'Recursion refinement and prior-MC uncertainty are separate checks.',
                      'This pilot is conditional and is not manuscript evidence.']}


def mocktest():
    r = np.array([1., 2., 4.]); gb = np.array([3., 2., 1.])
    vo = np.array([20., 15., 9.]); ev = np.array([2., 2., 1.])
    a = _finite_support_probability(r, gb, vo, ev, 1., 1001, 0., posterior=False)
    b = _finite_support_probability(r, gb, vo, ev, 7., 1001, 0., posterior=False)
    if not math.isclose(a.log_probability, b.log_probability, rel_tol=0., abs_tol=1e-10):
        raise AssertionError('common Vmax scaling changed prior ordering probability')
    m1 = _prior_mc(r, gb, 1., 5000, 0., 1729, a.log_probability)
    m2 = _prior_mc(r, gb, 1., 5000, 0., 1729, a.log_probability)
    if m1.__dict__ != m2.__dict__:
        raise AssertionError('seeded synthetic MC is not reproducible')
    print('MOCKTEST PASS: synthetic scale invariance and seeded MC replay; no galaxy inputs used')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    mode = p.add_mutually_exclusive_group(required=True)
    mode.add_argument('--plan', action='store_true')
    mode.add_argument('--prepare-order-receipt', action='store_true')
    mode.add_argument('--run-denominator-only', action='store_true')
    mode.add_argument('--mocktest', action='store_true')
    p.add_argument('--galaxy', default='NGC4068')
    p.add_argument('--canonical-csv', type=Path, default=Path('outputs/rerun_boundary_20261006/constant/profile_grid.shard_02.csv'))
    p.add_argument('--run-metadata', type=Path, default=Path('outputs/rerun_boundary_20261006/constant/run_metadata.shard_02.json'))
    p.add_argument('--receipt-output', type=Path)
    p.add_argument('--receipt', type=Path)
    p.add_argument('--output', type=Path)
    p.add_argument('--data-root', type=Path)
    p.add_argument('--dens-root', type=Path)
    p.add_argument('--table-path', type=Path)
    p.add_argument('--seed', type=int, default=SEED)
    p.add_argument('--draws', type=int, default=DRAWS)
    p.add_argument('--grid-sizes', type=int, nargs='+', default=list(DEFAULT_GRIDS))
    a = p.parse_args()
    if a.mocktest:
        mocktest(); return
    if a.plan:
        print(json.dumps({'mode':'one-force receipt preparation, then denominator-only ratio',
            'prepare_command':'timeout 300s python3 -B bayes_factor_denominator_pilot_20261021.py --prepare-order-receipt --galaxy NGC4068 --canonical-csv outputs/rerun_boundary_20261006/constant/profile_grid.shard_02.csv --run-metadata outputs/rerun_boundary_20261006/constant/run_metadata.shard_02.json --receipt-output outputs/bayes_factor_pilot_20261021/NGC4068_receipt.json',
            'ratio_command':'timeout 300s python3 -B bayes_factor_denominator_pilot_20261021.py --run-denominator-only --receipt outputs/bayes_factor_pilot_20261021/NGC4068_receipt.json --output outputs/bayes_factor_pilot_20261021/NGC4068_denominator.json --seed 1729 --draws 100000 --grid-sizes 8001 16001 32001',
            'prepare_calls_one_constant_disk_force_solver_only':True,'ratio_reuses_saved_conditional_logP_and_order':True,
            'ratio_recomputes_Q_or_numerator':False,'no_real_calculation_executed':True}, indent=2)); return
    if a.prepare_order_receipt:
        if a.receipt_output is None:
            p.error('--prepare-order-receipt requires --receipt-output')
        _prepare_order_receipt(a.galaxy, a.canonical_csv, a.run_metadata, a.receipt_output,
                               a.data_root, a.dens_root, a.table_path)
        return
    if a.receipt is None or a.output is None:
        p.error('--run-denominator-only requires --receipt and --output')
    if a.output.exists():
        raise FileExistsError('refusing to overwrite ' + str(a.output))
    result = denominator_only(a.receipt, a.grid_sizes, a.draws, a.seed)
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
