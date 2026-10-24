"""Compare zero-start Casertano disk curves to raw SPARC rotmod data at disk M/L=1."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import mono_SPARC_RAR as rar

ML_DISK = 1.0
EPSILON_FRACTIONS = (1e-3, 1e-6, 1e-9)
ROT_COLUMNS = ['Rad', 'Vobs', 'errV', 'Vgas', 'Vdisk', 'Vbul', 'SBdisk', 'SBbul']
OUT_DIR = Path(__file__).resolve().parent / 'outputs' / 'zero_radius_rotmod_ml1'


def compare_case(radius, sb, z0, rot):
    velocity = rar.vcdisk(radius, ML_DISK * sb, z0=z0, rhoz='exp')
    term = rar.signed_disk_term(velocity)
    keep = np.isfinite(radius) & (radius > 0) & np.isfinite(term)
    if np.count_nonzero(keep) < 2:
        raise ValueError('fewer than two finite positive-radius solver points')
    use = (np.isfinite(rot.Rad) & np.isfinite(rot.Vdisk)
           & (rot.Rad >= np.min(radius[keep])) & (rot.Rad <= np.max(radius[keep])))
    r_obs = rot.Rad.to_numpy(float)[use]
    v_sparc = rot.Vdisk.to_numpy(float)[use]
    if r_obs.size == 0:
        raise ValueError('no rotmod radii in solver coverage')
    term_at_obs = rar.interpolate_disk_force(radius[keep], term[keep], r_obs)
    v_model = np.sign(term_at_obs) * np.sqrt(np.abs(term_at_obs))
    delta = v_model - v_sparc
    mse = float(np.mean(delta**2))
    return (len(delta), mse, float(np.sqrt(mse)), float(np.mean(np.abs(delta))),
            float(100*np.sqrt(np.mean((delta/v_sparc)**2))), float(np.max(np.abs(delta))))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--galaxy', action='append', help='limit to one or more SPARC names')
    parser.add_argument('--output-dir', type=Path, default=OUT_DIR)
    args = parser.parse_args()
    table = rar.load_sparc_table()
    selected = set(args.galaxy) if args.galaxy else None
    rows, errors = [], []
    for _, info in table.iterrows():
        galaxy = str(info.Galaxy)
        if selected and galaxy not in selected:
            continue
        try:
            radius, sb = rar.get_SBdisk(galaxy, table=table)
            if len(radius) < 2 or not np.isclose(radius[0], 0.0, atol=1e-12):
                continue
            rot_path = rar.DEFAULT_DATA_ROOT / f'{galaxy}_rotmod.dat'
            rot = pd.read_csv(rot_path, sep=r'\s+', comment='#', header=None, names=ROT_COLUMNS)
            z0 = 0.196 * float(info.Rdisk)**0.633
            cases = [('deleted', radius[1:], sb[1:])]
            cases += [(f'eps={f:g}*R1', np.r_[f*radius[1], radius[1:]],
                       np.r_[sb[0], sb[1:]]) for f in EPSILON_FRACTIONS]
            for label, r_case, sb_case in cases:
                rows.append((galaxy, label, *compare_case(r_case, sb_case, z0, rot)))
        except Exception as exc:
            errors.append((galaxy, f'{type(exc).__name__}: {exc}'))
    columns = ['galaxy', 'case', 'n', 'MSE_km2_s-2', 'RMSE_km_s', 'MAE_km_s', 'RMS_fraction_pct', 'max_abs_km_s']
    args.output_dir.mkdir(parents=True, exist_ok=True)
    results = pd.DataFrame(rows, columns=columns)
    results.to_csv(args.output_dir / 'zero_radius_rotmod_ml1_results.csv', index=False)
    print(f'zero-start galaxies compared: {results.galaxy.nunique()}; errors: {len(errors)}')
    print('case n_galaxies median_MSE_km2_s-2 median_RMSE_km_s median_MAE_km_s median_RMS_fraction_pct')
    for label in dict.fromkeys(results.case):
        group = results[results.case == label]
        print(f'{label} {len(group)} {group["MSE_km2_s-2"].median():.6g} '
              f'{group.RMSE_km_s.median():.6g} {group.MAE_km_s.median():.6g} {group.RMS_fraction_pct.median():.6g}')
    for galaxy, error in errors:
        print(f'ERROR {galaxy}: {error}')


if __name__ == '__main__':
    main()
