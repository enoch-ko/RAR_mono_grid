#!/usr/bin/env python3
"""Plot saved GP-search frontiers; never rerun the probability/search calculation.

Three panels: signed acceleration components, fixed M/L profiles, and RAR.
RAR intervals propagate Gaussian velocity errors exactly through V**2/R and
sample the declared independent lognormal gas/bulge factors. They are pointwise
68.27% intervals, not GP-profile posteriors or simultaneous confidence bands.
"""
import argparse
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import foldnorm

import inverse_rar_diagnostic as inverse
import mono_SPARC_RAR as mono
from utils import DEFAULT_DENS_ROOT, DEFAULT_TABLE_PATH, vel2acc

Q = np.array([0.15865525393145707, 0.8413447460685429])


def observed_interval(v, e, r):
    """Exact marginal interval of squared Gaussian velocity, including both signs."""
    speeds = foldnorm.ppf(Q[:, None], np.abs(v) / e, scale=e)
    return speeds**2 * vel2acc(np.ones_like(r), r)


def nuisance_draws(settings, draws, seed):
    rng = np.random.default_rng(seed)
    z = rng.standard_normal((2, draws))
    gas = settings["gas_ml"] * 10.0**(settings["gas_ml_dex"] * z[0])
    bulge = settings["bulge_ml"] * 10.0**(settings["bulge_ml_dex"] * z[1])
    return gas, bulge


def baryonic_interval(disc, gas_coefficient, bulge_coefficient, factors):
    gas, bulge = factors
    values = (disc[None, :] + gas[:, None] * gas_coefficient[None, :]
              + bulge[:, None] * bulge_coefficient[None, :])
    return np.quantile(values, Q, axis=0), np.mean(values <= 0, axis=0)


def check_provenance(data, data_root, dens_root):
    recorded = data["provenance"]
    root = Path(mono.__file__).parent
    for name, expected in recorded["source_sha256"].items():
        if inverse.sha256(root / name) != expected:
            raise ValueError(f"source mismatch: {name}")
    script_file = recorded.get("script_file", "gp_rar_search.py")
    if script_file not in ("gp_rar_search.py", "gp_rar_search_retained.py"):
        raise ValueError("unsupported search script provenance")
    if inverse.sha256(root / script_file) != recorded["script_sha256"]:
        raise ValueError("search script differs from saved result")
    paths = {"rotmod": data_root / f'{data["galaxy"]}_rotmod.dat',
             "dens": dens_root / f'{data["galaxy"]}.dens',
             "table": Path(DEFAULT_TABLE_PATH)}
    for name, path in paths.items():
        if inverse.sha256(path) != recorded["input_sha256"][name]:
            raise ValueError(f"input mismatch: {path}")


def direct_disc(rad, brightness, r, zd, knots, xi):
    density = brightness * (inverse.hat_matrix(rad, knots) @ xi)
    force, seen, error = mono.solve_disk_basis(rad, density, zd, r)
    if error or seen or force is None or not np.all(np.isfinite(force)):
        raise ValueError(f"plot force solve rejected: {error}; {seen}")
    return vel2acc(force, r)


def targets(path, data, output_root):
    # Include every non-baseline frontier entry, not just the highest probability.
    pool = data.get("retained_profiles", data["frontier"])
    entries = [f for f in pool if f["record_index"] != 0]
    if not entries:
        entries = [f for f in pool if f["record_index"] == 0]
    folder = output_root / path.parent.name
    return [(entry, folder / f'{data["galaxy"]}_frontier_{entry["record_index"]:06d}.png')
            for entry in entries]


def plot_result(path, data, entry, target, args, loaded, force_cache):
    galaxy = data["galaxy"]
    settings = data["settings"]
    check_provenance(data, args.data_root, args.dens_root)
    if galaxy not in loaded:
        loaded[galaxy] = inverse.galaxy_data(galaxy, args.data_root, args.dens_root)
    rad, brightness, r, v, e, vg, vb, zd, excluded = loaded[galaxy]
    knots = np.asarray(data["knots_kpc"], float)
    if not np.array_equal(r, np.asarray(data["retained_radii_kpc"])):
        raise ValueError("retained radii differ from saved result")
    xi0 = np.full(len(knots), 0.5)
    xi1 = np.asarray(entry["xi"], float)
    discs = []
    for xi in (xi0, xi1):
        key = (galaxy, knots.tobytes(), xi.tobytes())
        if key not in force_cache:
            force_cache[key] = direct_disc(rad, brightness, r, zd, knots, xi)
        discs.append(force_cache[key])
    gas = vel2acc(vg * np.abs(vg), r)
    bulge = vel2acc(vb**2, r)
    gobs = vel2acc(v**2, r)
    yinterval = observed_interval(v, e, r)
    seed = args.seed + sum(map(ord, galaxy))
    factors = nuisance_draws(settings, args.draws, seed)
    centers = [disk + settings["gas_ml"] * gas + settings["bulge_ml"] * bulge
               for disk in discs]
    intervals, nonpositive = zip(*(baryonic_interval(d, gas, bulge, factors)
                                  for d in discs))
    record = data["records"][entry["record_index"]]
    baseline = data["records"][0]
    if not record["valid"] or not baseline["valid"]:
        raise ValueError("cannot plot invalid probability as a fitted result")
    fig = plt.figure(figsize=(13, 8), layout="constrained")
    grid = fig.add_gridspec(2, 2, width_ratios=(1, 1.15))
    acc = fig.add_subplot(grid[0, 0])
    ml = fig.add_subplot(grid[1, 0], sharex=acc)
    rar = fig.add_subplot(grid[:, 1])
    colors = ("#0077bb", "#cc3311")
    names = ("Initial constant 0.5", "Selected frontier")
    for j, (disk, center, band, color, name) in enumerate(
            zip(discs, centers, intervals, colors, names)):
        acc.plot(r, disk, "--", color=color, alpha=.65, label=f"Disc: {name}")
        acc.plot(r, center, color=color, label=f"Total: {name}")
        acc.fill_between(r, *band, color=color, alpha=.12)
        ml.plot(knots, xi0 if j == 0 else xi1, color=color, label=name)
        use = ((center > 0) & (gobs > 0) & (band[0] > 0)
               & (yinterval[0] > 0))
        rar.plot(center[use], gobs[use], "o-", ms=3, color=color, alpha=.75,
                 label=name)
        rar.hlines(gobs[use], band[0, use], band[1, use], color=color, alpha=.28)
        rar.vlines(center[use], yinterval[0, use], yinterval[1, use],
                   color=color, alpha=.28)
    acc.plot(r, settings["gas_ml"] * gas, color="#228833", label="Gas")
    acc.plot(r, settings["bulge_ml"] * bulge, color="#aa3377", label="Bulge")
    acc.plot(r, gobs, "k.", label="Observed")
    acc.vlines(r, *yinterval, color="k", alpha=.3)
    acc.set_yscale("symlog", linthresh=1e-12)
    acc.set_ylabel("Signed acceleration [m s$^{-2}$]")
    acc.legend(fontsize=7, ncol=2)
    acc.tick_params(labelbottom=False)
    ml.set_xlabel("Radius [kpc]")
    ml.set_ylabel("Disc M/L [solar units]")
    ml.legend(fontsize=8)
    rar.set(xscale="log", yscale="log", xlabel="$g_{bar}$ [m s$^{-2}$]",
            ylabel="$g_{obs}$ [m s$^{-2}$]")
    rar.legend(fontsize=8, loc="lower right")
    annotation = (f'Initial: P={math.exp(baseline["marginal_logp"]):.6g}, '
                  f'ln P={baseline["marginal_logp"]:.5f}, J=0\n'
                  f'Final: P={math.exp(record["marginal_logp"]):.6g}, '
                  f'ln P={record["marginal_logp"]:.5f}, J={record["J_ref"]:.6g}\n'
                  'P marginalized over gas/bulge; profile held fixed.\n'
                  'Bars: pointwise 68.27% propagated intervals.\n'
                  'No disc-profile, distance or inclination uncertainty.')
    rar.text(.03, .97, annotation, transform=rar.transAxes, va="top", fontsize=8,
             bbox=dict(facecolor="white", alpha=.88, edgecolor=".8"))
    fig.suptitle(f'{galaxy}: {settings["directions_per_ell"]} directions/lengthscale, '
                 f'{settings["amplitude_steps"]} amplitude positions\n'
                 f'Selected record {entry["record_index"]}: {record["source"]}; '
                 'finite-search diagnostic, not a global optimum', fontsize=11)
    for axis in (acc, ml, rar):
        axis.grid(alpha=.15)
    metadata = {"source_result": str(path.resolve()),
                "source_result_sha256": inverse.sha256(path),
                "plot_script_sha256": inverse.sha256(__file__),
                "search_provenance": data["provenance"], "settings": settings,
                "record_index": entry["record_index"], "record": record,
                "initial_record": baseline, "radii_kpc": r.tolist(),
                "initial_disc_ml": xi0.tolist(), "final_disc_ml": xi1.tolist(),
                "observed_interval": yinterval.tolist(),
                "initial_gbar_interval": intervals[0].tolist(),
                "final_gbar_interval": intervals[1].tolist(),
                "nonpositive_gbar_draw_fractions": [x.tolist() for x in nonpositive],
                "uncertainty_method": "Exact squared-Gaussian velocity quantiles; independent global gas/bulge lognormal Monte Carlo factors, shared across radii and initial/final profiles",
                "interval_quantiles": Q.tolist(), "mc_draws": args.draws,
                "mc_seed": seed, "distance_inclination_fixed": True,
                "disc_profile_uncertainty_included": False,
                "probabilities_recomputed": False, "excluded_kinematic_rows": excluded}
    check_provenance(data, args.data_root, args.dens_root)
    target.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(target, dpi=160)
    plt.close(fig)
    with target.with_suffix(".plotmeta.json").open("x") as stream:
        json.dump(metadata, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(target, flush=True)


def self_test():
    r = np.array([1., 2.])
    v = np.array([0., 20.])
    e = np.array([2., 1.])
    out = observed_interval(v, e, r)
    speeds = np.sqrt(out / vel2acc(np.ones_like(r), r))
    np.testing.assert_allclose(foldnorm.cdf(speeds, np.abs(v) / e, scale=e),
                               np.broadcast_to(Q[:, None], (2, 2)), atol=1e-10)
    settings = dict(gas_ml=1., bulge_ml=.7, gas_ml_dex=0., bulge_ml_dex=0.)
    factors = nuisance_draws(settings, 100, 3)
    bands, negative = baryonic_interval(np.ones(2), r, r, factors)
    np.testing.assert_allclose(bands, np.tile(1 + 1.7 * r, (2, 1)))
    assert not np.any(negative)
    np.testing.assert_array_equal(factors[0], nuisance_draws(settings, 100, 3)[0])
    print("Three uncertainty checks passed; no force/probability calculations.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", nargs="*", type=Path)
    parser.add_argument("--output-root", type=Path,
                        default=Path(mono.__file__).parent / "outputs/gp_rar_diagnostics")
    parser.add_argument("--data-root", type=Path, default=Path(mono.DEFAULT_DATA_ROOT))
    parser.add_argument("--dens-root", type=Path, default=Path(DEFAULT_DENS_ROOT))
    parser.add_argument("--draws", type=int, default=20000)
    parser.add_argument("--seed", type=int, default=20261013)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    if not args.results or args.draws < 1000:
        parser.error("provide saved results and at least 1000 uncertainty draws")
    jobs = []
    for path in args.results:
        data = json.loads(path.read_text())
        if data["status"] != "finite_search_completed":
            raise ValueError(f"incomplete result: {path}")
        jobs.extend((path, data, entry, target)
                    for entry, target in targets(path, data, args.output_root))
    destinations = [p for _, _, _, t in jobs for p in (t, t.with_suffix(".plotmeta.json"))]
    if len(set(destinations)) != len(destinations) or any(p.exists() for p in destinations):
        raise FileExistsError("duplicate or existing plot destination; use a new output root")
    loaded, force_cache = {}, {}
    for path, data, entry, target in jobs:
        plot_result(path, data, entry, target, args, loaded, force_cache)


if __name__ == "__main__":
    main()
