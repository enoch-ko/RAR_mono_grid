#!/usr/bin/env python3
"""Plot the RAR, velocity, and acceleration components for a SPARC galaxy.

Run from the project directory with:
    python3 plot_RAR_components.py --galaxy UGC07399
Use --help for profile, probability, marginalisation, and path options.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from vcdisk import vcdisk

from utils import (
    compute_A_B,
    get_SBdisk,
    load_sparc_table,
    log_monotonic_probability_for_gbar,
    log_weighted_average,
    lognormal_quadrature,
    vel2acc,
)

ROOT = Path(__file__).resolve().parent
DATA_ROOT = Path("/mnt/users/koe/data")
DENS_ROOT = Path("/mnt/users/koe/SPARC_RAR/BulgeDiskDec_LTG")
TABLE_PATH = ROOT / "SPARC_table.txt"


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--galaxy", default="UGC07399",
                        help="SPARC galaxy name (must occur in SPARC_table.txt)")
    parser.add_argument("--profile", choices=("constant", "linear", "both"), default="both")
    parser.add_argument("--prob-grid-size", type=int, default=800)
    parser.add_argument("--ml-samples", type=int, default=11,
                        help="Gauss-Hermite nodes per gas/bulge M/L term")
    parser.add_argument("--gas-ml", type=float, default=1.0)
    parser.add_argument("--bulge-ml", type=float, default=0.7)
    parser.add_argument("--gas-ml-dex", type=float, default=0.04)
    parser.add_argument("--bulge-ml-dex", type=float, default=0.1)
    parser.add_argument("--output-root", type=Path, default=ROOT / "outputs")
    parser.add_argument("--data-root", type=Path, default=DATA_ROOT)
    parser.add_argument("--dens-root", type=Path, default=DENS_ROOT)
    parser.add_argument("--table-path", type=Path, default=TABLE_PATH)
    parser.add_argument("--tie-atol", type=float, default=0.0)
    parser.add_argument("--zpadding", type=float, default=8.0)
    parser.add_argument("--max-tie-permutations", type=int, default=720)
    return parser.parse_args()


def probability(gbar, obs, err, radius, args):
    return float(log_monotonic_probability_for_gbar(
        gbar, obs, err, radius, ngrid=args.prob_grid_size,
        zpadding=args.zpadding, tie_atol=args.tie_atol,
        max_tie_permutations=args.max_tie_permutations,
    ))


def profile_components(galaxy, profile, setting, data, table, args):
    radius = data["Rad"].to_numpy(float)
    vobs = data["Vobs"].to_numpy(float)
    err = data["errV"].to_numpy(float)
    vgas = data["Vgas"].to_numpy(float)
    vbul = data["Vbul"].to_numpy(float)

    if profile == "constant":
        vdisk = data["Vdisk"].to_numpy(float) * np.sqrt(setting)
    else:
        match = table.loc[table["Galaxy"].astype(str) == galaxy]
        if len(match) != 1:
            raise ValueError(f"Expected one SPARC table row for {galaxy}; found {len(match)}")
        rdisk = float(match.iloc[0]["Rdisk"])
        z0 = 0.196 * rdisk**0.633
        r_dens, sb_disk = get_SBdisk(galaxy, dens_root=args.dens_root, table=table)
        a, b = compute_A_B(r_dens, sb_disk)
        # setting is the linear slope; normalise so the luminosity-weighted mean is 0.5.
        ml_profile = 0.5 + setting * (r_dens - b / a)
        if np.any(~np.isfinite(ml_profile) | (ml_profile <= 0.0)):
            raise ValueError(f"Linear M/L profile is not positive for slope {setting}")
        vdisk_density = vcdisk(r_dens, sb_disk * ml_profile, z0=z0, rhoz="exp")
        vdisk = np.interp(radius, r_dens, vdisk_density)

    disk_term = vdisk**2
    gas_term = vgas * np.abs(vgas) * args.gas_ml
    bulge_term = vbul**2 * args.bulge_ml
    vbar_squared = disk_term + gas_term + bulge_term
    if np.any(~np.isfinite(vbar_squared) | (vbar_squared < 0.0)):
        raise ValueError("Non-finite or negative total baryonic velocity squared")
    vbar = np.sqrt(vbar_squared)

    gdisk = vel2acc(disk_term, radius)
    ggas = vel2acc(gas_term, radius)
    gbul = vel2acc(bulge_term, radius)
    gbar = gdisk + ggas + gbul
    gobs = vel2acc(vobs**2, radius)

    conditional = probability(gbar, vobs, err, radius, args)
    gas_nodes, gas_weights = lognormal_quadrature(
        args.gas_ml, args.gas_ml_dex, args.ml_samples)
    bulge_nodes, bulge_weights = lognormal_quadrature(
        args.bulge_ml, args.bulge_ml_dex, args.ml_samples)
    logps, weights = [], []
    for gas_ml, gas_weight in zip(gas_nodes, gas_weights):
        for bulge_ml, bulge_weight in zip(bulge_nodes, bulge_weights):
            marginal_gbar = (
                gdisk
                + vel2acc(vgas * np.abs(vgas) * gas_ml, radius)
                + vel2acc(vbul**2 * bulge_ml, radius)
            )
            logps.append(probability(marginal_gbar, vobs, err, radius, args))
            weights.append(gas_weight * bulge_weight)
    marginal = float(log_weighted_average(logps, weights))

    return {
        "radius": radius, "vobs": vobs, "err": err, "vbar": vbar,
        "vdisk": vdisk, "vgas": vgas, "gobs": gobs, "gbar": gbar,
        "gdisk": gdisk, "ggas": ggas, "gbul": gbul,
        "conditional": conditional, "marginal": marginal,
    }


def make_figure(galaxy, profile, settings, labels, data, table, args):
    colors = {"obs": "black", "disk": "tab:blue", "gas": "tab:green", "bar": "tab:red"}
    fig, axes = plt.subplots(3, 3, figsize=(15, 12), sharex="col")
    fig.suptitle(f"{galaxy}: {profile} stellar-disk M/L profiles", y=0.995)

    for col, (setting, label) in enumerate(zip(settings, labels)):
        c = profile_components(galaxy, profile, setting, data, table, args)
        radius = c["radius"]
        order = np.argsort(radius)
        rr = radius[order]
        axes[0, col].set_title(label)

        ax = axes[0, col]
        ax.scatter(c["gbar"], c["gobs"], color=colors["obs"], s=22)
        low = min(float(np.min(c["gbar"])), float(np.min(c["gobs"])))
        high = max(float(np.max(c["gbar"])), float(np.max(c["gobs"])))
        ax.plot([low, high], [low, high], "--", color="0.55", lw=1)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_ylabel(r"$g_{\rm obs}$ [m s$^{-2}$]" if col == 0 else "")
        ax.set_xlabel(r"$g_{\rm bar}$ [m s$^{-2}$]")
        ax.text(0.03, 0.04,
                f"Conditional log P={c['conditional']:.3f}\nGas+bulge marginal log P={c['marginal']:.3f}",
                transform=ax.transAxes, fontsize=7, va="bottom",
                bbox={"facecolor": "white", "alpha": 0.82, "edgecolor": "0.7"})

        ax = axes[1, col]
        ax.errorbar(rr, c["vobs"][order], yerr=c["err"][order], color=colors["obs"],
                    marker="o", ms=2.5, lw=1, label=r"$V_{\rm obs}$")
        ax.plot(rr, c["vdisk"][order], color=colors["disk"], marker=".", label=r"$V_{\rm disk}$")
        ax.plot(rr, c["vgas"][order], color=colors["gas"], marker=".", label=r"$V_{\rm gas}$")
        ax.plot(rr, c["vbar"][order], color=colors["bar"], marker=".", label=r"$V_{\rm bar}$")
        ax.set_ylabel("Velocity [km/s]" if col == 0 else "")
        ax.set_xlabel("Radius [kpc]")

        ax = axes[2, col]
        ax.plot(rr, c["gobs"][order], color=colors["obs"], marker=".", label=r"$g_{\rm obs}$")
        ax.plot(rr, c["gdisk"][order], color=colors["disk"], marker=".", label=r"$g_{\rm disk}$")
        ax.plot(rr, c["ggas"][order], color=colors["gas"], marker=".", label=r"$g_{\rm gas}$")
        ax.plot(rr, c["gbar"][order], color=colors["bar"], marker=".", label=r"$g_{\rm bar}$")
        ax.set_yscale("log")
        ax.set_ylabel(r"Acceleration [m s$^{-2}$]" if col == 0 else "")
        ax.set_xlabel("Radius [kpc]")

    handles, legend_labels = axes[1, 0].get_legend_handles_labels()
    fig.legend(handles, legend_labels, loc="lower center", ncol=4, frameon=False)
    fig.tight_layout(rect=(0, 0.04, 1, 0.97))
    folder = args.output_root / f"{galaxy.lower()}_{profile}_rar_plots"
    folder.mkdir(parents=True, exist_ok=True)
    output = folder / f"{galaxy}_{profile}_RAR_components.png"
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {output}")


def main():
    args = parse_args()
    table = load_sparc_table(args.table_path)
    if not (table["Galaxy"].astype(str) == args.galaxy).any():
        raise ValueError(f"Galaxy {args.galaxy!r} is not in the SPARC table")
    columns = ["Rad", "Vobs", "errV", "Vgas", "Vdisk", "Vbul", "SBdisk", "SBbul"]
    raw = np.loadtxt(args.data_root / f"{args.galaxy}_rotmod.dat")
    data = pd.DataFrame(raw, columns=columns)

    if args.profile in ("constant", "both"):
        settings = [0.05, 0.5, 5.0]
        labels = [rf"Constant $\Upsilon_{{disk}}={x:g}$" for x in settings]
        make_figure(args.galaxy, "constant", settings, labels, data, table, args)
    if args.profile in ("linear", "both"):
        settings = [-0.1, 0.0, 0.2]
        labels = [rf"Linear $s={x:+.2f}$ kpc$^{{-1}}$; mean=0.5" for x in settings]
        make_figure(args.galaxy, "linear", settings, labels, data, table, args)


if __name__ == "__main__":
    main()
