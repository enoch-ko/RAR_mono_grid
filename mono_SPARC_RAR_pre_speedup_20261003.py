#!/usr/bin/env python
"""Grid and prior-marginalised SPARC RAR monotonicity analysis."""

import argparse
import csv
import hashlib
import json
import platform
import sys
import warnings
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from vcdisk import vcdisk
from utils import (
    compute_A_B,
    compute_A_beta,
    get_SBdisk,
    load_sparc_table,
    log_monotonic_probability_for_gbar,
    log_weighted_average,
    lognormal_quadrature,
    vel2acc,
)

VALID_PROFILES = ("constant", "linear", "log_linear")
DEFAULT_PROJECT_ROOT = Path("/mnt/users/koe/RAR_mono_grid")
DEFAULT_DATA_ROOT = Path("/mnt/users/koe/data")
DEFAULT_DENS_ROOT = Path("/mnt/users/koe/SPARC_RAR/BulgeDiskDec_LTG")
DEFAULT_TABLE_PATH = DEFAULT_PROJECT_ROOT / "SPARC_table.txt"



def parse_args():
    """Define the reproducible command-line controls for a grid run."""
    parser = argparse.ArgumentParser(description="Grid-search SPARC RAR monotonicity probabilities.")
    parser.add_argument("--profile", choices=VALID_PROFILES, default="constant")
    parser.add_argument("--grid-size", type=int, default=51)
    parser.add_argument("--prob-grid-size", type=int, default=400)
    parser.add_argument("--zpadding", type=float, default=8.0)
    parser.add_argument("--tie-atol", type=float, default=0.0)
    parser.add_argument("--max-tie-permutations", type=int, default=720)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_PROJECT_ROOT / "outputs")
    parser.add_argument("--data-root", type=Path, default=DEFAULT_DATA_ROOT)
    parser.add_argument("--dens-root", type=Path, default=DEFAULT_DENS_ROOT)
    parser.add_argument("--table-path", type=Path, default=DEFAULT_TABLE_PATH)
    parser.add_argument("--gas-ml", type=float, default=1.0,
                        help="Median gas M/L in the SPARC convention.")
    parser.add_argument("--disk-ml", type=float, default=0.5,
                        help="Median luminosity-weighted stellar-disk M/L.")
    parser.add_argument("--bulge-ml", type=float, default=0.7,
                        help="Median bulge M/L.")
    parser.add_argument("--gas-ml-dex", type=float, default=0.04,
                        help="SPARC 1-sigma uncertainty in log10 gas M/L.")
    parser.add_argument("--disk-ml-dex", type=float, default=0.1,
                        help="SPARC 1-sigma uncertainty in log10 disk M/L.")
    parser.add_argument("--bulge-ml-dex", type=float, default=0.1,
                        help="SPARC 1-sigma uncertainty in log10 bulge M/L.")
    parser.add_argument("--disk-grid-half-width-dex", type=float, default=1.0,
                        help="Half-width of the disk grid around its log10 median.")
    parser.add_argument("--slope-min", type=float, default=-5.0,
                        help="Lower bound of the uniform radial-slope prior.")
    parser.add_argument("--slope-max", type=float, default=5.0,
                        help="Upper bound of the uniform radial-slope prior.")
    parser.add_argument("--ml-samples", type=int, default=9,
                        help="Gauss-Hermite nodes per gas or bulge M/L term.")
    parser.add_argument("--marginalise-gas-bulge", action=argparse.BooleanOptionalAction,
                        default=True, help="Write grid probabilities marginalised over gas and bulge M/L.")
    parser.add_argument("--marginalise-disk", action=argparse.BooleanOptionalAction,
                        default=True, help="Integrate the gas/bulge-marginalised grid over the SPARC disk prior.")
    parser.add_argument("--testing", action="store_true", help="Process only the first SPARC galaxy.")
    parser.add_argument("--galaxy", action="append", help="Process only the named galaxy; may be repeated.")
    parser.add_argument("--galaxy-range", nargs=2, type=int, metavar=("START", "END"), help="Process the inclusive 1-based SPARC catalogue row range.")
    parser.add_argument("--shard-id", type=int, help="Write concurrency-safe shard artifacts for later merging.")
    parser.add_argument("--overwrite-shard", action="store_true", help="Replace existing artifacts for this shard ID.")
    args = parser.parse_args()
    if args.galaxy and args.galaxy_range:
        parser.error("--galaxy and --galaxy-range are mutually exclusive")
    if args.galaxy_range is not None:
        start, end = args.galaxy_range
        if start < 1 or end < start:
            parser.error("--galaxy-range requires 1 <= START <= END")
    if args.shard_id is not None and args.shard_id < 1:
        parser.error("shard-id must be a positive integer")
    if args.overwrite_shard and args.shard_id is None:
        parser.error("--overwrite-shard requires --shard-id")
    if args.grid_size < 2 or args.prob_grid_size < 2:
        parser.error("grid-size and prob-grid-size must both be at least 2")
    if args.gas_ml <= 0 or args.disk_ml <= 0 or args.bulge_ml <= 0:
        parser.error("gas, disk, and bulge M/L medians must be positive")
    if args.gas_ml_dex < 0 or args.disk_ml_dex <= 0 or args.bulge_ml_dex < 0:
        parser.error("M/L uncertainties must be non-negative, with disk-ml-dex positive")
    if args.disk_grid_half_width_dex <= 0:
        parser.error("disk-grid-half-width-dex must be positive")
    if args.slope_min >= args.slope_max:
        parser.error("slope-min must be smaller than slope-max")
    if args.ml_samples < 1:
        parser.error("ml-samples must be positive")
    if args.marginalise_disk and not args.marginalise_gas_bulge:
        parser.error("--marginalise-disk requires --marginalise-gas-bulge")
    return args


def source_sha256():
    """Hash this driver so every output records the exact implementation used."""
    digest = hashlib.sha256()
    with open(__file__, "rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def evaluate_probability(g_bar, Vobs, errV, r, args):
    """Evaluate the ordered-acceleration probability and capture non-fatal warnings."""
    caught = []
    try:
        with warnings.catch_warnings(record=True) as records:
            warnings.simplefilter("always")
            log_probability = log_monotonic_probability_for_gbar(
                g_bar, Vobs, errV, r,
                ngrid=args.prob_grid_size,
                zpadding=args.zpadding,
                tie_atol=args.tie_atol,
                max_tie_permutations=args.max_tie_permutations,
            )
            caught = [str(record.message) for record in records]
        if not np.isfinite(log_probability):
            return np.nan, False, "nonfinite_probability", caught
        return float(log_probability), True, "", caught
    except Exception as exc:
        return np.nan, False, f"{type(exc).__name__}: {exc}", caught


def evaluate_baryonic_probability(disk_term, Vgas, Vbul, r, Vobs, errV,
                                  gas_ml, bulge_ml, args):
    """Build g_bar for one pair of gas/bulge M/L values and evaluate it."""
    # Keep the SPARC signed-gas term: Vgas*abs(Vgas), rather than Vgas**2.
    Vbar_squared = disk_term + Vgas * np.abs(Vgas) * gas_ml + Vbul**2 * bulge_ml
    if not np.all(np.isfinite(Vbar_squared)):
        return np.nan, False, "nonfinite_baryonic_acceleration", []
    g_bar = vel2acc(Vbar_squared, r)
    return evaluate_probability(g_bar, Vobs, errV, r, args)


def evaluate_disk_term(disk_term, Vgas, Vbul, r, Vobs, errV, args,
                       gas_nodes, gas_weights, bulge_nodes, bulge_weights):
    """Return conditional and gas/bulge-marginalised probabilities for one disk model."""
    # The conditional result uses the SPARC central gas and bulge values.
    conditional = evaluate_baryonic_probability(
        disk_term, Vgas, Vbul, r, Vobs, errV, args.gas_ml, args.bulge_ml, args)
    if not args.marginalise_gas_bulge:
        return conditional, None

    # The nuisance-parameter integral is an expectation in probability space.
    # Work in log probability to avoid numerical underflow for long rotation curves.
    log_values = []
    weights = []
    marginal_warnings = []
    for gas_ml, gas_weight in zip(gas_nodes, gas_weights):
        for bulge_ml, bulge_weight in zip(bulge_nodes, bulge_weights):
            result = evaluate_baryonic_probability(
                disk_term, Vgas, Vbul, r, Vobs, errV, gas_ml, bulge_ml, args)
            logp, valid, reason, warnings_for_sample = result
            marginal_warnings.extend(warnings_for_sample)
            if not valid:
                # Numerical failures are not silently dropped from the integral.
                return conditional, (np.nan, False, f"marginal_sample_invalid: {reason}", marginal_warnings)
            log_values.append(logp)
            weights.append(gas_weight * bulge_weight)
    try:
        marginal_logp = log_weighted_average(log_values, weights)
    except Exception as exc:
        return conditional, (np.nan, False, f"{type(exc).__name__}: {exc}", marginal_warnings)
    return conditional, (marginal_logp, True, "", marginal_warnings)


def append_row(rows, galaxy, profile, grid_i, grid_j, log10_upsilon, upsilon,
               slope, log_probability, valid, reason, warning_messages, args,
               probability_kind):
    """Append one self-describing conditional or marginalised grid result."""
    rows.append({
        "galaxy": galaxy,
        "profile": profile,
        "probability_kind": probability_kind,
        "grid_i": grid_i,
        "grid_j": grid_j,
        "log10_upsilon_disk": log10_upsilon,
        "upsilon_disk_mean": upsilon,
        "disk_slope": slope,
        "log_probability": log_probability,
        "valid": valid,
        "invalid_reason": reason,
        "warnings": " | ".join(warning_messages),

    })


def append_probability_rows(rows, galaxy, profile, grid_i, grid_j, log10_upsilon,
                            upsilon, slope, results, args):
    """Write conditional and, when requested, marginalised versions of one grid cell."""
    kinds = [("conditional", results[0])]
    if results[1] is not None:
        kinds.append(("marginalised_gas_bulge", results[1]))
    for probability_kind, result in kinds:
        append_row(rows, galaxy, profile, grid_i, grid_j, log10_upsilon, upsilon,
                   slope, *result, args, probability_kind)


def profile_grids(args):
    """Construct the mean disk-M/L and slope grids used by both result families."""
    log10_upsilon_sparc = np.log10(args.disk_ml)
    half_width = args.disk_grid_half_width_dex
    log10_upsilon = np.linspace(
        log10_upsilon_sparc - half_width,
        log10_upsilon_sparc + half_width,
        args.grid_size,
    )
    upsilon = 10.0 ** log10_upsilon
    slopes = np.linspace(args.slope_min, args.slope_max, args.grid_size)
    return log10_upsilon_sparc, log10_upsilon, upsilon, slopes


def trapezoid_weights(coordinates):
    """Return integration-cell widths for a one-dimensional trapezoidal rule."""
    coordinates = np.asarray(coordinates, dtype=float)
    if coordinates.ndim != 1 or coordinates.size < 2:
        raise ValueError("quadrature coordinates must be a one-dimensional grid with at least two points")
    differences = np.diff(coordinates)
    if np.any(~np.isfinite(coordinates)) or np.any(differences <= 0):
        raise ValueError("quadrature coordinates must be finite and strictly increasing")
    weights = np.empty_like(coordinates)
    weights[0] = differences[0] / 2.0
    weights[-1] = differences[-1] / 2.0
    weights[1:-1] = (coordinates[2:] - coordinates[:-2]) / 2.0
    return weights


def grid_prior_weights(log10_upsilon, slopes, args):
    """Build normalised disk and slope prior weights on the declared grids."""
    # The SPARC prior is Gaussian in x=log10(Upsilon_disk), hence lognormal
    # in Upsilon_disk. Integrating in x avoids an explicit Jacobian.
    x0 = np.log10(args.disk_ml)
    sigma = args.disk_ml_dex
    disk_density = np.exp(-0.5 * ((log10_upsilon - x0) / sigma) ** 2)
    disk_density /= np.sqrt(2.0 * np.pi) * sigma
    disk_weights = disk_density * trapezoid_weights(log10_upsilon)
    covered_disk_prior_mass = float(np.sum(disk_weights))
    if not np.isfinite(covered_disk_prior_mass) or covered_disk_prior_mass <= 0:
        raise ValueError("disk prior has no finite mass on the declared grid")
    disk_weights /= covered_disk_prior_mass

    # For radial profiles, the current model assigns a continuous uniform prior
    # over the declared slope interval. Trapezoidal weights approximate its integral.
    slope_weights = trapezoid_weights(slopes)
    slope_weights /= np.sum(slope_weights)
    return disk_weights, slope_weights, covered_disk_prior_mass


def fully_marginalise_grid(galaxy, args, log10_upsilon, slopes,
                           marginal_log_probability_grid, marginal_valid_grid,
                           negative_grid):
    """Integrate a gas/bulge-marginalised grid over disk M/L and profile slope."""
    if not args.marginalise_disk:
        return None

    disk_weights, slope_weights, covered_disk_mass = grid_prior_weights(
        log10_upsilon, slopes, args)
    if args.profile == "constant":
        joint_weights = disk_weights[np.newaxis, :]
        slope_prior = "none"
    else:
        joint_weights = np.outer(slope_weights, disk_weights)
        slope_prior = f"uniform({args.slope_min},{args.slope_max})"

    # Non-positive or non-finite stellar profiles are outside the physical prior.
    # We condition on the remaining domain, but never hide failures in an allowed cell.
    physical_grid = ~negative_grid
    failed_allowed = physical_grid & ~marginal_valid_grid
    if np.any(failed_allowed):
        return {
            "galaxy": galaxy,
            "profile": args.profile,
            "probability_kind": "fully_marginalised_sparc_prior",
            "log_probability": np.nan,
            "valid": False,
            "invalid_reason": "one_or_more_physical_grid_cells_failed",
            "retained_profile_prior_mass": float(np.sum(joint_weights[physical_grid])),
            "covered_disk_prior_mass": covered_disk_mass,

            "slope_prior": slope_prior,

        }

    integration_mask = physical_grid & marginal_valid_grid & np.isfinite(marginal_log_probability_grid)
    retained_prior_mass = float(np.sum(joint_weights[integration_mask]))
    if not np.any(integration_mask) or retained_prior_mass <= 0:
        log_probability = np.nan
        valid = False
        reason = "no_valid_prior_mass"
    else:
        try:
            log_probability = log_weighted_average(
                marginal_log_probability_grid[integration_mask],
                joint_weights[integration_mask],
            )
            valid = True
            reason = ""
        except Exception as exc:
            log_probability = np.nan
            valid = False
            reason = f"{type(exc).__name__}: {exc}"

    return {
        "galaxy": galaxy,
        "profile": args.profile,
        "probability_kind": "fully_marginalised_sparc_prior",
        "log_probability": log_probability,
        "valid": valid,
        "invalid_reason": reason,
        "retained_profile_prior_mass": retained_prior_mass,
        "covered_disk_prior_mass": covered_disk_mass,

        "slope_prior": slope_prior,

    }


def merge_warnings(result, extra_warnings):
    """Attach warnings raised while constructing a disk model to a probability result."""
    if result is None or not extra_warnings:
        return result
    logp, valid, reason, warning_messages = result
    return logp, valid, reason, warning_messages + list(extra_warnings)


def process_galaxy(galaxy, table, args, rows, log10_upsilon, upsilon, slopes,
                   gas_nodes, gas_weights, bulge_nodes, bulge_weights):
    """Evaluate the profile grid and its SPARC-prior integral for one galaxy."""
    columns = ["Rad", "Vobs", "errV", "Vgas", "Vdisk", "Vbul", "SBdisk", "SBbul"]
    data_path = args.data_root / f"{galaxy}_rotmod.dat"
    rawdata = np.loadtxt(data_path)
    data = pd.DataFrame(rawdata, columns=columns)
    r = data["Rad"].to_numpy(float)
    Vobs = data["Vobs"].to_numpy(float)
    errV = data["errV"].to_numpy(float)
    Vgas = data["Vgas"].to_numpy(float)
    Vbul = data["Vbul"].to_numpy(float)

    shape = (1 if args.profile == "constant" else args.grid_size, args.grid_size)
    log_probability_grid = np.full(shape, np.nan)
    valid_grid = np.zeros(shape, dtype=bool)
    marginal_log_probability_grid = np.full(shape, np.nan)
    marginal_valid_grid = np.zeros(shape, dtype=bool)
    negative_grid = np.zeros(shape, dtype=bool)

    if args.profile == "constant":
        Vdisk = data["Vdisk"].to_numpy(float)
        for p, pdisk in enumerate(upsilon):
            # Constant-profile velocity-squared scales linearly with disk M/L.
            results = evaluate_disk_term(
                Vdisk**2 * pdisk, Vgas, Vbul, r, Vobs, errV, args,
                gas_nodes, gas_weights, bulge_nodes, bulge_weights)
            logp, valid, _, _ = results[0]
            log_probability_grid[0, p] = logp
            valid_grid[0, p] = valid
            if results[1] is not None:
                marginal_log_probability_grid[0, p] = results[1][0]
                marginal_valid_grid[0, p] = results[1][1]
            append_probability_rows(rows, galaxy, args.profile, 0, p, log10_upsilon[p],
                                    pdisk, np.nan, results, args)
    else:
        i_tab = np.flatnonzero(table["Galaxy"].to_numpy() == galaxy)
        if len(i_tab) != 1:
            raise ValueError(f"Expected one SPARC table row for {galaxy!r}, found {len(i_tab)}")
        rdisk = float(table.iloc[i_tab[0]]["Rdisk"])
        zd = 0.196 * rdisk**0.633
        rad, SBdisk = get_SBdisk(galaxy, dens_root=args.dens_root, table=table)
        A, B = compute_A_B(rad, SBdisk)
        for p, pdisk in enumerate(upsilon):
            for m, slope in enumerate(slopes):
                # pdisk is the luminosity-weighted mean. The intercept or
                # normalisation is chosen so every radial profile has that mean.
                if args.profile == "linear":
                    y_intercept = pdisk - slope * (B / A)
                    Upsilon_disk = slope * rad + y_intercept
                else:
                    _, beta = compute_A_beta(rad, SBdisk, slope)
                    norm_scale = pdisk / (beta / A)
                    Upsilon_disk = norm_scale * np.exp(slope * rad)
                negative = np.any(~np.isfinite(Upsilon_disk) | (Upsilon_disk <= 0))
                negative_grid[m, p] = negative
                if negative:
                    invalid = (np.nan, False, "nonpositive_or_nonfinite_disk_ml", [])
                    results = (invalid, invalid) if args.marginalise_gas_bulge else (invalid, None)
                    append_probability_rows(rows, galaxy, args.profile, m, p, log10_upsilon[p],
                                            pdisk, slope, results, args)
                    continue
                try:
                    with warnings.catch_warnings(record=True) as records:
                        warnings.simplefilter("always")
                        surface_density = SBdisk * Upsilon_disk
                        Vdisk_dens = vcdisk(rad, surface_density, z0=zd, rhoz="exp")
                        Vdisk = np.interp(r, rad, Vdisk_dens)
                        results = evaluate_disk_term(
                            Vdisk**2, Vgas, Vbul, r, Vobs, errV, args,
                            gas_nodes, gas_weights, bulge_nodes, bulge_weights)
                    construction_warnings = [str(record.message) for record in records]
                    results = (
                        merge_warnings(results[0], construction_warnings),
                        merge_warnings(results[1], construction_warnings),
                    )
                except Exception as exc:
                    invalid = (np.nan, False, f"{type(exc).__name__}: {exc}", [])
                    results = (invalid, invalid if args.marginalise_gas_bulge else None)
                logp, valid, _, _ = results[0]
                log_probability_grid[m, p] = logp
                valid_grid[m, p] = valid
                if results[1] is not None:
                    marginal_log_probability_grid[m, p] = results[1][0]
                    marginal_valid_grid[m, p] = results[1][1]
                append_probability_rows(rows, galaxy, args.profile, m, p, log10_upsilon[p],
                                        pdisk, slope, results, args)

    fully_marginalised = fully_marginalise_grid(
        galaxy, args, log10_upsilon, slopes,
        marginal_log_probability_grid, marginal_valid_grid, negative_grid,
    )
    return log_probability_grid, valid_grid, negative_grid, fully_marginalised


def finite_plot_grid(log_probability_grid, valid_grid):
    """Replace invalid cells for plotting while preserving their validity mask."""
    plotted = np.where(valid_grid, np.clip(log_probability_grid, -50.0, 0.0), np.nan)
    finite = plotted[np.isfinite(plotted)]
    if finite.size == 0:
        return np.full(plotted.shape, -50.0), np.linspace(-50.0, 0.0, 11)
    low = min(float(np.min(finite)), -1e-6)
    if low == 0.0:
        low = -1.0
    return np.where(np.isfinite(plotted), plotted, low), np.linspace(low, 0.0, 11)


def save_plot(galaxy, args, log10_upsilon_sparc, log10_upsilon, slopes,
              log_probability_grid, valid_grid, negative_grid):
    """Save a diagnostic plot of the conditional probability grid."""
    # The publication values are in the CSV files; this familiar plot remains
    # conditional at the central gas and bulge M/L values.
    plot_values, levels = finite_plot_grid(log_probability_grid, valid_grid)
    plot_root = args.output_root / "plots"
    if args.profile == "constant":
        output_dir = plot_root / "constant_ML_plots"
        output_dir.mkdir(parents=True, exist_ok=True)
        fig, ax = plt.subplots(figsize=(7, 5))
        ax.plot(log10_upsilon, plot_values[0], color="blue")
        ax.set_xlabel("log10(Upsilon_disk)")
        ax.set_ylabel("log(probability | gas, bulge M/L)")
        ax.set_title(f"Conditional log probability - {galaxy}")
    else:
        name = "linear_2D_grids" if args.profile == "linear" else "log_linear_2D_grids"
        output_dir = plot_root / name
        output_dir.mkdir(parents=True, exist_ok=True)
        fig, ax = plt.subplots(figsize=(10, 8))
        im = ax.contourf(log10_upsilon, slopes, plot_values, levels=levels, cmap="viridis")
        negative_plot = np.ma.masked_where(~negative_grid, np.zeros_like(negative_grid, dtype=float))
        ax.imshow(negative_plot, origin="lower",
                  extent=[log10_upsilon[0], log10_upsilon[-1], slopes[0], slopes[-1]],
                  cmap="gray", alpha=0.5, aspect="auto", interpolation="nearest", zorder=2)
        ax.plot(log10_upsilon_sparc, 0.0, marker="*", color="red",
                markersize=14, markeredgecolor="black")
        ax.set_xlabel("log10(Upsilon_disk)")
        ax.set_ylabel("disk slope" if args.profile == "linear" else "log-linear slope")
        ax.set_title(f"Conditional log probability matrix - {galaxy}")
        fig.colorbar(im, ax=ax, label="log(probability | gas, bulge M/L)")
    fig.tight_layout()
    fig.savefig(output_dir / f"{galaxy}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def write_csv(path, rows):
    """Write a list of self-describing result dictionaries to CSV."""
    fieldnames = list(rows[0]) if rows else []
    with path.open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=fieldnames)
        if fieldnames:
            writer.writeheader()
            writer.writerows(rows)


def output_paths(args):
    """Return per-shard paths, preserving legacy names for unsharded runs."""
    suffix = f".shard_{args.shard_id:02d}" if args.shard_id is not None else ""
    return {
        "profile_grid": args.output_root / f"profile_grid{suffix}.csv",
        "fully_marginalised": args.output_root / f"fully_marginalised{suffix}.csv",
        "metadata": args.output_root / f"run_metadata{suffix}.json",
    }


def main():
    """Run the selected profile grid and write results plus provenance metadata."""
    args = parse_args()
    args.output_root.mkdir(parents=True, exist_ok=True)
    paths = output_paths(args)
    if args.shard_id is not None and not args.overwrite_shard:
        existing = [str(path) for path in paths.values() if path.exists()]
        if existing:
            raise FileExistsError(
                "Shard artifacts already exist; choose a new shard ID or pass "
                f"--overwrite-shard: {', '.join(existing)}"
            )
    table = load_sparc_table(args.table_path)
    galaxies = table["Galaxy"].astype(str).tolist()
    if args.galaxy:
        requested = set(args.galaxy)
        missing = sorted(requested.difference(galaxies))
        if missing:
            raise ValueError(f"Galaxy names not found in SPARC table: {missing}")
        galaxies = [galaxy for galaxy in galaxies if galaxy in requested]
    if args.testing:
        galaxies = galaxies[:1]
    elif args.galaxy_range is not None:
        start, end = args.galaxy_range
        galaxies = galaxies[start - 1:end]

    # Gauss-Hermite quadrature integrates the independent lognormal gas and
    # bulge priors. Disk M/L and radial slope are integrated on the output grid.
    gas_nodes, gas_weights = lognormal_quadrature(args.gas_ml, args.gas_ml_dex, args.ml_samples)
    bulge_nodes, bulge_weights = lognormal_quadrature(args.bulge_ml, args.bulge_ml_dex, args.ml_samples)
    log10_upsilon_sparc, log10_upsilon, upsilon, slopes = profile_grids(args)
    rows = []
    fully_marginalised_rows = []
    run_metadata = {
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "profile": args.profile,
        "arguments": {key: str(value) for key, value in vars(args).items()},
        "python": sys.version,
        "platform": platform.platform(),
        "source_sha256": source_sha256(),
        "project_root": str(DEFAULT_PROJECT_ROOT),
        "conditional_probability": True,
        "gas_bulge_marginalisation": bool(args.marginalise_gas_bulge),
        "disk_marginalisation": bool(args.marginalise_disk),
        "gas_ml_median": args.gas_ml,
        "disk_ml_median": args.disk_ml,
        "bulge_ml_median": args.bulge_ml,
        "gas_ml_sigma_dex": args.gas_ml_dex,
        "disk_ml_sigma_dex": args.disk_ml_dex,
        "bulge_ml_sigma_dex": args.bulge_ml_dex,
        "gas_bulge_quadrature": "Gauss-Hermite in log10(M/L), normalised weights",
        "disk_quadrature": "trapezoidal in log10(M/L) with Gaussian SPARC prior",
        "slope_quadrature": "trapezoidal with uniform prior for radial profiles",
        "quadrature_samples_per_gas_or_bulge_term": args.ml_samples,
        "gas_term_convention": "Vgas*abs(Vgas)",
        "negative_disk_ml": "excluded; full integral is renormalised on the physical profile domain",
        "tie_policy": "arithmetic_mean_over_orderings_in_logspace",
        "disk_ml_reference": f"{args.disk_ml} with {args.disk_ml_dex} dex SPARC uncertainty",
        "galaxy_count": len(galaxies),
    }
    for index, galaxy in enumerate(galaxies, start=1):
        print(f"Processing galaxy {index}/{len(galaxies)}", flush=True)
        grid, valid, negative, fully_marginalised = process_galaxy(
            galaxy, table, args, rows, log10_upsilon, upsilon, slopes,
            gas_nodes, gas_weights, bulge_nodes, bulge_weights)
        if fully_marginalised is not None:
            fully_marginalised_rows.append(fully_marginalised)
        save_plot(galaxy, args, log10_upsilon_sparc, log10_upsilon, slopes,
                  grid, valid, negative)

    results_path = paths["profile_grid"]
    full_results_path = paths["fully_marginalised"]
    write_csv(results_path, rows)
    write_csv(full_results_path, fully_marginalised_rows)
    run_metadata["finished_utc"] = datetime.now(timezone.utc).isoformat()
    run_metadata["results_csv"] = str(results_path)
    run_metadata["fully_marginalised_csv"] = str(full_results_path)
    with paths["metadata"].open("w") as output:
        json.dump(run_metadata, output, indent=2)
    print(f"Wrote {len(rows)} grid rows to {results_path}")
    print(f"Wrote {len(fully_marginalised_rows)} fully marginalised rows to {full_results_path}")


if __name__ == "__main__":
    main()
