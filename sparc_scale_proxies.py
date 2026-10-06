#!/usr/bin/env python3
"""Post-process a SPARC profile grid into exploratory baryonic scale proxies.

The disk M/L parameter in mono_SPARC_RAR.py is the luminosity-weighted mean;
radial profiles are normalised to preserve that mean. This tool uses the same
profile definitions, estimates the bulge light fraction from the supplied
SBdisk/SBbul profiles, and reports mass-prior quadrature summaries. These
proxies are diagnostics, not galaxy cuts or theory tests.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from utils import (
    DEFAULT_DENS_ROOT,
    DEFAULT_TABLE_PATH,
    compute_A_B,
    compute_A_beta,
    load_sparc_table,
    lognormal_quadrature,
)

KPC_METRES = 3.0856775814913673e19
SOLAR_MASS_KG = 1.98847e30
G_SI = 6.67430e-11
A0_DEFAULT_M_S2 = 1.2e-10
G_KPC_KMS2_PER_MSUN = G_SI * SOLAR_MASS_KG * 1.0e-6 / KPC_METRES

GRID_COLUMNS = [
    "galaxy", "profile", "probability_kind", "grid_i", "grid_j",
    "log10_upsilon_disk", "upsilon_disk_mean", "disk_slope", "valid",
    "invalid_reason",
]
OUTPUT_FIELDS = [
    "galaxy", "grid_i", "grid_j", "upsilon_disk_mean", "disk_slope",
    "source_grid_valid", "source_invalid_reason", "mass_status",
    "M_disk_1e9_Msun", "M_bar_median_1e9_Msun", "M_bar_qmean_1e9_Msun",
    "M_bar_q16_1e9_Msun", "M_bar_q84_1e9_Msun", "r_M_median_kpc",
    "r_M_qmean_kpc", "r_M_q16_kpc", "r_M_q84_kpc",
    "eta_median", "eta_qmean", "eta_q16", "eta_q84",
]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def as_bool(value) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return False
    return str(value).strip().lower() in {"true", "1", "yes"}


def weighted_quantile(values, weights, quantile: float) -> float:
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    good = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    if not np.any(good):
        return float("nan")
    values, weights = values[good], weights[good]
    order = np.argsort(values)
    values, weights = values[order], weights[order]
    cdf = np.cumsum(weights) / np.sum(weights)
    return float(values[min(np.searchsorted(cdf, quantile, side="left"), len(values) - 1)])


def prepare_light_profile(path: Path, catalog_luminosity_1e9: float) -> dict:
    """Read the decomposition and scale its integrated component ratio to L3.6."""
    names = ["Rad", "SBdisk", "SBbulge"]
    frame = pd.read_csv(path, sep="\t", names=names, skiprows=1, comment="#")
    if not set(names).issubset(frame.columns) or frame.empty:
        raise ValueError("surface-brightness profile is empty or malformed")
    frame = frame[names].apply(pd.to_numeric, errors="coerce").dropna()
    frame = frame.sort_values("Rad")
    radius = frame["Rad"].to_numpy(dtype=float)
    disk = frame["SBdisk"].to_numpy(dtype=float) * 1.0e6  # Lsun/pc^2 -> Lsun/kpc^2
    bulge = frame["SBbulge"].to_numpy(dtype=float) * 1.0e6
    if len(radius) < 2 or np.any(~np.isfinite(radius)) or np.any(radius < 0):
        raise ValueError("radius grid must contain at least two finite non-negative values")
    if np.any(np.diff(radius) <= 0):
        raise ValueError("radius grid must be strictly increasing with no duplicates")
    if np.any(~np.isfinite(disk)) or np.any(~np.isfinite(bulge)) or np.any(disk < 0) or np.any(bulge < 0):
        raise ValueError("surface brightness must be finite and non-negative")
    raw_disk_lsun, disk_first_moment = compute_A_B(radius, disk)
    raw_bulge_lsun, _ = compute_A_B(radius, bulge)
    raw_total_lsun = raw_disk_lsun + raw_bulge_lsun
    if not np.isfinite(raw_total_lsun) or raw_total_lsun <= 0 or raw_disk_lsun <= 0:
        raise ValueError("integrated disk-plus-bulge light must be positive")
    if not np.isfinite(catalog_luminosity_1e9) or catalog_luminosity_1e9 <= 0:
        raise ValueError("catalogue L3.6 must be finite and positive")
    catalog_lsun = catalog_luminosity_1e9 * 1.0e9
    scale = catalog_lsun / raw_total_lsun
    disk_luminosity = raw_disk_lsun * scale
    bulge_luminosity = raw_bulge_lsun * scale
    return {
        "radius": radius,
        "disk_sb": disk * scale,
        "raw_disk_sb": disk,
        "raw_bulge_lsun": raw_bulge_lsun,
        "raw_disk_lsun": raw_disk_lsun,
        "raw_total_lsun": raw_total_lsun,
        "scale": scale,
        "disk_luminosity_1e9": disk_luminosity / 1.0e9,
        "bulge_luminosity_1e9": bulge_luminosity / 1.0e9,
        "bulge_fraction": bulge_luminosity / catalog_lsun,
        "rmax": float(radius[-1]),
        "n_radii": int(len(radius)),
        "bulge_method": "integrated SBbulge/(SBdisk+SBbulge) light fraction, jointly scaled to catalogue L3.6; no radial extrapolation",
    }


def disk_mass_for_profile(profile: str, mean_upsilon: float, slope: float,
                          light: dict) -> tuple[float, np.ndarray]:
    """Integrate the exact constant/linear/log-linear profile used by the grid."""
    radius = light["radius"]
    sb = light["raw_disk_sb"]
    A, B = compute_A_B(radius, sb)
    if not np.isfinite(A) or A <= 0 or not np.isfinite(mean_upsilon) or mean_upsilon <= 0:
        raise ValueError("non-positive or non-finite disk luminosity/M/L")
    if profile == "constant":
        upsilon = np.full_like(radius, mean_upsilon)
    elif profile == "linear":
        if not np.isfinite(slope):
            raise ValueError("linear profile slope is not finite")
        upsilon = slope * radius + (mean_upsilon - slope * B / A)
    elif profile == "log_linear":
        if not np.isfinite(slope):
            raise ValueError("log-linear profile slope is not finite")
        _, beta = compute_A_beta(radius, sb, slope)
        if not np.isfinite(beta) or beta <= 0:
            raise ValueError("log-linear normalization integral is not positive and finite")
        upsilon = (mean_upsilon / (beta / A)) * np.exp(slope * radius)
    else:
        raise ValueError(f"unsupported profile {profile!r}")
    if np.any(~np.isfinite(upsilon)) or np.any(upsilon <= 0):
        raise ValueError("radial disk M/L is non-positive or non-finite")
    # Restore the common catalogue luminosity scale after using the raw profile
    # for exactly the same luminosity-weighted normalization as the analysis.
    mass_lsun = compute_A_B(radius, sb * upsilon)[0] * light["scale"]
    return float(mass_lsun / 1.0e9), upsilon


def mass_distribution(m_disk_1e9: float, mhi_1e9: float,
                      lbulge_1e9: float, gas_nodes, gas_weights,
                      bulge_nodes, bulge_weights, gas_median: float,
                      bulge_median: float, helium_factor: float,
                      a0_m_s2: float, a_proxy_kpc: float) -> dict:
    gas_mass = np.asarray(gas_nodes, dtype=float)[:, None] * (helium_factor * mhi_1e9)
    bulge_mass = np.asarray(bulge_nodes, dtype=float)[None, :] * lbulge_1e9
    weights = np.asarray(gas_weights, dtype=float)[:, None] * np.asarray(bulge_weights, dtype=float)[None, :]
    mbar = m_disk_1e9 + gas_mass + bulge_mass
    if np.any(~np.isfinite(mbar)) or np.any(mbar <= 0):
        raise ValueError("baryonic mass quadrature contains a non-positive/non-finite value")
    if not np.isclose(np.sum(weights), 1.0, rtol=1e-10, atol=1e-12):
        raise ValueError("gas/bulge quadrature weights do not sum to one")
    a0_kms2_kpc = a0_m_s2 * KPC_METRES / 1.0e6
    r_m = np.sqrt(G_KPC_KMS2_PER_MSUN * mbar * 1.0e9 / a0_kms2_kpc)
    eta = a_proxy_kpc / r_m

    def stats(values, names):
        return {
            names[0]: float(np.sum(values * weights) / np.sum(weights)),
            names[1]: weighted_quantile(values.ravel(), weights.ravel(), 0.16),
            names[2]: weighted_quantile(values.ravel(), weights.ravel(), 0.84),
        }

    mbar_median = m_disk_1e9 + gas_median * helium_factor * mhi_1e9 + bulge_median * lbulge_1e9
    r_m_median = math.sqrt(G_KPC_KMS2_PER_MSUN * mbar_median * 1.0e9 / a0_kms2_kpc)
    result = {
        "M_bar_median_1e9_Msun": float(mbar_median),
        "r_M_median_kpc": float(r_m_median),
        "eta_median": float(a_proxy_kpc / r_m_median),
    }
    result.update(stats(mbar, ("M_bar_qmean_1e9_Msun", "M_bar_q16_1e9_Msun", "M_bar_q84_1e9_Msun")))
    result.update(stats(r_m, ("r_M_qmean_kpc", "r_M_q16_kpc", "r_M_q84_kpc")))
    result.update(stats(eta, ("eta_qmean", "eta_q16", "eta_q84")))
    return result


def get_run_settings(metadata_path: Path | None) -> tuple[dict, str]:
    defaults = {
        "gas_ml": 1.0, "bulge_ml": 0.7, "gas_ml_dex": 0.04,
        "bulge_ml_dex": 0.1, "quadrature_samples": 9,
    }
    if metadata_path is None or not metadata_path.exists():
        return defaults, "documented script defaults (no adjacent run metadata found)"
    metadata = json.loads(metadata_path.read_text())
    settings = metadata.get("settings", {})
    def pick(names, default):
        for name in names:
            if name in settings:
                return settings[name]
            if name in metadata:
                return metadata[name]
        return default
    values = {
        "gas_ml": float(pick(("gas_ml", "gas_ml_median"), defaults["gas_ml"])),
        "bulge_ml": float(pick(("bulge_ml", "bulge_ml_median"), defaults["bulge_ml"])),
        "gas_ml_dex": float(pick(("gas_ml_dex", "gas_ml_sigma_dex"), defaults["gas_ml_dex"])),
        "bulge_ml_dex": float(pick(("bulge_ml_dex", "bulge_ml_sigma_dex"), defaults["bulge_ml_dex"])),
        "quadrature_samples": int(pick(("ml_samples", "quadrature_samples_per_gas_or_bulge_term"), defaults["quadrature_samples"])),
    }
    return values, str(metadata_path)


def make_record(row: dict, catalogue_row, light: dict | None,
                settings: dict, gas_nodes, gas_weights,
                bulge_nodes, bulge_weights, args) -> dict:
    """Build a compact per-grid-cell row while preserving all invalid cells."""
    record = {field: "" for field in OUTPUT_FIELDS}
    for field in ("galaxy", "profile", "grid_i", "grid_j",
                  "upsilon_disk_mean", "disk_slope"):
        record[field] = row.get(field, "")
    record["source_grid_valid"] = as_bool(row.get("valid"))
    record["source_invalid_reason"] = row.get("invalid_reason", "")
    record["mass_status"] = "missing_or_invalid_light_profile"
    if light is None:
        return record
    if not record["source_grid_valid"]:
        record["mass_status"] = "invalid_source_grid_cell"
        return record
    try:
        m_disk, _ = disk_mass_for_profile(
            str(row["profile"]), float(row["upsilon_disk_mean"]),
            float(row["disk_slope"]) if pd.notna(row["disk_slope"]) else 0.0,
            light,
        )
        disk_check = float(row["upsilon_disk_mean"]) * light["disk_luminosity_1e9"]
        if not math.isclose(m_disk, disk_check, rel_tol=1e-10, abs_tol=1e-12):
            raise ValueError("profile-integrated disk mass fails luminosity-weighted mean check")
        proxy = mass_distribution(
            m_disk, float(catalogue_row.MHI), light["bulge_luminosity_1e9"],
            gas_nodes, gas_weights, bulge_nodes, bulge_weights,
            settings["gas_ml"], settings["bulge_ml"], args.helium_factor,
            args.a0, float(catalogue_row.Rdisk) / 0.85,
        )
        record.update({"mass_status": "ok", "M_disk_1e9_Msun": m_disk, **proxy})
    except (ArithmeticError, ValueError, FloatingPointError) as exc:
        record["mass_status"] = f"invalid_physical_profile: {exc}"
    return record


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile-grid", required=True, type=Path, help="Merged profile_grid.csv for one profile/settings batch")
    parser.add_argument("--table-path", type=Path, default=DEFAULT_TABLE_PATH)
    parser.add_argument("--dens-root", type=Path, default=DEFAULT_DENS_ROOT)
    parser.add_argument("--metadata", type=Path, default=None, help="Run metadata; defaults to profile_grid's sibling run_metadata.json")
    parser.add_argument("--output", type=Path, default=None, help="Output CSV; defaults to physical_scale_proxies.csv beside input")
    parser.add_argument("--galaxy", action="append", help="Optional galaxy filter; may be repeated")
    parser.add_argument("--chunk-size", type=int, default=100000)
    parser.add_argument("--helium-factor", type=float, default=1.33)
    parser.add_argument("--a0", type=float, default=A0_DEFAULT_M_S2, help="a0 in m/s^2")
    parser.add_argument("--replace-existing", action="store_true")
    args = parser.parse_args()
    if args.chunk_size < 1 or args.helium_factor <= 0 or args.a0 <= 0:
        parser.error("chunk size, helium factor, and a0 must be positive")
    return args


def main():
    args = parse_args()
    input_path = args.profile_grid.expanduser().resolve()
    if not input_path.is_file():
        raise FileNotFoundError(input_path)
    metadata_path = args.metadata or input_path.parent / "run_metadata.json"
    settings, settings_source = get_run_settings(metadata_path)
    gas_nodes, gas_weights = lognormal_quadrature(
        settings["gas_ml"], settings["gas_ml_dex"], settings["quadrature_samples"]
    )
    bulge_nodes, bulge_weights = lognormal_quadrature(
        settings["bulge_ml"], settings["bulge_ml_dex"], settings["quadrature_samples"]
    )
    output = args.output or input_path.parent / "physical_scale_proxies.csv"
    output = output.expanduser().resolve()
    out_metadata = output.with_name(output.stem + "_metadata.json")
    output.parent.mkdir(parents=True, exist_ok=True)
    if not args.replace_existing and (output.exists() or out_metadata.exists()):
        raise FileExistsError("output exists; choose another --output or explicitly pass --replace-existing")
    table = load_sparc_table(args.table_path).set_index("Galaxy")
    light_cache = {}
    galaxy_inputs = {}
    written = 0
    invalid = 0
    profiles = set()
    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=OUTPUT_FIELDS, extrasaction="ignore")
        writer.writeheader()
        for chunk in pd.read_csv(input_path, usecols=lambda name: name in GRID_COLUMNS, chunksize=args.chunk_size):
            if "probability_kind" not in chunk or "valid" not in chunk:
                raise ValueError("input grid is missing expected profile-grid columns")
            chunk = chunk[chunk["probability_kind"] == "conditional"]
            if args.galaxy:
                chunk = chunk[chunk["galaxy"].isin(args.galaxy)]
            for row in chunk.to_dict(orient="records"):
                galaxy = str(row["galaxy"])
                if galaxy not in table.index:
                    raise KeyError(f"{galaxy} is absent from SPARC table {args.table_path}")
                catalogue_row = table.loc[galaxy]
                if galaxy not in light_cache:
                    light_path = args.dens_root / f"{galaxy}.dens"
                    try:
                        light_cache[galaxy] = (
                            prepare_light_profile(light_path, float(catalogue_row.L)),
                            sha256_file(light_path),
                        )
                    except (OSError, ValueError) as exc:
                        light_cache[galaxy] = (None, "")
                        light_cache[(galaxy, "error")] = str(exc)
                light, light_hash = light_cache[galaxy]
                if galaxy not in galaxy_inputs:
                    galaxy_inputs[galaxy] = {
                        "catalogue": {
                            "L_total_1e9_Lsun": float(catalogue_row.L),
                            "Rdisk_kpc": float(catalogue_row.Rdisk),
                            "MHI_raw_1e9_Msun": float(catalogue_row.MHI),
                        },
                        "a_proxy_kpc": float(catalogue_row.Rdisk) / 0.85,
                        "light_profile": None if light is None else {
                            "sha256": light_hash,
                            "L_disk_1e9_Lsun": light["disk_luminosity_1e9"],
                            "L_bulge_1e9_Lsun": light["bulge_luminosity_1e9"],
                            "bulge_light_fraction": light["bulge_fraction"],
                            "raw_light_scale_to_catalog": light["scale"],
                            "rmax_kpc": light["rmax"],
                            "n_radii": light["n_radii"],
                            "bulge_light_method": light["bulge_method"],
                        },
                        "median_nuisance_masses_1e9_Msun": {
                            "gas": args.helium_factor * float(catalogue_row.MHI) * settings["gas_ml"],
                            "bulge": None if light is None else light["bulge_luminosity_1e9"] * settings["bulge_ml"],
                        },
                        "light_profile_error": light_cache.get((galaxy, "error")),
                    }
                rec = make_record(
                    row, catalogue_row, light, settings, gas_nodes, gas_weights,
                    bulge_nodes, bulge_weights, args,
                )
                if light is None:
                    rec["mass_status"] = "missing_or_invalid_light_profile: " + light_cache[(galaxy, "error")]
                writer.writerow(rec)
                written += 1
                invalid += rec["mass_status"] != "ok"
                profiles.add(str(row["profile"]))
    metadata = {
        "schema_version": 2,
        "tool": "sparc_scale_proxies.py",
        "input_profile_grid": str(input_path),
        "input_profile_grid_sha256": sha256_file(input_path),
        "run_metadata": str(metadata_path) if metadata_path.exists() else None,
        "run_metadata_sha256": sha256_file(metadata_path) if metadata_path.exists() else None,
        "settings_source": settings_source,
        "nuisance_settings": settings,
        "table_path": str(args.table_path.resolve()),
        "table_sha256": sha256_file(args.table_path),
        "dens_root": str(args.dens_root.resolve()),
        "source_probability_kind": "conditional",
        "profile_types": sorted(profiles),
        "rows_written": written,
        "invalid_mass_rows": int(invalid),
        "galaxy_inputs": galaxy_inputs,
        "shared_physical_settings": {
            "helium_factor": args.helium_factor,
            "a0_m_s2": args.a0,
            "gas_mass_convention": "Mgas = gas_ML * helium_factor * catalogue MHI; helium included exactly once; molecular gas omitted",
            "disk_profile_normalisation": "exact mono_SPARC_RAR.py luminosity-weighted mean; integrated disk mass checked against mean_ML * Ldisk",
            "a_proxy_definition": "Rdisk/0.85 kpc; exploratory Miyamoto-Nagai geometry correspondence",
            "r_M_definition": "sqrt(G*Mbar/a0); G=6.67430e-11 m^3 kg^-1 s^-2, Msun=1.98847e30 kg, kpc=3.0856775814913673e19 m",
            "mass_quadrature": "independent Gauss-Hermite lognormal gas and bulge nuisance priors; reported mean and inverse-weighted-CDF 16/84 percentiles are prior summaries, not probability-conditioned posteriors",
            "selection_policy": "invalid source cells and invalid light/profile cells are retained and flagged; this diagnostic does not exclude galaxies or alter monotonicity probabilities",
            "limitations": "surface-brightness profiles are integrated only over tabulated radial coverage; disk/bulge component ratio is jointly scaled to catalogue total L3.6; screen is not a theory test",
        },
        "output_csv": str(output),
        "output_schema": OUTPUT_FIELDS,
    }
    out_metadata.write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"Wrote {written} conditional grid rows ({invalid} flagged) to {output}")
    print(f"Provenance: {out_metadata}")


if __name__ == "__main__":
    main()