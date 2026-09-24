#!/usr/bin/env python3
"""Fit constant, linear, and log-linear K-band M/L profiles."""
from __future__ import annotations
import argparse, hashlib, json
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd

REQUIRED_COLUMNS = ["Galaxy", "Sample", "Radius_kpc", "ml_K"]

def rms(values):
    return float(np.sqrt(np.mean(values ** 2))) if len(values) else np.nan

def fit_group(group):
    galaxy, sample = str(group["Galaxy"].iloc[0]), str(group["Sample"].iloc[0])
    radius = pd.to_numeric(group["Radius_kpc"], errors="coerce")
    ml = pd.to_numeric(group["ml_K"], errors="coerce")
    finite = np.isfinite(radius) & np.isfinite(ml)
    usable = finite & (ml > 0)
    exclusions = []
    for idx, row in group.iterrows():
        flags = []
        r0 = pd.to_numeric(row["Radius_kpc"], errors="coerce")
        y0 = pd.to_numeric(row["ml_K"], errors="coerce")
        if not np.isfinite(r0):
            flags.append("missing_or_nonfinite_radius")
        if not np.isfinite(y0):
            flags.append("missing_or_nonfinite_ml_K")
        elif y0 <= 0:
            flags.append("nonpositive_ml_K")
        if flags:
            exclusions.append({"input_row": int(idx) + 2, "Galaxy": galaxy, "Sample": sample,
                               "exclusion_reason": ";".join(flags)})
    r = radius[usable].to_numpy(float)
    y = ml[usable].to_numpy(float)
    rp = float(np.mean(r)) if len(r) else np.nan
    duplicate_count = int(len(r) - len(np.unique(r))) if len(r) else 0
    span = float(np.ptp(r)) if len(r) else np.nan
    record = {
        "galaxy": galaxy,
        "sample": sample,
        "sample_group": "COSMOS" if sample.startswith("COSMOS") else sample,
        "cosmos_kin_flag": sample == "COSMOS_kin",
        "n_raw": int(len(group)),
        "n_usable": int(len(y)),
        "n_excluded": int(len(group) - len(y)),
        "n_duplicate_radius": duplicate_count,
        "duplicate_radius_flag": duplicate_count > 0,
        "profile_count_mismatch_flag": bool(len(group) != group["Radius_kpc"].nunique(dropna=True)),
        "mean_radius_kpc": rp,
        "radius_span_kpc": span,
        "constant_c": np.nan,
        "constant_rms": np.nan,
        "constant_success": bool(len(y) >= 1),
        "linear_intercept_c0": np.nan,
        "linear_c_at_Rp": np.nan,
        "linear_slope": np.nan,
        "linear_rms": np.nan,
        "linear_success": False,
        "log_linear_log_intercept_a0": np.nan,
        "log_linear_intercept_c0": np.nan,
        "log_linear_c_at_Rp": np.nan,
        "log_linear_slope": np.nan,
        "log_linear_rms_native": np.nan,
        "log_linear_rms_log": np.nan,
        "log_linear_success": False,
        "fit_status": "ok" if len(y) >= 2 else "insufficient_usable_points",
        "exclusion_notes": "",
    }
    if len(y):
        c = float(np.mean(y))
        record.update(constant_c=c, constant_rms=rms(y - c))
    if len(y) >= 2 and span > 0:
        slope, c0 = np.polyfit(r, y, 1)
        c_rp = c0 + slope * rp
        pred = c0 + slope * r
        record.update(linear_intercept_c0=float(c0), linear_c_at_Rp=float(c_rp),
                      linear_slope=float(slope), linear_rms=rms(y - pred), linear_success=True)
    elif len(y) >= 2:
        record["fit_status"] = "duplicate_radius_only"
    if len(y) >= 2:
        log_y = np.log(y)
        slope, log_c0 = np.polyfit(r, log_y, 1)
        log_c_rp = log_c0 + slope * rp
        pred_log = log_c0 + slope * r
        record.update(log_linear_log_intercept_a0=float(log_c0),
                      log_linear_intercept_c0=float(np.exp(log_c0)),
                      log_linear_c_at_Rp=float(np.exp(log_c_rp)),
                      log_linear_slope=float(slope),
                      log_linear_rms_native=rms(y - np.exp(pred_log)),
                      log_linear_rms_log=rms(log_y - pred_log),
                      log_linear_success=True)
    if exclusions:
        record["exclusion_notes"] = ";".join(sorted({e["exclusion_reason"] for e in exclusions}))
    return record, exclusions

def make_summary(records, label):
    subset = records if label == "ALL" else records[records.sample_group == label]
    return {
        "summary_group": label,
        "n_galaxy_sample_groups": int(len(subset)),
        "n_constant_success": int(subset.constant_success.sum()),
        "n_linear_success": int(subset.linear_success.sum()),
        "n_log_linear_success": int(subset.log_linear_success.sum()),
        "n_duplicate_radius_groups": int(subset.duplicate_radius_flag.sum()),
        "n_profile_count_mismatch_groups": int(subset.profile_count_mismatch_flag.sum()),
        "n_cosmos_kin_groups": int(subset.cosmos_kin_flag.sum()),
        "median_constant_c": float(subset.constant_c.median()) if len(subset) else np.nan,
        "median_linear_intercept_c0": float(subset.linear_intercept_c0.median()) if subset.linear_success.any() else np.nan,
        "median_linear_slope": float(subset.linear_slope.median()) if subset.linear_success.any() else np.nan,
        "median_log_linear_intercept_c0": float(subset.log_linear_intercept_c0.median()) if subset.log_linear_success.any() else np.nan,
        "median_log_linear_slope": float(subset.log_linear_slope.median()) if subset.log_linear_success.any() else np.nan,
    }

def run(input_csv, output_dir):
    data = pd.read_csv(input_csv)
    missing = [c for c in REQUIRED_COLUMNS if c not in data.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    output_dir.mkdir(parents=True, exist_ok=True)
    records, exclusions = [], []
    for _, group in data.groupby(["Sample", "Galaxy"], sort=True, dropna=False):
        record, excluded = fit_group(group)
        records.append(record)
        exclusions.extend(excluded)
    fits = pd.DataFrame(records)
    fits.to_csv(output_dir / "ml_profile_fits.csv", index=False)
    pd.DataFrame(exclusions, columns=["input_row", "Galaxy", "Sample", "exclusion_reason"]).to_csv(
        output_dir / "ml_profile_exclusions.csv", index=False)
    pd.DataFrame([make_summary(fits, x) for x in ["COSMOS", "XMM", "LADUMA", "ALL"]]).to_csv(
        output_dir / "ml_profile_summary.csv", index=False)
    metadata = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "input_csv": str(input_csv.resolve()),
        "input_sha256": hashlib.sha256(input_csv.read_bytes()).hexdigest(),
        "required_columns": REQUIRED_COLUMNS,
        "n_input_rows": int(len(data)),
        "n_unique_galaxies": int(data.Galaxy.nunique()),
        "n_unique_sample_galaxy_groups": int(data[["Sample", "Galaxy"]].drop_duplicates().shape[0]),
        "expected_129_130_count_note": "The supplied file contains 129 unique Galaxy and Sample/Galaxy groups; the 130-object reference count is not silently assumed.",
        "sample_counts": {str(k): int(v) for k, v in data.Sample.value_counts(dropna=False).items()},
        "unique_galaxies_by_sample": {str(k): int(v) for k, v in data.groupby("Sample").Galaxy.nunique().items()},
        "rows_by_sample": {str(k): int(v) for k, v in data.groupby("Sample").size().items()},
        "ml_uncertainties_used": False,
        "fit_definition": "unweighted least squares; linear model Upsilon=a+bR; log-linear model log(Upsilon)=a+bR",
        "normalization_convention": "linear_intercept_c0 and log_linear_intercept_c0 are the R=0 intercepts; *_c_at_Rp are derived values at the mean radius",
        "notes": ["Duplicate radii are retained and flagged.", "Missing/nonfinite radius or non-positive/nonfinite ml_K rows are excluded and listed.", "COSMOS_kin remains identifiable; COSMOS summaries combine COSMOS_DR1 and COSMOS_kin."]
    }
    (output_dir / "ml_profile_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path("ML_Andreea.csv"))
    parser.add_argument("--output-dir", type=Path, default=Path("outputs/ml_profile_fits"))
    args = parser.parse_args()
    run(args.input, args.output_dir)

if __name__ == "__main__":
    main()

