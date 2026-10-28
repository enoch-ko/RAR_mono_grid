#!/usr/bin/env python3
"""Build pinned profile inputs for exploratory force/RAR figures; no search calls."""
from pathlib import Path
import csv, hashlib, json, math, sys
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT))
from utils import load_sparc_table, get_SBdisk, compute_A_B, compute_A_beta

OUT = HERE
BF_ROOT = ROOT / "outputs/bayes_factors/bf_grid_rollout_20261023_v1"
GP_ROOT = ROOT / "outputs/GP_searches/science_20261016/gp_d5000_a16"
GALAXIES = ("UGC04325", "KK98-251")
FAMILIES = {"constant": "constant_01", "linear": "linear_04", "log_linear": "log_linear_04"}
PRESET_CUTS = (.15865525393145707, .022750131948179195, .5, .0013498980316300957)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def add_pin(p, pins):
    p = Path(p).resolve()
    if not p.is_file():
        raise FileNotFoundError(p)
    pins[str(p)] = sha(p)


def reconstructed_profile(row, family, rad, sb):
    mean = float(row["upsilon_disk_mean"])
    if family == "constant":
        return np.full(rad.shape, mean)
    slope = float(row["disk_slope"])
    A, B = compute_A_B(rad, sb)
    if family == "linear":
        return mean - slope * B / A + slope * rad
    _, beta = compute_A_beta(rad, sb, slope)
    return mean * A / beta * np.exp(slope * rad)


def finite_float(value):
    try:
        value = float(value)
        return value if math.isfinite(value) else None
    except (TypeError, ValueError):
        return None


def grid_selection(galaxy, family, run_name, rad, sb, pins):
    folder = BF_ROOT / run_name / galaxy
    summary_path = folder / "summary.json"
    context_path = folder / "context.json"
    cells_path = folder / "cells.csv"
    summary = json.loads(summary_path.read_text())
    if summary.get("status") != "complete_saved_grid_accuracy_unverified":
        raise RuntimeError(f"Unexpected grid status: {summary_path}: {summary.get('status')}")
    rows = list(csv.DictReader(cells_path.open(newline="")))
    valid = []
    for row in rows:
        logp = finite_float(row.get("saved_logP"))
        jrad = finite_float(row.get("J_rad"))
        if logp is not None and jrad is not None:
            valid.append(row)
    if len(valid) != int(summary["valid_cells"]):
        raise RuntimeError(f"Valid cell mismatch for {cells_path}: {len(valid)} vs {summary['valid_cells']}")
    if not valid:
        raise RuntimeError(f"No valid saved cells in {cells_path}")
    selected = min(valid, key=lambda r: (-float(r["saved_logP"]), float(r["J_rad"]), int(r["grid_i"]), int(r["grid_j"])))
    upsilon = reconstructed_profile(selected, family, rad, sb)
    if not np.all(np.isfinite(upsilon)) or np.any(upsilon <= 0):
        raise RuntimeError(f"Nonpositive or non-finite reconstructed profile: {cells_path}")
    span = float(rad[-1] - rad[0])
    slope = 0.0 if family == "constant" else float(selected["disk_slope"])
    if family == "log_linear" and slope != 0:
        x = abs(slope) * span
        if x < 1e-3:
            nodes, weights = np.polynomial.legendre.leggauss(8)
            z = 2 * upsilon[0] * np.expm1(slope * (nodes + 1) * span / 2) + (2 * upsilon[0] - 1)
            j_reconstructed = float(np.dot(weights, z * z) / 2)
        else:
            top = float(upsilon.max())
            m1 = top * (-math.expm1(-x)) / x
            m2 = top * top * (-math.expm1(-2 * x)) / (2 * x)
            j_reconstructed = max(0.0, 4 * m2 - 4 * m1 + 1)
    else:
        z = (upsilon - .5) / .5
        j_reconstructed = float(np.sum(np.diff(rad) * (z[:-1]**2 + z[:-1] * z[1:] + z[1:]**2) / 3) / span)
    j_saved = float(selected["J_rad"])
    if not math.isclose(j_reconstructed, j_saved, rel_tol=1e-6, abs_tol=1e-10):
        raise RuntimeError(f"J_rad mismatch for {cells_path}: {j_reconstructed} vs {j_saved}")
    for p in (summary_path, context_path, cells_path):
        add_pin(p, pins)
    return {
        "source": str(cells_path.resolve()),
        "status": summary["status"],
        "scientific_validated": summary.get("scientific_validated"),
        "valid_cells": len(valid),
        "max_saved_probability": math.exp(float(selected["saved_logP"])),
        "selection_rule": "highest saved_logP; ties broken by lower saved J_rad then grid indices",
        "selected": {"row": selected, "J_rad": j_saved, "upsilon": upsilon.tolist()},
    }


def gp_selection(galaxy, rad, pins):
    gp_path = GP_ROOT / f"{galaxy}_gp_search.json"
    checkpoint_dir = GP_ROOT / f"{galaxy}_gp_checkpoints"
    progress_path = checkpoint_dir / "progress.json"
    checkpoint_manifest_path = checkpoint_dir / "manifest.json"
    gp = json.loads(gp_path.read_text())
    progress = json.loads(progress_path.read_text())
    if gp.get("status") != "finite_search_completed" or not gp.get("diagnostic_only"):
        raise RuntimeError(f"Unexpected GP status: {gp_path}")
    records = gp["records"]
    valid = [r for r in records if r.get("valid") and finite_float(r.get("marginal_logp")) is not None and finite_float(r.get("J_ref")) is not None]
    if not valid:
        raise RuntimeError(f"No valid GP candidates in {gp_path}")
    maxp = max(valid, key=lambda r: (float(r["marginal_logp"]), -int(r["original_candidate_index"])))
    minj = min(valid, key=lambda r: (float(r["J_ref"]), int(r["original_candidate_index"])))
    frontier_ids = {int(records[int(f["record_index"])]["original_candidate_index"]) for f in gp["frontier"]}
    if int(maxp["original_candidate_index"]) not in frontier_ids or int(minj["original_candidate_index"]) not in frontier_ids:
        raise RuntimeError(f"Selected GP candidates are absent from the saved frontier: {gp_path}")
    for entry in gp["frontier"]:
        profile = entry.get("profile_on_photometric_grid", [])
        if len(profile) != len(rad):
            raise RuntimeError(f"GP profile/radius length mismatch for {gp_path}")
    baseline = next(r for r in valid if r.get("source") == "SPARC_constant_0.5")
    preset_winners = {}
    for cut in PRESET_CUTS:
        qualifying = [r for r in valid if float(r["marginal_logp"]) >= math.log(cut)]
        preset_winners[str(cut)] = min(qualifying, key=lambda r: (float(r["J_ref"]), int(r["original_candidate_index"]))) if qualifying else None
    display_cut = math.exp(float(maxp["marginal_logp"]))
    preset_winners[str(display_cut)] = maxp
    for p in (gp_path, progress_path, checkpoint_manifest_path):
        add_pin(p, pins)
    return {
        "cut": display_cut,
        "baseline_GP_marginal_P": math.exp(float(baseline["marginal_logp"])),
        "GP_domain": gp.get("physical_domain", {}),
        "normalization": gp.get("normalization", {}),
        "GP_cuts": preset_winners,
        "GP_max_saved_P": display_cut,
        "GP_retained_profiles": gp["frontier"],
        "GP_records": records,
        "GP_knots": gp["knots_kpc"],
        "GP_status": gp["status"],
        "GP_diagnostic_only": gp["diagnostic_only"],
        "GP_checkpoint_status": progress.get("status"),
        "GP_checkpoint_completed_ell_over_rd": progress.get("completed_ell_over_rd"),
        "GP_checkpoint_records_so_far": progress.get("records_so_far"),
        "GP_selected_candidates": {
            "max_saved_probability": int(maxp["original_candidate_index"]),
            "minimum_J_ref": int(minj["original_candidate_index"]),
        },
    }


def main():
    table = load_sparc_table(ROOT / "SPARC_table.txt")
    pins = {}
    for p in (ROOT / "SPARC_table.txt", ROOT / "utils.py", ROOT / "mono_SPARC_RAR.py",
              ROOT / "miscellaneous/support_code/prepare_paired_examples_20261021.py",
              ROOT / "miscellaneous/support_code/plot_paired_forces_20261021.py"):
        add_pin(p, pins)
    results = {}
    for galaxy in GALAXIES:
        rad, sb = get_SBdisk(galaxy, dens_root=ROOT / "../SPARC_RAR/BulgeDiskDec_LTG", table=table)
        rad, sb = np.asarray(rad, float), np.asarray(sb, float)
        report = {"families": {}, "rad_kpc": rad.tolist(), "sb": sb.tolist()}
        for family, run_name in FAMILIES.items():
            report["families"][family] = grid_selection(galaxy, family, run_name, rad, sb, pins)
        report.update(gp_selection(galaxy, rad, pins))
        results[galaxy] = report
        for p in (ROOT / "../data" / f"{galaxy}_rotmod.dat",
                  ROOT / "../SPARC_RAR/BulgeDiskDec_LTG" / f"{galaxy}.dens"):
            add_pin(p, pins)
    manifest = {
        "schema_version": 1,
        "purpose": "Exploratory disc-force and RAR figures for UGC04325 and KK98-251",
        "results": results,
        "inputs_sha256": pins,
        "driver_sha256": pins[str((ROOT / "miscellaneous/support_code/plot_paired_forces_20261021.py").resolve())],
        "preparer_sha256": sha(__file__),
        "no_force_probability_search_calls": True,
        "selection_note": "Each grid family uses its highest saved_logP profile because no preset probability cut has candidates across these galaxies and families. GP panels show the maximum saved marginal probability and minimum J_ref candidates from the completed finite search. These cut/profile choices are post-hoc and illustrative, not calibrated significance or optimizer claims. Grid summaries say accuracy_unverified/scientific_validated=false. Both GP files are diagnostic_only; checkpoint progress is partial_search at ell/Rd=4. The newer UGC04325 production_05 folder has no final result/completion artifact and is not used here.",
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + chr(10))
    print(json.dumps({g: {"cut_display": results[g]["cut"], "grid": {f: {"selected_cell": {k: results[g]["families"][f]["selected"]["row"][k] for k in ("grid_i", "grid_j", "saved_logP", "J_rad")}, "valid_cells": results[g]["families"][f]["valid_cells"]} for f in FAMILIES}, "GP": results[g]["GP_selected_candidates"], "checkpoint": results[g]["GP_checkpoint_status"]} for g in GALAXIES}, indent=2))

if __name__ == "__main__":
    main()
