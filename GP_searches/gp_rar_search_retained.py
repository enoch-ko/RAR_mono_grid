#!/usr/bin/env python3
"""Exploratory uncertainty-aware GP-direction search for disc M/L profiles.

The proposal GP is a search device, not a prior or posterior. All lengthscales
are scored with the same reference covariance metric. A finite frontier is a
best-found diagnostic, not a certified optimum or a probability threshold.
"""
import argparse
import hashlib
import json
import math
import sys
import warnings
from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import scipy
from scipy.linalg import solve_triangular

from gp_candidate_retention import CandidateRetention, DEFAULT_THRESHOLDS, limited_json
import inverse_rar_diagnostic as inverse
import mono_SPARC_RAR as mono
import probability_recursion as boundary
from utils import DEFAULT_DENS_ROOT, DEFAULT_PROJECT_ROOT, DEFAULT_TABLE_PATH
from utils import _tie_groups, load_sparc_table, lognormal_quadrature, vel2acc

PILOTS = ("NGC3521", "UGC02916", "UGC07399", "D512-2")
SOURCE_FILES = ("inverse_rar_diagnostic.py", "mono_SPARC_RAR.py",
                "utils.py", "probability_recursion.py", "gp_candidate_retention.py")


def nested_direction_rng(master_seed, galaxy, ell_over_rd, replicate=0):
    """Independent per-galaxy/ell stream; larger budgets retain every prefix."""
    payload = json.dumps(["rar-gp-directions-v1", int(master_seed), str(galaxy),
                          float(ell_over_rd).hex(), int(replicate)], separators=(",", ":"))
    digest = hashlib.sha256(payload.encode("utf-8")).digest()
    return np.random.default_rng(int.from_bytes(digest[:16], "big"))


def source_hashes():
    root = Path(__file__).resolve().parent
    return {name: inverse.sha256(root / name) for name in SOURCE_FILES}


def proposal_cholesky(knots, rd_kpc, ell_over_rd, nugget=1e-8):
    """The proposal ell/Rd is physical; inverse.kernel uses ell/knot span."""
    if not (np.isfinite(rd_kpc) and rd_kpc > 0 and
            np.isfinite(ell_over_rd) and ell_over_rd > 0):
        raise ValueError("positive finite Rd and proposal lengthscale required")
    ell_kpc = rd_kpc * ell_over_rd
    distance = (knots[:, None] - knots[None, :]) / ell_kpc
    covariance = np.exp(-0.5 * distance**2) + nugget * np.eye(len(knots))
    return np.linalg.cholesky(covariance), ell_kpc


def physical_interval(mean, direction, lower, upper, disk_response,
                      nuisance_terms, floor):
    """Exact nonnegative-amplitude interval from knot bounds and signed gbar."""
    if not (0 < lower < 0.5 < upper and np.isfinite(floor) and floor > 0):
        raise ValueError("invalid bounds or signed acceleration floor")
    q0 = np.concatenate(((mean - lower), (upper - mean),
                         (nuisance_terms + disk_response @ mean - floor).ravel()))
    q1 = np.concatenate((direction, -direction,
                         np.broadcast_to(disk_response @ direction,
                                         nuisance_terms.shape).ravel()))
    lo, hi = 0.0, np.inf
    for intercept, slope in zip(q0, q1):
        if not np.isfinite(intercept) or not np.isfinite(slope):
            return None
        if slope > 0:
            lo = max(lo, -intercept / slope)
        elif slope < 0:
            hi = min(hi, -intercept / slope)
        elif intercept < 0:
            return None
    if not np.isfinite(hi) or hi < lo or hi < 0:
        return None
    return max(0.0, lo), hi


def amplitude_grid(interval, count):
    if interval is None:
        return []
    lo, hi = interval
    if count < 2:
        raise ValueError("at least two amplitude steps required")
    # Include the admissible edges and resolve near-reference amplitudes.
    fractions = np.geomspace(1e-3, 1.0, count - 1)
    values = np.r_[lo, lo + (hi - lo) * fractions]
    return np.unique(values).tolist()


def frontier_indices(records):
    """Every retained point raises logP as fixed J increases."""
    eligible = [i for i, row in enumerate(records) if row["valid"]]
    eligible.sort(key=lambda i: (records[i]["J_ref"], -records[i]["marginal_logp"]))
    best, kept = -np.inf, []
    for i in eligible:
        if records[i]["marginal_logp"] > best + 1e-12:
            kept.append(i)
            best = records[i]["marginal_logp"]
    return kept


class ProbabilityEvaluator:
    """Reuse the existing uncertainty model and memoise identical orderings."""
    def __init__(self, D, r, vobs, errv, vgas, vbul, args, gas_nodes,
                 gas_weights, bulge_nodes, bulge_weights, nuisance_terms, floor):
        self.D, self.r = D, r
        self.vobs, self.errv, self.vgas, self.vbul = vobs, errv, vgas, vbul
        self.args = args
        self.gas_nodes, self.gas_weights = gas_nodes, gas_weights
        self.bulge_nodes, self.bulge_weights = bulge_nodes, bulge_weights
        self.nuisance_terms, self.floor = nuisance_terms, floor
        self.memo = {}
        self.cache_hits = 0

    def _ordered_logp(self, gbar, *values, **options):
        groups = _tie_groups(np.asarray(gbar), options.get("tie_atol", 0.0))
        multiplicity = math.prod(math.factorial(len(group)) for group in groups)
        if multiplicity > options.get("max_tie_permutations", 720):
            raise ValueError("tied-order cap exceeded; no stable-order fallback")
        key = tuple(tuple(int(j) for j in group) for group in groups)
        if key in self.memo:
            self.cache_hits += 1
            return self.memo[key]
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            result = boundary.boundary_log_probability_for_gbar(
                gbar, *values, log_tolerance=self.args.boundary_log_tolerance,
                max_refinements=self.args.boundary_max_refinements, **options)
        if caught:
            for warning in caught:
                warnings.warn(str(warning.message), warning.category)
        else:
            self.memo[key] = result
        return result

    def evaluate(self, xi):
        signed = self.D @ xi
        all_gbar = self.nuisance_terms + signed
        if not np.all(np.isfinite(all_gbar)):
            return {"valid": False, "reason": "nonfinite_signed_gbar"}
        if np.min(all_gbar) < self.floor:
            node, radius = np.unravel_index(np.argmin(all_gbar), all_gbar.shape)
            return {"valid": False, "reason": "signed_gbar_below_floor",
                    "worst_nuisance_node_index": int(node),
                    "worst_radius_index": int(radius),
                    "minimum_nuisance_gbar_m_s2": float(all_gbar[node, radius]),
                    "floor_m_s2": float(self.floor)}
        # evaluate_disk_term expects signed W=Rg, in (km/s)^2 units.
        disk_term = signed / (1e6 / (self.r * 3.086e19))
        original = mono.log_monotonic_probability_for_gbar
        mono.log_monotonic_probability_for_gbar = self._ordered_logp
        try:
            conditional, marginal = mono.evaluate_disk_term(
                disk_term, self.vgas, self.vbul, self.r, self.vobs, self.errv,
                self.args, self.gas_nodes, self.gas_weights,
                self.bulge_nodes, self.bulge_weights)
        finally:
            mono.log_monotonic_probability_for_gbar = original
        if not conditional[1] or not marginal[1]:
            return {"valid": False, "reason": "conditional=" + str(conditional[2]) +
                    "; marginal=" + str(marginal[2])}
        messages = list(conditional[3]) + list(marginal[3])
        if messages:
            return {"valid": False, "reason": "probability_warning: " + " | ".join(messages)}
        return {"valid": True, "conditional_logp": float(conditional[0]),
                "marginal_logp": float(marginal[0]),
                "minimum_nuisance_gbar_m_s2": float(np.min(all_gbar))}


def warm_reference_profile(path, data, galaxy, knots, lower, upper, args, input_hashes):
    """Import only the accepted warm fit from the targeted OSQP comparison."""
    if (data.get("galaxy") != galaxy or data.get("status") != "completed_fit_comparison"
            or data.get("input_sha256") != input_hashes or
            args.length_fraction != 0.3 or args.gas_ml != 1.0 or args.bulge_ml != 0.7 or
            args.gbar_floor_relative != 1e-8):
        raise ValueError("warm reference data or QP model/settings differ")
    source = data.get("source_sha256", {})
    root = Path(__file__).resolve().parent
    if (not source or any(not (root / name).is_file() or
                          inverse.sha256(root / name) != digest
                          for name, digest in source.items())):
        raise ValueError("warm reference source hashes differ")
    accepted = [run for run in data.get("runs", [])
                if run.get("mode") == "warm_slsqp" and
                run.get("fit", {}).get("status") == "feasible_weak_order"]
    if len(accepted) != 1:
        raise ValueError("warm reference has no unique accepted fit")
    run = accepted[0]
    fit = run["fit"]
    if (fit.get("qp_backend") != "osqp" or
            fit.get("constraint_version") != "positive_net_gbar_v1" or
            fit.get("gbar_floor_relative") != args.gbar_floor_relative or
            not fit.get("optimality", {}).get("passed", False) or
            not np.isfinite(fit.get("max_scaled_violation", np.nan)) or
            fit["max_scaled_violation"] > 1e-8 or
            not np.isfinite(fit.get("max_bound_violation", np.nan)) or
            fit["max_bound_violation"] > 1e-8 or
            not np.isfinite(run.get("direct_cache_qc", {}).get("relative_to_peak", np.nan)) or
            run["direct_cache_qc"]["relative_to_peak"] > 1e-8):
        raise ValueError("warm reference saved acceptance audit failed")
    xi = np.asarray(fit.get("xi"), float)
    if (xi.shape != knots.shape or not np.all(np.isfinite(xi)) or
            np.min(xi) < lower - 1e-8 or np.max(xi) > upper + 1e-8):
        raise ValueError("warm reference profile outside current knot geometry/bounds")
    return xi, {"status": "saved_warm_OSQP_audits_passed_current_sources_checked",
                "path": str(path), "sha256": inverse.sha256(path),
                "format": "targeted_OSQP_warm_comparison", "mode": run["mode"],
                "saved_objective": fit.get("objective"),
                "saved_response_sha256": data.get("response_sha256")}


def audited_finite_le(value, limit):
    try:
        return bool(np.isfinite(value) and value <= limit)
    except (TypeError, ValueError):
        return False


def reference_profile(path, galaxy, knots, lower, upper, args, input_hashes):
    """Import only a saved, audited central-order fit with matching model/data."""
    if path is None:
        return None, {"status": "not_supplied"}
    data = json.loads(Path(path).read_text())
    if isinstance(data.get("runs"), list):
        return warm_reference_profile(path, data, galaxy, knots, lower, upper, args, input_hashes)
    fit = data.get("optimization", {})
    provenance = data.get("provenance", {})
    audit = fit.get("optimality", {})
    if data.get("galaxy") != galaxy or fit.get("status") != "feasible_weak_order":
        raise ValueError("reference JSON is not an accepted fit for this galaxy")
    if data.get("input_sha256") != input_hashes:
        raise ValueError("reference inputs differ from current SPARC data")
    if not np.array_equal(np.asarray(data.get("knots_kpc"), float), knots):
        raise ValueError("reference knot geometry differs")
    if (data.get("bounds_ml") != [lower, upper] or
            data.get("kernel", {}).get("length_fraction") != args.length_fraction or
            data.get("kernel", {}).get("nugget") != 1e-8 or
            data.get("gas_ml_qp") != args.gas_ml or
            data.get("bulge_ml_qp") != args.bulge_ml or
            fit.get("constraint_version") != "positive_net_gbar_v1" or
            fit.get("gbar_floor_relative") != args.gbar_floor_relative):
        raise ValueError("reference QP model/settings differ")
    if (not audit.get("passed", False) or
            not audited_finite_le(fit.get("max_scaled_violation"), 1e-8) or
            not audited_finite_le(fit.get("max_bound_violation"), 1e-8) or
            not audited_finite_le(data.get("optimized_direct_cache_qc", {}).get("relative_to_peak"), 1e-8)):
        raise ValueError("reference QP saved physical/optimality audit failed")
    if (provenance.get("driver_sha256") != inverse.sha256(mono.__file__) or
            provenance.get("utils_sha256") != inverse.sha256(Path(DEFAULT_PROJECT_ROOT) / "utils.py")):
        raise ValueError("reference physics/probability source hashes differ")
    xi = np.asarray(fit.get("xi"), float)
    if (xi.shape != knots.shape or not np.all(np.isfinite(xi)) or
            np.min(xi) < lower - 1e-8 or np.max(xi) > upper + 1e-8):
        raise ValueError("reference profile is outside current physical bounds")
    return xi, {"status": "saved_audits_passed_current_model_checked",
                "path": str(path), "sha256": inverse.sha256(path),
                "saved_inverse_script_sha256": provenance.get("script_sha256"),
                "current_inverse_script_sha256": inverse.sha256(inverse.__file__)}


def galaxy_run(galaxy, args, reference_path=None):
    sampling_scheme = getattr(args, "sampling_scheme", "legacy-pilot")
    replicate = getattr(args, "replicate", 0)
    checkpoint_mode = getattr(args, "checkpoint_mode", "per-candidate")
    initial_hashes = source_hashes()
    initial_script_hash = inverse.sha256(__file__)
    gate = inverse.cache_gate(galaxy, args.validation_root)
    input_hashes = {"rotmod": inverse.sha256(Path(args.data_root) / f"{galaxy}_rotmod.dat"),
                    "dens": inverse.sha256(Path(args.dens_root) / f"{galaxy}.dens"),
                    "table": inverse.sha256(DEFAULT_TABLE_PATH)}
    rad, brightness, r, vobs, errv, vgas, vbul, zd, excluded = inverse.galaxy_data(
        galaxy, args.data_root, args.dens_root)
    row = load_sparc_table().loc[lambda t: t.Galaxy == galaxy]
    if len(row) != 1:
        raise ValueError("catalogue Rd unavailable")
    rd = float(row.iloc[0]["Rdisk"])
    knots = np.asarray(rad, float)
    H, D, force_warnings = inverse.build_response(rad, brightness, r, zd, knots)
    mean = np.full(len(knots), 0.5)
    baseline_qc = inverse.direct_check(rad, brightness, r, zd, H, D, mean)
    Lref = np.linalg.cholesky(inverse.kernel(knots, args.length_fraction))
    gas_nodes, gas_weights = lognormal_quadrature(
        args.gas_ml, args.gas_ml_dex, args.ml_samples)
    bulge_nodes, bulge_weights = lognormal_quadrature(
        args.bulge_ml, args.bulge_ml_dex, args.ml_samples)
    nuisance = np.array([
        vel2acc(vgas * np.abs(vgas) * gas + vbul**2 * bulge, r)
        for gas in gas_nodes for bulge in bulge_nodes])
    central = vel2acc(vgas * np.abs(vgas) * args.gas_ml +
                      vbul**2 * args.bulge_ml, r)
    scale = max(float(np.max(np.abs(central + D @ mean))), 1e-30)
    floor = args.gbar_floor_relative * scale
    pargs = SimpleNamespace(gas_ml=args.gas_ml, bulge_ml=args.bulge_ml,
                            marginalise_gas_bulge=True,
                            prob_grid_size=args.prob_grid_size,
                            zpadding=args.zpadding, tie_atol=0.0,
                            max_tie_permutations=720,
                            boundary_log_tolerance=args.boundary_log_tolerance,
                            boundary_max_refinements=args.boundary_max_refinements)
    evaluator = ProbabilityEvaluator(D, r, vobs, errv, vgas, vbul, pargs,
                                     gas_nodes, gas_weights, bulge_nodes,
                                     bulge_weights, nuisance, floor)
    ref, reference_info = reference_profile(
        reference_path, galaxy, knots, args.min_ml, args.max_ml, args, input_hashes)
    thresholds = dict(DEFAULT_THRESHOLDS)
    thresholds.update({f"P_{p:.17g}": p for p in getattr(args, "retain_probability", [])})
    retention = CandidateRetention(thresholds)
    byte_limit = getattr(args, "max_retained_bytes", 8 * 1024**2)

    def retained_snapshot():
        rows, vectors, _, summary = retention.export()
        return {"records": rows, "profiles": [x.tolist() for x in vectors], "summary": summary}
    checkpoint_root = args.output_root / f"{galaxy}_gp_checkpoints"
    checkpoint_root.mkdir(parents=True, exist_ok=False)
    with (checkpoint_root / "manifest.json").open("x") as stream:
        json.dump({"galaxy": galaxy, "seed": (args.seed + PILOTS.index(galaxy)
                   if sampling_scheme == "legacy-pilot" else args.seed),
                   "sampling_scheme": sampling_scheme, "replicate": replicate,
                   "knots_kpc": knots.tolist(), "rd_kpc": rd,
                   "script_sha256": initial_script_hash,
                   "source_sha256": initial_hashes, "input_sha256": input_hashes,
                   "reference": reference_info,
                   "args": {key: str(value) for key, value in vars(args).items()}},
                  stream, indent=2, allow_nan=False)
        stream.write("\n")

    def add(xi, source, ell=None, direction=None, amplitude=None):
        xi = np.asarray(xi, float)
        z = solve_triangular(Lref, xi - mean, lower=True)
        entry = {"source": source, "ell_over_rd": ell, "direction": direction,
                 "amplitude_ml": amplitude, "J_ref": float(z @ z)}
        entry.update(evaluator.evaluate(xi))
        retention.add(entry, xi)
        if retention.total % 1000 == 0:
            limited_json(inverse.json_safe(retained_snapshot()), byte_limit)

    add(mean, "SPARC_constant_0.5")
    reference_qc = None
    if ref is not None:
        reference_qc = inverse.direct_check(rad, brightness, r, zd, H, D, ref)
        add(ref, "saved_audited_central_order_QP")
    galaxy_seed = args.seed + PILOTS.index(galaxy) if sampling_scheme == "legacy-pilot" else args.seed
    rng = np.random.default_rng(galaxy_seed)
    directions_rejected = 0
    expected_candidate_count = 1 + int(ref is not None)
    for ell in args.ell_over_rd:
        if sampling_scheme == "nested-v1":
            rng = nested_direction_rng(args.seed, galaxy, ell, replicate)
        Lprop, _ = proposal_cholesky(knots, rd, ell)
        for direction_id in range(args.directions):
            direction = Lprop @ rng.standard_normal(len(knots))
            direction /= np.max(np.abs(direction))
            interval = physical_interval(mean, direction, args.min_ml, args.max_ml,
                                         D, nuisance, floor)
            if interval is None:
                directions_rejected += 1
            amplitudes = amplitude_grid(interval, args.amplitude_steps)
            expected_candidate_count += sum(amplitude != 0 for amplitude in amplitudes)
            for amplitude in amplitudes:
                if amplitude == 0:
                    continue  # the exact reference was evaluated once above
                add(mean + amplitude * direction, "GP_direction",
                    ell, direction_id, amplitude)
        progress = {"galaxy": galaxy, "completed_ell_over_rd": ell,
                    "records_so_far": retention.total, "status": "partial_search",
                    "retained_snapshot": retained_snapshot(),
                    "created_utc": datetime.now(timezone.utc).isoformat()}
        text = limited_json(inverse.json_safe(progress), byte_limit)
        with (checkpoint_root / "progress.json").open("w") as stream:
            stream.write(text)
    if retention.total != expected_candidate_count:
        raise RuntimeError("evaluated candidate count differs from actual amplitude schedule")
    records, profiles, frontier_ids, retention_summary = retention.export()
    retained_profiles = []
    for index, xi in enumerate(profiles):
        qc = inverse.direct_check(rad, brightness, r, zd, H, D, xi)
        retained_profiles.append({"record_index": index, "xi": xi.tolist(),
                                  "profile_on_photometric_grid": (H @ xi).tolist(),
                                  "direct_cache_qc": qc})
    frontier = [retained_profiles[index] for index in frontier_ids]
    storage = {"mode": "retained-only", "schema_version": 1,
               "profiles": "inline retained_profiles indexed by saved records",
               "all_candidates_saved": False, "summary": retention_summary,
               "max_retained_json_bytes": byte_limit, "selection_probability": "marginal_logp",
               "sigma_labels_are_not_calibrated_significance": True,
               "discard_policy": "exact dominance; baseline, threshold minima and maximum P pinned",
               "partial_progress": "progress.json contains per-ell retained snapshots; not automatically resumable"}
    current_inputs = {"rotmod": inverse.sha256(Path(args.data_root) / f"{galaxy}_rotmod.dat"),
                      "dens": inverse.sha256(Path(args.dens_root) / f"{galaxy}.dens"),
                      "table": inverse.sha256(DEFAULT_TABLE_PATH)}
    if (source_hashes() != initial_hashes or
            inverse.sha256(__file__) != initial_script_hash or
            current_inputs != input_hashes or
            (reference_path is not None and
             inverse.sha256(reference_path) != reference_info["sha256"])):
        raise RuntimeError("source, input, or reference changed during galaxy search")
    return {"galaxy": galaxy, "diagnostic_only": True,
            "status": "finite_search_completed" if frontier else "no_valid_candidate",
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "settings": {"seed": args.seed, "directions_per_ell": args.directions,
                         "amplitude_steps": args.amplitude_steps,
                         "proposal_ell_over_rd": args.ell_over_rd,
                         "reference_length_fraction_of_span": args.length_fraction,
                         "reference_nugget": 1e-8,
                         "bounds_ml": [args.min_ml, args.max_ml],
                         "signed_gbar_floor_relative": args.gbar_floor_relative,
                         "gas_ml": args.gas_ml, "gas_ml_dex": args.gas_ml_dex,
                         "bulge_ml": args.bulge_ml, "bulge_ml_dex": args.bulge_ml_dex,
                         "quadrature_nodes_per_term": args.ml_samples,
                         "probability_starting_grid": args.prob_grid_size,
                         "zpadding": args.zpadding,
                         "boundary_log_tolerance": args.boundary_log_tolerance,
                         "boundary_max_refinements": args.boundary_max_refinements},
            "sampling": {"scheme": sampling_scheme, "master_seed": args.seed,
                         "galaxy_seed": galaxy_seed, "replicate": replicate,
                         "nested_across_direction_budgets": sampling_scheme == "nested-v1",
                         "seed_recipe": "sha256(version,master,galaxy,float.hex(ell),replicate) first128bits"
                         if sampling_scheme == "nested-v1" else "legacy seed + pilot index"},
            "checkpoint_storage": storage,
            "retention_schema_version": 1, "retained_profiles": retained_profiles,
            "rd_kpc": rd, "photometric_span_kpc": float(knots[-1] - knots[0]),
            "knots_kpc": knots.tolist(), "retained_radii_kpc": r.tolist(),
            "excluded_kinematic_rows": excluded, "cache_gate": gate,
            "basis_warnings": force_warnings,
            "baseline_direct_cache_qc": baseline_qc,
            "reference_direct_cache_qc": reference_qc,
            "physical_domain": {"policy": "signed net gbar at every gas/bulge quadrature node must meet the existing numerical floor",
                                "floor_m_s2": floor,
                                "excluded_directions": directions_rejected,
                                "attempted_directions": len(args.ell_over_rd) * args.directions,
                                "expected_evaluated_candidates": expected_candidate_count,
                                "strictly_more_restrictive_than_fixed_nuisance_QP": True},
            "probability_cache": {"distinct_warning_free_orders": len(evaluator.memo),
                                  "hits": evaluator.cache_hits},
            "references": {"central_QP": reference_info},
            "checkpoint_directory": str(checkpoint_root),
            "records": records, "frontier": frontier,
            "provenance": {"script_sha256": initial_script_hash,
                           "script_file": "gp_rar_search_retained.py",
                           "numerical_parent_script_sha256": "50ae153f2689a511c62bfa89f33410b41ff1bbab6a9b491491f2a8084bfa8934",
                           "source_sha256": initial_hashes,
                           "input_sha256": input_hashes,
                           "python": sys.version, "python_executable": sys.executable,
                           "numpy": np.__version__, "scipy": scipy.__version__}}


def main():
    parser = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    parser.add_argument("--galaxy", action="append",
                        help="Repeat for a subset; defaults to all four pilots")
    parser.add_argument("--output-root", type=Path,
                        default=Path(DEFAULT_PROJECT_ROOT) / "outputs" / "gp_rar_retained_manual")
    parser.add_argument("--data-root", type=Path, default=mono.DEFAULT_DATA_ROOT)
    parser.add_argument("--dens-root", type=Path, default=Path(DEFAULT_DENS_ROOT))
    parser.add_argument("--validation-root", type=Path,
                        default=Path(DEFAULT_PROJECT_ROOT) / "outputs")
    parser.add_argument("--reference-json", action="append", default=[],
                        help="Accepted inverse JSON as GALAXY=PATH")
    parser.add_argument("--length-fraction", type=float, default=0.3)
    parser.add_argument("--ell-over-rd", type=float, nargs="+", default=[0.5, 1.0, 2.0, 4.0])
    parser.add_argument("--min-ml", type=float, default=0.01)
    parser.add_argument("--max-ml", type=float, default=10.0)
    parser.add_argument("--gbar-floor-relative", type=float, default=1e-8)
    parser.add_argument("--seed", type=int, default=20261012)
    parser.add_argument("--sampling-scheme", choices=["legacy-pilot", "nested-v1"], default="legacy-pilot")
    parser.add_argument("--replicate", type=int, default=0)
    parser.add_argument("--checkpoint-mode", choices=["retained-only"], default="retained-only")
    parser.add_argument("--max-retained-bytes", type=int, default=8 * 1024**2,
                        help="Per-galaxy retained snapshot/result JSON cap; fail rather than truncate frontier")
    parser.add_argument("--retain-probability", action="append", type=float, default=[],
                        help="Additional probability cuts to pin minimum-J candidates; defaults include one-sided 2/5 labels")
    parser.add_argument("--directions", type=int, default=8)
    parser.add_argument("--amplitude-steps", type=int, default=6)
    parser.add_argument("--gas-ml", type=float, default=1.0)
    parser.add_argument("--bulge-ml", type=float, default=0.7)
    parser.add_argument("--gas-ml-dex", type=float, default=0.04)
    parser.add_argument("--bulge-ml-dex", type=float, default=0.1)
    parser.add_argument("--ml-samples", type=int, default=9)
    parser.add_argument("--prob-grid-size", type=int, default=400)
    parser.add_argument("--zpadding", type=float, default=8.0)
    parser.add_argument("--boundary-log-tolerance", type=float, default=2e-5)
    parser.add_argument("--boundary-max-refinements", type=int, default=8)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    if args.max_retained_bytes < 1024:
        parser.error("retained byte limit must be at least 1024")
    if any(not math.isfinite(p) or not 0 < p <= 1 for p in args.retain_probability):
        parser.error("retention probabilities must be finite in (0,1]")
    galaxies = tuple(args.galaxy or PILOTS)
    catalogue = set(load_sparc_table()["Galaxy"])
    if any(g not in catalogue for g in galaxies):
        parser.error("unknown SPARC galaxy")
    if args.sampling_scheme == "legacy-pilot" and any(g not in PILOTS for g in galaxies):
        parser.error("non-pilot galaxies require --sampling-scheme nested-v1")
    if args.seed < 0 or args.replicate < 0 or len(set(args.ell_over_rd)) != len(args.ell_over_rd):
        parser.error("nonnegative seeds/replicates and distinct proposal lengthscales required")
    if (len(set(galaxies)) != len(galaxies) or args.directions < 1 or
            args.amplitude_steps < 2 or args.ml_samples < 1 or
            args.prob_grid_size < 2 or args.boundary_max_refinements < 1 or
            not 0 < args.min_ml < 0.5 < args.max_ml or
            not 0 < args.length_fraction <= 2 or
            not np.all(np.isfinite(args.ell_over_rd)) or
            min(args.ell_over_rd) <= 0 or
            not np.isfinite(args.gbar_floor_relative) or
            args.gbar_floor_relative <= 0 or
            not np.isfinite(args.zpadding) or args.zpadding <= 0 or
            not np.isfinite(args.boundary_log_tolerance) or
            args.boundary_log_tolerance <= 0):
        parser.error("invalid search, metric, bounds or probability settings")
    references = {}
    for item in args.reference_json:
        galaxy, sep, path = item.partition("=")
        if not sep or galaxy not in galaxies or galaxy in references or not path:
            parser.error("reference must be unique GALAXY=PATH for a selected galaxy")
        references[galaxy] = Path(path)
    paths = [args.output_root / f"{galaxy}_gp_search.json" for galaxy in galaxies]
    checkpoint_paths = [args.output_root / f"{galaxy}_gp_checkpoints"
                        for galaxy in galaxies]
    if any(path.exists() for path in paths + checkpoint_paths):
        raise FileExistsError("one or more pilot outputs/checkpoint roots already exist")
    if args.dry_run:
        print(json.dumps({"galaxies": galaxies, "outputs": [str(p) for p in paths],
                          "directions_per_ell": args.directions,
                          "amplitude_steps": args.amplitude_steps,
                          "ell_over_rd": args.ell_over_rd,
                          "reference_metric_length_fraction": args.length_fraction,
                          "max_candidate_count_per_galaxy": 1 + len(args.ell_over_rd) *
                          args.directions * args.amplitude_steps + 1,
                          "sampling_scheme": args.sampling_scheme, "replicate": args.replicate,
                          "checkpoint_mode": args.checkpoint_mode, "queued_jobs": 0,
                          "references": {galaxy: str(path) for galaxy, path in references.items()}},
                         indent=2))
        return
    for galaxy, path in zip(galaxies, paths):
        result = galaxy_run(galaxy, args, references.get(galaxy))
        path.parent.mkdir(parents=True, exist_ok=True)
        text = limited_json(inverse.json_safe(result), args.max_retained_bytes, indent=2)
        with path.open("x") as stream:
            stream.write(text)
        print(galaxy, result["checkpoint_storage"]["summary"]["evaluated"], "evaluated;",
              len(result["records"]), "retained;",
              len(result["frontier"]), "frontier entries;", path, flush=True)


if __name__ == "__main__":
    main()
