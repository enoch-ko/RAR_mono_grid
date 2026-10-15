#!/usr/bin/env python3
"""Audit saved GP pilot searches and report minimum FOUND J at explicit P levels.
No new probabilities, fits, checkpoint traversal or global-optimum claims.
"""
import argparse
import collections
import gc
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
THRESHOLDS = {"P_4.5_percent": .045,
              "2sigma_one_sided": .5*math.erfc(2/math.sqrt(2)),
              "5sigma_one_sided": .5*math.erfc(5/math.sqrt(2))}


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024*1024), b""):
            digest.update(block)
    return digest.hexdigest()


def minimum(records, threshold):
    best = None
    for index, record in enumerate(records):
        if record.get("valid") is not True:
            continue
        j, lp = record.get("J_ref"), record.get("marginal_logp")
        if (j is None or lp is None or not math.isfinite(j) or j < 0
                or not math.isfinite(lp) or lp > 1e-12 or lp < math.log(threshold)):
            continue
        if best is None or (j, index) < (best["J_ref"], best["record_index"]):
            best = {"record_index": index, "J_ref": j, "P": math.exp(lp),
                    "record": record}
    return best


def audit(path):
    with path.open() as stream:
        data = json.load(stream)
    if data.get("status") != "finite_search_completed":
        raise ValueError(f"incomplete search: {path}")
    records, settings = data["records"], data["settings"]
    domain = data["physical_domain"]
    attempted = len(settings["proposal_ell_over_rd"])*settings["directions_per_ell"]
    if domain["attempted_directions"] != attempted:
        raise ValueError(f"incomplete direction count: {path}")
    fixed = sum(r["source"] in ("SPARC_constant_0.5", "saved_audited_central_order_QP")
                for r in records)
    expected = ((attempted-domain["excluded_directions"])*(settings["amplitude_steps"]-1)
                + fixed)
    retained = data.get("checkpoint_storage", {}).get("mode") == "retained-only"
    counts = data.get("checkpoint_storage", {}).get("summary", {}) if retained else {}
    evaluated = counts.get("evaluated", len(records))
    if retained:
        expected = domain["expected_evaluated_candidates"]
        if (counts["retained"] != len(records) or counts["valid"] + counts["invalid"] != evaluated
                or counts["discarded"] + len(records) != evaluated):
            raise ValueError(f"retention count mismatch: {path}")
    if evaluated != expected:
        raise ValueError(f"record count mismatch: {path}: {evaluated} != {expected}")
    hits = {key: minimum(records, p) for key, p in THRESHOLDS.items()}
    # Independently check the supplied Pareto frontier gives the same minima.
    frontier_records = [records[f["record_index"]] for f in data["frontier"]]
    for key, p in THRESHOLDS.items():
        check = minimum(frontier_records, p)
        if (None if check is None else check["J_ref"]) != (None if hits[key] is None else hits[key]["J_ref"]):
            raise ValueError(f"frontier/full-record disagreement: {path}: {key}")
    invalid = (collections.Counter(counts["failure_counts"]) if retained else
               collections.Counter(r.get("reason", "unspecified") for r in records if not r.get("valid")))
    invalid_numeric = sum(r.get("valid") is True and
                          (r.get("marginal_logp") is None or not math.isfinite(r["marginal_logp"])
                           or r.get("J_ref") is None or not math.isfinite(r["J_ref"]))
                          for r in records)
    qc = [data["baseline_direct_cache_qc"]]
    if data.get("reference_direct_cache_qc"):
        qc.append(data["reference_direct_cache_qc"])
    qc.extend(f["direct_cache_qc"] for f in data.get("retained_profiles", data["frontier"]))
    model = {k: v for k, v in settings.items()
             if k not in ("seed", "directions_per_ell", "amplitude_steps")}
    fingerprint = json.dumps([model, data["provenance"]["script_sha256"],
                              data["provenance"]["source_sha256"],
                              data["provenance"]["input_sha256"]], sort_keys=True)
    valid_probabilities = [r["marginal_logp"] for r in records
                          if r.get("valid") and r.get("marginal_logp") is not None
                          and math.isfinite(r["marginal_logp"])]
    result = {"path": str(path.resolve()), "sha256": sha(path), "galaxy": data["galaxy"],
              "settings": settings, "provenance": data["provenance"],
              "records": evaluated, "retained_records": len(records), "expected_records": expected,
              "retention_mode": data.get("checkpoint_storage", {}).get("mode"),
              "valid_records": counts["valid"] if retained else sum(r.get("valid") is True for r in records),
              "invalid_reasons": dict(invalid), "invalid_numeric_valid_records": invalid_numeric,
              "physical_domain": domain, "basis_warnings": data["basis_warnings"],
              "max_direct_cache_relative_to_peak": max(x["relative_to_peak"] for x in qc),
              "baseline": records[0], "maximum_found_P": math.exp(max(valid_probabilities)),
              "minimum_found_J": hits, "frontier_cross_check_passed": True}
    del data, records, frontier_records
    gc.collect()
    return result, fingerprint


def self_test():
    records = [{"valid": True, "J_ref": 3., "marginal_logp": math.log(.05)},
               {"valid": True, "J_ref": 0., "marginal_logp": math.log(.02)},
               {"valid": False, "J_ref": 0., "marginal_logp": 0.}]
    assert minimum(records, .045)["J_ref"] == 3.
    assert minimum(records, THRESHOLDS["5sigma_one_sided"])["J_ref"] == 0.
    assert minimum(records, .1) is None
    assert minimum(records, .05)["record_index"] == 0
    print("Four threshold-selection checks passed.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, action="append")
    parser.add_argument("--output-root", type=Path,
                        default=ROOT/"outputs/gp_rar_diagnostics/thresholds_20261013")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    if args.output_root.exists():
        raise FileExistsError("choose a new output root")
    roots = args.root or [ROOT/"outputs/gp_rar_budget_20261012",
                          ROOT/"outputs/gp_rar_budget_20261013"]
    results, fingerprints = [], {}
    for root in roots:
        files = sorted(root.glob("*/*_gp_search.json"))
        if not files:
            raise ValueError(f"no result files in {root}")
        for path in files:
            result, fingerprint = audit(path)
            name = result["galaxy"]
            if name in fingerprints and fingerprint != fingerprints[name]:
                raise ValueError(f"incompatible model/source/input for {name}")
            fingerprints[name] = fingerprint
            results.append(result)
            print(name, result["settings"]["directions_per_ell"],
                  {k: None if v is None else v["J_ref"]
                   for k, v in result["minimum_found_J"].items()}, flush=True)
    results.sort(key=lambda r: (r["galaxy"], r["settings"]["directions_per_ell"], r["settings"]["seed"]))
    report = {"thresholds_P_min": THRESHOLDS, "probability": "marginal_logp",
              "selection": "minimum found J among valid saved records, including baseline and saved QP",
              "sigma_labels_are_not_calibrated_significance": True,
              "different_seeds_and_nonnested_direction_sets": True,
              "not_global_optima_or_posterior_samples": True,
              "script_sha256": sha(__file__), "results": results}
    args.output_root.mkdir(parents=True)
    with (args.output_root/"threshold_summary.json").open("x") as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print("SAVED", args.output_root/"threshold_summary.json", "searches", len(results))


if __name__ == "__main__":
    main()
