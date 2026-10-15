#!/usr/bin/env python3
"""M/L clouds: audited central-order placeholders or threshold-selected GP profiles.
Raw external radial values are connected for display, never refitted/extrapolated.
Radii are kpc or R/R_d; SPARC and external K_s are not bandpass-converted.
"""
import argparse
import collections
import csv
import hashlib
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

ROOT = Path(__file__).resolve().parent
COLORS = {"XMM": "#0077bb", "COSMOS_kin": "#ee7733",
          "COSMOS_DR1": "#cc3311", "LADUMA": "#009988"}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def valid_curve(radius, ml):
    radius, ml = np.asarray(radius, float), np.asarray(ml, float)
    if (radius.ndim != 1 or radius.shape != ml.shape or len(radius) < 1
            or not np.all(np.isfinite(radius)) or np.any(radius < 0)
            or not np.all(np.isfinite(ml)) or np.any(ml <= 0)):
        raise ValueError("invalid radius/M-L profile")
    order = np.argsort(radius, kind="stable")
    return radius[order].tolist(), ml[order].tolist()


def external_profiles(path, sample_column):
    groups, excluded = collections.defaultdict(list), []
    pairs = collections.Counter()
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        required = {"Galaxy", sample_column, "Radius_kpc", "ml_K"}
        if not required <= set(reader.fieldnames or []):
            raise ValueError(f"missing CSV columns: {required-set(reader.fieldnames or [])}")
        rows = list(reader)
    for rownum, row in enumerate(rows, start=2):
        pairs[(row.get("Sample", ""), row.get("sample", ""))] += 1
        name, sample = row["Galaxy"].strip(), row[sample_column].strip()
        try:
            r, ml = float(row["Radius_kpc"]), float(row["ml_K"])
            if not name or not sample:
                raise ValueError("missing galaxy or sample")
            valid_curve([r], [ml])
        except (ValueError, TypeError) as error:
            excluded.append({"csv_line": rownum, "galaxy": name, "reason": str(error)})
            continue
        groups[(sample, name)].append((r, ml))
    curves = []
    for (sample, name), points in sorted(groups.items()):
        r, ml = valid_curve(*zip(*points))
        curves.append({"galaxy": name, "dataset": sample, "radius_kpc": r, "ml": ml,
                       "duplicate_radius_count": len(r)-len(set(r))})
    return curves, {"source": str(path.resolve()), "sha256": sha(path),
                    "raw_rows": len(rows), "raw_unique_galaxies": len({x["Galaxy"] for x in rows}),
                    "sample_column": sample_column, "excluded_rows": excluded,
                    "sample_column_pairs": [{"Sample": a, "sample": b, "rows": n}
                                            for (a, b), n in sorted(pairs.items())],
                    "field_galaxy_counts": dict(collections.Counter(c["dataset"] for c in curves)),
                    "no_refitting_or_bandpass_conversion": True}


def finite_le(value, limit):
    return value is not None and math.isfinite(float(value)) and float(value) <= limit


def central_profile(data):
    opt = data.get("optimization", {})
    if data.get("status") != "feasible_weak_order" or opt.get("status") != "feasible_weak_order":
        raise ValueError(data.get("status", "missing status"))
    if opt.get("optimality", {}).get("passed") is not True:
        raise ValueError("saved optimality audit did not pass")
    if not all(finite_le(opt.get(k), 1e-8)
               for k in ("max_scaled_violation", "max_bound_violation")):
        raise ValueError("saved feasibility audit missing or failed")
    if not finite_le(data.get("optimized_direct_cache_qc", {}).get("relative_to_peak"), 1e-8):
        raise ValueError("saved direct/cache audit missing or failed")
    r, ml = valid_curve(data["knots_kpc"], opt["xi"])
    if len(r) < 2 or len(set(r)) != len(r):
        raise ValueError("invalid knot grid")
    lo, hi = data["bounds_ml"]
    if min(ml) < lo-1e-8 or max(ml) > hi+1e-8:
        raise ValueError("saved knot values violate bounds")
    gbar = np.asarray(opt.get("gbar_m_s2", []), float)
    if gbar.size == 0 or not np.all(np.isfinite(gbar)) or np.any(gbar <= 0):
        raise ValueError("saved net acceleration missing or nonpositive")
    return {"galaxy": data["galaxy"], "radius_kpc": r, "ml": ml,
            "J_ref": opt["objective"], "active_ties": bool(opt.get("boundary_ties_present")),
            "selection": "audited weak-central-order placeholder; NOT probability-threshold selected"}


def load_central(folder):
    curves, rejected, provenance = [], [], []
    files = sorted(folder.glob("*_inverse.json"))
    if not files:
        raise ValueError("no central result files found")
    names = set()
    for path in files:
        data = json.loads(path.read_text())
        name = data.get("galaxy", path.stem)
        if name in names:
            raise ValueError(f"duplicate central galaxy: {name}")
        names.add(name)
        provenance.append({"path": str(path.resolve()), "sha256": sha(path),
                           "galaxy": name, "status": data.get("status"),
                           "saved_provenance": data.get("provenance"),
                           "input_sha256": data.get("input_sha256")})
        try:
            curves.append(central_profile(data))
        except (ValueError, KeyError, TypeError) as error:
            rejected.append({"galaxy": name, "reason": str(error)})
    return {"Central-order placeholder": curves}, {"source_files": provenance,
            "n_files": len(files), "accepted": len(curves), "excluded": rejected,
            "no_retry_substitution": True, "threshold_filter_applied": False}


def threshold_probability(sigma, convention):
    if not math.isfinite(sigma) or sigma <= 0:
        raise ValueError("sigma label must be positive")
    return math.erfc(sigma/math.sqrt(2)) * (.5 if convention == "one-sided" else 1.)


def choose_record(records, probability):
    eligible = [(float(r["J_ref"]), i) for i, r in enumerate(records)
                if r.get("valid") and r.get("marginal_logp") is not None
                and math.isfinite(float(r["marginal_logp"]))
                and r.get("J_ref") is not None and math.isfinite(float(r["J_ref"]))
                and float(r["J_ref"]) >= 0 and float(r["marginal_logp"]) >= math.log(probability)]
    return min(eligible) if eligible else None


def load_gp(roots, sigmas, convention, read_mode="full"):
    """Select minimum FOUND J; no fitting, probability calculation or global claim."""
    panels = {f"GP minimum found J: {s:g}-sigma label": {} for s in sigmas}
    fingerprints, audited, coverage = {}, [], collections.defaultdict(set)
    files = sorted({p.resolve() for root in roots for p in root.rglob("*_gp_search.json")})
    if not files:
        raise ValueError("no completed GP search files")
    for file_number, path in enumerate(files, start=1):
        reader_hash = None
        if read_mode == "frontier-stream":
            import gp_frontier_reader as reader
            reader_hash = sha(reader.__file__)
            try:
                data, result_hash = reader.read_frontier_result(path)
            except reader.UnsupportedFrontierLayout:
                data = json.loads(path.read_text())
                result_hash = sha(path)
            cuts = [math.log(threshold_probability(s, convention)) for s in sigmas]
            if data.get("selection_reader") and any(
                    abs(float(r["marginal_logp"])-cut) <= 1e-12
                    for r in data["records"] if r.get("valid") for cut in cuts):
                data = json.loads(path.read_text())
                result_hash = sha(path)
                data["selection_reader"] = {"mode": "full_near_cut_fallback"}
            print(f"Read saved GP result {file_number}/{len(files)}", flush=True)
        else:
            data = json.loads(path.read_text())
            result_hash = sha(path)

        if data.get("status") != "finite_search_completed":
            raise ValueError(f"incomplete GP result: {path}")
        galaxy = data["galaxy"]
        model = {k:v for k,v in data["settings"].items()
                 if k not in ("seed", "directions_per_ell", "amplitude_steps")}
        fingerprint = json.dumps([model, data["provenance"]["source_sha256"],
                                  data["provenance"]["script_sha256"],
                                  data["provenance"]["input_sha256"]], sort_keys=True)
        if galaxy in fingerprints and fingerprints[galaxy] != fingerprint:
            raise ValueError(f"incompatible model/source/input files for {galaxy}")
        fingerprints[galaxy] = fingerprint
        audited.append({"path": str(path), "sha256": result_hash, "settings": data["settings"],
                        "selection_reader": data.get("selection_reader", {"mode": "full"}),
                        "reader_source_sha256": reader_hash,
                        "physical_domain": data["physical_domain"], "provenance": data["provenance"]})
        for sigma, (label, selected) in zip(sigmas, panels.items()):
            coverage[label].add(galaxy)
            hit = choose_record(data["records"], threshold_probability(sigma, convention))
            if hit is None or (galaxy in selected and selected[galaxy]["J_ref"] <= hit[0]):
                continue
            objective, index = hit
            profiles = data["frontier"] + data.get("retained_profiles", [])
            entry = next((f for f in profiles if f["record_index"] == index), None)
            if entry is not None:
                xi = entry["xi"]
            elif index == 0:
                xi = [.5] * len(data["knots_kpc"])
            elif data.get("checkpoint_storage", {}).get("mode") == "compact":
                storage = data["checkpoint_storage"]
                profile_path = path.parent / f"{galaxy}_gp_checkpoints" / "profiles.npy"
                if sha(profile_path) != storage["sha256"]:
                    raise ValueError("compact candidate profile hash mismatch")
                array = np.load(profile_path, mmap_mode="r", allow_pickle=False)
                if (list(array.shape) != storage["shape"] or len(array) != len(data["records"])
                        or array.shape[1] != len(data["knots_kpc"])):
                    raise ValueError("compact candidate geometry/index mismatch")
                xi = array[index].copy().tolist()
            else:
                checkpoint = path.parent / f"{galaxy}_gp_checkpoints" / f"candidate_{index:05d}.json"
                saved = json.loads(checkpoint.read_text())
                if (saved["record"] != data["records"][index]
                        or saved["script_sha256"] != data["provenance"]["script_sha256"]
                        or saved["input_sha256"] != data["provenance"]["input_sha256"]):
                    raise ValueError("selected candidate checkpoint provenance mismatch")
                xi = saved["xi"]
            r, ml = valid_curve(data["knots_kpc"], xi)
            selected[galaxy] = {"galaxy": galaxy, "radius_kpc": r, "ml": ml,
                    "J_ref": objective, "source": str(path), "record_index": index,
                    "record": data["records"][index], "probability_threshold": threshold_probability(sigma, convention)}
    report = {"source_files": audited, "threshold_convention": convention,
              "sigma_is_reporting_label_not_calibrated_significance": True,
              "unreachable_in_supplied_searches": {k:sorted(coverage[k]-set(v)) for k,v in panels.items()},
              "available_galaxies": sorted(fingerprints), "selection": "minimum found J across explicitly supplied compatible searches"}
    return {k:list(v.values()) for k,v in panels.items()}, report



def rescale_radii(panels, external, csv_path, coordinate):
    """Keep original kpc values and store the display coordinate separately."""
    if coordinate == "kpc":
        return {"coordinate": "physical_kpc", "normalization_applied": False}
    from utils import load_sparc_table
    table = load_sparc_table()
    rd_sparc = dict(zip(table["Galaxy"], table["Rdisk"]))
    rd_external = collections.defaultdict(set)
    with csv_path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        if "R_bar" not in (reader.fieldnames or []):
            raise ValueError("R/R_d requires the external R_bar column")
        for row in reader:
            rd_external[row["Galaxy"].strip()].add(float(row["R_bar"]))
    for name, values in rd_external.items():
        if len(values) != 1 or not all(math.isfinite(v) and v > 0 for v in values):
            raise ValueError(f"invalid or inconsistent R_bar for {name}: {values}")
    scales = {"SPARC": {}, "external": {}}
    for curves in panels.values():
        for curve in curves:
            name = curve["galaxy"]
            rd = float(rd_sparc[name])
            if not math.isfinite(rd) or rd <= 0:
                raise ValueError(f"invalid SPARC Rdisk for {name}")
            curve["radius_over_rd"] = (np.asarray(curve["radius_kpc"])/rd).tolist()
            scales["SPARC"][name] = rd
    for curve in external:
        name = curve["galaxy"]
        rd = next(iter(rd_external[name]))
        curve["radius_over_rd"] = (np.asarray(curve["radius_kpc"])/rd).tolist()
        scales["external"][name] = rd
    return {"coordinate": "R_over_Rd", "normalization_applied": True,
            "SPARC_scale_source": "SPARC_table.txt: Rdisk [kpc]",
            "SPARC_catalogue_sha256": sha(ROOT/"SPARC_table.txt"),
            "external_scale_source": "ML_Andreea.csv: R_bar, treated as kpc",
            "external_R_bar_as_Rd_is_provisional": True, "scales_kpc": scales}

def draw_clouds(panels, external, output, coordinate="kpc", coverage_note=None):
    radius_key = "radius_over_rd" if coordinate == "rd" else "radius_kpc"
    fig, axes = plt.subplots(1, len(panels)+1, figsize=(6*(len(panels)+1), 5),
                             sharex=True, sharey=True, layout="constrained", squeeze=False)
    for ax, (label, curves) in zip(axes[0], panels.items()):
        for c in curves:
            ax.plot(c[radius_key], c["ml"], color="#663399", alpha=.16, lw=.8)
        ax.axhline(.5, color="black", ls="--", lw=.8, label="SPARC reference 0.5")
        ax.set_title(f"SPARC 3.6 micron: {label}\n{len(curves)} galaxies", fontsize=10)
        ax.legend(fontsize=8)
    ax = axes[0, -1]
    samples = sorted({c["dataset"] for c in external})
    palette = {s:COLORS.get(s, plt.get_cmap("tab10")(i % 10)) for i,s in enumerate(samples)}
    for c in external:
        ax.plot(c[radius_key], c["ml"], "o-", ms=2, lw=.7,
                color=palette[c["dataset"]], alpha=.3)
    handles = [Line2D([0], [0], color=palette[s], label=f"{s}: {sum(c['dataset']==s for c in external)} galaxies")
               for s in samples]
    ax.legend(handles=handles, fontsize=8)
    ax.set_title(f"MIGHTEE-HI / LADUMA: raw K_s profiles\n{len(external)} galaxy profiles", fontsize=10)
    for axis in axes[0]:
        axis.set_xlabel("R / R_d" if coordinate == "rd" else "Radius [kpc]")
        axis.grid(alpha=.15)
    axes[0, 0].set_ylabel("Stellar M/L [solar units in the indicated band]")
    fig.suptitle("Profile-cloud comparison: bandpass/IMF compatibility unresolved"
                 + ("\nExternal R_d = R_bar (provisional)" if coordinate == "rd" else "")
                 + ("\n" + coverage_note if coverage_note else ""), fontsize=11)
    fig.savefig(output / "ml_profile_clouds.png", dpi=160)
    fig.savefig(output / "ml_profile_clouds.pdf")
    plt.close(fig)


def self_test():
    records = [{"valid":True, "J_ref":4., "marginal_logp":math.log(.1)},
               {"valid":True, "J_ref":2., "marginal_logp":math.log(.03)},
               {"valid":False, "J_ref":0., "marginal_logp":0.}]
    assert choose_record(records, .02) == (2., 1)
    assert choose_record(records, .05) == (4., 0)
    assert choose_record(records, .2) is None
    assert abs(threshold_probability(2., "one-sided")-.02275013194817922) < 1e-15
    assert threshold_probability(2., "two-sided") == 2*threshold_probability(2., "one-sided")
    try:
        central_profile({"status":"qp_failure"})
    except ValueError:
        pass
    else:
        raise AssertionError("failure accepted")
    assert valid_curve([2, 1], [.3, .4]) == ([1., 2.], [.4, .3])
    print("Seven selection/validation checks passed; no scientific calculations.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=["central-placeholder", "gp"], default="central-placeholder")
    parser.add_argument("--central-root", type=Path, default=ROOT/"outputs/rerun_boundary_20261006/binned_ML_opt/galaxies")
    parser.add_argument("--external-csv", type=Path, default=ROOT/"ML_Andreea.csv")
    parser.add_argument("--sample-column", default="Sample")
    parser.add_argument("--radius-coordinate", choices=["kpc", "rd"], default="kpc")
    parser.add_argument("--gp-root", type=Path, action="append", default=[])
    parser.add_argument("--gp-read-mode", choices=["full", "frontier-stream"], default="full",
                        help="Stream legacy frontiers/selected records with near-cut full-read fallback")
    parser.add_argument("--sigma", type=float, nargs="+", default=[2., 5.])
    parser.add_argument("--threshold-convention", choices=["one-sided", "two-sided"], default="one-sided")
    parser.add_argument("--output-root", type=Path, default=ROOT/"outputs/ml_profile_clouds/setup_20261013")
    parser.add_argument("--audit-only", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    if args.output_root.exists():
        raise FileExistsError("output root already exists; choose a new root")
    external, csv_audit = external_profiles(args.external_csv, args.sample_column)
    if args.mode == "central-placeholder":
        panels, sparc_audit = load_central(args.central_root)
    else:
        if not args.gp_root or len(set(args.sigma)) != len(args.sigma):
            parser.error("GP mode requires explicit roots and distinct sigma labels")
        panels, sparc_audit = load_gp(args.gp_root, args.sigma, args.threshold_convention, args.gp_read_mode)
    coverage_note = (f"GP completed subset: {len(sparc_audit['available_galaxies'])} / 175 SPARC galaxies; "
                     "sigma cuts are reporting labels, not calibrated significances" if args.mode == "gp" else None)
    radius_audit = rescale_radii(panels, external, args.external_csv, args.radius_coordinate)
    manifest = {"mode":args.mode, "coverage_note":coverage_note, "radius_coordinate":radius_audit["coordinate"],
                "radius_normalization":radius_audit, "external":csv_audit,
                "sparc":sparc_audit, "panels":panels, "raw_external_profiles":external,
                "script_sha256":sha(__file__), "no_fitting_or_probability_recalculation":True,
                "bandpass_conversion_applied":False, "r_E_normalization_applied":False,
                "profile_uncertainty_bands":False, "counts_not_point_pooled_statistics":True}
    args.output_root.mkdir(parents=True)
    if not args.audit_only:
        draw_clouds(panels, external, args.output_root, args.radius_coordinate, coverage_note)
    with (args.output_root/"manifest.json").open("x") as stream:
        json.dump(manifest, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({"output":str(args.output_root), "SPARC_panels":{k:len(v) for k,v in panels.items()},
                      "external_profiles":len(external), "excluded_CSV_rows":len(csv_audit["excluded_rows"]),
                      "audit_only":args.audit_only}, indent=2))


if __name__ == "__main__":
    main()
