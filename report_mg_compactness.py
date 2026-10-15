#!/usr/bin/env python3
"""Descriptive mass/support diagnostics; never a categorical MG test.
Read-only analytic profile specifications or completed retained-only GP records.
Finite disc integrals plus catalogue gas are mixed-support proxies, NOT totals.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import numpy as np
from utils import (DEFAULT_DENS_ROOT, DEFAULT_TABLE_PATH, compute_A_B,
                   compute_A_beta, get_SBdisk, load_sparc_table, lognormal_quadrature)
from sparc_scale_proxies import (G_KPC_KMS2_PER_MSUN, KPC_METRES, A0_DEFAULT_M_S2,
                                disk_mass_for_profile, mass_distribution, sha256_file)
ROOT = Path(__file__).resolve().parent


def bulge_information(rotmod, median):
    """Zero requires an all-zero supplied template; missing is never zero."""
    empty = {"presence": "unknown", "total_mass_1e9": None,
             "total_status": "unavailable_no_documented_total_mass",
             "enclosed_proxy_1e9": None, "enclosed_radius_kpc": None}
    if rotmod is None:
        return empty
    values = np.loadtxt(rotmod, comments="#", ndmin=2)
    if values.shape[1] < 6 or not np.all(np.isfinite(values[:, [0, 5]])):
        raise ValueError("invalid rotmod radius/bulge columns")
    r, v = values[:, 0], values[:, 5]
    if np.any(r <= 0) or np.any(np.diff(r) <= 0) or np.any(v < 0):
        raise ValueError("rotmod radii must increase and be positive; bulge speeds nonnegative")
    if np.all(v == 0):
        return {**empty, "presence": "absent_in_supplied_rotmod_template",
                "total_mass_1e9": 0., "total_status": "absent_template"}
    return {**empty, "presence": "present_in_supplied_rotmod_template",
            "enclosed_proxy_1e9": float(median*r[-1]*v[-1]**2/G_KPC_KMS2_PER_MSUN/1e9),
            "enclosed_radius_kpc": float(r[-1])}


def light_inputs(path, luminosity):
    values = np.loadtxt(path, skiprows=1, ndmin=2)
    if values.shape[1] != 3 or len(values) < 2 or not np.all(np.isfinite(values)):
        raise ValueError("density table requires >=2 finite three-column rows")
    r, disk, bulge = values.T
    if np.any(r < 0) or np.any(np.diff(r) <= 0) or np.any(disk < 0) or np.any(bulge < 0):
        raise ValueError("invalid density radius or brightness; no dropped rows")
    disk, bulge = disk*1e6, bulge*1e6
    ad, ab = compute_A_B(r, disk)[0], compute_A_B(r, bulge)[0]
    if not math.isfinite(luminosity) or luminosity <= 0 or ad <= 0:
        raise ValueError("invalid catalogue/disk luminosity")
    return {"radius": r, "raw_disk_sb": disk, "raw_bulge_sb": bulge,
            "raw_disk_lsun": ad, "raw_total_lsun": ad+ab,
            "scale": luminosity*1e9/(ad+ab), "catalogue_lsun": luminosity*1e9}


def selected_gp(path, galaxy, light, input_hashes):
    """Retained records only; no old full-candidate assumption or auto-pooling."""
    data = json.loads(path.read_text())
    if (data.get("status") != "finite_search_completed" or data.get("galaxy") != galaxy
            or data.get("checkpoint_storage", {}).get("mode") != "retained-only"):
        raise ValueError("requires a completed retained-only result for this galaxy")
    if data["provenance"]["input_sha256"] != input_hashes:
        raise ValueError("GP input hashes differ from supplied current inputs")
    for name, expected in data["provenance"]["source_sha256"].items():
        if sha256_file(ROOT/name) != expected:
            raise ValueError("GP source mismatch: "+name)
    script = ROOT/data["provenance"]["script_file"]
    if sha256_file(script) != data["provenance"]["script_sha256"]:
        raise ValueError("GP entrypoint hash mismatch")
    positive = light["radius"][light["radius"] > 0]
    if not np.array_equal(np.asarray(data["knots_kpc"]), positive):
        raise ValueError("unsupported GP knot geometry; no clamping or extrapolation")
    result = []
    for saved in data["retained_profiles"]:
        i = int(saved["record_index"])
        rec = data["records"][i]
        qc = saved.get("direct_cache_qc", {})
        if not rec.get("valid"):
            continue
        rel = qc.get("relative_to_peak")
        if rel is None or not math.isfinite(float(rel)) or float(rel) > 1e-8:
            raise ValueError("selected GP force audit missing/failed")
        for key in ("J_ref", "marginal_logp"):
            if rec.get(key) is None or not math.isfinite(float(rec[key])):
                raise ValueError("saved valid record has non-finite probability/objective")
        xi = np.asarray(saved["xi"], float)
        if xi.shape != positive.shape or np.any(~np.isfinite(xi)) or np.any(xi <= 0):
            raise ValueError("invalid saved profile")
        lo, hi = data["settings"]["bounds_ml"]
        if np.any(xi < lo-1e-8) or np.any(xi > hi+1e-8):
            raise ValueError("saved GP bounds fail")
        result.append(("saved_record_"+str(i), xi, rec))
    return result, data


def mass_rows(light, profiles, mhi, rd, nuisance, bulge, helium, a0):
    if not all(math.isfinite(x) and x > 0 for x in (rd, helium, a0)) or not math.isfinite(mhi) or mhi < 0:
        raise ValueError("invalid scale, gas mass or physical constants")
    gas_nodes, gas_weights = lognormal_quadrature(nuisance["gas_ml"], nuisance["gas_ml_dex"], nuisance["ml_samples"])
    bul_nodes, bul_weights = lognormal_quadrature(nuisance["bulge_ml"], nuisance["bulge_ml_dex"], nuisance["ml_samples"])
    r = light["radius"]
    rows = []
    for label, xi, record in profiles:
        gp = xi.shape != r.shape
        for support, mask in (("original_photometric", np.ones(len(r), bool)),
                              ("positive_radius_force", r > 0)):
            if gp and support == "original_photometric":
                continue  # GP is only defined on the positive solver support.
            rr = r[mask]
            yy = xi if gp else xi[mask]
            if len(rr) < 2 or yy.shape != rr.shape or np.any(~np.isfinite(yy)) or np.any(yy <= 0):
                raise ValueError("invalid M/L support")
            for convention in ("solver_disk_normalization", "decomposed_light"):
                scale = (light["catalogue_lsun"]/light["raw_disk_lsun"]
                         if convention == "solver_disk_normalization" else light["scale"])
                md = float(compute_A_B(rr, light["raw_disk_sb"][mask]*yy)[0]*scale/1e9)
                lb = float(compute_A_B(rr, light["raw_bulge_sb"][mask])[0]*light["scale"]/1e9)
                mb = (bulge["total_mass_1e9"] if convention == "solver_disk_normalization"
                      else lb*nuisance["bulge_ml"])
                row = {"profile_id": label, "convention": convention, "support": support,
                       "M_disk_support_1e9": md, "M_bulge_estimate_1e9": mb,
                       "M_bar_total_1e9": None, "r_M_total_kpc": None, "eta_total": None,
                       "total_status": "unavailable_finite_disc_support_not_infinite_total",
                       "M_bar_mixed_proxy_1e9": None, "r_M_mixed_proxy_kpc": None,
                       "eta_mixed_proxy": None, "J_ref": record.get("J_ref"),
                       "marginal_logp": record.get("marginal_logp")}
                if convention == "solver_disk_normalization" and mb is None:
                    row["proxy_status"] = "unavailable_solver_bulge_total"
                else:
                    # Catalogue gas coverage differs from finite photometric coverage.
                    # Explicit mixed proxy only; not a supported total or an MG class.
                    unit_bulge = 0. if convention == "solver_disk_normalization" else lb
                    stats = mass_distribution(md, mhi, unit_bulge, gas_nodes, gas_weights,
                                              bul_nodes, bul_weights, nuisance["gas_ml"],
                                              nuisance["bulge_ml"], helium, a0, rd/.85)
                    row.update({"proxy_status": "mixed_support_descriptive_only",
                                "M_bar_mixed_proxy_1e9": stats["M_bar_median_1e9_Msun"],
                                "r_M_mixed_proxy_kpc": stats["r_M_median_kpc"],
                                "eta_mixed_proxy": stats["eta_median"]})
                rows.append(row)
    return rows


def apply_model_specs(rows, specification, rd, a0):
    """Totals of an EXPLICIT declared model, never an inferred force match.
    Entries identify a profile/convention/support. Each component must supply
    mass_1e9, units, normalization, source and inner/outer boundary prescriptions.
    Missing components or no full-model declaration leave model totals unset.
    """
    entries = specification.get("entries", []) if specification else []
    indexed = {}
    for entry in entries:
        key = tuple(entry[k] for k in ("profile_id", "convention", "support"))
        if key in indexed: raise ValueError("duplicate model mass specification")
        indexed[key] = entry
    known = {tuple(r[k] for k in ("profile_id", "convention", "support")) for r in rows}
    if set(indexed)-known: raise ValueError("model specification targets unknown row")
    for row in rows:
        row.update(Mbar_model_1e9=None, rM_model_kpc=None, eta_model=None,
                   model_mass_status="not_defined", component_completeness="not_supplied",
                   force_model_match="not_assessed")
        key = tuple(row[k] for k in ("profile_id", "convention", "support"))
        if key not in indexed: continue
        entry = indexed[key];components = entry.get("components", {})
        missing = set(("disk", "gas", "bulge"))-set(components)
        if missing:
            row.update(model_mass_status="incomplete_specification",
                       component_completeness="missing:"+",".join(sorted(missing)))
            continue
        if set(components) != {"disk", "gas", "bulge"}:
            raise ValueError("unknown mass component")
        masses = {}
        for name, component in components.items():
            value = float(component["mass_1e9"])
            if not math.isfinite(value) or value < 0 or component.get("units") != "1e9_Msun":
                raise ValueError("invalid explicit component mass/units")
            for field in ("normalization", "source", "inner_boundary", "outer_boundary"):
                if not isinstance(component.get(field), str) or not component[field].strip():
                    raise ValueError("component requires "+field)
            masses[name] = value
        row["component_completeness"] = "complete_explicit_specification"
        if entry.get("full_model_defined") is not True or not entry.get("model_convention"):
            row["model_mass_status"] = "not_defined_missing_model_declaration"
            continue
        if not math.isclose(masses["disk"], row["M_disk_support_1e9"], rel_tol=1e-8, abs_tol=1e-12):
            raise ValueError("model disc mass conflicts with labelled profile/support")
        total = sum(masses.values())
        if total <= 0: raise ValueError("nonpositive declared model mass")
        rm = math.sqrt(G_KPC_KMS2_PER_MSUN*total*1e9/(a0*KPC_METRES/1e6))
        row.update(Mbar_model_1e9=total, rM_model_kpc=rm, eta_model=(rd/.85)/rm,
                   model_mass_status="defined_explicit_model_not_astronomical_total",
                   force_model_match="unverified_supplied_model")
    return rows


def self_test():
    import tempfile
    with tempfile.TemporaryDirectory() as directory:
        p = Path(directory)/"rotmod.dat"
        assert bulge_information(None, .7)["total_mass_1e9"] is None
        np.savetxt(p, [[1, 0, 0, 0, 0, 0], [2, 0, 0, 0, 0, 0]])
        assert bulge_information(p, .7)["total_mass_1e9"] == 0
        np.savetxt(p, [[1, 0, 0, 0, 0, 10], [2, 0, 0, 0, 0, 20]])
        bul = bulge_information(p, .7)
        assert bul["total_mass_1e9"] is None and bul["enclosed_proxy_1e9"] > 0
        assert bul["enclosed_radius_kpc"] == 2
        r = np.array([0., 1., 2.]);d = np.ones(3);ad = compute_A_B(r, d)[0]
        light = {"radius": r, "raw_disk_sb": d, "raw_bulge_sb": d,
                 "raw_disk_lsun": ad, "raw_total_lsun": 2*ad,
                 "scale": 1e9/(2*ad), "catalogue_lsun": 1e9}
        nuisance = dict(gas_ml=1., gas_ml_dex=.04, bulge_ml=.7, bulge_ml_dex=.1, ml_samples=9)
        rows = mass_rows(light, [("reference", np.full(3,.5), {})], .1, 1., nuisance, bul, 1.33, A0_DEFAULT_M_S2)
        assert len(rows) == 4 and all(x["eta_total"] is None for x in rows)
        assert rows[0]["eta_mixed_proxy"] is None and rows[1]["eta_mixed_proxy"] > 0
        assert math.isclose(rows[0]["M_disk_support_1e9"], .5)
        assert math.isclose(rows[1]["M_disk_support_1e9"], .25)
        assert rows[2]["M_disk_support_1e9"] < rows[0]["M_disk_support_1e9"]
        gp_rows = mass_rows(light, [("gp", np.full(2,.5), {})], .1, 1., nuisance, bul, 1.33, A0_DEFAULT_M_S2)
        assert len(gp_rows) == 2 and all(x["support"] == "positive_radius_force" for x in gp_rows)
        result_path = Path(directory)/"gp.json"
        fixture = {"status":"finite_search_completed", "galaxy":"synthetic",
                   "checkpoint_storage":{"mode":"retained-only"},
                   "settings":{"bounds_ml":[.01,10.]}, "knots_kpc":[1.,2.],
                   "records":[{"valid":True,"J_ref":0.,"marginal_logp":-.1}],
                   "retained_profiles":[{"record_index":0,"xi":[.5,.5],
                                         "direct_cache_qc":{"relative_to_peak":0.}}],
                   "provenance":{"input_sha256":{"dens":"synthetic"},
                                 "source_sha256":{"utils.py":sha256_file(ROOT/"utils.py")},
                                 "script_file":Path(__file__).name,"script_sha256":sha256_file(Path(__file__))}}
        result_path.write_text(json.dumps(fixture))
        profiles, _ = selected_gp(result_path, "synthetic", light, {"dens":"synthetic"})
        assert len(profiles) == 1 and np.array_equal(profiles[0][1],[.5,.5])
        try:
            selected_gp(result_path, "synthetic", light, {"dens":"wrong"})
        except ValueError:
            pass
        else:
            raise AssertionError("input hash mismatch accepted")
        fixture["retained_profiles"][0]["xi"][0] = 11.
        result_path.write_text(json.dumps(fixture))
        try:
            selected_gp(result_path, "synthetic", light, {"dens":"synthetic"})
        except ValueError:
            pass
        else:
            raise AssertionError("out-of-bounds GP profile accepted")
    print("Reporter synthetic checks passed; no real galaxy calculation.")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--self-test", action="store_true")
    p.add_argument("--galaxy")
    p.add_argument("--profile", choices=["constant", "linear", "log_linear"], default="constant")
    p.add_argument("--mean-ml", type=float, default=.5)
    p.add_argument("--slope", type=float, default=0.)
    p.add_argument("--gp-result", type=Path)
    p.add_argument("--mass-spec", type=Path, help="Explicit component model definitions; never inferred from mixed proxies")
    p.add_argument("--rotmod", type=Path)
    p.add_argument("--dens-root", type=Path, default=Path(DEFAULT_DENS_ROOT))
    p.add_argument("--table", type=Path, default=Path(DEFAULT_TABLE_PATH))
    p.add_argument("--helium-factor", type=float, default=1.33)
    p.add_argument("--a0", type=float, default=A0_DEFAULT_M_S2)
    p.add_argument("--output-root", type=Path)
    args = p.parse_args()
    if args.self_test:
        self_test();return
    if not args.galaxy or args.output_root is None:
        p.error("--galaxy and a fresh --output-root required")
    if args.galaxy not in set(load_sparc_table(args.table).Galaxy):
        p.error("unknown SPARC galaxy")
    if args.output_root.exists():
        raise FileExistsError("output already exists; no overwrite")
    table = load_sparc_table(args.table);cat = table.loc[table.Galaxy == args.galaxy]
    if len(cat) != 1: raise ValueError("catalogue galaxy is not unique")
    cat = cat.iloc[0];dens = args.dens_root/(args.galaxy+".dens")
    files = {"table": args.table, "dens": dens}
    if args.rotmod is not None: files["rotmod"] = args.rotmod
    if args.gp_result is not None: files["gp_result"] = args.gp_result
    if args.mass_spec is not None: files["mass_spec"] = args.mass_spec
    sources = {name: ROOT/name for name in ("utils.py", "sparc_scale_proxies.py", Path(__file__).name)}
    hashes = {k: sha256_file(v) for k,v in files.items()}
    source_hashes = {k: sha256_file(v) for k,v in sources.items()}
    light = light_inputs(dens, float(cat.L))
    rr, sb = get_SBdisk(args.galaxy, args.dens_root, table)
    if not np.array_equal(rr, light["radius"]) or not np.allclose(sb, light["raw_disk_sb"]*light["catalogue_lsun"]/light["raw_disk_lsun"], rtol=1e-13, atol=0):
        raise ValueError("normalization does not match current solver utility")
    nuisance = dict(gas_ml=1., gas_ml_dex=.04, bulge_ml=.7, bulge_ml_dex=.1, ml_samples=9)
    gp_data = None
    if args.gp_result is None:
        _, xi = disk_mass_for_profile(args.profile, args.mean_ml, args.slope, light)
        profiles = [(args.profile, xi, {})]
        selection = "explicit_analytic_profile_not_threshold_selected"
    else:
        if args.rotmod is None: p.error("--gp-result requires matching --rotmod")
        profiles, gp_data = selected_gp(args.gp_result, args.galaxy, light, {k:hashes[k] for k in ("table","dens","rotmod")})
        for k in nuisance:
            nuisance[k] = gp_data["settings"]["quadrature_nodes_per_term" if k == "ml_samples" else k]
        selection = "retained_valid_profiles_not_population_or_posterior_samples"
    bulge = bulge_information(args.rotmod, nuisance["bulge_ml"])
    rows = mass_rows(light, profiles, float(cat.MHI), float(cat.Rdisk), nuisance, bulge, args.helium_factor, args.a0)
    if not rows: raise ValueError("no valid profiles to report")
    specification = json.loads(args.mass_spec.read_text()) if args.mass_spec else None
    if specification and specification.get("galaxy") != args.galaxy:
        raise ValueError("model specification galaxy mismatch")
    apply_model_specs(rows, specification, float(cat.Rdisk), args.a0)
    manifest = {"schema_version": 2, "galaxy": args.galaxy, "selection": selection,
                "model_mass_specification": specification,
                "model_total_policy": "only explicitly declared complete component model; not a verified force match or astronomical true total",
                "units": {"mass":"1e9_Msun","radius":"kpc","SB_input":"Lsun/pc2"},
                "input_paths": {k:str(v.resolve()) for k,v in files.items()}, "input_sha256": hashes,
                "source_sha256": source_hashes, "nuisance": nuisance, "bulge": bulge,
                "a_proxy_kpc": float(cat.Rdisk)/.85, "a0_m_s2": args.a0,
                "helium_factor": args.helium_factor, "gas_catalogue_1e9_Msun": float(cat.MHI)*args.helium_factor,
                "gas_support": "catalogue_HI_mass; molecular_gas_omitted; radial_distribution_unavailable",
                "bulge_enclosed_method": "median_ML*R_last*Vbul_template_last^2/G; spherical; not total",
                "conventions": {"solver_disk_normalization": "disc alone scaled to full catalogue L; bulge total unavailable unless template absent",
                                "decomposed_light": "disc+bulge light jointly scaled to catalogue L; photometric bulge unverified against rotmod"},
                "supports": {name:{"rmin_kpc":float(r[0]), "rmax_kpc":float(r[-1]), "bins":len(r)}
                             for name,r in (("original_photometric", light["radius"]), ("positive_radius_force", light["radius"][light["radius"]>0]))},
                "integration": "trapezoidal 2*pi*R*SB*ML over labelled tabulated support; not exact implicit vcdisk mass",
                "profile_parameters": None if gp_data is not None else {"family":args.profile,"mean_ml":args.mean_ml,"slope":args.slope},
                "force_or_probability_recomputed": False, "radial_acceleration_fields": "unavailable_no_audited_saved_field_adapter",
                "categorical_MG_classification": False, "rows":len(rows)}
    if hashes != {k:sha256_file(v) for k,v in files.items()} or source_hashes != {k:sha256_file(v) for k,v in sources.items()}:
        raise RuntimeError("input/source changed during report")
    args.output_root.mkdir(parents=True)
    with (args.output_root/"profiles.csv").open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    with (args.output_root/"manifest.json").open("x") as stream:
        json.dump(manifest, stream, indent=2, allow_nan=False);stream.write("\n")
    print("Wrote descriptive rows",len(rows),"to",args.output_root)


if __name__ == "__main__":
    main()
