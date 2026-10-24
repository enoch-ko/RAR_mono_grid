"""Export descriptive constant/linear products from frozen raw-P grid CSVs only."""
from pathlib import Path
import csv, json, math, hashlib, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("/mnt/users/koe/RAR_mono_grid")
OUT = ROOT / "outputs/saved_profile_comparison_20261019_exactcut_v1"
sys.path.insert(0, str(ROOT))
from utils import compute_A_B, get_SBdisk, load_sparc_table
TABLE = load_sparc_table(ROOT / "SPARC_table.txt")
DOC = json.loads((OUT / "comparison.json").read_text())
CUT = 0.0013498980316300957
LOGCUT = math.log(CUT)
assert DOC["raw_probability_cut"] == CUT

# All outputs are new names; refuse to replace anything.
def save_csv(path, fields, rows):
    if path.exists():
        raise FileExistsError(path)
    with path.open("x", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader(); w.writerows(rows)

def save_fig(fig, path):
    if path.exists():
        raise FileExistsError(path)
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)

def rows_csv(path):
    with path.open(newline="") as f:
        yield from csv.DictReader(f)

def qualifies_marginalised_raw_cut(row):
    if row.get("probability_kind") != "marginalised_gas_bulge":
        return False
    if row.get("valid", "").lower() != "true":
        return False
    try:
        logp = float(row["log_probability"])
    except (KeyError, TypeError, ValueError):
        return False
    return math.isfinite(logp) and logp >= LOGCUT

def run_selection_regression():
    fixture = [
        {"probability_kind": "marginalised_gas_bulge", "valid": "True", "log_probability": str(LOGCUT), "id": "inclusive-boundary"},
        {"probability_kind": "marginalised_gas_bulge", "valid": "True", "log_probability": str(math.nextafter(LOGCUT, -math.inf)), "id": "below-cut"},
        {"probability_kind": "marginalised_gas_bulge", "valid": "False", "log_probability": "0", "id": "invalid"},
        {"probability_kind": "conditional", "valid": "True", "log_probability": "0", "id": "wrong-probability-kind"},
        {"probability_kind": "marginalised_gas_bulge", "valid": "True", "log_probability": "nan", "id": "nonfinite"},
    ]
    kept = [row["id"] for row in fixture if qualifies_marginalised_raw_cut(row)]
    assert kept == ["inclusive-boundary"], kept
    return True

SELECTION_REGRESSION_PASSED = run_selection_regression()

ledger = DOC["ledger"]
for fam in ("constant", "linear"):
    assert sum(x["family"] == fam for x in ledger) == 175

# Constant: preserve every contiguous run of qualifying sampled grid nodes.
constants = [x for x in ledger if x["family"] == "constant"]
paths = [Path(p) for p in sorted({x["source"]["grid"] for x in constants if x.get("source")})]
bygal = {}
for x in constants:
    bygal[x["galaxy"]] = {"required_source": x["required_source"],
        "status": x["status"], "mass_rank_value": 0.5 * float(x["normalization"]["Ltotal_1e9"]),
        "source_grid": x["source"]["grid"]}
accepted = {g: [] for g in bygal}
cell_audit = {}
for galaxy, meta in bygal.items():
    p = Path(meta["source_grid"])
    if not p.exists():
        raise FileNotFoundError(p)
    galaxy_rows = [r for r in rows_csv(p) if r["galaxy"] == galaxy
                   and r["profile"] == "constant"
                   and r["probability_kind"] == "marginalised_gas_bulge"]
    grid_values = sorted({float(r["log10_upsilon_disk"]) for r in galaxy_rows})
    grid_position = {value: i for i, value in enumerate(grid_values)}
    counts = {"marginalised_rows": len(galaxy_rows), "invalid": 0,
              "nonfinite_log_probability": 0, "malformed_log_probability": 0,
              "below_cut": 0, "qualifying": 0}
    cell_ids = [(r["grid_i"], r["grid_j"]) for r in galaxy_rows]
    if len(cell_ids) != len(set(cell_ids)):
        raise ValueError(f"{galaxy}: duplicate marginalised grid cells in source CSV")
    for r in galaxy_rows:
        if r["valid"].lower() != "true":
            counts["invalid"] += 1
            continue
        try:
            logp = float(r["log_probability"])
        except (TypeError, ValueError):
            counts["malformed_log_probability"] += 1
            continue
        if not math.isfinite(logp):
            counts["nonfinite_log_probability"] += 1
            continue
        if logp < LOGCUT:
            counts["below_cut"] += 1
            continue
        if not qualifies_marginalised_raw_cut(r):
            raise AssertionError(f"{galaxy}: qualifying-cell predicate disagrees")
        value = float(r["log10_upsilon_disk"])
        accepted[galaxy].append((grid_position[value], value))
        counts["qualifying"] += 1
    if counts["marginalised_rows"] != sum(counts[k] for k in
            ("invalid", "nonfinite_log_probability", "malformed_log_probability", "below_cut", "qualifying")):
        raise AssertionError(f"{galaxy}: cell audit does not reconcile")
    cell_audit[galaxy] = counts

interval_rows = []
for galaxy, vals in accepted.items():
    vals = sorted(set(vals))
    chunks = []
    for item in vals:
        if not chunks or item[0] != chunks[-1][-1][0] + 1:
            chunks.append([item])
        else:
            chunks[-1].append(item)
    for k, chunk in enumerate(chunks, 1):
        meta = bygal[galaxy]
        interval_rows.append({"galaxy": galaxy, "required_source": meta["required_source"],
            "mass_order_value_0p5Ltotal_1e9": meta["mass_rank_value"], "status": meta["status"],
            "interval_id": k, "first_grid_position": chunk[0][0], "last_grid_position": chunk[-1][0],
            "log10_ml_min_sample": chunk[0][1], "log10_ml_max_sample": chunk[-1][1],
            "n_qualifying_grid_nodes": len(chunk), "probability_cut_raw": CUT,
            "note": "sampled grid-node run; no interpolation across gaps"})
interval_rows.sort(key=lambda r: (r["mass_order_value_0p5Ltotal_1e9"], r["galaxy"], r["interval_id"]))
interval_path = OUT / "constant_rawP_accepted_grid_runs_source_scoped_v7.csv"
planned_outputs = [interval_path, OUT / "constant_rawP_accepted_grid_runs_source_scoped_v7.png",
                   OUT / "linear_marginalised_rawP_leastJrad_selected_profiles_v7.csv",
                   OUT / "linear_marginalised_rawP_profiles_kpc_v7.png",
                   OUT / "linear_marginalised_rawP_profiles_Rd_v7.png",
                   OUT / "grid_figure_manifest_20261020_source_scoped_v7.json"]
existing_outputs = [str(p) for p in planned_outputs if p.exists()]
if existing_outputs:
    raise FileExistsError("Refusing to reuse existing v7 products: " + ", ".join(existing_outputs))
save_csv(interval_path, list(interval_rows[0]), interval_rows)

fig, ax = plt.subplots(figsize=(9, max(12, .115 * len(bygal))))
ordered = sorted(bygal, key=lambda g: (bygal[g]["mass_rank_value"], g))
ypos = {g: i for i, g in enumerate(ordered)}
colors = {"original": "#276FBF", "corrected": "#D97706"}
for row in interval_rows:
    y = ypos[row["galaxy"]]
    ax.plot([row["log10_ml_min_sample"], row["log10_ml_max_sample"]], [y, y],
        color=colors[row["required_source"]], lw=2.0, solid_capstyle="butt")
prior_lo, prior_hi = math.log10(.5)-.1, math.log10(.5)+.1
ax.axvspan(prior_lo, prior_hi, color="grey", alpha=.16, label="SPARC 0.5 ± 0.1 dex prior")
ax.axvline(math.log10(.5), color="black", ls="--", lw=.8)
ax.set_yticks(range(len(ordered))); ax.set_yticklabels(ordered, fontsize=4.2)
ax.set_xlabel(r"$\log_{10}(\Upsilon_{\rm disk})$  [sampled grid nodes]")
ax.set_ylabel("Galaxies (ordered by fixed 0.5 × total stellar luminosity; display only)")
ax.set_title("Constant profiles passing raw P ≥ Φ(−3); gaps retained")
from matplotlib.lines import Line2D
handles = [Line2D([0],[0],color=colors[k],lw=2,label=k) for k in colors]
handles.append(plt.Rectangle((0,0),1,1,color="grey",alpha=.16,label="SPARC prior"))
ax.legend(handles=handles, fontsize=7, loc="best")
fig.tight_layout(); save_fig(fig, OUT / "constant_rawP_accepted_grid_runs_source_scoped_v7.png")

# Linear: rescan every valid nuisance-marginalised raw-P grid cell. The
# profile parameterisation and full-density support are taken from the actual
# production source and get_SBdisk/compute_A_B implementation.
def linear_jrad(mean_ml, slope_ml_per_kpc, pivot_kpc, rmin_kpc, rmax_kpc):
    """Exact uniform-dR J_rad for the producer's linear profile on full support."""
    if not (math.isfinite(rmin_kpc) and math.isfinite(rmax_kpc) and rmax_kpc > rmin_kpc):
        raise ValueError("linear J_rad requires finite increasing full .dens support")
    u0 = mean_ml + slope_ml_per_kpc * (rmin_kpc - pivot_kpc)
    u1 = mean_ml + slope_ml_per_kpc * (rmax_kpc - pivot_kpc)
    if not (math.isfinite(u0) and math.isfinite(u1) and u0 > 0 and u1 > 0):
        raise ValueError("linear profile must remain positive on full .dens support")
    d0 = (u0 - 0.5) / 0.5
    d1 = (u1 - 0.5) / 0.5
    return (d0 * d0 + d0 * d1 + d1 * d1) / 3.0, u0, u1

def run_linear_metric_regression():
    # For a linear profile, squared fractional departure is quadratic; this
    # numerical integral checks the exact endpoint formula on a known fixture.
    rr = np.linspace(0.0, 2.0, 2001)
    profile = 0.4 + 0.1 * rr
    exact, _, _ = linear_jrad(0.4, 0.1, 0.0, 0.0, 2.0)
    squared_departure = ((profile - 0.5) / 0.5) ** 2
    numerical = np.sum(squared_departure[0:-1:2] + 4.0 * squared_departure[1::2] + squared_departure[2::2]) * (rr[1] - rr[0]) / 3.0 / (rr[-1] - rr[0])
    assert np.isclose(exact, numerical, rtol=1e-12, atol=1e-14), (exact, numerical)
    zero, _, _ = linear_jrad(0.5, 0.0, 4.0, 0.0, 2.0)
    assert abs(zero) < 1e-15, zero
    return True

LINEAR_METRIC_REGRESSION_PASSED = run_linear_metric_regression()
linear_entries = [x for x in ledger if x["family"] == "linear"]
linear_rows = []
linear_profiles = {}
linear_cell_audit = {}
linear_skips = {}
source_cells = {}
linear_paths = sorted({Path(x["source"]["grid"]) for x in linear_entries
                       if x.get("source") and x["source"].get("grid")})
for source_path in linear_paths:
    if not source_path.exists():
        raise FileNotFoundError(source_path)
    for row in rows_csv(source_path):
        if row["profile"] == "linear" and row["probability_kind"] == "marginalised_gas_bulge":
            source_cells.setdefault((str(source_path), row["galaxy"]), []).append(row)
rdisk_by_galaxy = {str(row.Galaxy): float(row.Rdisk) for row in TABLE.itertuples()}
for entry in linear_entries:
    galaxy = entry["galaxy"]
    source = entry.get("source") or {}
    grid_path = source.get("grid")
    if not grid_path:
        linear_skips[galaxy] = "no ledger-linked linear grid CSV"
        continue
    grid_path_obj = Path(grid_path)
    cell_rows = source_cells.get((str(grid_path_obj), galaxy), [])
    audit = {"marginalised_rows": len(cell_rows), "invalid": 0,
             "nonfinite_log_probability": 0, "malformed_log_probability": 0,
             "below_cut": 0, "qualifying": 0, "nonpositive_full_support": 0}
    ids = [(r["grid_i"], r["grid_j"]) for r in cell_rows]
    if len(ids) != len(set(ids)):
        raise ValueError(f"{galaxy}: duplicate marginalised linear grid cells")
    if not cell_rows:
        linear_cell_audit[galaxy] = audit
        linear_skips[galaxy] = "no marginalised linear cells in linked CSV"
        continue
    rad, sb_disk = get_SBdisk(galaxy, table=TABLE)
    if (rad.ndim != 1 or rad.size < 2 or np.any(~np.isfinite(rad))
            or np.any(np.diff(rad) <= 0)):
        raise ValueError(f"{galaxy}: invalid original .dens radius grid")
    A, B = compute_A_B(rad, sb_disk)
    if not (math.isfinite(A) and A > 0 and math.isfinite(B)):
        raise ValueError(f"{galaxy}: invalid producer luminosity integrals")
    pivot = B / A
    candidates = []
    for row in cell_rows:
        if row["valid"].lower() != "true":
            audit["invalid"] += 1
            continue
        try:
            logp = float(row["log_probability"])
        except (TypeError, ValueError):
            audit["malformed_log_probability"] += 1
            continue
        if not math.isfinite(logp):
            audit["nonfinite_log_probability"] += 1
            continue
        if logp < LOGCUT:
            audit["below_cut"] += 1
            continue
        if not qualifies_marginalised_raw_cut(row):
            raise AssertionError(f"{galaxy}: linear selection predicate disagrees")
        try:
            mean_ml = float(row["upsilon_disk_mean"])
            slope = float(row["disk_slope"])
            jrad, u0, u1 = linear_jrad(mean_ml, slope, pivot, float(rad[0]), float(rad[-1]))
        except ValueError:
            audit["nonpositive_full_support"] += 1
            continue
        audit["qualifying"] += 1
        candidates.append({"source": entry.get("required_source"), "source_grid": str(grid_path_obj),
            "galaxy": galaxy, "probability_kind": "marginalised_gas_bulge",
            "raw_probability": math.exp(logp), "log_probability": logp,
            "grid_i": row["grid_i"], "grid_j": row["grid_j"],
            "mean_ml_parameter": mean_ml, "disk_slope": slope, "pivot_kpc": pivot,
            "a_ml_at_R0": mean_ml - slope * pivot, "b_ml_per_kpc": slope,
            "J_rad": jrad, "u_at_support_min": u0, "u_at_support_max": u1,
            "mean_ml_uniform_dR": (u0 + u1) / 2.0,
            "minimum_ml": min(u0, u1), "maximum_ml": max(u0, u1),
            "log10_range_dex": math.log10(max(u0, u1) / min(u0, u1)),
            "rms_departure_ml": 0.5 * math.sqrt(jrad),
            "support_min_kpc": float(rad[0]), "support_max_kpc": float(rad[-1]),
            "support_includes_R0": bool(rad[0] == 0.0),
            "rd_kpc": rdisk_by_galaxy.get(galaxy),
            "selection": "least J_rad among all valid marginalised cells with inclusive raw P >= Phi(-3); full original .dens support"})
    if audit["marginalised_rows"] != sum(audit[k] for k in
            ("invalid", "nonfinite_log_probability", "malformed_log_probability",
             "below_cut", "qualifying", "nonpositive_full_support")):
        raise AssertionError(f"{galaxy}: linear cell audit does not reconcile")
    linear_cell_audit[galaxy] = audit
    if not candidates:
        linear_skips[galaxy] = "no qualifying valid cell positive over full original .dens support"
        continue
    winner = min(candidates, key=lambda r: (r["J_rad"], -r["log_probability"],
                    int(r["grid_i"]), int(r["grid_j"])))
    linear_rows.append(winner)
    linear_profiles[galaxy] = (rad.copy(), winner["mean_ml_parameter"]
        + winner["disk_slope"] * (rad - winner["pivot_kpc"]))
linear_rows.sort(key=lambda r: r["galaxy"])
linear_path = OUT / "linear_marginalised_rawP_leastJrad_selected_profiles_v7.csv"
if not linear_rows:
    raise RuntimeError("No linear galaxy has a qualifying physical marginalised grid cell")
save_csv(linear_path, list(linear_rows[0]), linear_rows)
for coord, xlabel, filename in [("kpc", "Radius [kpc]", "linear_marginalised_rawP_profiles_kpc_v7.png"),
                                ("rd", r"Radius / $R_d$", "linear_marginalised_rawP_profiles_Rd_v7.png")]:
    fig, ax = plt.subplots(figsize=(9, 5.5))
    means = [r["mean_ml_uniform_dR"] for r in linear_rows]
    norm = plt.Normalize(min(means), max(means)) if max(means) > min(means) else plt.Normalize(min(means)-.01, max(means)+.01)
    cmap = plt.get_cmap("viridis")
    for row in linear_rows:
        xx, yy = linear_profiles[row["galaxy"]]
        if coord == "rd":
            rd = row.get("rd_kpc")
            if rd is None or not math.isfinite(rd) or rd <= 0:
                continue
            xx = xx / rd
        ax.plot(xx, yy, color=cmap(norm(row["mean_ml_uniform_dR"])), alpha=.28, lw=.7)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
    fig.colorbar(sm, ax=ax, label=r"Uniform-$dR$ mean $\Upsilon_{\rm disk}$")
    ax.axhline(.5, color="black", ls="--", lw=.8, label="SPARC reference 0.5")
    ax.set_yscale("log"); ax.set_xlabel(xlabel); ax.set_ylabel(r"$\Upsilon_{\rm disk}$")
    ax.set_title(f"Linear profiles: least full-support $J_{{rad}}$ among raw-P accepted cells ({len(linear_rows)}/175)")
    ax.legend(fontsize=7); fig.tight_layout(); save_fig(fig, OUT / filename)
linear_status = {"status": "complete", "selected_profiles": len(linear_rows),
    "galaxies_without_selected_profile": linear_skips,
    "probability_kind": "marginalised_gas_bulge", "probability_cut_inclusive": CUT,
    "jrad_support": "full original .dens endpoints, uniform continuous dR, including R=0",
    "profile_parameterisation": "Upsilon(R)=mean_ml+slope*(R-B/A), with A,B from utils.compute_A_B on the full original .dens profile"}

# Reproducibility manifest for both source-scoped profile selections.
constant_png = OUT / "constant_rawP_accepted_grid_runs_source_scoped_v7.png"
linear_products = [OUT / "linear_marginalised_rawP_profiles_kpc_v7.png",
                   OUT / "linear_marginalised_rawP_profiles_Rd_v7.png"]
products = [interval_path, constant_png, linear_path, *linear_products]
source_grid_records = []
for source_path in sorted(set(paths + linear_paths), key=str):
    if not source_path.exists():
        raise FileNotFoundError(source_path)
    raw = source_path.read_bytes()
    source_grid_records.append({"path": str(source_path), "bytes": len(raw),
        "sha256": hashlib.sha256(raw).hexdigest()})
calculation_sources = []
for source_path in (ROOT / "mono_SPARC_RAR.py", ROOT / "utils.py"):
    raw = source_path.read_bytes()
    calculation_sources.append({"path": str(source_path), "bytes": len(raw),
        "sha256": hashlib.sha256(raw).hexdigest()})
manifest = {
    "diagnostic_only": True,
    "probability_cut_inclusive": CUT,
    "selection_probability_kind": "marginalised_gas_bulge",
    "constant_selection": "all valid marginalised_gas_bulge cells with log_probability >= log(Phi(-3)); source scoped by ledger",
    "synthetic_selection_regression_passed": SELECTION_REGRESSION_PASSED,
    "constant_galaxies_with_accepted_nodes": sum(bool(v) for v in accepted.values()),
    "constant_accepted_grid_nodes": sum(len(v) for v in accepted.values()),
    "constant_cell_audit_by_galaxy": cell_audit,
    "constant_cell_audit_totals": {key: sum(row[key] for row in cell_audit.values())
        for key in next(iter(cell_audit.values()))},
    "constant_accepted_grid_runs": len(interval_rows),
    "linear_selection": linear_status,
    "linear_cell_audit_by_galaxy": linear_cell_audit,
    "linear_cell_audit_totals": {key: sum(row[key] for row in linear_cell_audit.values())
        for key in next(iter(linear_cell_audit.values()))},
    "linear_metric_regression_passed": LINEAR_METRIC_REGRESSION_PASSED,
    "linear_profile_parameterisation": "Upsilon(R)=upsilon_disk_mean + disk_slope*(R-B/A); producer computes A,B=compute_A_B(R,SBdisk) on full .dens arrays",
    "linear_jrad": "exact uniform-dR integral of ((Upsilon-0.5)/0.5)^2 over [R_dens_min,R_dens_max], using endpoint deviations; full original .dens support, R=0 included when present",
    "source_comparison_json": {"path": str(OUT / "comparison.json"),
        "sha256": hashlib.sha256((OUT / "comparison.json").read_bytes()).hexdigest()},
    "source_grid_csvs": source_grid_records,
    "calculation_sources": calculation_sources,
    "common_q": "not available; not recomputed",
    "fully_marginalised_table2_support": "not present in source fully_marginalised CSV schema",
    "products": [{"path": str(p), "bytes": p.stat().st_size,
        "sha256": hashlib.sha256(p.read_bytes()).hexdigest()} for p in products],
    "builder_script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
}
mp = OUT / "grid_figure_manifest_20261020_source_scoped_v7.json"
if mp.exists(): raise FileExistsError(mp)
mp.write_text(json.dumps(manifest, indent=2) + "\n")
print("SAVED", sum(bool(v) for v in accepted.values()), "constant galaxies /",
      sum(len(v) for v in accepted.values()), "constant qualifying cells /", len(interval_rows), "contiguous intervals")
print("LINEAR", linear_status["status"], linear_status["selected_profiles"], "selected profiles")
print("MANIFEST", mp)
