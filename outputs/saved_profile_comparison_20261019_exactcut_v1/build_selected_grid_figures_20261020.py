"""Export descriptive constant/linear products from frozen raw-P grid CSVs only."""
from pathlib import Path
import csv, json, math, hashlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("/mnt/users/koe/RAR_mono_grid")
OUT = ROOT / "outputs/saved_profile_comparison_20261019_exactcut_v1"
DOC = json.loads((OUT / "comparison.json").read_text())
CUT = 0.0013498980316300957
LOGCUT = math.log(CUT)
assert DOC["raw_probability_cut"] == CUT

# All outputs are new names; refuse to replace anything.
def save_csv(path, fields, rows):
    if path.exists():
        return
    with path.open("x", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader(); w.writerows(rows)

def save_fig(fig, path):
    if path.exists():
        return
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)

def rows_csv(path):
    with path.open(newline="") as f:
        yield from csv.DictReader(f)

ledger = DOC["ledger"]
for fam in ("constant", "linear"):
    assert sum(x["family"] == fam for x in ledger) == 175

# Constant: preserve every contiguous run of qualifying sampled grid nodes.
constants = [x for x in ledger if x["family"] == "constant"]
paths = sorted({x["source"]["grid"] for x in constants if x.get("source")})
bygal = {}
for x in constants:
    bygal[x["galaxy"]] = {"required_source": x["required_source"],
        "status": x["status"], "mass_rank_value": 0.5 * float(x["normalization"]["Ltotal_1e9"]),
        "source_grid": x["source"]["grid"]}
accepted = {g: [] for g in bygal}
for galaxy, meta in bygal.items():
    p = Path(meta["source_grid"])
    if not p.exists():
        raise FileNotFoundError(p)
    galaxy_rows = [r for r in rows_csv(p) if r["galaxy"] == galaxy
                   and r["profile"] == "constant"
                   and r["probability_kind"] == "conditional"]
    grid_values = sorted({float(r["log10_upsilon_disk"]) for r in galaxy_rows})
    grid_position = {value: i for i, value in enumerate(grid_values)}
    for r in galaxy_rows:
        if r["valid"].lower() != "true" or float(r["log_probability"]) < LOGCUT:
            continue
        value = float(r["log10_upsilon_disk"])
        accepted[galaxy].append((grid_position[value], value))

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
interval_path = OUT / "constant_rawP_accepted_grid_runs_source_scoped_v4.csv"
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
fig.tight_layout(); save_fig(fig, OUT / "constant_rawP_accepted_grid_runs_source_scoped_v4.png")

# Linear: one least-J_rad qualifying saved grid point per galaxy; export physical support summaries.
selected = [r for r in DOC["grid_selected"] if r["family"] == "linear"]
assert len(selected) == 77
linear_rows = []
for r in selected:
    x = np.asarray(r["radii_kpc"], float); y = np.asarray(r["profile_ml"], float)
    slope, intercept = np.polyfit(x, y, 1)

    ps = r["physical_summary"]
    linear_rows.append({"galaxy": r["galaxy"], "required_source": next(z["required_source"] for z in ledger if z["galaxy"] == r["galaxy"] and z["family"] == "linear"),
        "raw_probability": math.exp(r["log_probability"]), "log_probability": r["log_probability"],
        "common_q": "unavailable", "grid_i": r["grid_i"], "grid_j": r["grid_j"],
        "a_ml_at_R0": intercept, "b_ml_per_kpc": slope,
        "mean_ml_uniform_dR": ps["uniform_radial_mean_ml"], "min_ml": ps["minimum_ml"],
        "max_ml": ps["maximum_ml"], "log10_range_dex": ps["log10_range"],
        "rms_departure_from_0p5": ps["rms_departure_ml"], "J_rad": ps["J_rad"],
        "support_min_kpc": ps["support_kpc"][0], "support_max_kpc": ps["support_kpc"][1],
        "rd_kpc": r.get("rd_kpc"), "selection": r["selection"]})
linear_rows.sort(key=lambda r: r["galaxy"])
linear_path = OUT / "linear_rawP_leastJrad_selected_profiles.csv"
save_csv(linear_path, list(linear_rows[0]), linear_rows)

for coord, xlabel, filename in [("kpc", "Radius [kpc]", "linear_rawP_selected_profiles_kpc.png"),
                                ("rd", r"Radius / $R_d$", "linear_rawP_selected_profiles_Rd.png")]:
    fig, ax = plt.subplots(figsize=(9, 5.5))
    cmap = plt.get_cmap("viridis")
    vals = [r["mean_ml_uniform_dR"] for r in linear_rows]
    norm = plt.Normalize(min(vals), max(vals)) if max(vals) > min(vals) else plt.Normalize(min(vals)-.01, max(vals)+.01)
    for r in selected:
        xx = np.asarray(r["radii_kpc"], float)
        if coord == "rd":
            rd = r.get("rd_kpc")
            if not rd or rd <= 0: continue
            xx = xx / rd
        yy = np.asarray(r["profile_ml"], float)
        mr = next(q for q in linear_rows if q["galaxy"] == r["galaxy"])
        ax.plot(xx, yy, color=cmap(norm(mr["mean_ml_uniform_dR"])), alpha=.28, lw=.7)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
    fig.colorbar(sm, ax=ax, label=r"Uniform-$dR$ mean $\Upsilon_{\rm disk}$")
    ax.axhline(.5, color="black", ls="--", lw=.8, label="SPARC reference 0.5")
    ax.set_yscale("log"); ax.set_xlabel(xlabel); ax.set_ylabel(r"$\Upsilon_{\rm disk}$")
    ax.set_title("Linear least-$J_{rad}$ saved profile per galaxy; raw P cut only (77/175)")
    ax.legend(fontsize=7); fig.tight_layout(); save_fig(fig, OUT / filename)

# Reproducibility manifest for the corrected constant product and saved linear products.
constant_png = OUT / "constant_rawP_accepted_grid_runs_source_scoped_v4.png"
products = [interval_path, constant_png, linear_path,
            OUT / "linear_rawP_selected_profiles_kpc.png", OUT / "linear_rawP_selected_profiles_Rd.png"]
manifest = {"diagnostic_only": True, "probability_cut_inclusive": CUT,
    "grid_selections": {"constant": len([r for r in DOC["grid_selected"] if r["family"] == "constant"]), "linear": len(selected)},
    "constant_source_scope": "each galaxy read only from its ledger-recorded source grid CSV",
    "constant_galaxies_with_accepted_nodes": sum(bool(v) for v in accepted.values()),
    "constant_accepted_grid_runs": len(interval_rows),
    "linear_selected_profiles": len(linear_rows),
    "source_comparison_json": str(OUT / "comparison.json"),
    "source_grid_csvs": paths, "common_q": "not available; not recomputed",
    "fully_marginalised_table2_support": "not present in source fully_marginalised CSV schema",
    "products": [{"path": str(p), "bytes": p.stat().st_size,
        "sha256": hashlib.sha256(p.read_bytes()).hexdigest()} for p in products]}
mp = OUT / "grid_figure_manifest_20261020_source_scoped_v4.json"
if mp.exists(): raise FileExistsError(mp)
mp.write_text(json.dumps(manifest, indent=2) + "\n")
print("SAVED", sum(bool(v) for v in accepted.values()), "constant galaxies /", len(interval_rows), "contiguous intervals;", len(linear_rows), "linear profiles")
print("MANIFEST", mp)
