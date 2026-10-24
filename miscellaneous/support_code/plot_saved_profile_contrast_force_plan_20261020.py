#!/usr/bin/env python3
"""Saved-data plot and force-only plan for the D564-8 comparison; never runs forces."""
import argparse, csv, json, math
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "outputs/saved_profile_comparison_20261019_exactcut_v1/plot_verification_20261020_v1"
CUT = 0.0013498980316300957

def mock_test():
    calls = {("constant_0p5", "drop_exact_R0"), ("linear_v7", "drop_exact_R0"), ("gp_frontier", "drop_exact_R0")}
    assert len(calls) == 3
    print("PASS: four plot labels; three unique force profiles/policy pairs. No solver invoked.")

def saved_inputs():
    lin = ROOT / "outputs/saved_profile_comparison_20261019_exactcut_v1/linear_marginalised_rawP_leastJrad_selected_profiles_v7.csv"
    gp = ROOT / "outputs/science_20261016/gp_d5000_a16/D564-8_gp_search.json"
    dens = ROOT / "../SPARC_RAR/BulgeDiskDec_LTG/D564-8.dens"
    ledger = ROOT / "outputs/rerun_boundary_20261006/log_linear"
    return lin, gp, dens, ledger

def read_loglinear(ledger):
    rows=[]
    for p in sorted(ledger.glob("profile_grid.shard_*.csv")):
        with p.open(newline="") as f:
            for r in csv.DictReader(f):
                try:
                    if r["galaxy"] == "D564-8" and r["probability_kind"] == "marginalised_gas_bulge" and r["valid"].lower() == "true" and math.exp(float(r["log_probability"])) >= CUT:
                        rows.append(r)
                except (KeyError, ValueError, OverflowError):
                    pass
    return rows

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--mock-test", action="store_true")
    ap.add_argument("--plan", action="store_true", help="print proposed calls only; does not run a solver")
    ap.add_argument("--plot", action="store_true", help="plot saved M/L curves only")
    a=ap.parse_args()
    if a.mock_test: mock_test(); return
    linp,gpp,densp,ledger=saved_inputs()
    assert linp.exists() and gpp.exists() and densp.exists(), "Required saved input missing"
    lin=next(r for r in csv.DictReader(linp.open()) if r["galaxy"] == "D564-8")
    logrows=read_loglinear(ledger)
    if not logrows: raise SystemExit("No valid qualifying marginalised log-linear rows")
    d=np.loadtxt(densp, comments="#")
    rad,sb=d[:,0],d[:,1]
    A=np.trapz(sb*rad,rad)
    def jrad(r):
        s=float(r["disk_slope"]); p=float(r["upsilon_disk_mean"]); span=rad[-1]-rad[0]; x=abs(s)*span
        beta=np.trapz(sb*np.exp(s*rad)*rad,rad); top=max(p*A/beta*np.exp(s*rad[0]),p*A/beta*np.exp(s*rad[-1]))
        m1=top*(-math.expm1(-x)/x) if x else top
        m2=top*top*(-math.expm1(-2*x)/(2*x)) if x else top*top
        return max(0.,4*m2-4*m1+1)
    best=min(logrows,key=jrad)
    gp=json.load(gpp.open())
    rec=next(x for x in gp["records"] if x.get("original_candidate_index")==185803 and x.get("source")=="GP_direction")
    print("Saved-data summary only; no solver/probability/search call.")
    print(f"loglinear: {len(logrows)} qualifying valid cells; min Jrad={jrad(best):.12g}; grid=({best['grid_i']},{best['grid_j']}), p={best['upsilon_disk_mean']}, slope={best['disk_slope']}/kpc")
    print(f"linear v7: rawP={lin.get('raw_probability','see CSV')}; grid=({lin.get('grid_i')},{lin.get('grid_j')}), a={lin.get('a_ml_at_R0')}, b={lin.get('b_ml_per_kpc')}/kpc")
    print(f"GP retained: original index={rec['original_candidate_index']}; source={rec['source']}; saved candidate record available")
    if a.plan:
        print("PLAN ONLY: 3 unique profiles × drop_exact_R0 policy: constant 0.5, linear v7, and GP retained candidate. Log-linear is included in the saved-profile plot only.")
        print("Proposed solver calls: one vcdisk basis solve per unique profile (3 total); no RAR probability integration. Runtime not benchmarked; run only after approval.")
        print(f"Output location: {OUT}; reserve <1 MB for one PNG plus JSON manifest (estimate).")
    if a.plot:
        import matplotlib.pyplot as plt
        gpfront=next(x for x in gp["frontier"] if x.get("record_index")==gp["records"].index(rec))
        gr=np.asarray(gp["knots_kpc"],float); gv=np.asarray(gpfront["profile_on_photometric_grid"],float)
        p=float(best["upsilon_disk_mean"]); s=float(best["disk_slope"]); beta=np.trapz(sb*np.exp(s*rad)*rad,rad); logy=p*A/beta*np.exp(s*rad)
        rlin=np.linspace(rad[0],rad[-1],200); lin_y=float(lin["a_ml_at_R0"])+float(lin["b_ml_per_kpc"])*rlin
        fig,ax=plt.subplots(figsize=(7,4.5)); ax.axhline(.5,color="black",ls="--",lw=1,label="constant 0.5")
        ax.plot(rlin,lin_y,label="linear v7"); ax.plot(rad,logy,label="log-linear least-Jrad")
        ax.plot(gr,gv,label="GP retained",alpha=.8); ax.set(xlabel="R [kpc]",ylabel=r"$\Upsilon_{disk}$",title="D564-8 saved-profile contrast (not force-recomputed)")
        ax.legend(); fig.tight_layout(); OUT.mkdir(parents=True,exist_ok=True); fig.savefig(OUT/"D564-8_saved_profiles.png",dpi=180); plt.close(fig)
        print(f"Saved profile-only figure: {OUT/'D564-8_saved_profiles.png'}")
    if not a.plan and not a.plot: ap.print_help()

if __name__=="__main__": main()
