"""Saved-data example preparation. No force/probability/search calls."""
import argparse, csv, json, math, hashlib
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from utils import load_sparc_table, get_SBdisk, compute_A_B, compute_A_beta
ROOT=Path(__file__).resolve().parent
OUT=ROOT/"outputs/paired_examples_20261021_v1"
NAMES=("NGC4068","UGC08550")
CUTS=(.15865525393145707,.022750131948179195,.5,.0013498980316300957)
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rows(p):
    with Path(p).open(newline="") as f: return list(csv.DictReader(f))
def good(r):
    try: return r["valid"].lower()=="true" and math.isfinite(float(r["log_probability"]))
    except (ValueError,KeyError): return False
def jrad(rad,u): return float(np.trapz((u/.5-1)**2,rad)/(rad[-1]-rad[0]))
def profile(row,family,rad,sb):
    mean=float(row["upsilon_disk_mean"])
    if family=="constant": return np.full(rad.shape,mean)
    slope=float(row["disk_slope"])
    A,B=compute_A_B(rad,sb)
    if family=="linear": return mean-slope*B/A+slope*rad
    _,beta=compute_A_beta(rad,sb,slope)
    return mean*A/beta*np.exp(slope*rad)
def collect():
    data={n:{} for n in NAMES}; inputs={}
    for fam in ("constant","linear","log_linear"):
        paths=sorted((ROOT/"outputs/rerun_boundary_20261006"/fam).glob("profile_grid.shard_*.csv"))
        if fam=="log_linear": paths+=sorted((ROOT/"outputs/grid_recovery_20261019").glob("log_linear*/profile_grid.shard_*.csv"))
        for p in paths:
            rr=rows(p)
            for n in NAMES:
                subset=[r for r in rr if r["galaxy"]==n]
                if not subset: continue
                inputs[str(p)]=sha(p)
                key=(fam,str(p)); data[n][key]=subset
    return data,inputs
def prepare():
    table=load_sparc_table(ROOT/"SPARC_table.txt")
    allrows,inputs=collect(); result={}; grids={}
    for n in NAMES:
        rad,sb=get_SBdisk(n,dens_root=ROOT/"../SPARC_RAR/BulgeDiskDec_LTG",table=table)
        rad,sb=np.asarray(rad,float),np.asarray(sb,float)
        gpfile=ROOT/"outputs/science_20261016/gp_d5000_a16"/(n+"_gp_search.json")
        gp=json.loads(gpfile.read_text()); inputs[str(gpfile)]=sha(gpfile)
        baseline=next(r for r in gp["records"] if r.get("source")=="SPARC_constant_0.5")
        cut=CUTS[0] if n=="NGC4068" else CUTS[3]
        report={"cut":cut,"baseline_GP_marginal_P":math.exp(baseline["marginal_logp"]),"GP_domain":gp["physical_domain"],"normalization":gp["normalization"],"families":{},"GP_cuts":{}}
        for c in CUTS:
            valid=[r for r in gp["records"] if r.get("valid") and r["marginal_logp"]>=math.log(c)]
            win=min(valid,key=lambda r:(r["J_ref"],r["original_candidate_index"])) if valid else None
            report["GP_cuts"][str(c)]=win
        for fam in ("constant","linear","log_linear"):
            sources=[(p,rr) for (f,p),rr in allrows[n].items() if f==fam]
            if not sources:
                report["families"][fam]={"status":"no saved galaxy rows found in original or recovery shards"}; continue
            # Prefer complete original grid; list all recovery coverage rather than silently pooling.
            sources.sort(key=lambda z:("rerun_boundary_20261006" not in z[0],-len(z[1]),z[0]))
            p,rr=sources[0]; grids[(n,fam)]=rr
            marg=[r for r in rr if r["probability_kind"]=="marginalised_gas_bulge"]
            valid=[r for r in marg if good(r)]
            profiles=[(r,profile(r,fam,rad,sb)) for r in valid]
            profiles=[(r,u) for r,u in profiles if np.all(np.isfinite(u)) and np.all(u>0)]
            qualifying=[(r,u) for r,u in profiles if float(r["log_probability"])>=math.log(cut)]
            selected=min(qualifying,key=lambda ru:(jrad(rad,ru[1]),int(ru[0]["grid_i"]),int(ru[0]["grid_j"]))) if qualifying else None
            base=[r for r in valid if abs(float(r["upsilon_disk_mean"])-.5)<1e-12 and (fam=="constant" or abs(float(r["disk_slope"]))<1e-12)]
            report["families"][fam]={"source":p,"coverage":[{"path":x,"rows":len(y)} for x,y in sources],"conditional_rows":sum(r["probability_kind"]=="conditional" for r in rr),"marginal_rows":len(marg),"unique_marginal_cells":len({(r["grid_i"],r["grid_j"]) for r in marg}),"valid_marginal_rows":len(valid),"baseline_P":[math.exp(float(r["log_probability"])) for r in base],"max_saved_P":max((math.exp(float(r["log_probability"])) for r in valid),default=None),"qualifying":len(qualifying),"selected":None if selected is None else {"row":selected[0],"J_rad":jrad(rad,selected[1]),"upsilon":selected[1].tolist()}}
        report["rad_kpc"]=rad.tolist(); report["sb"]=sb.tolist()
        report["GP_max_saved_P"]=max(math.exp(r["marginal_logp"]) for r in gp["records"] if r.get("valid"))
        report["GP_retained_profiles"]=gp["frontier"]; report["GP_records"]=gp["records"]; report["GP_knots"]=gp["knots_kpc"]
        result[n]=report
    return result,grids,inputs
def render(result,grids):
    for n,d in result.items():
        fig,axes=plt.subplots(3,2,figsize=(11,12),constrained_layout=True)
        for i,fam in enumerate(("constant","linear","log_linear")):
            rr=grids.get((n,fam),[])
            for j,kind in enumerate(("conditional","marginalised_gas_bulge")):
                ax=axes[i,j]; a=[r for r in rr if r["probability_kind"]==kind]; ax.set_title(f"{fam}: {kind.replace(chr(95),chr(32))}")
                ax.set_xlabel("log10 mean disk M/L")
                if not a: ax.text(.5,.5,"Saved grid unavailable",ha="center",transform=ax.transAxes); continue
                xs=sorted({float(r["log10_upsilon_disk"]) for r in a})
                selected=d["families"][fam].get("selected")
                if fam=="constant":
                    a=sorted(a,key=lambda r:float(r["log10_upsilon_disk"]))
                    ax.plot([float(r["log10_upsilon_disk"]) for r in a],[float(r["log_probability"]) if good(r) else np.nan for r in a])
                    ax.axhline(math.log(d["cut"]),ls=":",color="red"); ax.axvline(math.log10(.5),ls="--",color="black")
                    if selected: ax.axvline(float(selected["row"]["log10_upsilon_disk"]),color="orange",lw=2)
                else:
                    ys=sorted({float(r["disk_slope"]) for r in a}); z=np.full((len(ys),len(xs)),np.nan)
                    for r in a:
                        if good(r): z[ys.index(float(r["disk_slope"])),xs.index(float(r["log10_upsilon_disk"]))]=float(r["log_probability"])
                    im=ax.pcolormesh(xs,ys,np.ma.masked_invalid(z),shading="nearest",vmin=-10,vmax=0); fig.colorbar(im,ax=ax,label="saved log P; colour floor -10")
                    if np.nanmin(z)<math.log(d["cut"])<np.nanmax(z): ax.contour(xs,ys,np.ma.masked_invalid(z),levels=[math.log(d["cut"])],colors="white")
                    ax.plot(math.log10(.5),0,"kx",ms=8)
                    if selected: ax.plot(float(selected["row"]["log10_upsilon_disk"]),float(selected["row"]["disk_slope"]),"*",color="orange",mec="black",ms=12)
                    ax.set_ylabel("slope [1/kpc]" if fam=="log_linear" else "slope [M/L/kpc]")
                ax.grid(alpha=.15)
        fig.suptitle(f"{n}: raw P cut {d[chr(99)+chr(117)+chr(116)]:.6g}; illustrative, not calibrated significance")
        for ext in ("png","pdf"): fig.savefig(OUT/(n+"_saved_scans."+ext),dpi=160)
        plt.close(fig)
        fig,axes=plt.subplots(1,2,figsize=(11,4),constrained_layout=True)
        rec=d["GP_records"]; winner=d["GP_cuts"][str(d["cut"])]; winid=None if winner is None else winner["original_candidate_index"]
        for f in d["GP_retained_profiles"]:
            r=rec[f["record_index"]]; selected=r["original_candidate_index"]==winid
            axes[0].plot(d["GP_knots"],f["profile_on_photometric_grid"],color="red" if selected else "gray",lw=2.5 if selected else .7,alpha=1 if selected else .35)
        axes[0].axhline(.5,color="black",ls="--",label="baseline .5"); axes[0].set(xlabel="R [kpc]",ylabel="disk M/L",title="Saved frontier profiles; red = least qualifying Jref found")
        axes[1].scatter([r["J_ref"] for r in rec if r.get("valid")],[math.exp(r["marginal_logp"]) for r in rec if r.get("valid")],s=15)
        axes[1].axhline(d["cut"],color="red",ls=":"); axes[1].set(xlabel="fixed Jref",ylabel="saved nuisance-marginalised raw P",title="Retained candidates only; older positive-floor domain")
        fig.suptitle(n+" finite saved GP search; post-hoc illustration")
        for ext in ("png","pdf"): fig.savefig(OUT/(n+"_saved_GP."+ext),dpi=160)
        plt.close(fig)
def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--saved-plots",action="store_true"); ap.add_argument("--plan",action="store_true"); args=ap.parse_args()
    result,grids,inputs=prepare()
    for n,d in result.items():
        print(n,"cut",d["cut"],"GP baseline",d["baseline_GP_marginal_P"],"GP max",d["GP_max_saved_P"])
        for fam,v in d["families"].items(): print(fam,{k:v.get(k) for k in ("status","marginal_rows","baseline_P","max_saved_P","qualifying")},"selection",None if not v.get("selected") else {k:v["selected"]["row"].get(k) for k in ("grid_i","grid_j","log_probability","upsilon_disk_mean","disk_slope")})
    if args.saved_plots:
        OUT.mkdir(exist_ok=False); render(result,grids)
        manifest={"results":result,"inputs_sha256":inputs,"driver_sha256":sha(__file__),"no_force_probability_search_calls":True,"selection_note":"NGC4068 cut is post-hoc illustrative; neither cut is calibrated significance; grid Jrad differs from GP Jref"}
        (OUT/"manifest.json").write_text(json.dumps(manifest,indent=2)+chr(10))
if __name__=="__main__": main()
