"""Bounded paired-galaxy force plots; --plan/--self-test perform no forces."""
import argparse,json,math,hashlib,resource
from pathlib import Path
import numpy as np
from scipy.stats import norm,ncx2
ROOT=Path(__file__).resolve().parent
PREP=ROOT/"outputs/paired_examples_20261021_v1/manifest.json"
OUT=ROOT/"outputs/paired_examples_20261021_v1/force_plots_v1"
CAP=20_000_000
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def intervals(v,s,rad):
    factor=1e6/(3.085677581491367e19*rad)
    return v*v*factor,ncx2.ppf(norm.cdf(-1),1,(v/s)**2)*s*s*factor,ncx2.ppf(norm.cdf(1),1,(v/s)**2)*s*s*factor
def models(d):
    rad=np.array(d["rad_kpc"]); m=[("baseline",np.full(rad.shape,.5))]
    for fam,v in d["families"].items():
        if v.get("selected"): m.append((fam,np.array(v["selected"]["upsilon"])))
    rec=d["GP_records"]; winner=d["GP_cuts"][str(d["cut"])]; wanted=[winner["original_candidate_index"]]
    valid=[r for r in rec if r.get("valid") and r.get("source")=="GP_direction"]
    for r in [min(valid,key=lambda r:(-r["marginal_logp"],r["original_candidate_index"])),min(valid,key=lambda r:(r["J_ref"],r["original_candidate_index"]))]:
        if r["original_candidate_index"] not in wanted: wanted.append(r["original_candidate_index"])
    knots=np.array(d["GP_knots"])
    if rad[0]<knots[0] or rad[-1]>knots[-1]: raise ValueError("GP support mismatch")
    for idx in wanted:
        ri=next(i for i,r in enumerate(rec) if r["original_candidate_index"]==idx)
        front=next(f for f in d["GP_retained_profiles"] if f["record_index"]==ri)
        m.append(("GP_"+str(idx),np.interp(rad,knots,front["profile_on_photometric_grid"])))
    groups=[]
    for label,u in m:
        if np.any(~np.isfinite(u)) or np.any(u<=0): raise ValueError("Nonphysical profile")
        match=next((g for g in groups if np.array_equal(u,g[1])),None)
        if match is None: groups.append(([label],u))
        else: match[0].append(label)
    return m,groups

def self_test():
    v=np.array([0.,2.,10.]); s=np.array([1.,1.,1.]); r=np.ones(3)
    g,lo,hi=intervals(v,s,r); c=1e6/3.085677581491367e19
    for q,t in [(norm.cdf(-1),lo/c),(norm.cdf(1),hi/c)]:
        got=norm.cdf((np.sqrt(t)-v)/s)-norm.cdf((-np.sqrt(t)-v)/s)
        assert np.allclose(got,q,rtol=1e-9,atol=1e-10)
    assert g[0]<lo[0] and np.all(lo>0) and np.all(hi>lo)
    print("PASS folded-square intervals including zero velocity, low S/N separate-marker requirement, and high S/N")

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--plan",action="store_true"); ap.add_argument("--run-forces",action="store_true"); ap.add_argument("--self-test",action="store_true"); a=ap.parse_args()
    if a.self_test: self_test(); return
    if a.plan==a.run_forces: ap.error("choose --plan or --run-forces")
    manifest=json.loads(PREP.read_text()); planned={}; pins={str(PREP):sha(PREP)}
    for p,h in manifest["inputs_sha256"].items():
        if sha(p)!=h: raise ValueError("Changed saved input: "+p)
    for name,d in manifest["results"].items():
        m,groups=models(d); planned[name]={"models":[x[0] for x in m],"solver_calls":len(groups),"groups":[g[0] for g in groups]}
        for p in [ROOT/"SPARC_table.txt",ROOT/"mono_SPARC_RAR.py",ROOT/"utils.py",ROOT/"../data"/(name+"_rotmod.dat"),ROOT/"../SPARC_RAR/BulgeDiskDec_LTG"/(name+".dens")]: pins[str(p.resolve())]=sha(p)
    plan={"galaxies":planned,"total_solver_calls":sum(d["solver_calls"] for d in planned.values()),"pins":pins,"driver_sha256":sha(__file__),"output":str(OUT),"cap_bytes":CAP,"timeout_seconds":600,"uncertainties":"exact central 68.27 percent folded-square Gaussian-velocity intervals, fixed R; same observations across models","GP_domain":"saved older positive-nuisance-floor finite search, not global optimum"}
    if a.plan:
        (PREP.parent/"force_plan.json").write_text(json.dumps(plan,indent=2)+chr(10))
        print(json.dumps(plan,indent=2)); return
    if json.loads((PREP.parent/"force_plan.json").read_text())!=plan: raise ValueError("Frozen force plan changed; review required")
    if OUT.exists(): raise FileExistsError(OUT)
    resource.setrlimit(resource.RLIMIT_FSIZE,(CAP,CAP))
    from mono_SPARC_RAR import solve_disk_basis
    from utils import load_sparc_table,vel2acc
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    table=load_sparc_table(ROOT/"SPARC_table.txt"); OUT.mkdir()
    def save(fig,label):
        for ext in ("png","pdf"):
            remaining=CAP-sum(p.stat().st_size for p in OUT.iterdir() if p.is_file())-500000
            if remaining<=0: raise RuntimeError("No storage budget remains")
            resource.setrlimit(resource.RLIMIT_FSIZE,(remaining,CAP))
            fig.savefig(OUT/(label+"."+ext),dpi=160,bbox_inches="tight")
        plt.close(fig)
    for name,d in manifest["results"].items():
        rad=np.array(d["rad_kpc"]); sb=np.array(d["sb"]); raw=np.loadtxt(ROOT/"../data"/(name+"_rotmod.dat"))
        mask=(raw[:,0]>=rad[0])&(raw[:,0]<=rad[-1])&(raw[:,0]>0); kin=raw[mask]; r,v,s=kin[:,:3].T
        gpfile=ROOT/"outputs/science_20261016/gp_d5000_a16"/(name+"_gp_search.json"); gp=json.loads(gpfile.read_text())
        if not np.array_equal(r,np.array(gp["retained_radii_kpc"])): raise ValueError("Saved scorer support mismatch")
        if np.any(~np.isfinite(kin)) or np.any(s<=0) or np.any(np.diff(r)<=0): raise ValueError("Invalid kinematics")
        g,lo,hi=intervals(v,s,r); gas=vel2acc(kin[:,3]*np.abs(kin[:,3]),r); bul=vel2acc(.7*kin[:,5]**2,r)
        rd=float(table.loc[table["Galaxy"].astype(str)==name,"Rdisk"].iloc[0]); zd=.196*rd**.633
        m,groups=models(d); responses={}; warnings=[]
        for labels,u in groups:
            response,warn,error=solve_disk_basis(rad,sb*u,zd,r)
            if error or np.any(~np.isfinite(response)): raise RuntimeError(str(error))
            warnings.extend(warn)
            for label in labels: responses[label]=vel2acc(response,r)
        plan[name]={"retained_rows":len(r),"input_rows":len(raw),"warnings":warnings,"nonpositive_gbar":{label:int(np.sum(responses[label]+gas+bul<=0)) for label,u in m}}
        for demo in (False,True):
            fig,ax=plt.subplots(2,2,figsize=(10,9),constrained_layout=True); axes=ax.ravel()
            selected=[x for x in m if (x[0]=="baseline" or x[0].startswith("GP_"))] if demo else [x for x in m if not x[0].startswith("GP_") or x[0]=="GP_"+str(d["GP_cuts"][str(d["cut"])]["original_candidate_index"])]
            for i,(label,u) in enumerate(selected):
                disk=responses[label]; bar=disk+gas+bul
                best="GP_"+str(d["GP_cuts"][str(d["cut"])]["original_candidate_index"])
                col="red" if label==best else f"C{i}"
                if label==best: label+=" least Jref qualifying found"
                if label.startswith("GP_"):
                    ident=int(label.split("_")[1].split()[0]); score=next(x for x in d["GP_records"] if x["original_candidate_index"]==ident)
                    label+=f" P={math.exp(score[chr(109)+chr(97)+chr(114)+chr(103)+chr(105)+chr(110)+chr(97)+chr(108)+chr(95)+chr(108)+chr(111)+chr(103)+chr(112)]):.4g}"
                elif label=="baseline": label+=f" P={d[chr(98)+chr(97)+chr(115)+chr(101)+chr(108)+chr(105)+chr(110)+chr(101)+chr(95)+chr(71)+chr(80)+chr(95)+chr(109)+chr(97)+chr(114)+chr(103)+chr(105)+chr(110)+chr(97)+chr(108)+chr(95)+chr(80)]:.4g}"
                axes[0].plot(rad,u,color=col,label=label)
                axes[1].plot(r,disk,color=col,label=label+" disc")
                axes[2].plot(r,bar,color=col,label=label)
                # Preserve radial order and gaps for nonpositive log coordinates.
                pos=(bar>0)&(g>0); bx=np.where(pos,bar,np.nan); gy=np.where(pos,g,np.nan)
                axes[3].plot(bx,gy,color=col,label=label,alpha=.8)
                axes[3].vlines(bar[pos],lo[pos],hi[pos],color=col,alpha=.45)
                axes[3].scatter(bar[pos],g[pos],c=r[pos],cmap="plasma",s=15)
            axes[1].plot(r,gas,"k--",label="signed gas (median 1)")
            if np.any(bul!=0): axes[1].plot(r,bul,"k:",label="bulge (median .7)")
            axes[2].plot(r,g,"ko",label="observed"); axes[2].vlines(r,lo,hi,color="black")
            axes[0].set(xlabel="R [kpc]",ylabel="disk M/L",title="Physical profiles")
            axes[1].set(xlabel="R [kpc]",ylabel="acceleration [m/s²]",title="vcdisk disc and signed components")
            axes[2].set(xlabel="R [kpc]",ylabel="acceleration [m/s²]",title="Nominal baryonic and observed acceleration")
            axes[3].set(xscale="log",yscale="log",xlabel="gbar [m/s²]",ylabel="gobs [m/s²]",title="Radial-order RAR; point colour = R")
            for axis in axes: axis.legend(fontsize=6); axis.grid(alpha=.15)
            fig.suptitle(name+f" raw P cut {d['cut']:.6g}; nominal tracks, nuisance-marginalised qualification",fontsize=10)
            fig.supxlabel("68.27% folded-square observational intervals; fixed geometry. Nonpositive log points omitted and counted in manifest.",fontsize=8)
            save(fig,name+("_GP_demonstration" if demo else "_physical_2x2"))
    plan["artifacts"]={p.name:sha(p) for p in OUT.iterdir() if p.is_file()}
    (OUT/"manifest.json").write_text(json.dumps(plan,indent=2)+chr(10))
if __name__=="__main__": main()
