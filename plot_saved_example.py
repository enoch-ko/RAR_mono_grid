#!/usr/bin/env python3
"""Render saved RAR examples; no force, search or monotonicity-probability imports."""
import argparse,csv,json,math,hashlib,shutil
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
REPO=Path(__file__).resolve().parent
FAMILIES=("constant","linear","log_linear")
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def bracket(records,cut):
    valid=sorted((math.exp(r["marginal_logp"]),r["J_ref"],r["original_candidate_index"]) for r in records if r.get("valid"))
    return {"below":max((r for r in valid if r[0]<cut),default=None),"above":min((r for r in valid if r[0]>=cut),default=None),"valid_retained":len(valid),"within_0p02":sum(abs(r[0]-cut)<=.02 for r in valid)}
def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--galaxy",default="NGC4068",choices=("NGC4068","UGC08550"))
    ap.add_argument("--input-root",type=Path,default=REPO/"outputs/paired_examples_20261021_v1")
    ap.add_argument("--output-root",type=Path,required=True,help="New directory; existing paths refused")
    ap.add_argument("--repository-root",type=Path,default=REPO,help="Relocate grid dependencies after manual copying")
    ap.add_argument("--kind",choices=("saved","scans","gp","physical-assets"),default="saved")
    ap.add_argument("--font-size",type=float,default=10)
    ap.add_argument("--line-width",type=float,default=1.5)
    ap.add_argument("--colour-floor",type=float,default=-10)
    ap.add_argument("--cmap",default="viridis")
    ap.add_argument("--gp-probability-limits",type=float,nargs=2,help="Display zoom only; does not resample/reselect")
    ap.add_argument("--self-test",action="store_true")
    a=ap.parse_args()
    if a.self_test:
        b=bracket([{"valid":True,"marginal_logp":math.log(p),"J_ref":j,"original_candidate_index":i} for p,j,i in [(.1,1,1),(.2,2,2)]],.15)
        assert b["below"][2]==1 and b["above"][2]==2
        print("PASS saved bracket and inclusive cut comparison"); return
    if a.output_root.exists(): raise FileExistsError(a.output_root)
    mp=a.input_root/"manifest.json"; m=json.loads(mp.read_text()); d=m["results"][a.galaxy]; cut=d["cut"]
    receipt={"input_manifest":str(mp.resolve()),"input_sha256":sha(mp),"script_sha256":sha(__file__),"arguments":{k:str(v) for k,v in vars(a).items()},"galaxy":a.galaxy,"cut":cut,"audit":bracket(d["GP_records"],cut),"draw_mode":a.kind,"inputs":{},"new_force_search_probability_evaluations":0}
    def grid_path(original):
        prefix="/mnt/users/koe/RAR_mono_grid/"
        return a.repository_root/original[len(prefix):] if original.startswith(prefix) else Path(original)
    def save(fig,name):
        for ext in ("png","pdf"): fig.savefig(a.output_root/(name+"."+ext),dpi=180,bbox_inches="tight")
        plt.close(fig)
    plt.rcParams.update({"font.size":a.font_size,"lines.linewidth":a.line_width})
    a.output_root.mkdir(parents=True,exist_ok=False)
    if a.kind in ("saved","scans"):
        fig,axs=plt.subplots(3,2,figsize=(12,12),constrained_layout=True)
        for i,fam in enumerate(FAMILIES):
            v=d["families"][fam]; source=v["source"]; p=grid_path(source)
            if sha(p)!=m["inputs_sha256"][source]: raise ValueError("Grid hash mismatch: "+str(p))
            receipt["inputs"][str(p)]=sha(p)
            with p.open(newline="") as f: rr=[r for r in csv.DictReader(f) if r["galaxy"]==a.galaxy]
            for j,kind in enumerate(("conditional","marginalised_gas_bulge")):
                ax=axs[i,j]; rows=[r for r in rr if r["probability_kind"]==kind]; selected=v["selected"]["row"]
                xs=sorted({float(r["log10_upsilon_disk"]) for r in rows})
                def value(r):
                    x=float(r["log_probability"])
                    return x if r["valid"].lower()=="true" and math.isfinite(x) else np.nan
                if fam=="constant":
                    rows.sort(key=lambda r:float(r["log10_upsilon_disk"]))
                    ax.plot([float(r["log10_upsilon_disk"]) for r in rows],[value(r) for r in rows])
                    ax.axhline(math.log(cut),color="red",ls=":",label="raw P cut")
                    ax.axvline(math.log10(.5),color="black",ls="--",label="baseline .5")
                    ax.axvline(float(selected["log10_upsilon_disk"]),color="orange",label="selected minimum Jrad")
                    ax.set_ylabel("saved ln P"); ax.legend(fontsize=7)
                else:
                    ys=sorted({float(r["disk_slope"]) for r in rows}); z=np.full((len(ys),len(xs)),np.nan)
                    for r in rows: z[ys.index(float(r["disk_slope"])),xs.index(float(r["log10_upsilon_disk"]))]=value(r)
                    cmap=plt.get_cmap(a.cmap).copy(); cmap.set_bad("lightgray")
                    im=ax.pcolormesh(xs,ys,np.ma.masked_invalid(z),shading="nearest",vmin=a.colour_floor,vmax=0,cmap=cmap)
                    ax.contour(xs,ys,np.ma.masked_invalid(z),levels=[math.log(cut)],colors="white")
                    ax.plot(math.log10(.5),0,"kx",label="baseline .5")
                    ax.plot(float(selected["log10_upsilon_disk"]),float(selected["disk_slope"]),"*",color="orange",mec="black",ms=10,label="selected min Jrad")
                    ax.legend(fontsize=7); ax.set_ylabel("slope [1/kpc]" if fam=="log_linear" else "slope [M/L/kpc]")
                    fig.colorbar(im,ax=ax,label="saved ln P (colour clipped only)")
                ax.set_title(fam+": "+kind.replace("_"," ")); ax.set_xlabel("log10 mean disk M/L")
        fig.suptitle(a.galaxy+f" raw P cut={cut:.8g}; post-hoc illustration, not calibrated significance")
        save(fig,a.galaxy+"_saved_scans")
    if a.kind in ("saved","gp"):
        records=d["GP_records"]; valid=[r for r in records if r.get("valid")]; winner=min((r for r in valid if r["marginal_logp"]>=math.log(cut)),key=lambda r:(r["J_ref"],r["original_candidate_index"]))
        fig,axs=plt.subplots(1,2,figsize=(12,5),constrained_layout=True)
        for f in d["GP_retained_profiles"]:
            r=records[f["record_index"]]; best=r["original_candidate_index"]==winner["original_candidate_index"]
            axs[0].plot(d["GP_knots"],f["profile_on_photometric_grid"],color="red" if best else "gray",alpha=1 if best else .25,lw=2.5 if best else .7,label="least qualifying Jref found" if best else None)
        axs[0].axhline(.5,color="black",ls="--",label="baseline .5"); axs[0].legend(fontsize=8)
        axs[0].set(xlabel="R [kpc]",ylabel="disk M/L",title="All saved frontier profiles")
        axs[1].scatter([r["J_ref"] for r in valid],[math.exp(r["marginal_logp"]) for r in valid],s=20,label="all retained records")
        axs[1].scatter(winner["J_ref"],math.exp(winner["marginal_logp"]),color="red",marker="*",s=100,label="qualifying winner")
        axs[1].axhline(cut,color="red",ls=":",label="raw P cut")
        for key,b in receipt["audit"].items():
            if key in ("below","above") and b: axs[1].annotate(str(b[2]),(b[1],b[0]),xytext=(4,5),textcoords="offset points",fontsize=8)
        if a.gp_probability_limits: axs[1].set_ylim(*a.gp_probability_limits)
        axs[1].set(xlabel="fixed Jref",ylabel="saved nuisance-marginalised raw P",title="Retained frontier, not sampling density"); axs[1].legend(fontsize=8)
        fig.suptitle(a.galaxy+f" cut={cut:.8g}; finite older positive-floor search")
        save(fig,a.galaxy+"_saved_GP")
    if a.kind=="physical-assets":
        fm=a.input_root/"force_plots_v1/manifest.json"; fd=json.loads(fm.read_text()); receipt["force_manifest_sha256"]=sha(fm)
        receipt["draw_mode"]="byte-identical reuse of completed figure assets; no numerical redraw"
        for name,h in fd["artifacts"].items():
            if name.startswith(a.galaxy+"_"):
                p=fm.parent/name
                if sha(p)!=h: raise ValueError("Figure hash mismatch")
                shutil.copy2(p,a.output_root/name)
        receipt["limitation"]="Original force run did not persist numerical arrays. Internal physical/RAR styling requires separately approved force regeneration with array persistence."
    receipt["artifacts"]={p.name:sha(p) for p in a.output_root.iterdir() if p.is_file()}
    (a.output_root/"render_manifest.json").write_text(json.dumps(receipt,indent=2)+chr(10))
    print(json.dumps(receipt["audit"])); print(a.output_root)
if __name__=="__main__": main()
