#!/usr/bin/env python3
"""Six threshold summaries from completed saved grid shards; no RAR reruns.
Reuse the established contour extraction and shard-argument normalization.
Each figure separates fixed gas/bulge and nuisance-marginalized probabilities.
Constant profiles are 1D: display threshold-allowed intervals per galaxy row.
"""
import argparse
import collections
import hashlib
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import plot_linear_grid_contours as old
from merge_run_shards import shared_arguments
from utils import load_sparc_table

ROOT = Path(__file__).resolve().parent
PROFILES = ("constant", "linear", "log_linear")
KINDS = old.KINDS


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def threshold_intervals(x, z, valid, level):
    """Linear interpolation in unclipped logP; never bridge invalid/missing cells."""
    pieces = [(float(a), float(a)) for a, b, ok in zip(x, z, valid) if ok and b >= level]
    for i in range(len(x)-1):
        if not (valid[i] and valid[i+1]):
            continue
        a, b, u, v = x[i], x[i+1], z[i], z[i+1]
        if u >= level and v >= level:
            pieces.append((float(a), float(b)))
        elif (u < level <= v) or (v < level <= u):
            crossing = a + (b-a)*(level-u)/(v-u)
            pieces.append((float(crossing if u < level else a),
                           float(b if u < level else crossing)))
    merged = []
    for a, b in sorted(pieces):
        if merged and a <= merged[-1][1]+1e-12:
            merged[-1][1] = max(b, merged[-1][1])
        else:
            merged.append([a, b])
    return merged


def shade_threshold(ax, x, y, z, valid, level, color, alpha):
    """Shade P>=cut from the masked logP field, never by closing open outlines.

    Use the same Matplotlib corner-mask convention as the existing contours:
    only valid three-corner triangles may survive next to a masked grid node.
    Masked nodes, gaps and values below the cut are not filled.
    """
    valid = np.asarray(valid, bool) & np.isfinite(z)
    finite = np.asarray(z)[valid]
    if not len(finite) or np.max(finite) < level or alpha == 0:
        return []
    masked = np.ma.masked_where(~valid, z)
    upper = max(float(np.max(finite)), float(level)) + 1.
    cs = ax.contourf(x, y, masked, levels=[level, upper], colors=[color],
                     alpha=alpha, corner_mask=True, antialiased=False, zorder=1)
    kinds = cs.allkinds[0] if cs.allkinds is not None else [None]*len(cs.allsegs[0])
    return [{"vertices": points.tolist(), "codes": None if codes is None else codes.tolist()}
            for points, codes in zip(cs.allsegs[0], kinds) if len(points) >= 3]


def read_profile(folder, profile, recorded_inputs=None):
    metas = sorted(folder.glob("run_metadata.shard_*.json"))
    if recorded_inputs is not None:
        metas = []
        for item in recorded_inputs:
            path = Path(item["metadata_path"])
            grid = Path(item["grid_csv"])
            if path.parent.resolve() != folder.resolve() or grid.parent.resolve() != folder.resolve():
                raise ValueError("snapshot input root differs from requested root")
            if sha(path) != item["metadata_sha256"] or sha(grid) != item["grid_sha256"]:
                raise RuntimeError("snapshot source changed; do not silently mix input versions")
            metas.append(path)
        metas.sort()
    if not metas:
        raise ValueError(f"no finished shard metadata: {folder}")
    chunks, inputs, seen = [], [], set()
    reference = None
    for path in metas:
        metadata = json.loads(path.read_text())
        shard = int(metadata["arguments"]["shard_id"])
        grid_path = folder / f"profile_grid.shard_{shard:02d}.csv"
        full_path = folder / f"fully_marginalised.shard_{shard:02d}.csv"
        if not grid_path.is_file() or not full_path.is_file() or not metadata.get("finished_utc"):
            raise ValueError(f"incomplete finished shard {shard}")
        signature = (shared_arguments(metadata), metadata["source_sha256"])
        if reference is None:
            reference = signature
        elif signature != reference:
            raise ValueError(f"incompatible settings/source in {path}")
        if metadata["profile"] != profile:
            raise ValueError("profile metadata mismatch")
        start_hash = sha(grid_path)
        columns = list(pd.read_csv(grid_path, nrows=0).columns)
        use = list(old.REQUIRED) + [c for c in ("invalid_reason", "warnings") if c in columns]
        frame = pd.read_csv(grid_path, usecols=use, low_memory=False)
        if set(frame["profile"]) != {profile}:
            raise ValueError("CSV profile mismatch")
        galaxies = set(frame["galaxy"])
        if seen & galaxies or len(galaxies) != int(metadata["galaxy_count"]):
            raise ValueError("duplicate coverage or metadata galaxy count mismatch")
        seen.update(galaxies)
        inputs.append({"shard":shard, "grid_csv":str(grid_path.resolve()),
                       "grid_sha256":start_hash, "metadata_path":str(path.resolve()),
                       "metadata_sha256":sha(path), "metadata":metadata,
                       "galaxies":sorted(galaxies)})
        chunks.append(frame)
        if sha(grid_path) != start_hash:
            raise RuntimeError("grid changed while being read")
    return pd.concat(chunks, ignore_index=True), inputs, reference[0]


def classify(values, level):
    if not len(values):
        return "no_valid_cells"
    if values.min() > level:
        return "all_valid_cells_above_threshold"
    if values.max() < level:
        return "all_valid_cells_below_threshold"
    return "threshold_bracketed"


def summarize(profile, frame, inputs, args, catalog, shade=False, fill_alpha=.025):
    known = set(catalog)
    present = set(frame["galaxy"])
    if present-known:
        raise ValueError("grid contains galaxies outside the reference catalogue")
    ngrid = int(args["grid_size"])
    groups = {(g,k): part for (g,k),part in frame.groupby(["galaxy", "probability_kind"], sort=False)}
    statuses, geometry, filled_geometry = [], [], []
    xlimits = (math.log10(float(args["disk_ml"]))-float(args["disk_grid_half_width_dex"]),
               math.log10(float(args["disk_ml"]))+float(args["disk_grid_half_width_dex"]))
    for sigma in (2,5):
        level = old.log_tail(sigma)
        fig, axes = plt.subplots(1,2,figsize=(14,15 if profile=="constant" else 6),
                                 sharex=True, sharey=True, layout="constrained")
        counts = []
        for ax, kind in zip(axes, KINDS):
            local = []
            for index, galaxy in enumerate(catalog, start=1):
                row = {"profile":profile,"sigma":sigma,"probability_kind":kind,
                       "galaxy":galaxy,"catalogue_index":index}
                group = groups.get((galaxy,kind))
                if group is None:
                    row["status"] = "missing_saved_grid"
                    statuses.append(row)
                    local.append(row)
                    continue
                group = group.sort_values(["grid_i","grid_j"])
                expected = ngrid if profile=="constant" else ngrid**2
                if (len(group)!=expected or group.duplicated(["grid_i","grid_j"]).any()
                        or set(group["grid_j"])!=set(range(ngrid))
                        or set(group["grid_i"])!=set(range(1 if profile=="constant" else ngrid))):
                    row["status"]="incomplete_or_duplicate_grid"
                    statuses.append(row)
                    local.append(row)
                    continue
                ok = old.valid_mask(group["valid"]).to_numpy() & np.isfinite(group["log_probability"].to_numpy(float))
                reason = group.get("invalid_reason",pd.Series("",index=group.index)).fillna("").astype(str)
                row.update(valid_cells=int(ok.sum()), physical_excluded_cells=int(reason.eq("nonpositive_or_nonfinite_disk_ml").sum()),
                           other_invalid_cells=int((~ok & ~reason.eq("nonpositive_or_nonfinite_disk_ml").to_numpy()).sum()))
                row["threshold_relation"] = classify(group.loc[ok,"log_probability"].to_numpy(float),level)
                color = plt.get_cmap("turbo")((index-1)/max(len(catalog)-1,1))
                if profile=="constant":
                    x=np.log10(group["upsilon_disk_mean"].to_numpy(float))
                    z=group["log_probability"].to_numpy(float)
                    segments=threshold_intervals(x,z,ok,level)
                    for lo,hi in segments:
                        ax.plot([lo,hi],[index,index],color=color,lw=2,marker="|",ms=3)
                    touches=any(abs(lo-x[0])<1e-10 or abs(hi-x[-1])<1e-10 for lo,hi in segments)
                    row["status"]=("allowed_intervals_touch_grid_boundary" if touches else "allowed_intervals_interior") if segments else row["threshold_relation"]
                else:
                    x,y,z,mask=old.make_grid(group)
                    if shade:
                        filled = shade_threshold(ax,x,y,z,mask,level,color,fill_alpha)
                        row["filled_region_count"] = len(filled)
                        row["threshold_satisfying_nodes"] = int(np.sum(mask & (z >= level)))
                        if filled:
                            filled_geometry.append({"profile":profile,"sigma":sigma,
                                "probability_kind":kind,"galaxy":galaxy,"regions":filled})
                    segments,status=old.contour_segments(ax,x,y,z,mask,level,color,"solid",linewidth=.55)
                    row["status"]=status if segments else row["threshold_relation"] if row["threshold_relation"]!="threshold_bracketed" else "no_drawable_segment"
                    segments=[s.tolist() for s in segments]
                row["segment_count"]=len(segments)
                statuses.append(row)
                local.append(row)
                if segments:
                    geometry.append({"profile":profile,"sigma":sigma,"probability_kind":kind,"galaxy":galaxy,"segments":segments})
            summary=dict(collections.Counter(r["status"] for r in local))
            counts.append({"probability_kind":kind,"statuses":summary})
            drawable=sum(r.get("segment_count",0)>0 for r in local)
            available=sum(r["status"]!="missing_saved_grid" for r in local)
            filled_count=sum(r.get("filled_region_count",0)>0 for r in local)
            title=f"{kind}: {available}/{len(catalog)} saved, {drawable} contours"
            if shade and profile!="constant":
                title+=f", {filled_count} shaded"
            ax.set_title(title,fontsize=10)
            ax.set_xlim(*xlimits)
            ax.set_xlabel("log10 luminosity-weighted mean disc M/L")
            ax.grid(alpha=.15)
        if profile=="constant":
            axes[0].set_ylim(len(catalog)+1,0)
            tick=np.arange(1,len(catalog)+1,10)
            axes[0].set_yticks(tick,[f"{i}: {catalog[i-1]}" for i in tick],fontsize=7)
            axes[0].set_ylabel("Catalogue galaxy row; intervals satisfy threshold")
        else:
            axes[0].set_ylim(float(args["slope_min"]),float(args["slope_max"]))
            axes[0].set_ylabel("Linear slope [M/L per kpc]" if profile=="linear" else "d ln(M/L) / dR [kpc^-1]")
        fig.suptitle(f"{profile}: {sigma}-sigma reporting level, P(monotonic) >= {math.exp(level):.10g}; ln P = {level:.8f}\n"
                     "One-sided Gaussian-tail label, NOT calibrated theory significance; current boundary batch"
                     + (f"\nFill: P >= cut (alpha={fill_alpha:g}); outlines: P = cut; colours identify galaxies."
                        "\nUnshaded may be below cut or masked. Overlap is not posterior density."
                        if shade and profile!="constant" else ""),fontsize=10)
        filename=f"{profile}_{sigma}sigma_summary.png"
        fig.savefig(args_output/filename,dpi=160)
        plt.close(fig)
        print(profile,sigma,json.dumps(counts),flush=True)
    report={"profile":profile,"available_galaxies":sorted(present),"missing_galaxies":sorted(known-present),
            "inputs":inputs,"statuses":statuses,"geometry":geometry,"filled_geometry":filled_geometry,
            "merged_csv_present":(Path(inputs[0]["grid_csv"]).parent/"profile_grid.csv").exists()}
    return report


def shading_self_test():
    """Analytic masks and rendered pixels test direction, edges and invalid holes."""
    x=np.array([0.,.5,1.])
    X,Y=np.meshgrid(x,x)
    valid=np.ones_like(X,dtype=bool)
    def panel():
        fig=plt.figure(figsize=(2,2),dpi=80)
        ax=fig.add_axes([0,0,1,1])
        ax.set_axis_off()
        ax.set_xlim(0,1)
        ax.set_ylim(0,1)
        return fig,ax
    def is_filled(fig,ax,point):
        fig.canvas.draw()
        rgba=np.asarray(fig.canvas.buffer_rgba())
        px,py=ax.transData.transform(point)
        rgb=rgba[rgba.shape[0]-1-int(py),int(px),:3]
        return bool(np.any(rgb < 240))
    fig,ax=panel()
    regions=shade_threshold(ax,x,x,X-.5,valid,0.,"blue",1.)
    assert regions and is_filled(fig,ax,(.8,.3)) and not is_filled(fig,ax,(.2,.3))
    plt.close(fig)
    fig,ax=panel()
    regions=shade_threshold(ax,x,x,np.ones_like(X),valid,0.,"blue",1.)
    assert regions and is_filled(fig,ax,(.1,.1)) and is_filled(fig,ax,(.9,.9))
    plt.close(fig)
    fig,ax=panel()
    assert not shade_threshold(ax,x,x,-np.ones_like(X),valid,0.,"blue",1.)
    assert not is_filled(fig,ax,(.5,.5))
    plt.close(fig)
    fig,ax=panel()
    hole=valid.copy()
    hole[1,1]=False
    regions=shade_threshold(ax,x,x,np.ones_like(X),hole,0.,"blue",1.)
    assert regions and not is_filled(fig,ax,(.5,.5)) and is_filled(fig,ax,(.1,.1))
    plt.close(fig)
    fig,ax=panel()
    assert not shade_threshold(ax,x,x,np.ones_like(X),np.zeros_like(valid),0.,"blue",1.)
    assert not shade_threshold(ax,x,x,np.ones_like(X),valid,0.,"blue",0.)
    plt.close(fig)
    fig,ax=panel()
    regions=shade_threshold(ax,x,x,np.zeros_like(X),valid,0.,"blue",1.)
    assert regions and is_filled(fig,ax,(.5,.5))
    plt.close(fig)
    fig,ax=panel()
    five=np.linspace(0,1,5)
    z=np.tile(np.array([1.,-1.,-1.,-1.,1.]),(5,1))
    regions=shade_threshold(ax,five,five,z,np.ones_like(z,bool),0.,"blue",1.)
    assert len(regions)==2 and is_filled(fig,ax,(.05,.5)) and is_filled(fig,ax,(.95,.5))
    assert not is_filled(fig,ax,(.5,.5))
    plt.close(fig)
    print("Eight shading/render checks passed, including threshold equality and masked holes.")


def self_test():
    x=np.array([-1.,0.,1.])
    assert threshold_intervals(x,np.array([-2.,0.,-2.]),np.ones(3,bool),-1.)==[[-.5,.5]]
    assert threshold_intervals(x,np.ones(3),np.array([True,False,True]),0.)==[[-1.,-1.],[1.,1.]]
    assert threshold_intervals(x,np.ones(3),np.ones(3,bool),0.)==[[-1.,1.]]
    assert not threshold_intervals(x,-np.ones(3),np.ones(3,bool),0.)
    assert classify(np.array([]),0.)=="no_valid_cells"
    assert abs(math.exp(old.log_tail(2))-.02275013194817922)<1e-15
    print("Six threshold/interval checks passed; no grid/probability calculations.")
    shading_self_test()


def main():
    global args_output
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root",type=Path,default=ROOT/"outputs/rerun_boundary_20261006")
    parser.add_argument("--output-root",type=Path,default=ROOT/"outputs/profile_grid_summaries/boundary_20261013")
    parser.add_argument("--self-test",action="store_true")
    parser.add_argument("--profiles",choices=PROFILES,nargs="+",default=list(PROFILES))
    parser.add_argument("--shade-allowed",action="store_true",help="Very faint P>=threshold fills for 2D profiles")
    parser.add_argument("--fill-alpha",type=float,default=.025)
    parser.add_argument("--reuse-input-report",type=Path,help="Reuse exactly the source files/hashes in a previous summary report")
    options=parser.parse_args()
    if options.self_test:
        self_test()
        return
    if not np.isfinite(options.fill_alpha) or not 0 <= options.fill_alpha <= 1:
        parser.error("fill alpha must be finite and between zero and one")
    if len(set(options.profiles)) != len(options.profiles):
        parser.error("profile selection must be unique")
    snapshot = json.loads(options.reuse_input_report.read_text()) if options.reuse_input_report else None
    snapshots = {p["profile"]:p for p in snapshot["profiles"]} if snapshot else {}
    if snapshot and any(p not in snapshots for p in options.profiles):
        parser.error("reference report does not contain every selected profile")
    if options.output_root.exists():
        raise FileExistsError("output root exists; choose a new destination")
    args_output=options.output_root
    args_output.mkdir(parents=True)
    catalog=load_sparc_table()["Galaxy"].astype(str).tolist()
    results=[]
    for profile in options.profiles:
        recorded=snapshots[profile]["inputs"] if snapshot else None
        frame,inputs,args=read_profile(options.input_root/profile,profile,recorded)
        results.append(summarize(profile,frame,inputs,args,catalog,options.shade_allowed,options.fill_alpha))
    for result in results:
        for item in result["inputs"]:
            if sha(Path(item["grid_csv"]))!=item["grid_sha256"] or sha(Path(item["metadata_path"]))!=item["metadata_sha256"]:
                raise RuntimeError("source shard changed during plotting")
    report={"script_sha256":sha(Path(__file__)),"contour_helper_sha256":sha(Path(old.__file__)),
            "catalogue_order":catalog,"profiles":results,"probability_kinds":list(KINDS),
            "thresholds":{str(s):{"P_min":math.exp(old.log_tail(s)),"lnP":old.log_tail(s)} for s in (2,5)},
            "direction":"P >= threshold", "no_heatmap_clipping_applied":True,
            "contours_may_end_at_invalid_or_physical_masks":True,"CSV_merger_invoked":False,
            "no_RAR_recalculation":True,
            "shading":{"enabled":options.shade_allowed,"alpha":options.fill_alpha,
                       "direction":"unclipped logP >= threshold", "corner_mask":True,
                       "invalid_nodes_masked":True,"zorder":1,"outlines_preserved":True,
                       "overlap_is_not_posterior_density":True},
            "input_snapshot_report":None if not snapshot else {"path":str(options.reuse_input_report.resolve()),
                                        "sha256":sha(options.reuse_input_report)},
            "legacy_linear_cache":"outputs/production_linear_51/contour_plots preserved; different input batch, geometry not mixed"}
    with (args_output/"coverage_and_geometry.json").open("x") as stream:
        json.dump(report,stream,indent=2,allow_nan=False)
        stream.write("\n")


if __name__=="__main__":
    main()
