#!/usr/bin/env python3
"""One frozen common SPARC baseline per galaxy; no production submissions.
Build requires explicit B, seed, normalization manifest and NEW output path.
Dry-run reads/validates inputs only; no Casertano solve, draws or output writes.
Legacy GP/QP compatibility is deliberately blocked pending input remediation.
"""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import signal
import numpy as np
import common_isotonic_reference as common
import utils

ROOT=Path(__file__).resolve().parent
BAD_INVERSE_HASH="c2d4cae176bc05ae718da306548add5a4976a623a6357ccf222f1d8e58835a47"
NUISANCE={"gas_median":1.,"gas_dex":.04,"bulge_median":.7,"bulge_dex":.1,"samples":9}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_inputs(galaxy, normalization_path, data_root, dens_root, table_path):
    table=utils.load_sparc_table(table_path)
    matched=table.loc[table.Galaxy==galaxy]
    if len(matched)!=1:raise ValueError("galaxy catalogue join absent or nonunique")
    cat=matched.iloc[0]
    density_path=Path(dens_root)/(galaxy+".dens")
    density=np.loadtxt(density_path,ndmin=2)
    if density.shape[1]!=3 or len(density)<3 or not np.isfinite(density).all():
        raise ValueError("invalid density schema")
    rad=density[:,0]
    if np.any(rad<0) or np.any(np.diff(rad)<=0) or np.any(density[:,1:]<0):
        raise ValueError("invalid full photometric profile")
    norm=json.loads(Path(normalization_path).read_text())
    if norm.get("policy")!="disc_to_catalogue_total_minus_published_bulge_v1":
        raise ValueError("explicit supported normalization policy required")
    rec=norm["galaxies"][galaxy]
    if (not np.isclose(float(cat.L),rec["Ltotal_1e9"],rtol=0,atol=1e-10) or
        not np.isclose(float(cat.D),rec["reference_distance_Mpc"],rtol=0,atol=1e-10)):
        raise ValueError("normalization catalogue reference mismatch")
    ld=float(rec["Ldisk_1e9"]); lb=float(rec["Lbulge_1e9"])
    if not (np.isfinite([ld,lb,rec["kappa"]]).all() and ld>0 and lb>=0 and
            np.isclose(ld+lb,float(cat.L),rtol=1e-12) and
            np.isclose(ld/float(cat.L),rec["kappa"],rtol=1e-12)):
        raise ValueError("inconsistent disc/bulge normalization record")
    rr,base=utils.get_SBdisk(galaxy,dens_root=dens_root,table=table)
    if not np.array_equal(rr,rad):raise ValueError("density reader radius alignment mismatch")
    brightness=base*rec["kappa"]
    if not np.isclose(utils.compute_A_B(rad,brightness)[0],ld*1e9,rtol=1e-12):
        raise ValueError("disc luminosity integral mismatch")
    rotmod=Path(data_root)/(galaxy+"_rotmod.dat")
    data,parse=common.read_rotmod(rotmod,(rad[0],rad[-1]))
    positive=rad[rad>0]
    if len(positive)<2 or np.any(data["Rad"]<positive[0]):
        raise ValueError("retained kinematics outside positive-radius force support; no clamping")
    distance=float(rotmod.read_text().splitlines()[0].split("=")[1].split()[0])
    if not np.isclose(distance,float(cat.D),rtol=0,atol=.005):
        raise ValueError("rotmod/catalogue distance mismatch")
    zd=.196*float(cat.Rdisk)**.633
    if not np.isfinite(zd) or zd<=0:raise ValueError("invalid scale height")
    hashes={"rotmod":sha(rotmod),"dens":sha(density_path),"table":sha(table_path),
            "normalization":sha(normalization_path)}
    import inspect
    from mono_SPARC_RAR import vcdisk
    solver_path=Path(inspect.getfile(vcdisk))
    force_model={"geometry":{"distance_Mpc":float(cat.D),"inclination_deg":float(cat.Inc),
        "zd_kpc":zd,"vertical_profile":"exp"},
        "force_source_sha256":{"mono_SPARC_RAR.py":sha(ROOT/"mono_SPARC_RAR.py"),
                               "utils.py":sha(ROOT/"utils.py"),"vcdisk":sha(solver_path)},
        "density_support":"full supplied profile; remove R=0 before Casertano",
        "kinematic_support":"overlap only; refuse outside positive force support",
        "interpolation":"signed W/R force interpolation then multiply observed R",
        "signed_W":"V*abs(V)","units":{"radius":"kpc","velocity":"km/s","acceleration":"m/s^2"},
        "kpc_metres":3.086e19}
    invariants={"force_model":force_model,"galaxy":galaxy,"parser_id":common.PARSER_ID,"input_sha256":hashes,
        "normalization_id":norm["policy"],"nuisance_model":NUISANCE,
        "retained_mask":parse["mask"],"kinematics_sha256":common.digest({
        k:data[k].tolist() for k in ("Rad","Vobs","errV","Vgas","Vbul")})}
    return data,rad,brightness,zd,invariants,rec


def build(args, solve=None):
    data,rad,brightness,zd,invariants,normalization=load_inputs(
        args.galaxy,args.normalization,args.data_root,args.dens_root,args.table)
    settings=common.Settings(args.prob_grid_size,args.zpadding,args.log_tolerance,args.max_refinements)
    invariants["probability_settings"]=asdict(settings)
    invariants["probability_sha256"]=sha(ROOT/"probability_recursion.py")
    context={**invariants,"invariants":invariants}
    plan={"galaxy":args.galaxy,"B":args.B,"seed":args.seed,"baseline_disc_ml":.5,
          "nominal_nuisance":NUISANCE,"probability_settings":asdict(settings),
          "retained_rows":len(data["Rad"]),"full_density_rows":len(rad),
          "R0_nodes_removed_by_force_solver":int(np.count_nonzero(rad==0)),
          "normalization":normalization,"invariants":invariants,
          "maximum_output_bytes":8*1024*1024,"force_solves":1,"mock_force_solves":0}
    if args.dry_run:return {"status":"prepared_no_force_no_draws",**plan}
    if args.output.exists():raise FileExistsError(args.output)
    if solve is None:
        from mono_SPARC_RAR import solve_disk_basis
        solve=solve_disk_basis
    sources={n:sha(ROOT/n) for n in ("common_reference_driver.py","common_isotonic_reference.py",
        "mono_SPARC_RAR.py","probability_recursion.py","utils.py")}
    from mono_SPARC_RAR import vcdisk
    import inspect
    solver_path=Path(inspect.getfile(vcdisk))
    sources[str(solver_path)]=sha(solver_path)
    unit_w,warnings,error=solve(rad,brightness,zd,data["Rad"])
    if error or warnings or unit_w is None:
        raise ValueError("constant Casertano basis unavailable: "+str((error,warnings)))
    scorer=common.prepare_sparc_baseline(data,unit_w,settings)
    result=common.run_reference(data["Vobs"],scorer,B=args.B,seed=args.seed,
                                context=context,save_vectors=args.save_vectors)
    result["build_plan"]=plan
    result["driver_source_sha256"]=sources
    result["unit_disk_W_kms2"]=np.asarray(unit_w).tolist()
    result["unit_disk_W_sha256"]=common.digest(np.asarray(unit_w).tolist())
    if sources!={n:sha(Path(n) if Path(n).is_absolute() else ROOT/n) for n in sources}:
        raise ValueError("source changed during reference construction")
    _,_,_,_,after,_=load_inputs(args.galaxy,args.normalization,args.data_root,args.dens_root,args.table)
    if after!={k:v for k,v in invariants.items() if k not in ("probability_settings","probability_sha256")}:
        raise ValueError("input changed during reference construction")
    common.save_reference(result,args.output)
    return {"status":result["status"],"output":str(args.output),"B":args.B,
            "failed_evaluations":len(result["failures"]),"reference_sha256":sha(args.output)}


def annotate_document(reference, document):
    """Only explicit new evidence schema. Never promote a legacy hash to parser proof.

    A converter must independently construct invariant evidence, not copy the
    reference context. Candidate profile parameters are intentionally NOT invariants.
    """
    provenance=document.get("provenance",{})
    if BAD_INVERSE_HASH in json.dumps(provenance):
        raise ValueError("confirmed wrong-gas-column production source; rescoring/selection review required")
    if document.get("schema")!="common-reference-candidates-v1":
        raise ValueError("legacy product lacks independently verified invariant evidence; no automatic conversion")
    evidence=document.get("invariants")
    if not evidence or evidence!=reference["context"].get("invariants"):
        raise ValueError("candidate invariant provenance mismatch")
    expected=reference["context"]["source_sha256"]["common_isotonic_reference.py"]
    if document.get("parser_source_sha256")!=expected or evidence.get("parser_id")!=common.PARSER_ID:
        raise ValueError("unverified candidate input parser")
    if reference.get("status")!="complete":raise ValueError("reference CDF unavailable")
    results=[]
    for record in document["records"]:
        if not record.get("valid") or record.get("warnings"):
            raise ValueError("invalid or warned candidate cannot be annotated")
        if record.get("probability_kind")!="marginalised_gas_bulge":
            raise ValueError("candidate must use original nuisance-marginalised probability")
        results.append({**record,**common.annotate(reference,float(record["log_probability"]),
                                                   reference["context"])})
    metric=document.get("deviation_metric")
    if any(r.get("J_ref") is not None for r in results) and not (isinstance(metric,dict) and metric.get("id")):
        raise ValueError("minimum-J selection requires separate declared deviation-metric provenance")
    candidates=[(i,r) for i,r in enumerate(results) if r["q_common_ref"]>=.05 and
                r.get("J_ref") is not None and np.isfinite(r["J_ref"]) and r["J_ref"]>=0]
    selected=min(candidates,key=lambda item:item[1]["J_ref"]) if candidates else None
    return {"schema":"common-reference-annotations-v1","records":results,
            "deviation_metric":metric,"q_threshold":.05,"minimum_found_J_record_index":None if selected is None else selected[0],
            "selection_scope":"supplied compatible finite candidate pool; unchanged J_ref",
            "reference_context_sha256":reference["context_sha256"],
            "candidate_document_sha256":common.digest(document),
            "interpretation":"operational common benchmark, not calibrated significance"}


def parser():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest="command",required=True)
    b=sub.add_parser("build")
    b.add_argument("--galaxy",required=True);b.add_argument("--normalization",type=Path,required=True)
    b.add_argument("--B",type=int,required=True);b.add_argument("--seed",type=int,required=True)
    b.add_argument("--output",type=Path,required=True)
    b.add_argument("--data-root",type=Path,default=Path("/mnt/users/koe/data"))
    b.add_argument("--dens-root",type=Path,default=Path(utils.DEFAULT_DENS_ROOT))
    b.add_argument("--table",type=Path,default=Path(utils.DEFAULT_TABLE_PATH))
    b.add_argument("--prob-grid-size",type=int,default=400);b.add_argument("--zpadding",type=float,default=8.)
    b.add_argument("--log-tolerance",type=float,default=2e-5);b.add_argument("--max-refinements",type=int,default=8)
    b.add_argument("--seconds",type=int,default=60);b.add_argument("--save-vectors",action="store_true")
    b.add_argument("--dry-run",action="store_true")
    a=sub.add_parser("annotate");a.add_argument("--reference",type=Path,required=True)
    a.add_argument("--candidates",type=Path,required=True);a.add_argument("--output",type=Path,required=True)
    return p


def main():
    args=parser().parse_args()
    if args.command=="build":
        if args.B<1 or args.seed<0 or args.seconds<1:raise ValueError("positive B/time and nonnegative seed required")
        def timeout(*unused):raise TimeoutError("declared reference budget exhausted")
        signal.signal(signal.SIGALRM,timeout);signal.alarm(args.seconds)
        try:result=build(args)
        finally:signal.alarm(0)
    else:
        reference=json.loads(args.reference.read_text());document=json.loads(args.candidates.read_text())
        result=annotate_document(reference,document)
        common.save_reference(result,args.output)
        result={"status":"annotated","output":str(args.output),"records":len(result["records"])}
    print(json.dumps(result,indent=2,allow_nan=False))

if __name__=="__main__":main()
