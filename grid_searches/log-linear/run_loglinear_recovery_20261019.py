#!/mnt/users/koe/RAR_mono_grid/.venv-osqp/bin/python
"""Isolated 59-galaxy log-linear recovery; positional shard CLI, --plan is read-only.
No queue submission or automatic retry. Frozen science sources are unchanged.
"""
import argparse, hashlib, json, os, re, sys, traceback
from datetime import datetime, timezone
from pathlib import Path
ROOT = Path(__file__).resolve().parent
PREP = ROOT / "outputs/grid_recovery_20261019/preparation_v1"
MANIFEST = PREP / "manifest.json"
OUTPUT_PARENT = ROOT / "outputs/grid_recovery_20261019"
SETTINGS = {"profile":"log_linear", "grid_size":51, "prob_grid_size":400,
    "ml_samples":9, "zpadding":8, "boundary_log_tolerance":2e-5,
    "boundary_max_refinements":8, "disk_ml":.5, "disk_ml_dex":.1,
    "disk_grid_half_width_dex":1, "slope_min":-5, "slope_max":5,
    "gas_ml":1, "gas_ml_dex":.04, "bulge_ml":.7, "bulge_ml_dex":.1,
    "tie_atol":0, "max_tie_permutations":720}

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def stamp():
    return datetime.now(timezone.utc).isoformat()

def write_new(path, value):
    with Path(path).open("x") as stream:
        json.dump(value,stream,indent=2,allow_nan=False); stream.write("\n")

def validate(m):
    if m.get("schema") != "loglinear-recovery-v1" or m.get("settings") != SETTINGS:
        raise ValueError("unexpected recovery schema/settings")
    out = Path(m["output_root"]).resolve()
    if out.parent != OUTPUT_PARENT.resolve() or out.name != "log_linear_v1":
        raise ValueError("new isolated log_linear_v1 output namespace required")
    groups = m["shards"]
    names = [g for group in groups.values() for g in group]
    if (len(groups)!=22 or sorted(map(len,groups.values()))!=[2]*7+[3]*15 or
            len(names)!=59 or len(set(names))!=59 or any(not re.fullmatch(r"[A-Za-z0-9_-]+",g) for g in names)):
        raise ValueError("22 distinct shards must allocate 59 unique galaxies")
    norm = json.loads(Path(m["normalization"]).read_text())
    bulges = {g for g,r in norm["galaxies"].items() if r["Lbulge_1e9"]>0}
    controls = set(m["failed_bulgeless_galaxies"])
    if len(bulges)!=32 or len(controls)!=27 or bulges & controls or set(names)!=bulges|controls:
        raise ValueError("exact 32 bulge plus 27 failed bulgeless allocation required")
    if any(norm["galaxies"][g]["Lbulge_1e9"] != 0 for g in controls):
        raise ValueError("failed controls must be bulgeless")
    for kind in ("source_sha256","input_sha256"):
        if not m.get(kind): raise ValueError("missing frozen dependency pins")
        for path,h in m[kind].items():
            if not Path(path).is_file() or sha(path)!=h:
                raise ValueError("frozen dependency changed: " + path)
    return out

def load_context():
    pin = json.loads((PREP/"launch_pin.json").read_text())
    if pin["manifest_sha256"] != sha(MANIFEST) or pin["worker_sha256"] != sha(__file__):
        raise ValueError("manifest/worker differs from reviewed launch pin")
    m = json.loads(MANIFEST.read_text()); validate(m)
    return m

def build_plan(m, shard):
    if shard not in m["shards"]: raise ValueError("unknown fixed shard; no refill")
    argv = [str(ROOT/"mono_SPARC_RAR_boundary.py")]
    for name,value in SETTINGS.items(): argv += ["--"+name.replace("_","-"),str(value)]
    argv += ["--marginalise-gas-bulge","--marginalise-disk", "--output-root",m["output_root"],
        "--shard-id",str(int(shard)),"--data-root",m["data_root"],"--dens-root",m["dens_root"],
        "--table-path",m["table"]]
    for galaxy in m["shards"][shard]: argv += ["--galaxy",galaxy]
    return {"shard":shard,"galaxies":m["shards"][shard],"argv":argv,
        "output_root":m["output_root"],"manifest_sha256":sha(MANIFEST)}

def worker(m, plan, event):
    # Import scientific modules only for an approved queued execution.
    import run_corrected_recovery_pilot as recovery
    import mono_SPARC_RAR as mono
    import mono_SPARC_RAR_boundary as boundary
    proof = recovery.normalization_proof(m,MANIFEST)
    recovery.install(m,proof)  # Exactly one Ldisk/Ltotal scaling in this new process.
    process, hashes, previous_argv = mono.process_galaxy, boundary.hashes, sys.argv
    def monitored(galaxy,*args,**kwargs):
        load_context()
        event({"state":"galaxy_start","galaxy":galaxy,"utc":stamp()})
        result = process(galaxy,*args,**kwargs)
        event({"state":"galaxy_finished","galaxy":galaxy,"utc":stamp()})
        return result
    mono.process_galaxy = monitored
    boundary.hashes = lambda:dict(hashes(),new_grid_launcher=sha(__file__),
        normalization_input=sha(m["normalization"]),recovery_manifest=sha(MANIFEST))
    try:
        sys.argv = plan["argv"]
        boundary.main()
        metadata = Path(plan["output_root"])/("run_metadata.shard_"+plan["shard"]+".json")
        data = json.loads(metadata.read_text())
        data["normalization"] = proof
        data["recovery_scope"] = {"manifest_sha256":sha(MANIFEST),"galaxies":plan["galaxies"],
            "normalization_applied_once":True,"original_sources_modified":False}
        metadata.write_text(json.dumps(data,indent=2,allow_nan=False)+"\n")
        load_context()
    finally:
        mono.process_galaxy, boundary.hashes, sys.argv = process, hashes, previous_argv

def run(m,plan,job_id,work=None):
    if not job_id: raise ValueError("scheduler job ID required before science/output creation")
    root = Path(plan["output_root"])
    for prefix,suffix in (("profile_grid","csv"),("fully_marginalised","csv"),("run_metadata","json")):
        path = root/(prefix+".shard_"+plan["shard"]+"."+suffix)
        if path.exists(): raise FileExistsError("refusing existing shard output: "+str(path))
    root.mkdir(parents=True,exist_ok=True)
    claims = root/"claims"; claims.mkdir(exist_ok=True)
    write_new(claims/("shard_"+plan["shard"]+".json"),dict(plan,job_id=str(job_id),utc=stamp()))
    journal = root/("shard_"+plan["shard"]+".events.jsonl")
    with journal.open("x") as stream:
        def event(value):
            stream.write(json.dumps(value,allow_nan=False)+"\n"); stream.flush()
            print(json.dumps(value,allow_nan=False),flush=True)
        event({"state":"shard_start","job_id":str(job_id),"utc":stamp()})
        try:
            (worker if work is None else work)(m,plan,event)
            required = [root/(p+".shard_"+plan["shard"]+"."+s) for p,s in
                (("profile_grid","csv"),("fully_marginalised","csv"),("run_metadata","json"))]
            if not all(p.is_file() for p in required): raise RuntimeError("missing final shard products")
            event({"state":"shard_complete","utc":stamp(),"files":[str(p) for p in required]})
        except BaseException as error:
            event({"state":"shard_failed","utc":stamp(),"error":repr(error),"traceback":traceback.format_exc()})
            raise

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("shard"); parser.add_argument("--plan",action="store_true")
    args = parser.parse_args()
    if not args.plan and not os.environ.get("SLURM_JOB_ID"):
        parser.error("queued execution requires SLURM_JOB_ID; use --plan for read-only inspection")
    m = load_context(); plan = build_plan(m,args.shard)
    if args.plan: print(json.dumps(plan,indent=2)); return
    for key in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"):
        os.environ[key] = "1"
    run(m,plan,os.environ["SLURM_JOB_ID"])
if __name__ == "__main__": main()
