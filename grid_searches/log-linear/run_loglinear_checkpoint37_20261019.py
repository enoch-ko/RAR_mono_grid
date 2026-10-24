#!/mnt/users/koe/RAR_mono_grid/.venv-osqp/bin/python
"""37 useful bulgeless replacements with numeric checkpoints. --plan is read-only.
Keep935526; no submission/cancellation/automatic retry/resume implementation.
"""
import importlib.util, json, re, sys
from pathlib import Path
import grid_numeric_checkpoint_20261019 as checkpoint
ROOT=Path(__file__).resolve().parent
PREP=ROOT/"outputs/grid_recovery_20261019/checkpoint37_preparation_v1"
MANIFEST=PREP/"manifest.json"
OUTPUT=ROOT/"outputs/grid_recovery_20261019/log_linear_checkpoint37_v1"
ENGINE=ROOT/"run_loglinear_recovery_20261019.py"
spec=importlib.util.spec_from_file_location("_private_checkpoint37_engine",ENGINE)
engine=importlib.util.module_from_spec(spec);spec.loader.exec_module(engine)
SETTINGS=dict(engine.SETTINGS)
sha,write_new=engine.sha,engine.write_new
JOBS={935528:7,935529:8,935531:10}
KEEP=935526

def validate(m):
    if m.get("schema")!="loglinear-checkpoint37-v1" or m.get("settings")!=SETTINGS:
        raise ValueError("unexpected37 scope/settings")
    if Path(m["output_root"]).resolve()!=OUTPUT.resolve():raise ValueError("isolated37 namespace required")
    for kind in ("source_sha256","input_sha256"):
        if not m.get(kind):raise ValueError("missing pins")
        for path,digest in m[kind].items():
            if not Path(path).is_file() or sha(path)!=digest:raise ValueError("changed frozen dependency: "+path)
    audit=Path(m["allocation_audit"])
    if sha(audit)!=m["allocation_audit_sha256"]:raise ValueError("allocation audit changed")
    jobs={int(j["job_id"]):j for j in json.loads(audit.read_text())["jobs"]}
    if {j:int(jobs[j]["shard"]) for j in JOBS}!=JOBS:raise ValueError("old-job scope mismatch")
    expected={g for job in JOBS for g in jobs[job]["unique_bulgeless"]}
    keep=set(jobs[KEEP]["unique_bulgeless"])
    names=[g for group in m["shards"].values() for g in group]
    if (len(expected)!=37 or len(keep)!=10 or expected&keep or
            set(m["shards"])!={f"{i:02d}" for i in range(1,38)} or
            any(len(group)!=1 for group in m["shards"].values()) or
            len(names)!=37 or len(set(names))!=37 or set(names)!=expected or set(names)&keep):
        raise ValueError("exact37 onegalaxy allocation excluding ALL10 kept galaxies required")
    if m.get("keep_job")!=KEEP or set(m.get("replacing_jobs",[]))!=set(JOBS):
        raise ValueError("keep/cancel declaration mismatch")
    if m.get("verified_completed_galaxies"):raise ValueError("new completed coverage needs reviewed version")
    norm=json.loads(Path(m["normalization"]).read_text())["galaxies"]
    if any(norm[g]["Lbulge_1e9"]!=0 or norm[g]["kappa"]!=1 for g in names):
        raise ValueError("all37 must have zero bulge and kappa1")
    if any(not re.fullmatch(r"[A-Za-z0-9_-]+",g) for g in names):raise ValueError("bad galaxy name")
    for path in (ENGINE,Path(__file__).resolve(),Path(checkpoint.__file__).resolve()):
        if str(path) not in m["source_sha256"]:raise ValueError("launcher/checkpoint pin missing")
    if m["normalization"] not in m["input_sha256"]:raise ValueError("normalization pin missing")
    return OUTPUT

def load_context():
    pin=json.loads((PREP/"launch_pin.json").read_text())
    if any(pin[k]!=sha(p) for k,p in (("manifest_sha256",MANIFEST),("worker_sha256",__file__),
            ("engine_sha256",ENGINE),("checkpoint_sha256",checkpoint.__file__))):
        raise ValueError("reviewed37 launch pin changed")
    m=json.loads(MANIFEST.read_text());validate(m);return m

def worker(m,plan,event):
    # Same scientific implementation; only process-local I/O hooks are new.
    import run_corrected_recovery_pilot as recovery
    import mono_SPARC_RAR as mono
    import mono_SPARC_RAR_boundary as boundary
    proof=recovery.normalization_proof(m,MANIFEST)
    recovery.install(m,proof)
    original_process,original_append,original_full=mono.process_galaxy,mono.append_probability_rows,mono.fully_marginalise_grid
    original_hashes,previous_argv=boundary.hashes,sys.argv
    def monitored(galaxy,table,args,rows,x,means,slopes,gas_nodes,gas_weights,bulge_nodes,bulge_weights):
        load_context()
        identity=dict(galaxy=galaxy,profile="log_linear",manifest_path=str(MANIFEST),
            manifest_sha256=sha(MANIFEST),source_sha256=m["source_sha256"],input_sha256=m["input_sha256"],
            settings=m["settings"],normalization=proof,shard=plan["shard"])
        axes=dict(log10_mean_ml=list(x),mean_ml=list(means),slopes=list(slopes),
            gas_nodes=list(gas_nodes),gas_weights=list(gas_weights),bulge_nodes=list(bulge_nodes),bulge_weights=list(bulge_weights))
        cp=checkpoint.GridCheckpoint(Path(plan["output_root"])/"checkpoints"/galaxy,identity,axes)
        def append(*positional,**keywords):
            target=positional[0] if positional else keywords["rows"]
            before=len(target)
            result=original_append(*positional,**keywords)
            cp.append_cell(target[before:])
            return result
        def full(*positional,**keywords):
            negative=positional[6] if len(positional)>6 else keywords["negative_grid"]
            cp.finish_grids(checkpoint.json_value(negative))
            event(dict(state="numeric_grids_checkpointed",galaxy=galaxy,cells=len(cp.seen),utc=engine.stamp()))
            result=original_full(*positional,**keywords)
            cp.save_fully(result)
            return result
        mono.append_probability_rows,mono.fully_marginalise_grid=append,full
        event(dict(state="galaxy_start",galaxy=galaxy,utc=engine.stamp()))
        try:
            result=original_process(galaxy,table,args,rows,x,means,slopes,gas_nodes,gas_weights,bulge_nodes,bulge_weights)
            load_context()
            cp.finish_galaxy(result[4])
            checkpoint.verify_complete(cp.path,identity)
            event(dict(state="galaxy_finished",galaxy=galaxy,numeric_checkpoint=str(cp.path),utc=engine.stamp()))
            return result
        except BaseException as error:
            try:cp.fail(error)
            except BaseException as checkpoint_error:
                event(dict(state="checkpoint_failure_receipt_unavailable",galaxy=galaxy,error=repr(checkpoint_error)))
            raise
        finally:
            mono.append_probability_rows,mono.fully_marginalise_grid=original_append,original_full
    mono.process_galaxy=monitored
    boundary.hashes=lambda:dict(original_hashes(),new_grid_launcher=sha(__file__),
        numeric_checkpoint=sha(checkpoint.__file__),normalization_input=sha(m["normalization"]),recovery_manifest=sha(MANIFEST))
    try:
        sys.argv=plan["argv"];boundary.main()
        meta=Path(plan["output_root"])/("run_metadata.shard_"+plan["shard"]+".json")
        data=json.loads(meta.read_text());data["normalization"]=proof
        data["recovery_scope"]=dict(manifest_sha256=sha(MANIFEST),galaxies=plan["galaxies"],
            normalization_applied_once=True,original_sources_modified=False,numeric_checkpoint_schema="grid-numeric-checkpoint-v1",
            resume_implemented=False)
        # A separate atomic authoritative sidecar avoids rewriting the original final metadata.
        checkpoint.publish_json(meta.with_name("checkpoint_provenance.shard_"+plan["shard"]+".json"),data)
        load_context()
    finally:
        mono.process_galaxy,mono.append_probability_rows,mono.fully_marginalise_grid=original_process,original_append,original_full
        boundary.hashes,sys.argv=original_hashes,previous_argv

engine.PREP,engine.MANIFEST=PREP,MANIFEST
engine.load_context,engine.validate,engine.worker=load_context,validate,worker
build_plan,run=engine.build_plan,engine.run

def main():engine.main()
if __name__=="__main__":main()
