#!/usr/bin/env python3
"""Execute only Enoch-approved3cancellations and15checkpoint37submissions.
Attempt receipts prevent ambiguous duplicate submission. No retry/refill/cleanup.
"""
import argparse,json,os,re,shlex,subprocess,sys
from pathlib import Path
ROOT=Path("/mnt/users/koe/RAR_mono_grid")
sys.path.insert(0,str(ROOT))
import run_loglinear_checkpoint37_20261019 as r
PREP=r.PREP
RECEIPTS=PREP/"execution_receipts_v1"
CANCEL=["935528","935529","935531"]
APPROVED=[f"{i:02d}" for i in range(1,16)]
PIN="86b4254555095f63e5a94fe012c3105a382a0e2573a1bbe825eaea595688c4ea"
WORKER="2413da69d5e7e406e0662aca03db8a6c1c6e4afe11ecf341f3184947ffec361e"
def command(sid):
    return ["addqueue","-q","cmb","-n","1","-m","8","-c","2 days",str(Path(r.__file__)),sid]
def call(argv,cwd=ROOT):
    p=subprocess.run(argv,cwd=cwd,capture_output=True,text=True,timeout=45)
    return dict(argv=argv,cwd=str(cwd),utc=r.engine.stamp(),returncode=p.returncode,stdout=p.stdout,stderr=p.stderr)
def write(name,value):r.write_new(RECEIPTS/name,value)
def queue():
    result=call(["squeue","-u","koe","-h","-o","%i|%T|%M"])
    if result["returncode"]:raise RuntimeError(result)
    jobs={line.split("|")[0]:line.split("|")[1:] for line in result["stdout"].splitlines() if line.strip()}
    return jobs,result

def guard():
    m=r.load_context()
    if r.sha(r.MANIFEST)!=PIN or r.sha(r.__file__)!=WORKER:raise ValueError("reviewed37pins differ")
    proposed=[line for line in (PREP/"commands_proposed.sh").read_text().splitlines() if line.startswith("addqueue ")]
    if [shlex.split(line) for line in proposed[:15]]!=[command(sid) for sid in APPROVED]:
        raise ValueError("approved exact command list differs")
    if not os.access(r.__file__,os.X_OK):raise ValueError("payload not executable")
    old=ROOT/"outputs/rerun_boundary_20261006/log_linear"
    for sid in ("07","08","10"):
        for p,s in (("profile_grid","csv"),("fully_marginalised","csv"),("run_metadata","json")):
            path=old/(p+".shard_"+sid+"."+s)
            if path.exists():raise RuntimeError("old numerical products arrived; HOLD: "+str(path))
    return m

def approved_context():
    record=json.loads((RECEIPTS/"preflight.json").read_text())
    if record["operator_sha256"]!=r.sha(__file__):raise RuntimeError("operator changed after preflight")
    return record

def preflight():
    m=guard();jobs,q=queue()
    if any(j not in jobs or jobs[j][0]!="RUNNING" for j in CANCEL):raise RuntimeError("old job states changed; HOLD")
    if len(jobs)-3+15>40:raise RuntimeError("approved15would exceed40cap; HOLD")
    if r.OUTPUT.exists():raise RuntimeError("new output namespace already exists; reconcile, do not repeat")
    write("preflight.json",dict(status="explicit_human_approval_verified_ready",utc=r.engine.stamp(),
        authorization_thread="01a0c386-c9b2-70e2-89be-34027b06c4c8",
        authorization_turn="01a1560c-afe2-77a3-a34f-108cb468a551",
        authorization_user_message="01a1560c-b097-7471-855b-2ac6fc4cbf32",
        operator_sha256=r.sha(__file__),manifest_sha256=PIN,worker_sha256=WORKER,
        source_pins=len(m["source_sha256"]),input_pins=len(m["input_sha256"]),
        queue=q,active_before=len(jobs),projected_after=len(jobs)-3+15,
        cancellation_argv=["scancel",*CANCEL],approved_submissions={sid:command(sid) for sid in APPROVED},
        kept_job="935526",unapproved_hold=list(range(16,38))))
    print("PREFLIGHT_OK",len(jobs),"active; projected",len(jobs)-3+15,"cap40; allpinsmatch",flush=True)

def cancel():
    approved_context();guard();jobs,q=queue()
    if any(j not in jobs or jobs[j][0]!="RUNNING" for j in CANCEL):raise RuntimeError("old state changed before cancellation; HOLD")
    if len(jobs)-3+15>40:raise RuntimeError("capacity changed; HOLD")
    write("cancellation_attempt.json",dict(utc=r.engine.stamp(),argv=["scancel",*CANCEL],queue=q))
    result=call(["scancel",*CANCEL]);write("cancellation_receipt.json",result)
    if result["returncode"]:raise RuntimeError(result)
    print("CANCELLATION_REQUEST_ACCEPTED",','.join(CANCEL),flush=True)

def submit(sid):
    if sid not in APPROVED:raise ValueError("only01–15explicitly approved;16–37HOLD")
    approved_context();guard()
    cancel_receipt=json.loads((RECEIPTS/"cancellation_receipt.json").read_text())
    if cancel_receipt["returncode"]:raise RuntimeError("cancellation not accepted")
    jobs,q=queue()
    if set(CANCEL)&set(jobs):raise RuntimeError("old cancellation not complete; wait/recheck, never expand")
    if len(jobs)>=40:raise RuntimeError("account cap reached; HOLD")
    if sid!="01" and not (RECEIPTS/"first_startup_gate.json").exists():
        raise RuntimeError("first actual startup gate not verified")
    if (r.OUTPUT/"claims"/("shard_"+sid+".json")).exists():raise RuntimeError("existing worker claim")
    # This exclusive attempt is written BEFORE scheduler call. No automatic retry.
    write("attempt_"+sid+".json",dict(utc=r.engine.stamp(),argv=command(sid),queue=q))
    try:result=call(command(sid),cwd=PREP/"submission_logs")
    except BaseException as error:
        write("ambiguous_"+sid+".json",dict(utc=r.engine.stamp(),error=repr(error),automatic_retry=False));raise
    matched=re.findall(r"Submitted batch job\s+(\d+)",result["stdout"]+"\n"+result["stderr"])
    result["shard"]=sid;result["job_id"]=matched[-1] if len(matched)==1 else None
    write("receipt_"+sid+".json",result)
    if result["returncode"] or result["job_id"] is None:raise RuntimeError("unconfirmed submission; reconcile, do not retry: "+str(result))
    print("SUBMITTED",sid,result["job_id"],flush=True)

def first_gate():
    m=guard();approved_context()
    receipt=json.loads((RECEIPTS/"receipt_01.json").read_text());job=receipt["job_id"]
    jobs,q=queue()
    if job not in jobs or jobs[job][0]!="RUNNING":raise RuntimeError("firstjob not RUNNING")
    claim=json.loads((r.OUTPUT/"claims/shard_01.json").read_text())
    if claim["job_id"]!=job or claim["manifest_sha256"]!=PIN or claim["galaxies"]!=m["shards"]["01"]:
        raise RuntimeError("firstclaim identity differs")
    context=json.loads((r.OUTPUT/"checkpoints"/m["shards"]["01"][0]/"context.json").read_text())
    if context["identity"]["manifest_sha256"]!=PIN:raise RuntimeError("firstcheckpointmanifest differs")
    if context["identity"]["source_sha256"]!=m["source_sha256"]:raise RuntimeError("firstcheckpointsource differs")
    status=call(["sacct","-n","-P","-j",job,"-o","JobIDRaw,State,ExitCode"])
    if not any(line.startswith(job+".0|RUNNING|") for line in status["stdout"].splitlines()):
        raise RuntimeError("actualfirst.0step notRUNNING: "+status["stdout"])
    write("first_startup_gate.json",dict(utc=r.engine.stamp(),job_id=job,shard="01",claim=claim,
        context_manifest=PIN,actual_step=status,queue=q,status="startup_identity_verified_not_scientific_validation"))
    print("FIRST_STARTUP_GATE_OK",job,"sourcepinnedcontext present; actual.0RUNNING",flush=True)

def audit(label):
    approved_context();m=guard();jobs,q=queue();records=[]
    for sid in APPROVED:
        path=RECEIPTS/("receipt_"+sid+".json")
        if not path.is_file():continue
        receipt=json.loads(path.read_text());job=receipt["job_id"]
        claim_path=r.OUTPUT/"claims"/("shard_"+sid+".json")
        claim=json.loads(claim_path.read_text()) if claim_path.is_file() else None
        galaxy=m["shards"][sid][0];context_path=r.OUTPUT/"checkpoints"/galaxy/"context.json"
        context=json.loads(context_path.read_text()) if context_path.is_file() else None
        if claim and (claim["job_id"]!=job or claim["manifest_sha256"]!=PIN):raise RuntimeError("claim differs")
        if context and (context["identity"]["manifest_sha256"]!=PIN or context["identity"]["source_sha256"]!=m["source_sha256"]):
            raise RuntimeError("checkpointcontext differs")
        journal=r.OUTPUT/("shard_"+sid+".events.jsonl")
        events=[]
        if journal.is_file():
            for line in journal.read_text().splitlines():
                try:events.append(json.loads(line))
                except json.JSONDecodeError:pass  # A concurrent last write is not a failure.
        partial=r.OUTPUT/"checkpoints"/galaxy/"conditional.partial.csv"
        records.append(dict(shard=sid,galaxy=galaxy,job_id=job,queue_state=jobs.get(job),
            claim_present=claim is not None,source_pinned_context=context is not None,
            events=[e.get("state") for e in events],failed=any(e.get("state")=="shard_failed" for e in events),
            partial_csv_bytes=partial.stat().st_size if partial.exists() else 0,
            checkpoint_complete=(r.OUTPUT/"checkpoints"/galaxy/"galaxy_complete.json").is_file()))
    cancel_status=call(["sacct","-n","-P","-j",','.join(CANCEL),"-o","JobIDRaw,State,ExitCode"])
    actual=call(["sacct","-n","-P","-j",','.join(x["job_id"] for x in records),"-o","JobIDRaw,State,ExitCode"]) if records else None
    write(label+".json",dict(utc=r.engine.stamp(),queue=q,active=len(jobs),old_absent=not(set(CANCEL)&set(jobs)),
        kept_state=jobs.get("935526"),cancellation_accounting=cancel_status,new_actual_steps=actual,records=records,
        submitted=len(records),unapproved16to37_submitted=False,source_pins_verified=True))
    print("AUDIT",label,"active",len(jobs),"submitted",len(records),"oldabsent",not(set(CANCEL)&set(jobs)),
        "kept935526",jobs.get("935526"),flush=True)
    for row in records:print(json.dumps(row),flush=True)
    print("CANCEL_ACCOUNTING",cancel_status["stdout"].strip(),flush=True)
    if actual:print("ACTUAL_NEW_STEPS",actual["stdout"].strip(),flush=True)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode",choices=("preflight","cancel","submit","gate","audit"))
    parser.add_argument("values",nargs="*");args=parser.parse_args()
    RECEIPTS.mkdir(exist_ok=True)
    if args.mode=="preflight":preflight()
    elif args.mode=="cancel":cancel()
    elif args.mode=="gate":first_gate()
    elif args.mode=="submit":
        for sid in args.values:submit(sid)
    elif args.mode=="audit":
        if len(args.values)!=1 or not re.fullmatch(r"[A-Za-z0-9_-]+",args.values[0]):raise ValueError("unique auditlabel needed")
        audit(args.values[0])
if __name__=="__main__":main()
