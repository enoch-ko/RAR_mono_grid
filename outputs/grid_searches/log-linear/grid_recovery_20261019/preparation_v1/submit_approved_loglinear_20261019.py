#!/usr/bin/env python3
"""Execute only Enoch's exact approved22commands, once; first-job startup gate.
No retries, cancellations, extra jobs, or changes to pinned scientific sources.
"""
import argparse, hashlib, json, os, re, subprocess, sys
from datetime import datetime, timezone
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT))
import run_loglinear_recovery_20261019 as worker
PREP=worker.PREP
LOGS=PREP/"submission_logs"
RECEIPTS=PREP/"submission_receipts"
APPROVAL="Alice direct human approval turn01a154c6-7a4b-7002-90a6-ac0c209a06cc: Yep, let's approve these"
EXPECTED={"manifest.json":"f23528f761118f6dd44ccb6063b29216cfff07317f060f587541b19d8d321651",
    "commands.sh":"1d53d236dbabff66bb0be915619ae16598a9657f0b7371f179cd1e3228d448ce"}
WORKER_SHA="7b6cb9f043ee983eee9aa5cb2f250b7e6e8033b47496d3a917a6c2a48dece048"
def stamp():return datetime.now(timezone.utc).isoformat()
def write(path,value):worker.write_new(path,value)
def call(argv,timeout=30):
    p=subprocess.run(argv,cwd=LOGS,capture_output=True,text=True,timeout=timeout)
    return dict(argv=argv,returncode=p.returncode,stdout=p.stdout,stderr=p.stderr,utc=stamp())
def queue():
    q=call(["squeue","-h","-u","koe","-o","%i|%T|%j"])
    if q["returncode"]:raise RuntimeError("cannot verify account occupancy: "+q["stderr"])
    q["rows"]=[s for s in q["stdout"].splitlines() if s.strip()]
    return q
def receipts():
    return [json.loads(p.read_text()) for p in sorted(RECEIPTS.glob("receipt_*.json"))] if RECEIPTS.exists() else []
def preflight(stage):
    for name,h in EXPECTED.items():
        if worker.sha(PREP/name)!=h:raise ValueError("approved artifact hash changed: "+name)
    if worker.sha(worker.__file__)!=WORKER_SHA:raise ValueError("approved worker changed")
    m=worker.load_context()
    approved=(PREP/"commands.sh").read_text().splitlines()[2:]
    expected=['addqueue -q cmb -n 1 -m 8 -c "2 days" '+str(ROOT/"run_loglinear_recovery_20261019.py")+" "+f"{i:02d}" for i in range(1,23)]
    if approved!=expected:raise ValueError("exact approved22command list differs")
    existing=receipts();q=queue()
    if len(q["rows"])+(22-len(existing))>30:
        raise RuntimeError("account occupancy would exceed30softcap: "+str(len(q["rows"])))
    root=Path(m["output_root"])
    if stage=="first" and (existing or root.exists()):raise FileExistsError("first submission requires fresh receipts and science root")
    for i in (range(1,23) if stage=="first" else range(2,23)):
        sid=f"{i:02d}"
        if (RECEIPTS/("attempt_"+sid+".json")).exists():raise FileExistsError("already attempted; no retry: "+sid)
        for path in [root/"claims"/("shard_"+sid+".json"),root/("profile_grid.shard_"+sid+".csv"),
                     root/("fully_marginalised.shard_"+sid+".csv"),root/("run_metadata.shard_"+sid+".json")]:
            if path.exists():raise FileExistsError("approved target not fresh: "+str(path))
    return m,q

def submit(sid):
    worker.load_context()
    if len(queue()["rows"])>=30:raise RuntimeError("soft account cap reached before "+sid)
    argv=["addqueue","-q","cmb","-n","1","-m","8","-c","2 days",str(ROOT/"run_loglinear_recovery_20261019.py"),sid]
    write(RECEIPTS/("attempt_"+sid+".json"),dict(approval=APPROVAL,argv=argv,utc=stamp()))
    result=call(argv,timeout=60)
    ids=set(re.findall(r"run_loglinear_recovery_20261019\.py-(\d+)\.out",result["stdout"]))
    result.update(shard=sid,approval=APPROVAL,job_id=next(iter(ids)) if len(ids)==1 else None,
        accepted=result["returncode"]==0 and len(ids)==1)
    generated=LOGS/"run_loglinear_recovery_20261019.py.sh"
    if generated.is_file():
        text=generated.read_text()
        (RECEIPTS/("generated_wrapper_"+sid+".txt")).write_text(text)
        result["generated_wrapper_sha256"]=hashlib.sha256(text.encode()).hexdigest()
        result["generated_wrapper_correct_shard"]=bool(re.search(re.escape(str(ROOT/"run_loglinear_recovery_20261019.py"))+r"\s+"+sid+r"(?:\s|$)",text))
    write(RECEIPTS/("receipt_"+sid+".json"),result)
    print(json.dumps(dict(shard=sid,job_id=result["job_id"],accepted=result["accepted"],wrapper_correct=result.get("generated_wrapper_correct_shard"))),flush=True)
    if not result["accepted"] or not result.get("generated_wrapper_correct_shard"):
        raise RuntimeError("submission/wrapper failed; stop without retry: "+sid)
    return result

def startup_gate():
    first=json.loads((RECEIPTS/"receipt_01.json").read_text())
    m=worker.load_context();root=Path(m["output_root"])
    claim=root/"claims/shard_01.json"
    journal=root/"shard_01.events.jsonl"
    if not claim.is_file() or not journal.is_file():raise RuntimeError("first job has not produced a startup claim/journal")
    c=json.loads(claim.read_text())
    events=[json.loads(s) for s in journal.read_text().splitlines() if s.strip()]
    if any(e["state"]=="shard_failed" for e in events):raise RuntimeError("first worker failed; stop without retry")
    if c["job_id"]!=first["job_id"] or c["manifest_sha256"]!=EXPECTED["manifest.json"] or c["argv"]!=worker.build_plan(m,"01")["argv"]:
        raise ValueError("actual first-worker claim differs from approved manifest/argv")
    if not any(e["state"]=="galaxy_start" for e in events):raise RuntimeError("first job has not entered the normalized grid calculation")
    # This callback is emitted only after recovery.install() and boundary.main().
    return dict(job_id=first["job_id"],manifest_sha256=c["manifest_sha256"],events=events,
        normalization_sha256=worker.sha(m["normalization"]),output_root=m["output_root"],
        observed="actual normalized-adapter installation and grid-entry callback; numerical results pending",utc=stamp())

def inspect():
    q=queue();out={"queue":q,"receipts":[],"utc":stamp()}
    m=worker.load_context();root=Path(m["output_root"])
    for r in receipts():
        sid=r["shard"];p=root/("shard_"+sid+".events.jsonl")
        events=[json.loads(s) for s in p.read_text().splitlines() if s.strip()] if p.is_file() else []
        out["receipts"].append(dict(shard=sid,job_id=r["job_id"],accepted=r["accepted"],events=events))
    ids=[r["job_id"] for r in receipts() if r["job_id"]]
    out["accounting"]=call(["sacct","-j",",".join(ids),"-n","-P","-o","JobID,State,Elapsed,ExitCode"]) if ids else None
    print(json.dumps(out,indent=2))
    return out

def main():
    p=argparse.ArgumentParser();p.add_argument("stage",choices=("preflight","first","gate","remaining","inspect"));a=p.parse_args()
    if a.stage=="preflight":
        m,q=preflight("first");print(json.dumps(dict(pins="passed",occupancy=len(q["rows"]),new_jobs=22,predicted_total=len(q["rows"])+22,output_root=m["output_root"])))
    elif a.stage=="first":
        m,q=preflight("first");RECEIPTS.mkdir(exist_ok=False);write(RECEIPTS/"preflight.json",q);submit("01")
    elif a.stage=="gate":print(json.dumps(startup_gate(),indent=2))
    elif a.stage=="remaining":
        gate=startup_gate();m,q=preflight("remaining");write(RECEIPTS/"first_startup_gate.json",gate)
        for i in range(2,23):submit(f"{i:02d}")
        result=inspect();write(RECEIPTS/"initial_reconciliation.json",result)
    else:inspect()
if __name__=="__main__":main()
