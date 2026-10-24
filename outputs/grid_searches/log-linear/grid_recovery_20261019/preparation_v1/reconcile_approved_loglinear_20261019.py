"""Read-only startup reconciliation of the exact22 accepted jobs; writes one receipt."""
import contextlib, io, json
from pathlib import Path
import submit_approved_loglinear_20261019 as control

def main():
    with contextlib.redirect_stdout(io.StringIO()):snapshot=control.inspect()
    m=control.worker.load_context();root=Path(m["output_root"])
    rs=control.receipts()
    if len(rs)!=22 or len({r["job_id"] for r in rs})!=22 or not all(r["accepted"] and r.get("generated_wrapper_correct_shard") for r in rs):
        raise ValueError("expected22distinct accepted submissions and matching wrappers")
    queue={s.split("|",2)[0]:s.split("|",2)[1] for s in snapshot["queue"]["rows"]}
    accounting={}
    for line in snapshot["accounting"]["stdout"].splitlines():
        f=line.split("|")
        if len(f)>=4:accounting[f[0]]=dict(state=f[1],elapsed=f[2],exitcode=f[3])
    # Enumerate directory entries before point lookups on the shared filesystem.
    claims={p.name:p for p in (root/"claims").glob("*.json")}
    jobs=[]
    for receipt in rs:
        sid=receipt["shard"];job=receipt["job_id"]
        claimpath=claims.get("shard_"+sid+".json")
        claim=json.loads(claimpath.read_text()) if claimpath is not None else None
        expected=control.worker.build_plan(m,sid)
        if claim is not None and (claim["job_id"]!=job or claim["manifest_sha256"]!=control.EXPECTED["manifest.json"] or claim["argv"]!=expected["argv"]):
            raise ValueError("startup claim mismatch: "+sid)
        journal=root/("shard_"+sid+".events.jsonl")
        events=[json.loads(s) for s in journal.read_text().splitlines() if s.strip()] if journal.is_file() else []
        logfile=control.LOGS/("run_loglinear_recovery_20261019.py-"+job+".out")
        log=logfile.read_text() if logfile.is_file() else ""
        failures=[e for e in events if e["state"]=="shard_failed"]
        error=bool(failures or "Traceback (most recent call last)" in log)
        jobs.append(dict(shard=sid,job_id=job,scheduler_state=queue.get(job),
            actual_step=accounting.get(job+".0"),claim_matches=claim is not None,
            normalized_grid_entry=any(e["state"]=="galaxy_start" for e in events),
            galaxies_started=[e["galaxy"] for e in events if e["state"]=="galaxy_start"],
            last_event=events[-1] if events else None,startup_error=error,
            log_path=str(logfile),log_size=logfile.stat().st_size if logfile.is_file() else None))
    result=dict(approval=control.APPROVAL,utc=control.stamp(),accepted=22,
        active_pending_total=len(queue),account_queue=snapshot["queue"],jobs=jobs,
        matched_claims=sum(j["claim_matches"] for j in jobs),normalized_grid_entries=sum(j["normalized_grid_entry"] for j in jobs),
        startup_errors=sum(j["startup_error"] for j in jobs),manifest_sha256=control.EXPECTED["manifest.json"],
        worker_sha256=control.WORKER_SHA,output_root=str(root),
        old_loglinear_states={str(j):queue.get(str(j)) for j in (935526,935528,935529,935531)},
        no_cancellations=True,GP_reference_submissions=0,accounting=snapshot["accounting"])
    control.write(control.RECEIPTS/"startup_reconciliation.json",result)
    print(json.dumps({k:v for k,v in result.items() if k not in ("account_queue","accounting")},indent=2))
if __name__=="__main__":main()
