#!/usr/bin/env python3
"""Prepare metadata and read-only job audit for the approved 59-galaxy grid scope.
Never submit/cancel jobs or execute galaxy calculations.
"""
import csv, hashlib, json, re, shlex, subprocess, sys
from datetime import datetime, timezone
from pathlib import Path
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
import run_loglinear_recovery_20261019 as runner
PREP=runner.PREP
PARENT=ROOT/"outputs/science_20261016/manifest.json"
NORM=ROOT/"outputs/bulge_normalized_20261015/inputs/normalization.json"
OLD=ROOT/"outputs/rerun_boundary_20261006/log_linear"
JOBS={935526:5,935528:7,935529:8,935531:10}
FIRST={"utc":"2026-10-19T15:10:39+00:00","origin":"read-only terminal sstat observation before preparation",
    "cpu":{"935526":"13-02:44:59","935528":"13-02:44:51","935529":"13-02:44:45","935531":"13-02:44:44"}}

def sha(path): return runner.sha(path)
def write(path,value): runner.write_new(path,value)
def command(argv,timeout=30):
    p=subprocess.run(argv,cwd=ROOT,capture_output=True,text=True,timeout=timeout)
    return dict(argv=argv,returncode=p.returncode,stdout=p.stdout,stderr=p.stderr)
def seconds(value):
    day,clock=(value.split("-",1) if "-" in value else ("0",value))
    fields=[int(x) for x in clock.split(":")]
    return int(day)*86400+sum(v*60**i for i,v in enumerate(reversed(fields)))
def roster(table,names):
    values=[]
    for line in Path(table).read_text().splitlines():
        f=line.split()
        if f and f[0] in names: values.append(f[0])
    if len(values)!=175 or len(set(values))!=175: raise ValueError("catalogue roster/order mismatch")
    return values

def audit(galaxies,allocations,targets,bulges):
    now=datetime.now(timezone.utc).isoformat()
    usage=command(["sstat","-a","-j",",".join(map(str,JOBS)),"-n","-P","-o","JobID,AveCPU,AveRSS,MaxRSS,NTasks"])
    rows={}
    for line in usage["stdout"].splitlines():
        fields=line.split("|")
        if fields[0].endswith(".0"): rows[fields[0].split(".")[0]]=fields
    plots={p.stem:p for p in OLD.rglob("*.png")}
    records=[]
    for job,shard in JOBS.items():
        names=allocations[shard]
        logfile=ROOT/("mono_SPARC_RAR_boundary.py-"+str(job)+".out")
        log=logfile.read_text() if logfile.is_file() else ""
        found=re.findall(r"Processing galaxy (\d+)/(\d+)",log)
        index,total=map(int,found[-1]) if found else (0,len(names))
        if total!=len(names) or index>len(names): raise ValueError("active log allocation mismatch")
        current=names[index-1] if index else None
        before=names[:max(index-1,0)]
        active=names[max(index-1,0):]
        evidence=[]
        for g in names:
            p=plots.get(g)
            if p is not None:
                stat=p.stat(); evidence.append(dict(galaxy=g,path=str(p),size=stat.st_size,mtime_unix=stat.st_mtime))
        cpu=rows.get(str(job))
        delta=seconds(cpu[1])-seconds(FIRST["cpu"][str(job)]) if cpu else None
        records.append(dict(job_id=job,shard=shard,allocation=names,
            overlap_with_new=sorted(set(names)&set(targets)),unique_bulgeless=sorted(set(names)-bulges),
            logged_current_index=index,current_galaxy=current,current_has_bulge=current in bulges,
            preceding_galaxies_in_memory=before,current_and_remaining=active,
            overlap_current_and_remaining=sorted(set(active)&set(targets)),
            unique_bulgeless_current_and_remaining=sorted(set(active)-bulges),
            log=dict(path=str(logfile),exists=logfile.is_file(),size=logfile.stat().st_size if logfile.is_file() else None,
                mtime_unix=logfile.stat().st_mtime if logfile.is_file() else None,tail=log[-3000:],
                error_flag=any(x in log.lower() for x in ("traceback","disk quota exceeded","no space left","killed"))),
            per_galaxy_plots=evidence,
            final_csv_exists=(OLD/("profile_grid.shard_"+f"{shard:02d}"+".csv")).is_file(),
            final_marginal_csv_exists=(OLD/("fully_marginalised.shard_"+f"{shard:02d}"+".csv")).is_file(),
            cpu_second_observation=cpu,cpu_increase_seconds=delta,
            disposition="retain_mixed_shard_unique_bulgeless_work; CPU activity alone is not proof of scientific progress"))
    processes=command(["ssh","-o","BatchMode=yes","-o","ConnectTimeout=10","comp48",
        "ps -u koe -o pid=,etimes=,cputimes=,pcpu=,rss=,stat=,wchan:20=,args="],timeout=25)
    processes["stdout"]="\n".join(line for line in processes["stdout"].splitlines()
        if "python " in line and "mono_SPARC_RAR_boundary.py" in line)
    return dict(first=FIRST,second_utc=now,sstat=usage,worker_processes=processes,
        queue=command(["squeue","-h","-u","koe","-o","%i|%T|%M|%o"]),jobs=records,
        new_unique_galaxies=len(targets),active_overlap_galaxies=sum(len(r["overlap_with_new"]) for r in records),
        active_unique_bulgeless=sum(len(r["unique_bulgeless"]) for r in records),
        cancellations=0,reason="No whole job redundant; all mix useful bulgeless work with superseded bulge work. Do not discard in-memory rows without approved replacement scope.")

def main():
    if PREP.exists() or (runner.OUTPUT_PARENT/"log_linear_v1").exists():
        raise FileExistsError("new versioned preparation/output namespace required")
    parent=json.loads(PARENT.read_text())
    for key in ("source_sha256","input_sha256"):
        bad=[p for p,h in parent[key].items() if not Path(p).is_file() or sha(p)!=h]
        if bad: raise ValueError("frozen parent changed: "+repr(bad))
    norm=json.loads(NORM.read_text())["galaxies"]
    names=roster(parent["table"],set(norm))
    review=ROOT/"submission_review_20261006.md"
    allocations={int(shard):names[int(start)-1:int(end)] for start,end,shard in
        re.findall(r"--profile log_linear.*?--galaxy-range (\d+) (\d+).*?--shard-id (\d+)",review.read_text())}
    if set(allocations)!=set(range(1,11)) or sum(map(len,allocations.values()))!=175:
        raise ValueError("historical log-linear allocation could not be verified")
    bulges={g for g,r in norm.items() if r["Lbulge_1e9"]>0}
    failed=set(allocations[3]+allocations[6])
    controls=failed-bulges
    targets=sorted(bulges|controls)
    if len(bulges)!=32 or len(failed)!=35 or len(controls)!=27 or len(targets)!=59:
        raise ValueError("approved target allocation mismatch")
    if controls & set(g for s in JOBS.values() for g in allocations[s]):
        raise ValueError("new bulgeless targets overlap active useful allocations")
    weights={}
    for g in targets:
        dens=[]
        for line in (Path(parent["dens_root"])/(g+".dens")).read_text().splitlines():
            try: values=list(map(float,line.split()))
            except ValueError: continue
            if len(values)==3: dens.append(values)
        raw=[]
        for line in (Path(parent["data_root"])/(g+"_rotmod.dat")).read_text().splitlines():
            if line.startswith("#"):continue
            try: values=list(map(float,line.split()))
            except ValueError:continue
            if len(values)==8:raw.append(values)
        count=sum(dens[0][0]<=r[0]<=dens[-1][0] for r in raw)
        weights[g]=sum(r[0]>0 for r in dens)*count**2
    groups=[dict(id=f"{i+1:02d}",capacity=3 if i<15 else 2,galaxies=[],weight=0) for i in range(22)]
    for g in sorted(targets,key=lambda g:(-weights[g],g)):
        group=min((r for r in groups if len(r["galaxies"])<r["capacity"]),key=lambda r:(r["weight"],r["id"]))
        group["galaxies"].append(g);group["weight"]+=weights[g]
    m=dict(schema="loglinear-recovery-v1",status="prepared_pending_explicit_queue_approval",
        created_utc=datetime.now(timezone.utc).isoformat(),settings=runner.SETTINGS,
        normalization=str(NORM),data_root=parent["data_root"],dens_root=parent["dens_root"],table=parent["table"],
        output_root=str(runner.OUTPUT_PARENT/"log_linear_v1"),
        source_sha256=dict(parent["source_sha256"]),input_sha256=dict(parent["input_sha256"]),
        parent_manifest=str(PARENT),parent_manifest_sha256=sha(PARENT),
        allocation_source=str(review),allocation_source_sha256=sha(review),
        failed_bulgeless_galaxies=sorted(controls),bulge_galaxies=sorted(bulges),
        shards={r["id"]:sorted(r["galaxies"]) for r in groups},allocation_weights=weights,
        weight_rule="positive .dens row count * retained kinematic count squared; heuristic only, no force/probability solve",
        resources={"cpus_per_job":1,"memory_gib_per_job":8,"queue":"cmb","comment":"2 days","hard_time_limit":None},
        scope_note="32 corrected bulges plus27 failed bulgeless, no duplicate8bulges. Existing active useful47 bulgeless galaxies retained.")
    for path in (Path(runner.__file__),Path(__file__)):
        m["source_sha256"][str(path.resolve())]=sha(path)
    runner.validate(m)
    observation=audit(names,allocations,targets,bulges)
    PREP.mkdir(parents=True)
    (PREP/"submission_logs").mkdir()
    write(runner.MANIFEST,m)
    write(PREP/"launch_pin.json",dict(manifest_sha256=sha(runner.MANIFEST),worker_sha256=sha(runner.__file__)))
    checked=runner.load_context()
    plans={sid:runner.build_plan(checked,sid) for sid in checked["shards"]}
    write(PREP/"resolved_plans.json",plans)
    write(PREP/"overlap_progress_audit.json",observation)
    commands=["# PREPARED ONLY.22 jobs require explicit approval; never run on the login node.",
        "cd "+shlex.quote(str(PREP/"submission_logs"))]
    commands += ['addqueue -q cmb -n 1 -m 8 -c "2 days" '+shlex.quote(str(ROOT/"run_loglinear_recovery_20261019.py"))+" "+sid for sid in plans]
    with (PREP/"commands.sh").open("x") as stream:stream.write("\n".join(commands)+"\n")
    frozen={p.name:sha(p) for p in PREP.iterdir() if p.is_file()}
    receipt=dict(status="prepared_all22_plans_validated_no_science",source_pins=len(m["source_sha256"]),
        input_pins=len(m["input_sha256"]),jobs=22,galaxies=59,bulges=32,failed_bulgeless=27,
        capacity_counts=sorted(map(len,m["shards"].values())),output_root_created=Path(m["output_root"]).exists(),
        science_evaluations=0,queue_submissions=0,cancellations=0,hashes=frozen)
    write(PREP/"preparation_receipt.json",receipt)
    print(json.dumps(dict(receipt,active_jobs=[dict(job_id=r["job_id"],current=r["current_galaxy"],overlap=len(r["overlap_with_new"]),
        useful_bulgeless=len(r["unique_bulgeless"]),cpu_delta=r["cpu_increase_seconds"]) for r in observation["jobs"]]),indent=2))
if __name__=="__main__":main()
