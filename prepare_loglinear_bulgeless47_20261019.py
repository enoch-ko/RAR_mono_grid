#!/usr/bin/env python3
"""Metadata-only preparation. No galaxy solves, submission, cancellation or merging."""
import csv, json, re, shlex, subprocess
from datetime import datetime, timezone
from pathlib import Path
import run_loglinear_bulgeless47_20261019 as runner
ROOT, PREP = runner.ROOT, runner.PREP
BASE = ROOT / "outputs/grid_recovery_20261019/preparation_v1/manifest.json"
AUDIT = BASE.parent / "overlap_progress_audit.json"
OLD = ROOT / "outputs/rerun_boundary_20261006/log_linear"

def command(argv):
    p = subprocess.run(argv, cwd=ROOT, capture_output=True, text=True, timeout=30)
    if p.returncode:
        raise RuntimeError(str(argv) + ": " + p.stderr)
    return dict(argv=argv, stdout=p.stdout, stderr=p.stderr, returncode=p.returncode)

def read_rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))

def main():
    if PREP.exists() or runner.OUTPUT.exists():
        raise FileExistsError("new preparation and science namespaces required")
    base = json.loads(BASE.read_text())
    for kind in ("source_sha256", "input_sha256"):
        for p, digest in base[kind].items():
            if runner.sha(p) != digest:
                raise ValueError("running-batch dependency changed: " + p)
    audit = json.loads(AUDIT.read_text())
    norm = json.loads(Path(base["normalization"]).read_text())["galaxies"]
    queue = command(["squeue", "-u", "koe", "-h", "-o", "%i|%T|%M"])
    active = [line.split("|")[0] for line in queue["stdout"].splitlines() if line.strip()]
    usage = command(["sstat", "-a", "-j", ",".join(map(str,runner.JOBS)),
        "-n", "-P", "-o", "JobID,AveCPU,AveRSS,MaxRSS,NTasks"])
    targets, previous, observations = [], set(), []
    for job in audit["jobs"]:
        sid = f"{int(job['shard']):02d}"
        products = [OLD / (p + ".shard_" + sid + "." + s) for p,s in
            (("profile_grid","csv"),("fully_marginalised","csv"),("run_metadata","json"))]
        # Refuse even partial arrivals: never race a newly finished old shard.
        if any(p.exists() for p in products):
            raise RuntimeError("old numeric products appeared; review compatible completed coverage before preparing: " + sid)
        logpath = ROOT / ("mono_SPARC_RAR_boundary.py-" + str(job["job_id"]) + ".out")
        log = logpath.read_text()
        found = re.findall(r"Processing galaxy (\d+)/(\d+)", log)
        index, total = map(int, found[-1])
        if total != len(job["allocation"]):
            raise ValueError("old allocation/log mismatch")
        useful = job["unique_bulgeless"]
        previous.update(set(job["allocation"][:index-1]) & set(useful))
        targets.extend(useful)
        observations.append(dict(job_id=job["job_id"], old_shard=sid,
            allocation=job["allocation"], useful_bulgeless=useful, current_index=index,
            current_galaxy=job["allocation"][index-1], log_sha256=runner.sha(logpath),
            log_size=logpath.stat().st_size, final_numeric_products_present=False))
    if len(targets) != 47 or len(set(targets)) != 47:
        raise ValueError("unexpected target coverage")
    ongoing59 = {g for group in base["shards"].values() for g in group}
    if set(targets) & ongoing59:
        raise ValueError("replacement targets overlap corrected running 59")
    weights = {}
    for g in targets:
        dens, raw = [], []
        for path, rows, width in ((Path(base["dens_root"])/(g+".dens"),dens,3),
                (Path(base["data_root"])/(g+"_rotmod.dat"),raw,8)):
            for line in path.read_text().splitlines():
                try:
                    values = list(map(float,line.split()))
                except ValueError:
                    continue
                if len(values) == width:
                    rows.append(values)
        retained = sum(dens[0][0] <= r[0] <= dens[-1][0] for r in raw)
        weights[g] = sum(r[0] > 0 for r in dens) * retained**2
    # Stable longest-proxy-first ordering, not a measured runtime prediction.
    ordered = sorted(targets, key=lambda g:(-weights[g],g))
    m = dict(base)
    m.update(schema="loglinear-bulgeless47-v1", status="prepared_not_approved_for_execution",
        created_utc=datetime.now(timezone.utc).isoformat(), settings=runner.SETTINGS,
        output_root=str(runner.OUTPUT), parent_manifest=str(BASE), parent_manifest_sha256=runner.sha(BASE),
        allocation_audit=str(AUDIT), allocation_audit_sha256=runner.sha(AUDIT),
        shards={f"{i:02d}":[g] for i,g in enumerate(ordered,1)},
        verified_completed_galaxies=[], replacement_galaxies=sorted(targets),
        replacing_jobs=list(runner.JOBS), previously_computed_ram_only=sorted(previous),
        allocation_weights=weights, scope_note="All47 useful zero-bulge kappa1 galaxies; scheduling replacement, not gas-bug correction.",
        weight_rule="positive dens row count * retained radius count squared, heuristic only",
        resources=dict(cpus_per_job=1,memory_gib_per_job=8,queue="cmb",comment="2 days",hard_time_limit=None))
    m.pop("failed_bulgeless_galaxies",None)
    m.pop("bulge_galaxies",None)
    m["source_sha256"] = dict(base["source_sha256"])
    m["input_sha256"] = dict(base["input_sha256"])
    for p in (Path(runner.__file__),Path(__file__),ROOT/"test_loglinear_bulgeless47_20261019.py"):
        m["source_sha256"][str(p.resolve())] = runner.sha(p)
    for p in (BASE,AUDIT):
        m["input_sha256"][str(p.resolve())] = runner.sha(p)
    runner.validate(m)
    # Historical same-batch finished metadata, not current hashes alone.
    historical = []
    per_galaxy_bytes = []
    for sid in ("01","02","04","09"):
        p = OLD/("run_metadata.shard_"+sid+".json")
        data = json.loads(p.read_text())
        pins = data["arguments"]["numerical_source_sha256"]
        pins = json.loads(pins) if isinstance(pins,str) else pins
        if not all(runner.sha(ROOT/name)==digest for name,digest in pins.items()):
            raise ValueError("historical numerical source mismatch")
        files = [OLD/(prefix+".shard_"+sid+".csv") for prefix in ("profile_grid","fully_marginalised")]
        n = int(data["galaxy_count"])
        size = sum(f.stat().st_size for f in files)
        per_galaxy_bytes.append(size/n)
        historical.append(dict(shard=sid,metadata_sha256=runner.sha(p),started_utc=data["started_utc"],
            finished_utc=data["finished_utc"],galaxies=n,csv_bytes=size,numerical_pins_match=True))
    PREP.mkdir(parents=True)
    (PREP/"submission_logs").mkdir()
    runner.write_new(runner.MANIFEST,m)
    runner.write_new(PREP/"launch_pin.json",dict(manifest_sha256=runner.sha(runner.MANIFEST),
        worker_sha256=runner.sha(runner.__file__),engine_sha256=runner.sha(runner.ENGINE)))
    checked = runner.load_context()
    plans = {sid:runner.build_plan(checked,sid) for sid in checked["shards"]}
    runner.write_new(PREP/"resolved_plans.json",plans)
    runner.write_new(PREP/"preparation_audit.json",dict(utc=m["created_utc"],queue=queue,usage=usage,
        jobs=observations,active_count=len(active),historical=historical,
        ram_only_repeated=len(previous),science_evaluations=0,submissions=0,cancellations=0,
        projected_csv_bytes_47=47*max(per_galaxy_bytes),
        storage_caveat="CSV-only observed maximum completed-shard bytes per galaxy; PNGs/events/provenance extra; not a hard cap"))
    commands = ['addqueue -q cmb -n 1 -m 8 -c "2 days" ' + shlex.quote(str(Path(runner.__file__))) + " " + sid for sid in plans]
    first = max(0,30-(len(active)-sum(str(j) in active for j in runner.JOBS)))
    first = min(first,len(commands))
    if len(active) > 30 or first != 6:
        raise RuntimeError("capacity changed; review initial tranche, never auto-submit")
    with (PREP/"commands_proposed.sh").open("x") as stream:
        stream.write("# PROPOSAL ONLY. DO NOT EXECUTE THIS FILE. All47 need staged approvals.\n")
        stream.write("# FIRST obtain explicit cancellation and initial-six approval.\n")
        stream.write("scancel 935526 935528 935529 935531\n")
        stream.write("cd " + shlex.quote(str(PREP/"submission_logs")) + "\n")
        stream.write("\n".join(commands[:first])+"\n")
        stream.write("# HOLD: later41 commands below are NOT initial-wave approval; no automatic refill.\n")
        stream.write("\n".join(commands[first:])+"\n")
    review = ["# 47-galaxy log-linear scheduling replacement: PREPARED ONLY", "",
        "No science/cancellation/submission/merge performed. Old dependencies and active jobs untouched.",
        "47 zero-bulge kappa=1 controls; source column correction does not require their rerun.",
        f"Repeat {len(previous)} useful calculations held only in RAM; all47 are covered, not just remaining30.",
        "One galaxy per shard gives durable CSV/fully-marginalised/metadata after each galaxy, plus exclusive claim and events.",
        "Settings identical to approved51/400/9, full nuisance/disc marginalisation, signed gas, same boundary/tie/physical-prior policies.",
        "Reuse private instance of frozen worker implementation; no existing source file modified.",
        "",
        "28 active now; cancelling these4 would leave24, allowing6 replacements under30 cap.",
        "Recheck queue and old completed products immediately before any approved action; hold if changed.",
        "All22 corrected-grid jobs973851,973864..973884 and GP969420/969421 remain untouched.",
        "No automated submission/refill. Later41 require explicit exact-command approval in capacity-bounded tranches.",
        "1CPU/8GiB each; -q cmb; -c 2 days is a comment, not a hard limit. Initial6 request6CPU/48GiB total.",
        "Old completed mixed17/18-galaxy shards ran days; elapsed13+days on active jobs is not a pergalaxy ETA.",
        "Smaller shards do not accelerate one difficult galaxy, and proxy ordering is not measured runtime.",
        "CSV-only storage estimate is in preparation_audit.json; no full per-candidate JSON output.",
        "Merging is deferred: remind user to run merge_run_shards.py after all47 have verified complete compatible products.",
        "Do not combine with other roots blindly: manifest/source signatures intentionally differ.",
        "", "## Exact allocations", "", "shard | galaxy | oldjob | repeated RAM-only", "--- | --- | --- | ---"]
    jobs_by_galaxy = {g:j["job_id"] for j in audit["jobs"] for g in j["unique_bulgeless"]}
    review += [sid+" | "+g+" | "+str(jobs_by_galaxy[g])+" | "+str(g in previous) for sid,(g,) in m["shards"].items()]
    review += ["", "## Proposed initial cancellation and six-job tranche (approval required)", "", "~~~sh", "scancel 935526 935528 935529 935531", "cd "+shlex.quote(str(PREP/"submission_logs")), *commands[:6], "~~~", "", "## Later41 (HOLD; exact staged approval required)", "", "~~~sh", *commands[6:], "~~~"]
    with (PREP/"review.md").open("x") as stream:
        stream.write("\n".join(review)+"\n")
    runner.write_new(PREP/"preparation_receipt.json",dict(status="prepared_pending_tests_and_explicit_execution_approval",
        galaxies=47,shards=47,repeated_ram_only=len(previous),source_pins=len(m["source_sha256"]),
        input_pins=len(m["input_sha256"]),manifest_sha256=runner.sha(runner.MANIFEST),
        worker_sha256=runner.sha(runner.__file__),science_evaluations=0,submissions=0,cancellations=0,
        science_output_exists=runner.OUTPUT.exists(),initial_tranche=6,remaining_hold=41))
    print("PREPARED",len(plans),"one-galaxy plans; RAM-only repeat",len(previous),"initial6 HOLD41; no science or jobs executed")

if __name__ == "__main__":
    main()
