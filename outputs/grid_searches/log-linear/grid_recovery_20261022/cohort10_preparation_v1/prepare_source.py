import ast,copy,csv,datetime,glob,hashlib,importlib.util,json,os,re,subprocess
from pathlib import Path
ROOT=Path.cwd().resolve()
assert str(ROOT)=='/mnt/users/koe/RAR_mono_grid'
old=ROOT/'run_loglinear_checkpoint37_20261019.py'
new=ROOT/'run_loglinear_cohort10_checkpoint_20261022.py'
prep=ROOT/'outputs/grid_recovery_20261022/cohort10_preparation_v1'
out=ROOT/'outputs/grid_recovery_20261022/log_linear_cohort10_checkpoint_v1'
assert not new.exists() and not prep.exists() and not out.exists(), 'new namespace must be unused'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
basepath=ROOT/'outputs/grid_recovery_20261019/checkpoint37_preparation_v1/manifest.json'
base=json.loads(basepath.read_text())
for field in ['source_sha256','input_sha256']:
 for p,h in base[field].items(): assert sha(p)==h, 'frozen dependency changed '+p
names='NGC4183 NGC4214 NGC4389 NGC4559 NGC5055 NGC5371 NGC5585 NGC5907 NGC6015 NGC6503'.split()
shards={f'{i:02d}':[g] for i,g in enumerate(names,1)}
done={r['galaxy'] for f in glob.glob('outputs/rerun_boundary_20261006/log_linear/fully_marginalised.shard_*.csv')+glob.glob('outputs/grid_recovery_20261019/log_linear*/fully_marginalised.shard_*.csv') for r in csv.DictReader(open(f)) if r['valid']=='True'}
ids=subprocess.check_output(['squeue','-u','koe','-h','-o','%i']).decode().split(); jobs=[]; active=set()
for i in ids:
 b=subprocess.check_output(['scontrol','write','batch_script',i,'-']).decode(); match=re.search(r'run_loglinear_(recovery|checkpoint37)_20261019.py\s+(\d+)',b)
 if not match: raise ValueError('unmapped active job '+i)
 kind,s=match.groups(); p='checkpoint37_preparation_v1' if kind=='checkpoint37' else 'preparation_v1'; d=json.load(open('outputs/grid_recovery_20261019/'+p+'/manifest.json')); ns=d['shards'][s.zfill(2)]; active.update(ns);jobs.append(dict(job_id=i,kind=kind,shard=s,galaxies=ns,batch_script_sha256=hashlib.sha256(b.encode()).hexdigest()))
assert not (set(names)&done or set(names)&active)
audit=dict(schema='cohort10-coverage-v1',utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),cancelled_job_id=935526,galaxies=names,completed_galaxies=sorted(done),active_jobs=jobs,queue_count=len(ids),status='all10_uncovered_no_submission_authorized')
source=old.read_text(); tree=ast.parse(source); node=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='validate')
validation='''def validate(m):
    # Fixed replacement scope; no refill or expansion through manifest changes.
    expected={f"{i:02d}":[g] for i,g in enumerate("NGC4183 NGC4214 NGC4389 NGC4559 NGC5055 NGC5371 NGC5585 NGC5907 NGC6015 NGC6503".split(),1)}
    if m.get("schema")!="loglinear-cohort10-checkpoint-v1" or m.get("settings")!=SETTINGS:
        raise ValueError("unexpected cohort10 scope/settings")
    if m.get("shards")!=expected:raise ValueError("exact10 singleton shards required; no refill")
    if Path(m["output_root"]).resolve()!=OUTPUT.resolve():raise ValueError("isolated cohort10 namespace required")
    for kind in ("source_sha256","input_sha256"):
        if not m.get(kind):raise ValueError("missing pins")
        for path,digest in m[kind].items():
            if not Path(path).is_file() or sha(path)!=digest:raise ValueError("changed frozen dependency: "+path)
    for key in ("data_root","dens_root","table","normalization"):
        if m.get(key)!=FIXED_INPUTS[key]:raise ValueError("changed input root: "+key)
    required=[m["table"],m["normalization"]]
    for group in expected.values():
        galaxy=group[0]
        required.extend([str(Path(m["data_root"])/(galaxy+"_rotmod.dat")),str(Path(m["dens_root"])/(galaxy+".dens"))])
    if not set(required).issubset(m["input_sha256"]):raise ValueError("required cohort input pins absent")
    audit_path=PREP/"coverage_audit.json"
    if m.get("coverage_audit")!=str(audit_path) or sha(audit_path)!=m.get("coverage_audit_sha256"):
        raise ValueError("coverage receipt changed")
    audit=json.loads(audit_path.read_text())
    names={g for group in expected.values() for g in group}
    active={g for job in audit["active_jobs"] for g in job["galaxies"]}
    if set(audit["galaxies"])!=names or names&set(audit["completed_galaxies"]) or names&active:
        raise ValueError("coverage receipt does not support exact uncovered10")
'''
lines=source.splitlines(True); source=''.join(lines[:node.lineno-1])+validation+''.join(lines[node.end_lineno:])
source=source.replace('outputs/grid_recovery_20261019/checkpoint37_preparation_v1','outputs/grid_recovery_20261022/cohort10_preparation_v1').replace('outputs/grid_recovery_20261019/log_linear_checkpoint37_v1','outputs/grid_recovery_20261022/log_linear_cohort10_checkpoint_v1').replace('JOBS={935528:7,935529:8,935531:10}\nKEEP=935526','FIXED_INPUTS='+repr({k:base[k] for k in ['data_root','dens_root','table','normalization']})).replace('reviewed37 launch pin changed','reviewed cohort10 launch pin changed')
assert 'JOBS=' not in source
compile(source,str(new),'exec')
oldworker=next(n for n in ast.parse(old.read_text()).body if isinstance(n,ast.FunctionDef) and n.name=='worker')
newworker=next(n for n in ast.parse(source).body if isinstance(n,ast.FunctionDef) and n.name=='worker')
assert ast.dump(oldworker,include_attributes=False)==ast.dump(newworker,include_attributes=False)
prep.mkdir(parents=True)
def write(p,d):
 with Path(p).open('x') as f:json.dump(d,f,indent=2,allow_nan=False);f.write(chr(10))
write(prep/'coverage_audit.json',audit)
with new.open('x') as f:f.write(source)
new.chmod(old.stat().st_mode & 0o777)
m={k:copy.deepcopy(base[k]) for k in ['settings','normalization','data_root','dens_root','table','source_sha256','input_sha256','resources']}
m.update(schema='loglinear-cohort10-checkpoint-v1',status='prepared_not_approved_for_submission',created_utc=audit['utc'],output_root=str(out),shards=shards,parent_manifest=str(basepath),parent_manifest_sha256=sha(basepath),coverage_audit=str(prep/'coverage_audit.json'),coverage_audit_sha256=sha(prep/'coverage_audit.json'),scope_note='Exact10 uncovered cancelled935526 galaxies; original37 worker body unchanged; no new scientific settings; no submission authorized')
m['source_sha256'][str(new)]=sha(new);m['source_sha256'][str(old)]=sha(old)
for g in names:
 for p in [Path(m['data_root'])/(g+'_rotmod.dat'),Path(m['dens_root'])/(g+'.dens')]:m['input_sha256'][str(p)]=sha(p)
write(prep/'manifest.json',m)
write(prep/'launch_pin.json',dict(manifest_sha256=sha(prep/'manifest.json'),worker_sha256=sha(new),engine_sha256=sha(ROOT/'run_loglinear_recovery_20261019.py'),checkpoint_sha256=sha(ROOT/'grid_numeric_checkpoint_20261019.py')))
spec=importlib.util.spec_from_file_location('cohort10',new);c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c);loaded=c.load_context(); plans={s:c.build_plan(loaded,s) for s in shards};write(prep/'resolved_plans.json',plans)
write(prep/'preparation_receipt.json',dict(status='prepared_validation_pending_no_submission',original_worker_body_ast_identical=True,old_worker_sha256=sha(old),new_worker_sha256=sha(new),manifest_sha256=sha(prep/'manifest.json'),plans=10,scientific_runs=0,queue_submissions=0,repository_git_status='not a Git repository; existing files preserved and hash-checked'))
print('PREPARED',new,'PLANS',len(plans),'PINS',len(m['source_sha256']),len(m['input_sha256']),'QUEUE',len(ids),'SCIENCE_RUNS',0,'SUBMISSIONS',0)
