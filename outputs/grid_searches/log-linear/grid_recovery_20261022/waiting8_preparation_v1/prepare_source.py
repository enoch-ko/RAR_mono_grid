import ast,copy,csv,datetime,glob,hashlib,importlib.util,json,re,subprocess
from pathlib import Path
ROOT=Path.cwd().resolve();assert str(ROOT)=='/mnt/users/koe/RAR_mono_grid'
parent=ROOT/'outputs/grid_recovery_20261022/cohort10_preparation_v1';PREP=ROOT/'outputs/grid_recovery_20261022/waiting8_preparation_v1';new=ROOT/'run_loglinear_waiting8_checkpoint_20261022.py';out=ROOT/'outputs/grid_recovery_20261022/log_linear_waiting8_checkpoint_v1'
assert not PREP.exists() and not new.exists() and not out.exists()
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
base=json.load(open(parent/'manifest.json'))
for kind in ['source_sha256','input_sha256']:
 for p,h in base[kind].items():assert sha(p)==h,'frozen dependency changed '+p
q=subprocess.check_output(['squeue','-u','koe','-h','-o','%i %T']).decode().splitlines();ids=[l.split()[0] for l in q]
done={r['galaxy'] for f in glob.glob('outputs/rerun_boundary_20261006/log_linear/fully_marginalised.shard_*.csv')+glob.glob('outputs/grid_recovery_20261019/log_linear*/fully_marginalised.shard_*.csv')+glob.glob('outputs/grid_recovery_20261022/log_linear*/fully_marginalised.shard_*.csv') for r in csv.DictReader(open(f)) if r.get('valid')=='True'}
jobs=[];other_assigned=set()
for job in ids:
 b=subprocess.check_output(['scontrol','write','batch_script',job,'-']).decode();match=re.search(r'run_loglinear_(recovery|checkpoint37|cohort10_checkpoint)_2026(?:1019|1022).py\s+(\d+)',b)
 if not match:raise ValueError('unmapped active job '+job)
 kind,s=match.groups();s=s.zfill(2)
 if kind=='recovery':pd=ROOT/'outputs/grid_recovery_20261019/preparation_v1'
 elif kind=='checkpoint37':pd=ROOT/'outputs/grid_recovery_20261019/checkpoint37_preparation_v1'
 else:pd=parent
 m=json.load(open(pd/'manifest.json'));ns=m['shards'][s]
 if len(ns)<=1:other_assigned.update(ns);continue
 journal=Path(m['output_root'])/('shard_'+s+'.events.jsonl');events=[json.loads(l) for l in journal.read_text().splitlines()];started=[e['galaxy'] for e in events if e.get('state')=='galaxy_start'];finished=[e['galaxy'] for e in events if e.get('state')=='galaxy_finished'];current=started[-1] if started else None
 jobs.append(dict(job_id=job,shard=s,galaxies=ns,started_galaxies=started,finished_galaxies=finished,current_galaxy=current,journal=str(journal),batch_script_sha256=hashlib.sha256(b.encode()).hexdigest()))
started_all={g for j in jobs for g in j['started_galaxies']}
waiting={g for j in jobs for g in j['galaxies'] if g not in done and g not in started_all and g not in other_assigned}
names='PGC51017 UGC00891 NGC4217 NGC1090 UGC00128 NGC2955 UGC01230 UGC02885'.split()
assert waiting==set(names),('roster changed; review',waiting)
shards={f'{i:02d}':[g] for i,g in enumerate(names,1)}
audit=dict(schema='waiting8-coverage-v1',utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),galaxies=names,unstarted_galaxies=sorted(waiting),completed_galaxies=sorted(done),active_jobs=jobs,other_active_singleton_galaxies=sorted(other_assigned),queue_snapshot=q,queue_count=len(q),prospective_queue_count=len(q)+8,duplicate_computation_accepted=True,prior40cap_relaxed_for_this_proposal=True,submission_authorized=False)
source=(ROOT/'run_loglinear_cohort10_checkpoint_20261022.py').read_text();source=source.replace('cohort10','waiting8').replace('exact10','exact8').replace('uncovered10','unstarted8').replace('NGC4183 NGC4214 NGC4389 NGC4559 NGC5055 NGC5371 NGC5585 NGC5907 NGC6015 NGC6503',' '.join(names))
assert 'names&active:' in source
source=source.replace('names&active:','not names.issubset(active):')
needle='        raise ValueError("coverage receipt does not support exact unstarted8")'
assert needle in source
source=source.replace(needle,needle+chr(10)+'    if set(audit.get("unstarted_galaxies",[]))!=names or any(names&set(job["started_galaxies"]) for job in audit["active_jobs"]):'+chr(10)+'        raise ValueError("waiting selection includes started or out-of-scope galaxies")')
compile(source,str(new),'exec')
getworker=lambda t:ast.dump(next(n for n in ast.parse(t).body if isinstance(n,ast.FunctionDef) and n.name=='worker'),include_attributes=False)
assert getworker(source)==getworker((ROOT/'run_loglinear_checkpoint37_20261019.py').read_text())
PREP.mkdir(parents=True);(PREP/'submission_logs').mkdir()
def write(p,v):
 with p.open('x') as f:json.dump(v,f,indent=2);f.write(chr(10))
write(PREP/'coverage_audit.json',audit)
with new.open('x') as f:f.write(source)
new.chmod((ROOT/'run_loglinear_cohort10_checkpoint_20261022.py').stat().st_mode & 0o777)
m={k:copy.deepcopy(base[k]) for k in ['settings','normalization','data_root','dens_root','table','source_sha256','input_sha256','resources']};m['resources']['comment']='5 days'
m.update(schema='loglinear-waiting8-checkpoint-v1',status='prepared_not_approved_for_submission',created_utc=audit['utc'],output_root=str(out),shards=shards,parent_manifest=str(parent/'manifest.json'),parent_manifest_sha256=sha(parent/'manifest.json'),coverage_audit=str(PREP/'coverage_audit.json'),coverage_audit_sha256=sha(PREP/'coverage_audit.json'),scope_note='All8 unstarted nonfinal galaxies across20active multi-galaxy tasks. Duplicated computation explicitly accepted; no original cancellation, no auto-refill;40cap relaxed for proposed44total, submission still requires approval.')
m['source_sha256'][str(new)]=sha(new)
for g in names:
 for p in [Path(m['data_root'])/(g+'_rotmod.dat'),Path(m['dens_root'])/(g+'.dens')]:m['input_sha256'][str(p)]=sha(p)
write(PREP/'manifest.json',m);write(PREP/'launch_pin.json',dict(manifest_sha256=sha(PREP/'manifest.json'),worker_sha256=sha(new),engine_sha256=sha(ROOT/'run_loglinear_recovery_20261019.py'),checkpoint_sha256=sha(ROOT/'grid_numeric_checkpoint_20261019.py')))
spec=importlib.util.spec_from_file_location('waiting8prep',new);c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c);loaded=c.load_context();plans={s:c.build_plan(loaded,s) for s in shards};write(PREP/'resolved_plans.json',plans)
print('PREPARED8 queue',len(q),'prospective',len(q)+8,'multi_jobs_checked',len(jobs),'finals',len(done),'source_pins',len(m['source_sha256']),'input_pins',len(m['input_sha256']))
for s,ns in shards.items():
 g=ns[0];j=next(j for j in jobs if g in j['galaxies']);print('PLAN',s,g,'parent',j['job_id'],'current',j['current_galaxy'])
