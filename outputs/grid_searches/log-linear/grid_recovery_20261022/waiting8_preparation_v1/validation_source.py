import ast,copy,datetime,hashlib,importlib.util,json,os,subprocess,sys
from pathlib import Path
ROOT=Path('/mnt/users/koe/RAR_mono_grid'); PREP=ROOT/'outputs/grid_recovery_20261022/waiting8_preparation_v1'
new=ROOT/'run_loglinear_waiting8_checkpoint_20261022.py';old=ROOT/'run_loglinear_checkpoint37_20261019.py'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
spec=importlib.util.spec_from_file_location('validated_waiting8',new);c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c);m=c.load_context(); checks=[]
def ok(name,test):
 assert test,name
 checks.append(name)
def refusal(name,fn,types=(ValueError,FileExistsError)):
 try:fn()
 except types:checks.append(name)
 else:raise AssertionError('expected refusal '+name)
def changed(name,change):
 d=copy.deepcopy(m);change(d);refusal(name,lambda:c.validate(d))
expected='PGC51017 UGC00891 NGC4217 NGC1090 UGC00128 NGC2955 UGC01230 UGC02885'.split()
ok('exact8_singletons',[m['shards'][f'{i:02d}'] for i in range(1,9)]==[[g] for g in expected])
getworker=lambda t:ast.dump(next(n for n in ast.parse(t).body if isinstance(n,ast.FunctionDef) and n.name=='worker'),include_attributes=False)
ok('scientific_worker_AST_unchanged',getworker(new.read_text())==getworker(old.read_text()))
base=json.load(open(m['parent_manifest']));ok('settings_identical',m['settings']==base['settings']);ok('parent_manifest_unchanged',sha(m['parent_manifest'])==m['parent_manifest_sha256'])
changed('wrong_schema',lambda d:d.update(schema='wrong'))
changed('wrong_settings',lambda d:d['settings'].update(grid_size=52))
changed('wrong_namespace',lambda d:d.update(output_root=base['output_root']))
changed('extra_shard',lambda d:d['shards'].update({'11':['NGC6503']}))
changed('wrong_galaxy',lambda d:d['shards'].update({'01':['NGC4138']}))
changed('missing_source_pins',lambda d:d.update(source_sha256={}))
changed('missing_input_pins',lambda d:d.update(input_sha256={}))
changed('changed_source_digest',lambda d:d['source_sha256'].update({str(old):'0'*64}))
changed('changed_input_digest',lambda d:d['input_sha256'].update({m['table']:'0'*64}))
changed('missing_galaxy_input_pin',lambda d:d['input_sha256'].pop(str(Path(m['data_root'])/'PGC51017_rotmod.dat')))
changed('changed_normalization',lambda d:d.update(normalization=m['table']))
changed('changed_coverage_receipt',lambda d:d.update(coverage_audit_sha256='0'*64))
refusal('unknown_shard',lambda:c.build_plan(m,'11'))
original_sha=c.sha
for key,path in [('manifest',c.MANIFEST),('worker',new),('engine',c.ENGINE),('checkpoint',c.checkpoint.__file__)]:
 c.sha=lambda p,target=str(path):'0'*64 if str(p)==target else original_sha(p)
 try:refusal('launch_pin_'+key,lambda:c.load_context())
 finally:c.sha=original_sha

original_loads=c.json.loads
coverage=(PREP/'coverage_audit.json').read_text();audit=original_loads(coverage)
for label,mutate in [('started_galaxy_in_receipt',lambda a:a['active_jobs'][0]['started_galaxies'].append('PGC51017')),('final_galaxy_in_receipt',lambda a:a['completed_galaxies'].append('PGC51017')),('wrong_unstarted_set',lambda a:a.update(unstarted_galaxies=[]))]:
 bad=copy.deepcopy(audit);mutate(bad);c.json.loads=lambda text,altered=bad:altered if text==coverage else original_loads(text)
 try:refusal(label,lambda:c.validate(m))
 finally:c.json.loads=original_loads

plans={}
env=dict(os.environ);env.pop('SLURM_JOB_ID',None);env['PYTHONDONTWRITEBYTECODE']='1'
for i,g in enumerate(expected,1):
 s=f'{i:02d}';p=subprocess.run([sys.executable,'-B',str(new),s,'--plan'],env=env,capture_output=True,text=True,timeout=30)
 assert p.returncode==0,p.stderr
 plan=json.loads(p.stdout);ok('cli_plan_'+s,plan['galaxies']==[g] and plan['output_root']==m['output_root'] and plan['argv'][plan['argv'].index('--galaxy')+1]==g);plans[s]=plan
p=subprocess.run([sys.executable,'-B',str(new),'01'],env=env,capture_output=True,text=True,timeout=30);ok('login_science_execution_refused',p.returncode==2 and 'SLURM_JOB_ID' in p.stderr)
fixture=PREP/'validation_fixture';fixture.mkdir()
plan=copy.deepcopy(plans['01']);plan['output_root']=str(fixture/'synthetic_output')
refusal('missing_scheduler_id',lambda:c.run(m,plan,''))
ok('no_output_before_scheduler',not Path(plan['output_root']).exists())
def fake_work(manifest,plan,event):
 root=Path(plan['output_root'])
 for stem,ext in [('profile_grid','csv'),('fully_marginalised','csv'),('run_metadata','json')]:
  (root/(stem+'.shard_'+plan['shard']+'.'+ext)).write_text('synthetic orchestration fixture only'+chr(10))
c.run(m,plan,'SYNTHETIC_VALIDATION_ONLY',work=fake_work);ok('synthetic_orchestration_success',Path(plan['output_root'],'shard_01.events.jsonl').is_file())
refusal('existing_final_output_refused',lambda:c.run(m,plan,'SYNTHETIC_DUPLICATE',work=fake_work))
# Direct exclusive claim writer test preserves the first claim file.
claim=fixture/'exclusive_claim.json';c.write_new(claim,{'fixture':True});refusal('duplicate_claim_refused',lambda:c.write_new(claim,{'fixture':False}));ok('claim_preserved',json.load(open(claim))=={'fixture':True})
ok('production_output_namespace_absent',not Path(m['output_root']).exists())
ids=subprocess.check_output(['squeue','-u','koe','-h','-o','%i %T']).decode().splitlines()
receipt=dict(status='passed_preparation_only_not_submitted',utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),checks=checks,check_count=len(checks),plans_validated=8,source_pins=len(m['source_sha256']),input_pins=len(m['input_sha256']),scientific_force_runs=0,probability_evaluations=0,queue_submissions=0,queue_snapshot=ids,queue_after8=len(ids)+8,worker_sha256=sha(new),manifest_sha256=sha(c.MANIFEST),original_frozen_dependencies_unchanged=True,science_output_namespace_created=False,limitation='Synthetic orchestration and read-only plans only; no real-galaxy numerical or checkpoint computation')
with (PREP/'validation_receipt.json').open('x') as f:json.dump(receipt,f,indent=2);f.write(chr(10))
commands=['addqueue -q cmb -n 1 -m 8 -c "5 days" '+str(new)+' '+f'{i:02d}' for i in range(1,9)]
with (PREP/'submission_commands.txt').open('x') as f:f.write(chr(10).join(commands)+chr(10))
print('VALIDATION_PASSED',len(checks),'PLANS',8,'PINS',len(m['source_sha256']),len(m['input_sha256']),'QUEUE',len(ids),'AFTER8',len(ids)+8,'SCIENCE_RUNS',0,'SUBMISSIONS',0)
print('WORKER_SHA256',sha(new));print('MANIFEST_SHA256',sha(c.MANIFEST));print('RECEIPT',PREP/'validation_receipt.json')
