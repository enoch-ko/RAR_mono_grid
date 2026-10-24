from pathlib import Path
import json,csv,hashlib,subprocess,datetime,os,importlib.util,time,collections
R=Path('/mnt/users/koe/RAR_mono_grid');P=R/'bf_grid_extension_20261024_v1';assert not P.exists();P.mkdir();(P/'frozen_native').mkdir();sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest();old=R/'bf_grid_rollout_20261023_v1';parent=json.loads((old/'manifest.json').read_text());G=['NGC1090','NGC2955','NGC5907','NGC7793','UGC02487'];assert all(g in parent['deferred_log_linear'] for g in G)
source=(old/'driver.py').read_bytes();assert sha(old/'driver.py')==parent['common_pins'][str(old/'driver.py')];(P/'driver.py').write_bytes(source)
wrapper='''#!/bin/bash
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
cd /mnt/users/koe/RAR_mono_grid
exec /mnt/users/koe/RAR_mono_grid/.venv-osqp/bin/python -B /mnt/users/koe/RAR_mono_grid/bf_grid_extension_20261024_v1/driver.py "$@"
''';(P/'worker.sh').write_text(wrapper);(P/'worker.sh').chmod(0o755)
normal=Path(parent['normalization']);norm=json.loads(normal.read_text())['galaxies'];roster={g:{f:dict(status='deferred_not_in_extension') for f in ['constant','linear','log_linear']} for g in G};frozen={}
def freeze(q):
 q=Path(q);data=q.read_bytes();h=hashlib.sha256(data).hexdigest();p=P/'frozen_native'/(h+'_'+q.name)
 if not p.exists():
  with p.open('xb') as f:f.write(data)
 assert q.read_bytes()==data,'source changed during freeze '+str(q);v=dict(original_path=str(q),sha256=h,snapshot=str(p));frozen[str(q)]=v;return v
def add(spec,q):
 v=freeze(q);spec['pins'][v['snapshot']]=v['sha256'];return v
roots=['grid_recovery_20261019/log_linear_v1','grid_recovery_20261019/log_linear_checkpoint37_v1','grid_recovery_20261022/log_linear_cohort10_checkpoint_v1','grid_recovery_20261022/log_linear_waiting8_checkpoint_v1']
for rel in roots:
 root=R/'outputs'/rel
 for q in sorted(root.glob('fully_marginalised.shard_*.csv')):
  suffix=q.name[len('fully_marginalised'):];cp=root/('profile_grid'+suffix);mp=root/('checkpoint_provenance'+suffix.replace('.csv','.json'));mp=mp if mp.exists() else root/('run_metadata'+suffix.replace('.csv','.json'));meta=json.loads(mp.read_text())
  with q.open() as f:gal=[v['galaxy'] for v in csv.DictReader(f) if v['valid']=='True' and v['profile']=='log_linear' and v['galaxy'] in G]
  for g in gal:
   assert norm[g]['kappa']==1. or meta.get('normalization',{}).get('normalization_sha256')==sha(normal),'uncorrected normalization'
   fc,fm,ff=freeze(cp),freeze(mp),freeze(q);spec=dict(status='ready',csv=fc['snapshot'],metadata=fm['snapshot'],fully_marginalised=ff['snapshot'],native_csv_provenance=fc,native_metadata_provenance=fm,native_final_provenance=ff,pins={v['snapshot']:v['sha256'] for v in [fc,fm,ff]})
   complete=root/'checkpoints'/g/'galaxy_complete.json'
   if complete.exists():
    d=json.loads(complete.read_text());ident=d['identity'];assert d['status']=='galaxy_calculation_complete' and ident['galaxy']==g and ident['profile']=='log_linear';assert sha(ident['manifest_path'])==ident['manifest_sha256']
    for name,h in d['files'].items():assert sha(complete.parent/name)==h,'checkpoint product drift'
    for role in ['source_sha256','input_sha256']:
     for name,h in ident[role].items():assert sha(name)==h,'checkpoint source/input drift';spec['pins'][name]=h
    with cp.open() as f:top=[v for v in csv.DictReader(f) if v['galaxy']==g]
    with (complete.parent/'profile_grid.csv').open() as f:part=list(csv.DictReader(f))
    assert top==part,'top/checkpoint mismatch'
    for path in [complete,complete.parent/'context.json',Path(ident['manifest_path'])]:add(spec,path)
    spec['checkpoint_proof']=dict(complete_sha256=sha(complete),all_saved_file_hashes_verified=True,top_level_rows_equal_checkpoint=True,normalization=ident['normalization'])
   roster[g]['log_linear']=spec
assert all(roster[g]['log_linear']['status']=='ready' for g in G)
common=dict(parent['common_pins']);common[str(P/'driver.py')]=sha(P/'driver.py');common[str(P/'worker.sh')]=sha(P/'worker.sh');deferred=[g for g in parent['deferred_log_linear'] if g not in G];queue=subprocess.run(['squeue','-h','-u','koe','-o','%i|%T|%o'],capture_output=True,text=True,check=True).stdout.splitlines()
m=dict(schema='frozen-grid-bf-extension-v1',status='prepared_not_submitted',created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),parent_manifest=str(old/'manifest.json'),parent_manifest_sha256=sha(old/'manifest.json'),parent_driver_byte_identical=True,galaxies=G,roster=roster,baselines={g:parent['baselines'][g] for g in G},normalization=str(normal),common_pins=common,settings=parent['settings'],batches=dict(log_linear_01=dict(family='log_linear',galaxies=G,force_solve_cap=250,process_wall_cap_seconds=14400,scheduler_comment='BF_extension_log_linear_01')),deferred_log_linear=deferred,reviewed_gate_receipt=str(P/'reviewed_gate.json'),output=str(R/'outputs/bf_grid_extension_20261024_v1'),resources=dict(jobs=1,cpus_per_job=1,memory_GiB_per_job=8,process_wall_cap_hours=4,per_force_timeout_seconds=900,queue_cap=40),queue_snapshot=queue,force_solve_upper_bounds=dict(log_linear=250),denominator_request_upper_bounds=dict(log_linear=5*2601*81),limits=['diagnostic accuracy_unverified','saved numerators unchanged; no numerator evaluator/mocks/GP search','exact queue approval and recorded limited-pilot review required','same full-cell/physical-exclusion/provenance/ambiguity/failure gates as completed driver','no overwrite/refill/retry/cancel','old historical input/per-node archive limitations remain explicit'])
(P/'manifest.json').write_text(json.dumps(m,indent=2)+chr(10));(P/'frozen_native_index.json').write_text(json.dumps(frozen,indent=2)+chr(10));command='addqueue -q cmb -n 1 -m 8 -c "BF_extension_log_linear_01" '+str(P/'worker.sh')+' log_linear_01';(P/'command.txt').write_text(command+chr(10));cwd=P/'submissions/log_linear_01';cwd.mkdir(parents=True);(P/'submission_plan.json').write_text(json.dumps(dict(status='pending_exact_approval',manifest_sha256=sha(P/'manifest.json'),exact_command=command,working_directory=str(cwd),galaxies=G,resources=m['resources']),indent=2)+chr(10))
spec=importlib.util.spec_from_file_location('extension_worker',P/'driver.py');w=importlib.util.module_from_spec(spec);spec.loader.exec_module(w);w.checkpins(common);checks=[];start=time.perf_counter()
for g in G:
 x=w.load(m,g,'log_linear');checks.append(dict(galaxy=g,valid_cells=len(x['valid']),physical_exclusions=len(x['invalid']),retained_bins=len(x['r']),photometric_knots=len(x['rad']),kappa=x['b']['kappa_applied_once'],checkpoint_proof='checkpoint_proof' in roster[g]['log_linear']))
w.write(P/'preflight.json',dict(status='passed_input_only',manifest_sha256=sha(P/'manifest.json'),checks=checks,seconds=time.perf_counter()-start,driver_byte_identical=True,force_solves=0,numerator_evaluations=0,denominator_evaluations=0,queue_count=len(queue)))
print('EXTENSION_READY',checks,'QUEUE',len(queue),'WITH_ONE',len(queue)+1);print('EXACT_COMMAND',command);print('MANIFEST_SHA',sha(P/'manifest.json'));print('NAMESPACE',P)
