#!/mnt/users/koe/RAR_mono_grid/.venv-osqp/bin/python
"""Queue-safe packaging ONLY: positional mode/shard; no queue submission code.

Use --plan for input/argument inspection without running scientific software.
Production execution requires a scheduler job ID. All outputs/claims are exclusive.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent
SCIENCE = ROOT/'outputs/science_20261016'
MANIFEST = SCIENCE/'manifest.json'
MANIFEST_SHA = '2dc247a520a6f9d84caab46f9ee3211dd0be22bec72988c652a4f1b0548086d2'
CONFIG = SCIENCE/'launch_worker_v2.json'
ALLOWED = {'gp':list(range(1,21)), 'reference':[1,6,8,9,15]}
LIMITS = {'gp':{'per_galaxy':172800,'wrapper':172740,'shard':1209600,'memory_gib':8},
          'reference':{'per_galaxy':21600,'wrapper':21540,'shard':172800,'memory_gib':4}}
ACTIVE = None


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def stamp():
    return datetime.now(timezone.utc).isoformat()


def write_new(path, document):
    with Path(path).open('x') as stream:
        json.dump(document,stream,indent=2,allow_nan=False)
        stream.write('\n')


def mode_shard(mode, shard):
    if mode not in ALLOWED:
        raise ValueError('unsupported mode')
    if not str(shard).isdigit():
        raise ValueError('shard must be a decimal integer')
    shard=int(shard)
    if shard not in ALLOWED[mode]:
        raise ValueError('mode/shard outside proposed first25 scope')
    return mode,shard


def load_context():
    config=json.loads(CONFIG.read_text())
    if config.get('schema')!='science-worker-launch-v1':
        raise ValueError('unsupported worker launch config')
    if config.get('worker_sha256')!=sha(__file__) or config.get('worker_path')!=str(Path(__file__).resolve()):
        raise ValueError('worker bytes/path differ from launch receipt')
    if config.get('production_manifest_sha256')!=MANIFEST_SHA or sha(MANIFEST)!=MANIFEST_SHA:
        raise ValueError('production manifest changed')
    if config.get('allowed_mode_shards')!=ALLOWED or config.get('limits')!=LIMITS:
        raise ValueError('launch scope or caps changed')
    manifest=json.loads(MANIFEST.read_text())
    shards={s['shard']:s['galaxies'] for s in manifest['shards']}
    roster=[g for galaxies in shards.values() for g in galaxies]
    if len(shards)!=25 or any(len(g)!=7 for g in shards.values()) or len(set(roster))!=175:
        raise ValueError('invalid or duplicated allocation')
    return config,manifest,shards


def build_argv(mode,galaxy,seed=None):
    common=[str(ROOT/'run_corrected_recovery_pilot.py'),mode,'--manifest',str(MANIFEST),
            '--seconds',str(LIMITS[mode]['wrapper']),'--']
    if mode=='gp':
        return common+['--galaxy',galaxy,'--output-root',str(SCIENCE/'gp_d5000_a16'),
            '--data-root','/mnt/users/koe/data','--dens-root','/mnt/users/koe/SPARC_RAR/BulgeDiskDec_LTG',
            '--validation-root',str(ROOT/'outputs'),'--directions','5000','--amplitude-steps','16',
            '--ell-over-rd','0.5','1','2','4','--length-fraction','0.3','--min-ml','0.01','--max-ml','10',
            '--gbar-floor-relative','1e-8','--seed','20261016','--sampling-scheme','nested-v1','--replicate','0',
            '--checkpoint-mode','retained-only','--max-retained-bytes','8388608',
            '--retain-probability','0.0013498980316300957','--gas-ml','1','--gas-ml-dex','0.04',
            '--bulge-ml','0.7','--bulge-ml-dex','0.1','--ml-samples','9','--prob-grid-size','400',
            '--zpadding','8','--boundary-log-tolerance','2e-5','--boundary-max-refinements','8']
    if mode!='reference' or seed is None:
        raise ValueError('reference mode needs a seed')
    return common+['build','--galaxy',galaxy,'--B','1024','--seed',str(seed),
        '--output',str(SCIENCE/'reference_B1024'/(galaxy+'_B1024.json')),
        '--normalization',str(ROOT/'outputs/bulge_normalized_20261015/inputs/normalization.json'),
        '--data-root','/mnt/users/koe/data','--dens-root','/mnt/users/koe/SPARC_RAR/BulgeDiskDec_LTG',
        '--table',str(ROOT/'SPARC_table.txt'),'--prob-grid-size','400','--zpadding','8',
        '--log-tolerance','2e-5','--max-refinements','8','--seconds','21540','--save-vectors']


def build_plan(mode,shard,galaxies):
    mode,shard=mode_shard(mode,shard)
    if len(galaxies)!=7 or len(set(galaxies))!=7:
        raise ValueError('worker must have exactly seven distinct galaxies')
    entries=[]
    for index,galaxy in enumerate(galaxies):
        if not galaxy or any(c not in 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-_' for c in galaxy):
            raise ValueError('unsafe galaxy token')
        seed=202610160000+7*(shard-1)+index+1 if mode=='reference' else None
        output=SCIENCE/('gp_d5000_a16' if mode=='gp' else 'reference_B1024')/(galaxy+('_gp_search.json' if mode=='gp' else '_B1024.json'))
        guards=[str(output)]
        if mode=='gp':guards.append(str(output.parent/(galaxy+'_gp_checkpoints')))
        entries.append({'galaxy':galaxy,'seed':seed,'argv':build_argv(mode,galaxy,seed),
            'output':str(output),'collision_guards':guards})
    return {'mode':mode,'shard':shard,'entries':entries,'limits':LIMITS[mode],
            'receipt_root':str(SCIENCE/'launch_receipts'), 'log_root':str(SCIENCE/'logs/galaxies')}


def stop_group(proc):
    try:os.killpg(proc.pid,signal.SIGTERM)
    except ProcessLookupError:return
    try:proc.wait(timeout=120)
    except subprocess.TimeoutExpired:
        try:os.killpg(proc.pid,signal.SIGKILL)
        except ProcessLookupError:pass
        proc.wait()


def interrupted(signum,frame):
    if ACTIVE is not None:stop_group(ACTIVE)
    raise SystemExit(128+signum)


def execute(argv,seconds,log_path):
    global ACTIVE
    started=time.monotonic()
    timed_out=False
    with Path(log_path).open('xb') as stream:
        proc=subprocess.Popen(argv,cwd=ROOT,stdin=subprocess.DEVNULL,stdout=stream,
            stderr=subprocess.STDOUT,start_new_session=True)
        ACTIVE=proc
        try:
            try:code=proc.wait(timeout=seconds)
            except subprocess.TimeoutExpired:
                timed_out=True;stop_group(proc);code=124
        finally:ACTIVE=None
    if code!=0:
        with Path(log_path).open('rb') as stream:
            stream.seek(0,2);size=stream.tell();stream.seek(max(0,size-65536))
            tail=stream.read().decode(errors='replace').lower()
        if 'disk quota exceeded' in tail or 'no space left on device' in tail:
            raise OSError('scientific child reported storage exhaustion; stopped remaining galaxies')
    return {'returncode':code,'timed_out':timed_out,'elapsed_seconds':time.monotonic()-started}


def run_plan(plan,job_id,runner=execute,clock=time.monotonic):
    """Runner injection is for unit tests only, never a queue CLI option."""
    if not job_id:
        raise ValueError('scheduler job ID required; scientific execution refused on login node')
    for item in plan['entries']:
        for name in item['collision_guards']:
            if Path(name).exists():raise FileExistsError('existing science/checkpoint output: '+name)
    receipts=Path(plan['receipt_root']);receipts.mkdir(exist_ok=True)
    logs=Path(plan['log_root']);logs.mkdir(exist_ok=True)
    prefix=plan['mode']+'_'+format(plan['shard'],'02d')
    claim=receipts/(prefix+'.claim.json')
    write_new(claim,{'schema':'science-shard-claim-v1','utc':stamp(),'job_id':str(job_id),
        'worker_sha256':sha(__file__),'manifest_sha256':MANIFEST_SHA,'plan':plan})
    journal=receipts/(prefix+'.jsonl')
    started=clock();failed=False;results=[]
    with journal.open('x') as stream:
        for item in plan['entries']:
            remaining=plan['limits']['shard']-(clock()-started)
            if remaining<=0:
                result={'galaxy':item['galaxy'],'status':'unattempted_shard_cap','returncode':124}
                failed=True
            else:
                limit=min(plan['limits']['per_galaxy'],remaining)
                log=logs/(prefix+'_'+item['galaxy']+'.log')
                try:
                    outcome=runner(item['argv'],limit,log)
                    exists=Path(item['output']).is_file()
                    good=outcome['returncode']==0 and exists
                    result=dict(outcome,galaxy=item['galaxy'],seed=item['seed'],argv=item['argv'],
                        log=str(log),output=item['output'],output_present=exists,
                        status='completed' if good else 'failed_or_missing_output',hard_timeout_seconds=limit)
                    if not good:failed=True
                except OSError as error:
                    # Storage/write/process-launch errors must stop, not flood the rest.
                    result={'galaxy':item['galaxy'],'status':'stopped_io_error','returncode':1,'reason':str(error)}
                    results.append(result);stream.write(json.dumps(result)+'\n');stream.flush()
                    raise
            results.append(result);stream.write(json.dumps(result,allow_nan=False)+'\n');stream.flush()
            print(json.dumps(result,allow_nan=False),flush=True)
    summary={'schema':'science-shard-summary-v1','mode':plan['mode'],'shard':plan['shard'],
        'job_id':str(job_id),'utc':stamp(),'manifest_sha256':MANIFEST_SHA,'worker_sha256':sha(__file__),
        'returncode':1 if failed else 0,'results':results,'elapsed_seconds':clock()-started}
    write_new(receipts/(prefix+'.summary.json'),summary)
    return summary['returncode']


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=ALLOWED)
    parser.add_argument('shard')
    parser.add_argument('--plan',action='store_true',help='read-only arguments; never run science')
    args=parser.parse_args()
    mode,shard=mode_shard(args.mode,args.shard)
    config,manifest,allocation=load_context()
    plan=build_plan(mode,shard,allocation[shard])
    if args.plan:
        print(json.dumps(plan,indent=2));return 0
    job_id=os.environ.get('SLURM_JOB_ID')
    if not job_id:raise ValueError('SLURM_JOB_ID absent: use --plan on login node')
    old_handlers={s:signal.signal(s,interrupted) for s in (signal.SIGTERM,signal.SIGINT,signal.SIGALRM)}
    signal.alarm(LIMITS[mode]['shard'])
    try:return run_plan(plan,job_id)
    finally:
        signal.alarm(0)
        for s,handler in old_handlers.items():signal.signal(s,handler)


if __name__=='__main__':
    sys.exit(main())
