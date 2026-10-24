from pathlib import Path
import os,json,re,subprocess,hashlib,datetime
ROOT=Path('/mnt/users/koe/RAR_mono_grid');AUDIT=ROOT/'cleanup_archive_20261024_v1';assert not AUDIT.exists();AUDIT.mkdir();(AUDIT/'active_job_scripts').mkdir();sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest();jobs=subprocess.run(['squeue','-h','-u','koe','-o','%i|%T|%j'],capture_output=True,text=True,check=True).stdout.splitlines();active=[];protected={};references={};texts=[];skipped=[]
def protect(p,reason):
    p=str(Path(p));protected.setdefault(p,[])
    if reason not in protected[p]:protected[p].append(reason)
def strings(v):
    if isinstance(v,str):yield v
    elif isinstance(v,dict):
        for k,z in v.items():yield from strings(k);yield from strings(z)
    elif isinstance(v,list):
        for z in v:yield from strings(z)
pattern=r'/mnt/users/koe/[^\s\"\x27;()<>\x60,\[\]{}]+'
for line in jobs:
    job,state,name=line.split('|',2);x=subprocess.run(['scontrol','show','job',job,'-o'],capture_output=True,text=True,check=True);desc=x.stdout;(AUDIT/'active_job_scripts'/(job+'.job.txt')).write_text(desc);fields={k:v for k,v in re.findall(r'(?:^|\s)(\w+)=([^\s]*)',desc)};script=AUDIT/'active_job_scripts'/(job+'.sh');x=subprocess.run(['scontrol','write','batch_script',job,str(script)],capture_output=True,text=True);assert x.returncode==0 and script.exists(),'active script capture failed '+job
    content=script.read_text();entry=dict(job_id=int(job),state=state,name=name,fields={k:fields.get(k) for k in ['Command','WorkDir','StdOut','StdErr','TimeLimit']},script=str(script),script_sha256=sha(script));active.append(entry)
    for k in ['Command','StdOut','StdErr']:
        if fields.get(k):protect(fields[k],'active job '+job+' '+k)
    for token in re.findall(pattern,content):
        if Path(token).exists():protect(token,'literal path in active job '+job)
manifest=ROOT/'bf_grid_rollout_20261023_v1/manifest.json';m=json.loads(manifest.read_text());protect(manifest,'approved active BF manifest');protect(ROOT/'bf_grid_rollout_20261023_v1','active BF worker/frozen inputs/approvals/submissions')
for v in [m['common_pins'],*[z['pins'] for z in m['baselines'].values()],*[z.get('pins',{}) for fs in m['roster'].values() for z in fs.values()]]:
    for p,h in v.items():assert sha(p)==h,'pin mismatch before cleanup '+p;protect(p,'active BF transitive pinned dependency')
for p in [ROOT/'.venv-osqp',ROOT/'bf_joint_adapter_20261023_v1',ROOT/'bf_three_force_gate_20261023_v1',AUDIT]:protect(p,'current runtime/tool/approval provenance')
for name in json.loads((ROOT/'bf_joint_adapter_20261023_v1/preserved_source_hashes.json').read_text()):protect(ROOT/name,'original source preservation hash record')
for rel in ['gp_rar_search.py','bayes_factor_vmax_scale_check_20261022.py','outputs/bayes_factor_vmax_scale_check_20261022/vmax_scale_check.png','outputs/rerun_boundary_20261006/log_linear/fully_marginalised.shard_01.csv','outputs/bayes_factor_vmax_scale_styled_20261022_v1/vmax_scale_check_styled.png','outputs/bayes_factor_prior_support_comparison_20261022_v1/prior_support_comparison.png','outputs/bayes_factor_ranked_20261022_v1/ranked_galaxies.csv']:protect(ROOT/rel,'open VS Code buffer/asset; tabs have no unsaved indicator')
files=[];dirs=[]
for base,sub,names in os.walk(ROOT,followlinks=False):
    bp=Path(base);sub[:]=[s for s in sub if s not in ['.venv-osqp','archive','__pycache__','.git'] and not (bp/s).is_symlink()]
    if AUDIT==bp or AUDIT in bp.parents:sub[:]=[];continue
    for s in sub:dirs.append(str((bp/s).relative_to(ROOT)))
    for name in names:
        p=bp/name;st=p.lstat();files.append(dict(path=str(p),relative=str(p.relative_to(ROOT)),size=st.st_size,mtime_ns=st.st_mtime_ns,symlink=p.is_symlink()))
        if p.is_symlink():protect(p,'symlink left intact');continue
        if p.suffix not in ['.py','.sh','.md','.txt','.json','.toml','.yaml','.yml','.ini']:continue
        if st.st_size>16*1024*1024:skipped.append(dict(path=str(p),reason='large reference material not read; parent held'));protect(p.parent,'large unexamined reference/provenance');continue
        try:t=p.read_text()
        except (UnicodeError,OSError):continue
        texts.append((str(p),t))
        if p.suffix=='.json':
            try:v=json.loads(t)
            except ValueError:continue
            for value in strings(v):
                if value.startswith('/mnt/users/koe/') and Path(value).exists():references.setdefault(value,[]).append(str(p))
for origin,t in texts:
    for token in re.findall(pattern,t):
        if Path(token).exists():references.setdefault(token,[]).append(origin)
    if origin.endswith(('.py','.sh','.md')):
        for rel in re.findall(r'outputs/[A-Za-z0-9_./-]+',t):
            p=ROOT/rel.rstrip('.')
            if p.exists():references.setdefault(str(p),[]).append(origin)
for p in (ROOT/'outputs').iterdir():
    if p.is_dir() and p.name!='archive':protect(p,'scientific/provenance parent held; only explicit obsolete child candidates considered')
for entry in active:
    for p in re.findall(pattern,Path(entry['script']).read_text()):
        q=Path(p)
        if q.suffix in ['.py','.sh'] and q.exists():
            try:content=q.read_text()
            except OSError:continue
            entry.setdefault('entrypoint_paths',[]).append(str(q))
            for token in re.findall(pattern,content):
                if Path(token).exists():protect(token,'active entrypoint literal dependency '+str(entry['job_id']))
record=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),active_jobs=active,protected_paths=protected,external_path_references=references,files=files,directories=dirs,skipped_large=skipped,exclusions=['.venv-osqp','archive','__pycache__','.git'],no_moves=True)
(AUDIT/'inventory.json').write_text(json.dumps(record,indent=2)+chr(10));logs=[f for f in files if (Path(f['path']).parent==ROOT or 'slurm_logs' in Path(f['path']).parts) and re.search(r'-(\d+)\.(?:out|err)$',f['path'])];generated=[f for f in files if Path(f['path']).parent==ROOT and (f['path'].endswith('.py.sh') or f['path'].endswith('.orig'))];small_dirs=[d for d in dirs if any(k in d for k in ['validation_00','saved_selection_v2_tests','saved_comparison_validation','definition_checks','example_renderer_check']) and d.count('/')<=2]
(AUDIT/'candidates.json').write_text(json.dumps(dict(logs=logs,generated=generated,small_output_dirs=small_dirs),indent=2)+chr(10));print('DEPENDENCY_INVENTORY',len(active),'activejobs',len(files),'files',len(protected),'protectedpaths',len(references),'referencedpaths',len(skipped),'largeheld');print('CANDIDATES',len(logs),'logs',len(generated),'generated',len(small_dirs),'small_output_dirs');print('INVENTORY_RECEIPT',AUDIT/'inventory.json')
