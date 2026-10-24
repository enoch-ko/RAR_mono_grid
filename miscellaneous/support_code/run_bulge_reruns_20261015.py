#!/usr/bin/env python
"""Frozen, isolated 32-galaxy bulge-normalized reruns; NEVER submits jobs.
Existing source files are read-only. Runtime injection applies only in this
process. Preparation and tests are not approval to submit the 36 commands.
"""
import os
for _key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_key] = '1'
import argparse
import functools
import hashlib
import io
import json
import shlex
import subprocess
import sys
import time
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT/'outputs/bulge_normalized_20261015'
GP_FOLDER = 'gp_d5000_a16'
COUNTS = {'constant':4, 'linear':8, 'qp':8, 'gp':16}
PUBLIC = {'bulges':'https://astroweb.case.edu/SPARC/Bulges.mrt',
          'catalogue':'https://astroweb.case.edu/SPARC/SPARC_Lelli2016c.mrt'}
GP = dict(directions=5000, amplitude_steps=16, ell_over_rd=[.5,1.,2.,4.],
          seed=20261013, replicate=0, sampling_scheme='nested-v1',
          checkpoint_mode='retained-only', max_retained_bytes=8*1024**2,
          min_ml=.01, max_ml=10., length_fraction=.3, gbar_floor_relative=1e-8,
          gas_ml=1., bulge_ml=.7, gas_ml_dex=.04, bulge_ml_dex=.1,
          ml_samples=9, prob_grid_size=400, zpadding=8.,
          boundary_log_tolerance=2e-5, boundary_max_refinements=8)
QP = dict(knots_from_dens=True, knots=8, min_ml=.01, max_ml=10., retry_max_ml=20.,
          length_fraction=.3, gas_ml=1., bulge_ml=.7, gas_ml_dex=.04, bulge_ml_dex=.1,
          ml_samples=9, prob_grid_size=400, zpadding=8., gbar_floor_relative=1e-8,
          qp_backend='osqp', osqp_eps=1e-9, osqp_max_iter=100000, osqp_time_limit=10.)
FILES = ('utils.py','mono_SPARC_RAR.py','mono_SPARC_RAR_boundary.py',
         'inverse_rar_diagnostic.py','probability_recursion.py','gp_candidate_retention.py',
         'gp_rar_search_retained.py','plot_inverse_rar_comparison.py',Path(__file__).name)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as stream:
        json.dump(data, stream, indent=2, allow_nan=False)
        stream.write('\n')


def load(output):
    manifest = json.loads((output/'batch_manifest.json').read_text())
    if manifest['output_root'] != str(output) or manifest['gp_settings'] != GP or manifest['qp_settings'] != QP:
        raise ValueError('manifest settings/root mismatch')
    for field in ('source_sha256','input_sha256'):
        for name, digest in manifest[field].items():
            if sha(name) != digest:
                raise RuntimeError('frozen dependency changed: '+name)
    return manifest


def allocate(items, count):
    capacity = 32//count
    result = [dict(shard=i+1, galaxies=[], weight=0) for i in range(count)]
    for item in sorted(items, key=lambda x:(-x['weight'],x['galaxy'])):
        row = min((r for r in result if len(r['galaxies'])<capacity), key=lambda r:(r['weight'],r['shard']))
        row['galaxies'].append(item['galaxy']);row['weight'] += item['weight']
    for row in result: row['galaxies'].sort()
    assert sorted(g for row in result for g in row['galaxies']) == sorted(x['galaxy'] for x in items)
    return result


def commands(output):
    result=[]
    for kind,count in COUNTS.items():
        comment = '2 days' if kind=='gp' else '5 hrs'
        for shard in range(1,count+1):
            command = ['addqueue','-q','cmb','-n','1','-m','8','-c',comment,
                       str(Path(__file__).resolve()),'--output-root',str(output),
                       '--kind',kind,'--shard',str(shard)]
            result.append((kind,shard,shlex.join(command)))
    return result


def prepare(output):
    import numpy as np
    import utils
    if output.exists(): raise FileExistsError('new output root required')
    downloaded={key:urllib.request.urlopen(url, timeout=30).read() for key,url in PUBLIC.items()}
    table=utils.load_sparc_table();names=set(table.Galaxy)
    bulges={}
    for line in downloaded['bulges'].decode().splitlines():
        fields=line.split()
        if fields and fields[0] in names:
            if fields[0] in bulges: raise ValueError('duplicate public bulge row')
            bulges[fields[0]]=float(fields[1])
    if set(bulges)!=names or len(names)!=175: raise ValueError('public/local galaxy coverage differs')
    official={}
    for line in downloaded['catalogue'].decode().splitlines():
        fields=line.split()
        if fields and fields[0] in names:
            if fields[0] in official:raise ValueError('duplicate official catalogue row')
            official[fields[0]]=dict(D=float(fields[2]),L=float(fields[7]))
    if set(official)!=names: raise ValueError('official reference catalogue parsing/coverage failed')
    records={};inputs={str(Path(utils.DEFAULT_TABLE_PATH)):sha(utils.DEFAULT_TABLE_PATH)}
    for _,row in table.iterrows():
        galaxy=row.Galaxy;lb=bulges[galaxy]
        if not np.isfinite(lb) or lb<0: raise ValueError('invalid public bulge luminosity')
        if not (np.isclose(row.D,official[galaxy]['D'],rtol=0,atol=1e-8) and
                np.isclose(row.L,official[galaxy]['L'],rtol=0,atol=1e-8)):
            raise ValueError('distance/luminosity differs from public reference: '+galaxy)
        if not 0<=lb<float(row.L): raise ValueError('nonpositive disc luminosity: '+galaxy)
        records[galaxy]=dict(Ltotal_1e9=float(row.L),Lbulge_1e9=lb,
                             Ldisk_1e9=float(row.L)-lb,kappa=1-lb/float(row.L),
                             reference_distance_Mpc=float(row.D))
    targets=sorted(g for g in names if bulges[g]>0)
    if len(targets)!=32: raise ValueError('expected 32 affected galaxies')
    oldroot=ROOT/'outputs/gp_rar_full_retained_20261014_d5000_a16'
    oldmanifest=json.loads((oldroot/'batch_manifest.json').read_text())
    weights={x['galaxy']:x for x in oldmanifest['galaxies']}
    items=[]
    import mono_SPARC_RAR as mono
    for galaxy in targets:
        dens=Path(utils.DEFAULT_DENS_ROOT)/(galaxy+'.dens');rot=mono.DEFAULT_DATA_ROOT/(galaxy+'_rotmod.dat')
        inputs.update({str(dens):sha(dens),str(rot):sha(rot)})
        a=np.loadtxt(dens);b=np.loadtxt(rot,ndmin=2)
        if not np.any(a[:,2]>0) or not np.any(b[:,5]!=0): raise ValueError('bulge classifications differ: '+galaxy)
        items.append(dict(galaxy=galaxy,weight=weights[galaxy]['weight'],
                          knots=weights[galaxy]['knots'],retained_points=weights[galaxy]['retained_points']))
    for galaxy in ['UGC07399']:
        for path in [Path(utils.DEFAULT_DENS_ROOT)/(galaxy+'.dens'),mono.DEFAULT_DATA_ROOT/(galaxy+'_rotmod.dat')]:
            inputs[str(path)]=sha(path)
    sources={str(ROOT/name):sha(ROOT/name) for name in FILES}
    for path in sorted((ROOT/'outputs').glob('real_disk_cache_all175_p101_*_20261005.json')):
        inputs[str(path)]=sha(path)
    if len(list((ROOT/'outputs').glob('real_disk_cache_all175_p101_*_20261005.json')))!=5:
        raise ValueError('five historical cache validation shards required')
    for path in [ROOT/'.venv-osqp/bin/python',ROOT/'.venv-osqp/pyvenv.cfg']:
        inputs[str(path)]=sha(path)
    output.mkdir(parents=True)
    parent=output/'inputs/allocation_parent_manifest.json'
    parent.parent.mkdir(exist_ok=True)
    with parent.open('xb') as stream:stream.write((oldroot/'batch_manifest.json').read_bytes())
    inputs[str(parent)]=sha(parent)
    for key,raw in downloaded.items():
        path=output/'inputs'/('Bulges.mrt' if key=='bulges' else 'SPARC_Lelli2016c.mrt')
        path.parent.mkdir(exist_ok=True)
        with path.open('xb') as stream:stream.write(raw)
        inputs[str(path)]=sha(path)
    write_json(output/'inputs/normalization.json',dict(policy='disc_to_catalogue_total_minus_published_bulge_v1',
        sources=PUBLIC,downloaded_utc=datetime.now(timezone.utc).isoformat(),galaxies=records,
        assumptions=['bulge table at same distance as matching public SPARC catalogue',
                     'shape, inclination, vertical geometry, radial support and zero-removal unchanged',
                     'gas and unit-M/L SPARC bulge force template unchanged',
                     'published nonparametric/extrapolated bulge luminosities are model dependent']))
    inputs[str(output/'inputs/normalization.json')]=sha(output/'inputs/normalization.json')
    manifest=dict(status='prepared_not_submitted',output_root=str(output),targets=targets,
        source_sha256=sources,input_sha256=inputs,gp_settings=GP,qp_settings=QP,
        plan={kind:allocate(items,count) for kind,count in COUNTS.items()},
        old_gp_root=str(oldroot),source_policy='original files read-only; per-worker runtime injection',
        stages={'A':'4 constant + 8 linear + 8 QP','B':'16 GP only after all corrected QP output JSONs exist'},
        note='No log-linear jobs. Finite-search results not calibrated significance/global optima.')
    write_json(output/'batch_manifest.json',manifest)
    review=['# Bulge-normalized rerun submission review','',
        'PREPARED ONLY; NO JOBS SUBMITTED. 32 galaxies per family; 36 jobs total.',
        '1 CPU and 8 GB each. Max36 CPU/288 GB if all concurrent; intended staged20+16.',
        'Comments are expectations, not time limits; GP uses 2 days instead of optimistic5hrs.',
        'Constant/linear:51/400/9, logmean .5+/-1dex, prior.5/.1dex, slope[-5,5], nuisance1/.04dex and.7/.1dex.',
        'Boundary tolerance2e-5,maxrefinements8,padding8; full marginalisation in same grid jobs.',
        'QP uses existing isolated .venv-osqp; OSQP100000iterations/10s,eps1e-9; .01..10 with proven-infeasible-only retry20.',
        'QP numerical failures remain failures; no automatic unapproved escalation or old warm-profile import.',
        'GP5000directions/ell,Rd scales.5/1/2/4,16amplitudes,seed20261013,retained-only8MiBcap.',
        'GP conditional on corrected QP artifacts; only accepted audited fits within GP bounds imported; others explicitly omitted.',
        'Shared output root: '+str(output),'Separate GP root: '+str(output/GP_FOLDER),
        'Old searches and source dependencies preserved; no automatic merge/replacement.',
        'Production is blocked until normalization_validation.json matches current manifest.',
        'No speed optimization is deployed in these commands.','', '## Exact commands',chr(96)*3+'bash']
    for kind,shard,command in commands(output):review.append(command)
    review += [chr(96)*3,'','## Allocation']
    for kind,rows in manifest['plan'].items():
        for row in rows:review.append(kind+' '+str(row['shard'])+': '+', '.join(row['galaxies']))
    review += ['','Preflight: python3 -B '+str(Path(__file__).resolve())+' --output-root '+str(output)+' --dry-run',
               'After each grid family completes: use merge_run_shards.py on its isolated family directory; do not merge old/new normalization shards.']
    with (output/'submission_review.md').open('x') as stream:stream.write('\n'.join(review)+'\n')
    load(output)
    print('PREPARED',output,'JOBS',len(commands(output)),'SUBMITTED',0)


def install(output, manifest):
    import utils
    original=utils.get_SBdisk
    normalization=json.loads((output/'inputs/normalization.json').read_text())
    records=normalization['galaxies']
    def corrected(galaxy, dens_root=utils.DEFAULT_DENS_ROOT, table=None):
        rad,sb=original(galaxy,dens_root=dens_root,table=table)
        return rad,sb*records[galaxy]['kappa']
    utils.get_SBdisk=corrected
    import mono_SPARC_RAR as mono
    import inverse_rar_diagnostic as inverse
    # These modules may already be imported by the preparation/test process.
    mono.get_SBdisk=corrected;inverse.get_SBdisk=corrected
    proof=dict(policy=normalization['policy'],manifest_sha256=sha(output/'batch_manifest.json'),
               normalization_sha256=sha(output/'inputs/normalization.json'),
               launcher_sha256=sha(__file__),sources=PUBLIC)
    historical=inverse.cache_gate
    def gate(galaxy, root):
        data=historical(galaxy,root)
        return dict(data,normalization=proof,
                    normalization_qc='historical basis-algorithm gate plus current scalar validation and current per-result direct checks')
    inverse.cache_gate=gate
    return utils,mono,inverse,records,proof,original


def validate(output):
    import numpy as np
    manifest=load(output);utils,mono,inverse,records,proof,original=install(output,manifest)
    table=utils.load_sparc_table();rows=[]
    for galaxy in manifest['targets']+['UGC07399']:
        rad,old=original(galaxy,table=table);rr,new=utils.get_SBdisk(galaxy,table=table)
        expected=records[galaxy]['Ldisk_1e9']*1e9
        assert np.array_equal(rad,rr)
        assert np.isclose(utils.compute_A_B(rad,new)[0],expected,rtol=1e-12)
        assert np.allclose(new,old*records[galaxy]['kappa'],rtol=0,atol=0)
        if galaxy=='UGC07399':assert np.array_equal(new,old)
    for galaxy in ['UGC07399','NGC0891','NGC7814']:
        rad,old=original(galaxy,table=table);_,new=utils.get_SBdisk(galaxy,table=table)
        rot=np.loadtxt(mono.DEFAULT_DATA_ROOT/(galaxy+'_rotmod.dat'),ndmin=2)
        mask=(rot[:,0]>=rad[0])&(rot[:,0]<=rad[-1]);r=rot[mask,0]
        zd=.196*float(table.loc[table.Galaxy==galaxy].iloc[0].Rdisk)**.633
        before,w1,e1=mono.solve_disk_basis(rad,old,zd,r);after,w2,e2=mono.solve_disk_basis(rad,new,zd,r)
        if e1 or e2 or w1 or w2:raise RuntimeError((galaxy,w1,w2,e1,e2))
        scale=max(float(np.max(np.abs(after))),1e-30)
        relative=float(np.max(np.abs(after-before*records[galaxy]['kappa']))/scale)
        if relative>1e-8:raise ValueError('scalar force check failed: '+galaxy)
        rows.append(dict(galaxy=galaxy,kappa=records[galaxy]['kappa'],force_relative_delta=relative))
    write_json(output/'normalization_validation.json',dict(status='passed_scalar_normalization_checks',
        **proof,luminosity_checks=33,force_checks=rows,
        scope='scalar normalization and common solver only; not probability/search convergence or physical calibration'))
    load(output);print('VALIDATED',rows)


def require_validation(output):
    data=json.loads((output/'normalization_validation.json').read_text())
    if (data.get('status')!='passed_scalar_normalization_checks' or
        data['manifest_sha256']!=sha(output/'batch_manifest.json') or data['launcher_sha256']!=sha(__file__)):
        raise ValueError('stale normalization validation')


def worker(output,kind,galaxies,shard):
    manifest=load(output);require_validation(output)
    utils,mono,inverse,records,proof,original=install(output,manifest)
    import probability_recursion as probability
    mono.log_monotonic_probability_for_gbar=functools.partial(
        probability.boundary_log_probability_for_gbar,log_tolerance=2e-5,max_refinements=8)
    destination=output/(GP_FOLDER if kind=='gp' else kind)
    destination.mkdir(parents=True,exist_ok=True)
    if kind in ('constant','linear'):
        import mono_SPARC_RAR_boundary as boundary
        historical_hashes=boundary.hashes
        boundary.hashes=lambda:dict(historical_hashes(),normalization_launcher=sha(__file__),
                                   normalization_input=sha(output/'inputs/normalization.json'))
        sys.argv=[str(ROOT/'mono_SPARC_RAR_boundary.py'),'--profile',kind,'--grid-size','51',
            '--prob-grid-size','400','--ml-samples','9','--marginalise-gas-bulge','--marginalise-disk',
            '--output-root',str(destination),'--shard-id',str(shard)]
        for galaxy in galaxies:sys.argv+=['--galaxy',galaxy]
        original_process=mono.process_galaxy
        def guarded(*args,**kwargs):
            load(output);return original_process(*args,**kwargs)
        mono.process_galaxy=guarded
        boundary.main()
        path=destination/f'run_metadata.shard_{shard:02d}.json'
        data=json.loads(path.read_text());data['normalization']=proof
        # Metadata is generated by this worker, not an existing production artifact.
        path.write_text(json.dumps(data,indent=2)+'\n')
        load(output);return
    for galaxy in galaxies:
        load(output);started=time.monotonic();target=destination/'galaxies'
        target.mkdir(exist_ok=True)
        if kind=='qp':
            args=SimpleNamespace(**QP,galaxy=galaxy,data_root=mono.DEFAULT_DATA_ROOT,
                dens_root=Path(utils.DEFAULT_DENS_ROOT),validation_root=ROOT/'outputs')
            data=inverse.analyze(args)
            data['normalization']=dict(proof,galaxy=galaxy,**records[galaxy])
            data['probability_numerical']=dict(method='boundary_log_cumulative_v1',log_tolerance=2e-5,max_refinements=8)
            path=target/(galaxy+'_inverse.json')
            import plot_inverse_rar_comparison as plotter
            data['provenance']['normalization_launcher_sha256']=sha(__file__)
            data['provenance']['plot_script_sha256']=sha(plotter.__file__)
            try:
                data['plot']=plotter.render_plot(data,path.with_suffix('.png'),args.data_root,
                                                  args.dens_root,inverse,False)
            except Exception as error:
                data['plot']=dict(status='failed',reason=repr(error))
            write_json(path,inverse.json_safe(data))
        else:
            import gp_rar_search_retained as gp
            gp.SOURCE_FILES=tuple(dict.fromkeys((*gp.SOURCE_FILES,Path(__file__).name,
                                                str(output/'inputs/normalization.json'))))
            qp_path=output/'qp/galaxies'/(galaxy+'_inverse.json')
            if not qp_path.exists():raise FileNotFoundError('run corrected QP stage first: '+str(qp_path))
            saved=json.loads(qp_path.read_text())
            if saved.get('normalization',{}).get('manifest_sha256')!=proof['manifest_sha256']:
                raise ValueError('reference normalization mismatch')
            args=SimpleNamespace(**GP,output_root=target,validation_root=ROOT/'outputs',
                data_root=mono.DEFAULT_DATA_ROOT,dens_root=Path(utils.DEFAULT_DENS_ROOT),retain_probability=[])
            reference=None;reason='not_an_accepted_current_domain_QP'
            rad=utils.get_SBdisk(galaxy)[0]
            hashes=dict(rotmod=sha(mono.DEFAULT_DATA_ROOT/(galaxy+'_rotmod.dat')),
                        dens=sha(Path(utils.DEFAULT_DENS_ROOT)/(galaxy+'.dens')),table=sha(utils.DEFAULT_TABLE_PATH))
            try:
                gp.reference_profile(qp_path,galaxy,rad,.01,10.,args,hashes)
                reference=qp_path;reason='accepted_corrected_QP'
            except (ValueError,KeyError,TypeError) as error:reason=str(error)
            data=gp.galaxy_run(galaxy,args,reference)
            data['normalization']=dict(proof,galaxy=galaxy,**records[galaxy])
            data['corrected_QP_reference_selection']=dict(path=str(qp_path),sha256=sha(qp_path),reason=reason)
            path=target/(galaxy+'_gp_search.json')
            text=gp.limited_json(inverse.json_safe(data),GP['max_retained_bytes'],indent=2)
            with path.open('x') as stream:stream.write(text)
        load(output)
        write_json(destination/(galaxy+'_completion.json'),dict(galaxy=galaxy,status=data['status'],
            seconds=time.monotonic()-started,path=str(path),manifest_sha256=proof['manifest_sha256']))
        print('FINISH',galaxy,data['status'],flush=True)


def run_shard(output,kind,shard):
    manifest=load(output);require_validation(output)
    plan=manifest['plan'][kind][shard-1]
    destination=output/(GP_FOLDER if kind=='gp' else kind)
    receipt=destination/f'shard_{shard:02d}.json'
    if receipt.exists():raise FileExistsError(receipt)
    galaxies=plan['galaxies']
    if kind=='gp':
        missing=[g for g in galaxies if not (output/'qp/galaxies'/(g+'_inverse.json')).exists()]
        if missing:raise ValueError('corrected QP artifacts required before GP submission: '+','.join(missing))
    logs=destination/'logs';logs.mkdir(parents=True,exist_ok=True)
    log=logs/f'shard_{shard:02d}.log'
    command=[sys.executable,'-B',str(Path(__file__).resolve()),'--output-root',str(output),
             '--kind',kind,'--shard',str(shard),'--worker']
    with log.open('x') as stream:
        result=subprocess.run(command,cwd=ROOT,stdout=stream,stderr=subprocess.STDOUT)
    load(output)
    statuses={};numerical=0
    if kind in ('qp','gp'):
        for galaxy in galaxies:
            path=destination/'galaxies'/(galaxy+('_inverse.json' if kind=='qp' else '_gp_search.json'))
            if not path.exists():
                statuses[galaxy]='missing';numerical+=1;continue
            saved=json.loads(path.read_text());statuses[galaxy]=saved['status']
            if kind=='qp':
                numerical+=int(saved['status'] not in ('feasible_weak_order','infeasible'))
                numerical+=int(saved.get('plot',{}).get('status')=='failed')
            else:
                numerical+=int(saved['status']!='finite_search_completed')
                numerical+=sum(n for reason,n in saved['checkpoint_storage']['summary']['failure_counts'].items() if reason!='signed_gbar_below_floor')
    else:
        import csv
        path=destination/f'fully_marginalised.shard_{shard:02d}.csv'
        if path.exists():
            for row in csv.DictReader(path.open()):
                statuses[row['galaxy']]=row['valid'];numerical+=int(row['valid'].lower()!='true')
        else:numerical+=1
    write_json(receipt,dict(kind=kind,shard=shard,galaxies=galaxies,returncode=result.returncode,
                           result_statuses=statuses,numerical_or_output_failures=numerical,
                           log=str(log),manifest_sha256=sha(output/'batch_manifest.json')))
    if result.returncode or numerical:raise SystemExit(result.returncode or 1)


def self_test():
    import numpy as np
    items=[dict(galaxy='G'+str(i),weight=i*i+1) for i in range(32)]
    for count in COUNTS.values():
        plan=allocate(items,count);assert all(len(x['galaxies'])==32//count for x in plan)
    for kappa in [1.,.8635435883,.199217]:
        x=np.array([.1,.5,2.]);D=np.array([[1.,2.,3.],[3.,2.,1.]])
        assert np.allclose((D*kappa)@(x/kappa),D@x)
    assert len(commands(DEFAULT_OUTPUT))==36
    print('PASS allocation, no duplicate targets, scalar physical-profile mapping,36commands; no force/search run')


def main():
    parser=argparse.ArgumentParser(description=__doc__,allow_abbrev=False)
    parser.add_argument('--output-root',type=Path,default=DEFAULT_OUTPUT)
    parser.add_argument('--prepare',action='store_true')
    parser.add_argument('--validate-normalization',action='store_true')
    parser.add_argument('--self-test',action='store_true')
    parser.add_argument('--dry-run',action='store_true')
    parser.add_argument('--kind',choices=COUNTS)
    parser.add_argument('--shard',type=int)
    parser.add_argument('--worker',action='store_true',help=argparse.SUPPRESS)
    args=parser.parse_args();output=args.output_root.resolve()
    if args.self_test:self_test();return
    if args.prepare:prepare(output);return
    if args.validate_normalization:validate(output);return
    manifest=load(output)
    if args.kind=='qp':
        import importlib.util
        if importlib.util.find_spec('osqp') is None:
            environment=ROOT/'.venv-osqp'
            if Path(sys.prefix).resolve()==environment.resolve():
                raise RuntimeError('OSQP unavailable in the isolated QP environment')
            interpreter=environment/'bin/python'
            os.execv(str(interpreter),[str(interpreter),'-B',str(Path(__file__).resolve()),*sys.argv[1:]])
    if args.dry_run:
        print(json.dumps(dict(jobs=36,submitted=0,plan=manifest['plan'],
            normalization_validation_exists=(output/'normalization_validation.json').exists(),
            gp_root=str(output/GP_FOLDER),settings=dict(gp=GP,qp=QP)),indent=2));return
    if args.kind is None or args.shard is None or not 1<=args.shard<=COUNTS[args.kind]:
        parser.error('valid --kind and --shard required')
    if args.worker:worker(output,args.kind,manifest['plan'][args.kind][args.shard-1]['galaxies'],args.shard)
    else:run_shard(output,args.kind,args.shard)


if __name__=='__main__':main()
