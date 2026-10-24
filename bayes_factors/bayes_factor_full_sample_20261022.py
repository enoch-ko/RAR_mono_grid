#!/usr/bin/env python3
"""Checkpointed conditional broad-limit Bayes-factor run for canonical SPARC rows.
Reuses saved conditional logP, reconstructs the production baseline ordering once,
and computes the unit-uniform signed-velocity prior denominator with refinement/MC."""
import argparse, csv, json, math, os, tempfile, time
from pathlib import Path
import numpy as np
from bounded_uniform_bayes_factor_20261021 import _finite_support_probability, _prior_mc
from bayes_factor_denominator_pilot_20261021 import _canonical_conditional_row, _sha256
AUDIT=Path('outputs/bayes_factor_reuse_audit_20261021.csv')
OUT=Path('outputs/bayes_factor_full_sample_20261022')
NORM=Path('outputs/bulge_normalized_20261015/inputs/normalization.json')
GRIDS=(8001,16001,32001); SEED=1729; DRAWS=100000

def atomic_json(path,obj):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    fd,tmp=tempfile.mkstemp(prefix=path.name+'.',suffix='.tmp',dir=path.parent)
    try:
        with os.fdopen(fd,'w') as f:
            json.dump(obj,f,indent=2,allow_nan=False); f.write('\n'); f.flush(); os.fsync(f.fileno())
        os.replace(tmp,path)
    except Exception:
        try: os.unlink(tmp)
        except OSError: pass
        raise

def roster(path):
    with Path(path).open(newline='') as f: rows=list(csv.DictReader(f))
    selected=[r for r in rows if r.get('probability_kind')=='conditional' and r.get('profile')=='constant']
    by={}
    for row in selected:
        name=row['galaxy']
        if name in by: raise ValueError('duplicate conditional row: '+name)
        by[name]=row
    if len(by)!=175: raise ValueError(f'expected 175 unique conditional galaxies; got {len(by)}')
    return by

def normalized_disk(sb,galaxy,record):
    total=float(record['Ltotal_1e9']); disk=float(record['Ldisk_1e9']); bulge=float(record['Lbulge_1e9']); kappa=float(record['kappa'])
    if not np.isclose(disk+bulge,total,rtol=1e-12,atol=1e-10) or not np.isclose(kappa,disk/total,rtol=1e-12,atol=1e-12):
        raise ValueError(f'{galaxy}: invalid once-only disc luminosity normalization')
    corrected=np.asarray(sb,float)*kappa
    if not np.all(np.isfinite(corrected)): raise ValueError(f'{galaxy}: normalized disc SB non-finite')
    return corrected,kappa

def calculate(galaxy,auditrow,dataroot,densroot,tablepath,norm,outroot,grids,draws,seed,normpath):
    from mono_SPARC_RAR import solve_disk_basis
    from utils import get_SBdisk,load_sparc_table,vel2acc,compute_A_B
    csvpath=Path(auditrow['source_csv']); metapath=Path(auditrow['metadata_file']); meta=json.loads(metapath.read_text())
    _,logp=_canonical_conditional_row(csvpath,galaxy)
    if meta.get('conditional_probability') is not True: raise ValueError('metadata does not pin conditional probability')
    if not math.isclose(float(meta.get('disk_ml_median',-1)),.5): raise ValueError('disk M/L baseline mismatch')
    if not math.isclose(float(meta.get('gas_ml_median',-1)),1.): raise ValueError('gas M/L baseline mismatch')
    bulge_ml=float(meta.get('bulge_ml_median',.7))
    if not math.isclose(bulge_ml,.7): raise ValueError('bulge M/L baseline mismatch')
    table=load_sparc_table(tablepath); hit=table.loc[table['Galaxy'].astype(str)==galaxy]
    if len(hit)!=1: raise ValueError(f'master rows={len(hit)}')
    master=hit.iloc[0]; rot=Path(dataroot)/f'{galaxy}_rotmod.dat'; raw=np.atleast_2d(np.loadtxt(rot))
    if raw.shape[1]!=8: raise ValueError(f'{galaxy}: rotmod columns={raw.shape[1]}')
    rad,sb=get_SBdisk(galaxy,dens_root=densroot,table=table)
    if galaxy not in norm: raise ValueError(f'{galaxy}: missing normalization record')
    record=norm[galaxy]; sb,kappa=normalized_disk(sb,galaxy,record)
    if not np.isclose(float(master['L']),float(record['Ltotal_1e9']),rtol=0,atol=1e-10): raise ValueError('normalization catalogue luminosity mismatch')
    if not np.isclose(compute_A_B(rad,sb)[0],float(record['Ldisk_1e9'])*1e9,rtol=1e-12): raise ValueError('normalized disk luminosity integral mismatch')
    dens=Path(densroot)/f'{galaxy}.dens'
    if len(rad)<2 or len(rad)!=len(sb) or not dens.is_file(): raise ValueError('invalid full .dens source/profile')
    rin=raw[:,0]; mask=np.isfinite(rin)&(rin>=float(rad[0]))&(rin<=float(rad[-1]))
    r=rin[mask]; vo=raw[mask,1]; ev=raw[mask,2]; vg=raw[mask,3]; vb=raw[mask,5]
    if len(r)<2 or not all(np.all(np.isfinite(x)) for x in (r,vo,ev,vg,vb)) or np.any(r<=0) or np.any(ev<=0): raise ValueError('production overlap mask has invalid scorer inputs')
    zd=.196*float(master['Rdisk'])**.633; basis,warnings,error=solve_disk_basis(rad,sb,zd,r)
    if error: raise RuntimeError('one baseline force solve failed: '+str(error))
    vbar2=.5*np.asarray(basis,float)+vg*np.abs(vg)+vb**2*bulge_ml
    if vbar2.shape!=r.shape or not np.all(np.isfinite(vbar2)): raise ValueError('non-finite/misaligned fixed baryonic field')
    gb=np.asarray(vel2acc(vbar2,r),float)
    if not np.all(np.isfinite(gb)): raise ValueError('non-finite gbar')
    order=np.argsort(-gb,kind='stable')
    receipt={'galaxy':galaxy,'saved_conditional_logP':logp,'r_kpc':r.tolist(),'gbar':gb.tolist(),'vobs':vo.tolist(),'sigma_v':ev.tolist(),'descending_gbar_order':order.tolist(),'tie_atol':0.,'bulge_ml':bulge_ml,'kappa_applied_once':kappa,'mask':'finite kinematic radius in inclusive full dens range; finite Vobs,errV,Vgas,Vbul then required, r>0,errV>0','source_csv':str(csvpath),'source_csv_sha256':_sha256(csvpath),'metadata':str(metapath),'metadata_sha256':_sha256(metapath),'rotmod_sha256':_sha256(rot),'dens_sha256':_sha256(dens),'normalization_sha256':_sha256(normpath),'table_sha256':_sha256(tablepath),'force_calls':1,'solver_warnings':list(warnings)}
    root=Path(outroot); rp=root/'receipts'/f'{galaxy}.json'; atomic_json(rp,receipt)
    refinement=[]
    for ngrid in grids:
        value=_finite_support_probability(r,gb,vo,ev,1.,ngrid,0.,posterior=False).log_probability
        refinement.append({'grid':ngrid,'log_c':value})
    logs=[x['log_c'] for x in refinement]
    if not all(math.isfinite(x) for x in logs):
        result={'galaxy':galaxy,'status':'unresolved_nonfinite_recursion','saved_conditional_logP':logp,'grid_refinement':refinement,'receipt':str(rp),'BF':None}
    else:
        logc=logs[-1]; mc=_prior_mc(r,gb,1.,draws,0.,seed,logc); logbf=logp-logc
        result={'galaxy':galaxy,'status':'complete_numerical_not_yet_claim_validated','saved_conditional_logP':logp,'log_c':logc,'log_BF':logbf,'BF':math.exp(logbf) if logbf<709 else None,'grid_refinement':refinement,'prior_mc':mc.__dict__,'seed':seed,'draws':draws,'receipt':str(rp),'interpretation':'conditional broad-prior-limit ratio; not nuisance-integrated evidence'}
    atomic_json(root/'results'/f'{galaxy}.json',result); return result

def self_test():
    sb,k=normalized_disk(np.array([2.,4.]),'mock',{'Ltotal_1e9':10.,'Ldisk_1e9':8.,'Lbulge_1e9':2.,'kappa':.8})
    assert k==.8 and np.array_equal(sb,np.array([1.6,3.2]))
    r=np.array([1.,2.,3.]); vo=np.array([8.,6.,4.]); ev=np.ones(3)
    a=_finite_support_probability(r,np.array([2.,-1.,-2.]),vo,ev,1.,1001,0.,posterior=False)
    assert math.isfinite(a.log_probability)
    tied=_finite_support_probability(r,np.ones(3),vo,ev,1.,1001,0.,posterior=False)
    assert math.isfinite(tied.log_probability) and tied.log_probability==_finite_support_probability(r,np.ones(3),vo,ev,1.,1001,0.,posterior=False).log_probability
    print('SELFTEST PASS: once-only kappa, signed negative gbar, deterministic exact-tie handling')

def main():
    p=argparse.ArgumentParser(); mode=p.add_mutually_exclusive_group(required=True); mode.add_argument('--plan',action='store_true'); mode.add_argument('--run',action='store_true'); mode.add_argument('--self-test',action='store_true')
    p.add_argument('--audit',type=Path,default=AUDIT); p.add_argument('--output',type=Path,default=OUT); p.add_argument('--data-root',type=Path); p.add_argument('--dens-root',type=Path); p.add_argument('--table',type=Path); p.add_argument('--normalization',type=Path,default=NORM); p.add_argument('--grids',nargs='+',type=int,default=GRIDS); p.add_argument('--draws',type=int,default=DRAWS); p.add_argument('--seed',type=int,default=SEED); p.add_argument('--limit',type=int); a=p.parse_args()
    if a.self_test: self_test(); return
    data=roster(a.audit)
    if a.plan: print(json.dumps({'galaxies':len(data),'one_force_solve_each':True,'saved_conditional_logP_reused':True,'grids':a.grids,'mc_draws_per_galaxy':a.draws,'normalization':str(a.normalization),'no_science_run':True},indent=2)); return
    if a.output.exists(): raise FileExistsError('refusing output collision: '+str(a.output))
    from mono_SPARC_RAR import DEFAULT_DATA_ROOT,DEFAULT_DENS_ROOT,DEFAULT_TABLE_PATH
    dr=a.data_root or Path(DEFAULT_DATA_ROOT); den=a.dens_root or Path(DEFAULT_DENS_ROOT); tab=a.table or Path(DEFAULT_TABLE_PATH)
    normdata=json.loads(a.normalization.read_text()); norm=normdata['galaxies']; a.output.mkdir(parents=True,exist_ok=False)
    start=time.monotonic(); ledger=[]; counts={'complete':0,'unresolved':0,'failed':0}
    for i,(galaxy,row) in enumerate(sorted(data.items())):
        if a.limit is not None and i>=a.limit: break
        try:
            result=calculate(galaxy,row,dr,den,tab,norm,a.output,a.grids,a.draws,a.seed+i,a.normalization); status=result['status']
            if status=='complete_numerical_not_yet_claim_validated': counts['complete']+=1
            else: counts['unresolved']+=1
            ledger.append({'galaxy':galaxy,'status':status,'log_BF':result.get('log_BF')})
        except Exception as e:
            counts['failed']+=1; ledger.append({'galaxy':galaxy,'status':'failed','error':f'{type(e).__name__}: {e}'})
        atomic_json(a.output/'checkpoint.json',{'counts':counts,'ledger':ledger,'elapsed_seconds':time.monotonic()-start,'next_index':i+1})
        if time.monotonic()-start>=590: break
    atomic_json(a.output/'completion.json',{'counts':counts,'not_started':max(0,len(data)-len(ledger)),'elapsed_seconds':time.monotonic()-start,'grids':a.grids,'draws':a.draws,'seed':a.seed,'ledger':ledger})
    print(json.dumps({'output':str(a.output),'counts':counts,'not_started':max(0,len(data)-len(ledger)),'elapsed_seconds':time.monotonic()-start},indent=2))
if __name__=='__main__': main()
