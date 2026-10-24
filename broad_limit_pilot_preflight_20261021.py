#!/usr/bin/env python3
"""Fixed-theta NGC4068 broad-limit BF pilot; --plan is input-only, --run computes."""
import argparse
import hashlib
import json
import platform
from pathlib import Path
import numpy as np
from bounded_uniform_bayes_factor_20261021 import estimate_broad_limit_bayes_factor
from mono_SPARC_RAR import solve_disk_basis
from utils import get_SBdisk, load_sparc_table, vel2acc
ROOT=Path('/mnt/users/koe/RAR_mono_grid')
DATA=Path('/mnt/users/koe/data')
DENS=Path('/mnt/users/koe/SPARC_RAR/BulgeDiskDec_LTG')
GALAXY='NGC4068'

def sha256(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''): h.update(block)
    return h.hexdigest()

def read_inputs(galaxy,data_root,dens_root,table_path):
    if galaxy!=GALAXY: raise ValueError('adapter is pinned to NGC4068 only')
    rotmod=data_root/f'{galaxy}_rotmod.dat'; dens=dens_root/f'{galaxy}.dens'
    raw=np.loadtxt(rotmod)
    if raw.ndim!=2 or raw.shape[1]!=8: raise ValueError('rotmod schema changed: expected verified 8 columns')
    r0,v0,e0,gas0,_,bul0,_,_=raw.T
    dr=np.loadtxt(dens)[:,0]
    mask=(r0>=np.min(dr))&(r0<=np.max(dr))
    if np.count_nonzero(mask)<2: raise ValueError('fewer than 2 radii overlap .dens support')
    table=load_sparc_table(file=table_path); rows=table.loc[table['Galaxy']==galaxy]
    if len(rows)!=1: raise ValueError('expected exactly one SPARC table row')
    return {'rotmod_path':rotmod,'dens_path':dens,'table_path':table_path,'mask':mask,
      'r':r0[mask],'vobs':v0[mask],'sigma_v':e0[mask],'vgas':gas0[mask],
      'vbul':bul0[mask],'table':table,'catalogue_row':rows.iloc[0]}

def make_receipt(galaxy,paths,v):
    return {'scope':'input/provenance receipt only; no force, P, MC, or search','galaxy':galaxy,
      'sources':{k:{'path':str(p),'sha256':sha256(p)} for k,p in paths.items()},
      'rotmod_columns':['Rad_kpc','Vobs_km_s','errV_km_s','Vgas_km_s','Vdisk_km_s','Vbul_km_s','SBdisk','SBbul'],
      'retained_rotmod_rows_zero_based':np.flatnonzero(v['mask']).tolist(),
      'excluded_rotmod_rows_zero_based':np.flatnonzero(~v['mask']).tolist(),
      'r_kpc':v['r'].tolist(),'Vobs_km_s':v['vobs'].tolist(),'sigmaV_km_s':v['sigma_v'].tolist(),
      'conditioning':{'disk_ML':0.5,'gas_ML_factor':1.0,'bulge':'absent; retained Vbul must be zero'},
      'units':{'R':'kpc','Vobs':'km/s','sigmaV':'km/s','gbar_order':'km^2 s^-2 kpc^-1'},
      'event':'sort decreasing gbar; require strictly decreasing gobs in that sequence; average exact ties only; tie_atol=0',
      'proposal':{'seed':1729,'prior_MC_draws':100000,'grid_sizes':[8001,16001,32001],
        'gaussian_tail_tolerance':1e-8,'wall_timeout_s':300}}

def fixed_gbar(v,dens_root):
    r=v['r']; rad,sb=get_SBdisk(GALAXY,dens_root=str(dens_root),table=v['table'])
    zd=0.196*float(v['catalogue_row']['Rdisk'])**0.633
    basis,warnings,error=solve_disk_basis(rad,sb,zd,r)
    if error: raise RuntimeError(f'Casertano solve failed: {error}; warnings={warnings}')
    if np.any(v['vbul']!=0): raise ValueError('NGC4068 no longer matches bulgeless pin')
    gbar=vel2acc(0.5*basis+v['vgas']*np.abs(v['vgas']),r)
    if not np.all(np.isfinite(gbar)): raise RuntimeError('non-finite fixed-theta gbar')
    return gbar,warnings

def main():
    p=argparse.ArgumentParser(description=__doc__); mode=p.add_mutually_exclusive_group(required=True)
    mode.add_argument('--plan',action='store_true'); mode.add_argument('--run',action='store_true')
    p.add_argument('--galaxy',choices=[GALAXY],required=True); p.add_argument('--seed',type=int,default=1729)
    p.add_argument('--draws',type=int,default=100000); p.add_argument('--grid-sizes',type=int,nargs='+',default=[8001,16001,32001])
    p.add_argument('--tail-tolerance',type=float,default=1e-8)
    p.add_argument('--output-root',type=Path,default=ROOT/'outputs/bayes_factor_pilot_20261021/NGC4068')
    p.add_argument('--data-root',type=Path,default=DATA); p.add_argument('--dens-root',type=Path,default=DENS)
    p.add_argument('--table',type=Path,default=ROOT/'SPARC_table.txt'); a=p.parse_args()
    paths={'rotmod':a.data_root/f'{a.galaxy}_rotmod.dat','dens':a.dens_root/f'{a.galaxy}.dens','table':a.table,
      'mono_SPARC_RAR.py':ROOT/'mono_SPARC_RAR.py','utils.py':ROOT/'utils.py',
      'probability_recursion.py':ROOT/'probability_recursion.py','broad_limit_module':ROOT/'bounded_uniform_bayes_factor_20261021.py'}
    missing=[str(x) for x in paths.values() if not x.is_file()]
    if missing: p.error('missing pinned file(s): '+', '.join(missing))
    v=read_inputs(a.galaxy,a.data_root,a.dens_root,a.table); rec=make_receipt(a.galaxy,paths,v)
    if a.plan: print(json.dumps(rec,indent=2)); return
    if a.seed<0 or a.draws<1 or any(n<32 for n in a.grid_sizes): p.error('seed>=0, draws>0, grids>=32 required')
    if a.output_root.exists(): raise FileExistsError(f'refusing overwrite: {a.output_root}')
    gbar,warnings=fixed_gbar(v,a.dens_root); order=np.argsort(-gbar,kind='stable'); sorted_g=gbar[order]; ties=[]; start=0
    while start<len(order):
        stop=start+1
        while stop<len(order) and sorted_g[stop]==sorted_g[start]: stop+=1
        if stop-start>1: ties.append(order[start:stop].tolist())
        start=stop
    results=[]
    for ngrid in a.grid_sizes:
        x=estimate_broad_limit_bayes_factor(v['r'],gbar,v['vobs'],v['sigma_v'],seed=a.seed,
          ngrid=ngrid,mc_draws=a.draws,tie_atol=0.0,tail_tolerance=a.tail_tolerance)
        results.append({'ngrid':ngrid,'prior_logp':x['prior'].log_probability,
          'logQ_lower':x['posterior_log_probability_lower'],'Q_lower':x['posterior_probability_lower'],
          'Q_upper':x['posterior_probability_upper'],'logBF_lower':x['log_bayes_factor_lower'],
          'logBF_upper':x['log_bayes_factor_upper'],'tail_mass_bound':x['posterior_probability_tail_upper_bound'],
          'prior_MC':x['prior_monte_carlo'].__dict__})
    out={'scope':'single fixed-theta exploratory conditional Q/c pilot; not manuscript evidence',
            'input_receipt':rec,'gbar_km2_s2_kpc':gbar.tolist(),'descending_gbar_order_zero_based':order.tolist(),
      'exact_tie_groups_zero_based':ties,'force_warnings':list(warnings),'results_by_grid':results,
      'python':platform.python_version(),'notes':['Gaussian tail bound is not a quadrature error bound.',
      'Report grid refinement and prior-MC uncertainty separately.','Fixed-theta Q/c is not nuisance-marginalised evidence.']}
    a.output_root.mkdir(parents=True,exist_ok=False)
    (a.output_root/'pilot_receipt.json').write_text(json.dumps(out,indent=2))
    print(json.dumps(out,indent=2))

if __name__=='__main__': main()
