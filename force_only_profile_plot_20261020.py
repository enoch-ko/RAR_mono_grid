#!/usr/bin/env python3
"""Bounded D564-8 force-only contrast; plan/mock never call Casertano or probabilities."""
import argparse, csv, hashlib, json, math
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent
SEL = ROOT / "outputs/saved_profile_comparison_20261019_exactcut_v1"
LEDGER = ROOT / "outputs/rerun_boundary_20261006"
GP_PATH = ROOT / "outputs/science_20261016/gp_d5000_a16/D564-8_gp_search.json"
DENS_PATH = ROOT / "../SPARC_RAR/BulgeDiskDec_LTG/D564-8.dens"
ROT_PATH = ROOT / "../data/D564-8_rotmod.dat"
TABLE_PATH = ROOT / "SPARC_table.txt"
CUT = 0.0013498980316300957
OUT = ROOT / "outputs/saved_profile_comparison_20261019_exactcut_v1/plot_verification_20261020_v1/force_only_run_20261020_v1"

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1<<20),b''): h.update(block)
    return h.hexdigest()

def read_csv(path):
    with Path(path).open(newline='') as f: return list(csv.DictReader(f))

def jrad(rad, sb, ups):
    q=np.asarray(ups,float)/0.5
    return float(np.trapz((q-1.0)**2,rad)/(rad[-1]-rad[0]))

def select_inputs():
    from utils import load_sparc_table, get_SBdisk, compute_A_B, compute_A_beta
    table=load_sparc_table(TABLE_PATH)
    ix=np.flatnonzero(table['Galaxy'].astype(str).to_numpy()=='D564-8')
    if len(ix)!=1: raise ValueError('D564-8 must occur exactly once in SPARC table')
    row=table.iloc[ix[0]]
    rad,sb=get_SBdisk('D564-8',dens_root=DENS_PATH.parent,table=table)
    rad=np.asarray(rad,float); sb=np.asarray(sb,float)
    if rad.ndim!=1 or sb.shape!=rad.shape or np.any(np.diff(rad)<=0): raise ValueError('Invalid .dens arrays')
    # Constant: least J_rad among valid nuisance-marginalised cells above the raw-P cut.
    cfiles=sorted((LEDGER/'constant').glob('profile_grid.shard_*.csv'))
    crows=[]
    for path in cfiles:
        crows += [r for r in read_csv(path) if r.get('galaxy')=='D564-8' and r.get('probability_kind')=='marginalised_gas_bulge' and r.get('valid','').lower()=='true' and math.exp(float(r['log_probability']))>=CUT]
    cbest=min(crows,key=lambda r:jrad(rad,sb,np.full(rad.size,float(r['upsilon_disk_mean'])))) if crows else None
    lfile=SEL/'linear_marginalised_rawP_leastJrad_selected_profiles_v7.csv'
    lrows=[r for r in read_csv(lfile) if r.get('galaxy')=='D564-8' and r.get('probability_kind')=='marginalised_gas_bulge']
    if len(lrows)!=1: raise ValueError('Expected one D564-8 linear v7 marginalised selection')
    lrow=lrows[0]
    logrows=[]
    for path in sorted((LEDGER/'log_linear').glob('profile_grid.shard_*.csv')):
        logrows += [r for r in read_csv(path) if r.get('galaxy')=='D564-8' and r.get('probability_kind')=='marginalised_gas_bulge' and r.get('valid','').lower()=='true' and math.exp(float(r['log_probability']))>=CUT]
    if not logrows: raise ValueError('No qualifying valid log-linear rows')
    A=float(np.trapz(sb*rad,rad))
    def log_profile(r):
        slope=float(r['disk_slope']); _,beta=compute_A_beta(rad,sb,slope)
        if not np.isfinite(beta) or beta<=0: raise ValueError('Invalid log-linear normalisation')
        return float(r['upsilon_disk_mean'])*A/float(beta)*np.exp(slope*rad)
    bestlog=min(logrows,key=lambda r:jrad(rad,sb,log_profile(r)))
    gp=json.loads(GP_PATH.read_text())
    rec=next(x for x in gp['records'] if x.get('original_candidate_index')==185803 and x.get('source')=='GP_direction')
    front=next(x for x in gp['frontier'] if x.get('record_index')==gp['records'].index(rec))
    gr=np.asarray(gp['knots_kpc'],float); gu=np.asarray(front['profile_on_photometric_grid'],float)
    if gr.ndim!=1 or gu.shape!=gr.shape or rad[0]<gr[0] or rad[-1]>gr[-1]: raise ValueError('GP profile does not cover full .dens support; refusing extrapolation')
    gpups=np.interp(rad,gr,gu)
    # The profile arrays are physical Upsilon(R); no probability or search is performed here.
    models=[{'id':'constant_baseline_0p5','source':'fixed baseline','ups':np.full(rad.size,0.5)}]
    if cbest is not None and not np.allclose(float(cbest['upsilon_disk_mean']),0.5,rtol=0,atol=1e-12):
        models.append({'id':'constant_selected_leastJrad','source':'constant saved grid','ups':np.full(rad.size,float(cbest['upsilon_disk_mean']))})
    a=float(lrow['a_ml_at_R0']); b=float(lrow['b_ml_per_kpc'])
    models.append({'id':'linear_v7','source':'linear saved v7','ups':a+b*rad})
    models.append({'id':'loglinear_leastJrad','source':'original log-linear grid','ups':log_profile(bestlog)})
    models.append({'id':'gp_retained_185803','source':'saved finite GP candidate','ups':gpups})
    for m in models:
        if np.any(~np.isfinite(m['ups'])) or np.any(m['ups']<=0): raise ValueError(f"Nonphysical selected profile: {m['id']}")
    raw=np.loadtxt(ROT_PATH)
    if raw.ndim!=2 or raw.shape[1]<8: raise ValueError('Expected 8-column SPARC rotmod file')
    r0=raw[:,0]; mask=np.isfinite(r0)&(r0>=rad[0])&(r0<=rad[-1])
    if mask.sum()<2: raise ValueError('Fewer than two rotmod rows in inclusive .dens overlap')
    kin={'r':r0[mask],'vobs':raw[mask,1],'errv':raw[mask,2],'vgas':raw[mask,3],'vbul':raw[mask,5]}
    meta={'galaxy':'D564-8','dens_support_kpc':[float(rad[0]),float(rad[-1])],'rotmod_rows_input':int(len(raw)),'rotmod_rows_retained':int(mask.sum()),'probability_cut_raw':CUT,'loglinear_qualifying_cells':len(logrows),'constant_qualifying_cells':len(crows),'constant_selection':None if cbest is None else {'upsilon':float(cbest['upsilon_disk_mean']),'log_probability':float(cbest['log_probability'])},'linear_selection':{k:lrow.get(k) for k in ['grid_i','grid_j','raw_probability','a_ml_at_R0','b_ml_per_kpc']},'loglinear_selection':{k:bestlog.get(k) for k in ['grid_i','grid_j','log_probability','upsilon_disk_mean','disk_slope']},'gp_selection':{'original_candidate_index':185803,'source':rec['source'],'J_ref':rec.get('J_ref'),'conditional_logp':rec.get('conditional_logp'),'marginal_logp':rec.get('marginal_logp')},'inputs':{str(p):sha(p) for p in [TABLE_PATH,DENS_PATH,ROT_PATH,lfile,GP_PATH]+cfiles+sorted((LEDGER/'log_linear').glob('profile_grid.shard_*.csv'))}}
    return rad,sb,float(0.196*float(row['Rdisk'])**0.633),kin,models,meta

def unique_models(models):
    groups=[]
    for m in models:
        found=next((g for g in groups if np.array_equal(g['ups'],m['ups'])),None)
        if found is None: groups.append({'ups':m['ups'],'ids':[m['id']]})
        else: found['ids'].append(m['id'])
    return groups

def mock_test():
    ms=[{'id':'base','ups':np.array([.5,.5])},{'id':'same','ups':np.array([.5,.5])},{'id':'different','ups':np.array([.4,.6])}]
    gs=unique_models(ms); calls=[]
    def fake_solver(g): calls.append(g['ids']); return np.ones(2)
    for g in gs: fake_solver(g)
    assert len(calls)==2 and calls[0]==['base','same']
    print('PASS mock: exact duplicate physical profiles share one solve; distinct profile gets one; probability/search functions are not imported or called.')

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--mock-test',action='store_true')
    ap.add_argument('--plan',action='store_true')
    ap.add_argument('--run-forces',action='store_true',help='Runs the explicitly listed force-only solves; requires separate user approval.')
    a=ap.parse_args()
    if a.mock_test: mock_test(); return
    if not (a.plan or a.run_forces) or (a.plan and a.run_forces): ap.error('choose exactly one of --plan or --run-forces')
    rad,sb,zd,kin,models,meta=select_inputs(); groups=unique_models(models)
    meta.update({'force_solver':'mono_SPARC_RAR.solve_disk_basis -> vcdisk(..., z0=zd, rhoz=exp)','zero_radius_policy':'drop exact R=0 before vcdisk','interpolation':'signed W/R over positive .dens nodes; multiply at retained rotmod radii','gas_term':'Vgas*abs(Vgas), median 1.0','bulge_term':'Vbul**2, median 0.7','ties':'not reordered or re-evaluated; saved probabilities only','solver_calls':len(groups),'candidate_profiles':[{'id':m['id'],'source':m['source'],'Jrad':jrad(rad,sb,m['ups']),'min_upsilon':float(m['ups'].min()),'max_upsilon':float(m['ups'].max())} for m in models],'output_directory':str(OUT),'storage_cap_bytes':5_000_000})
    if a.plan:
        print(json.dumps(meta,indent=2)); print('Exact mock: python3 force_only_profile_plot_20261020.py --mock-test'); print('Exact plan: python3 force_only_profile_plot_20261020.py --plan'); print('Execution command (DO NOT RUN before Enoch approval): python3 force_only_profile_plot_20261020.py --run-forces'); print('Time: unbenchmarked; bounded to listed one-galaxy solver calls. Proposed stop ceiling 10 min; output cap 5 MB. No solver or probability call made.')
        return
    if OUT.exists(): raise FileExistsError(f'Refusing to overwrite existing output directory: {OUT}')
    from mono_SPARC_RAR import solve_disk_basis
    from utils import vel2acc
    OUT.mkdir(parents=True,exist_ok=False)
    results={}
    for group in groups:
        response,warnings_list,error=solve_disk_basis(rad,sb*group['ups'],zd,kin['r'])
        if error: raise RuntimeError(f"{group['ids']}: {error}")
        results[tuple(group['ids'])]=(response,warnings_list)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axs=plt.subplots(1,3,figsize=(15,4.8))
    for m in models:
        axs[0].plot(rad,m['ups'],label=m['id'])
        group=next(g for g in groups if m['id'] in g['ids']); disk=results[tuple(group['ids'])][0]
        axs[1].plot(kin['r'],vel2acc(disk,kin['r']),'.-',label=m['id'])
        gbar=vel2acc(disk+kin['vgas']*np.abs(kin['vgas'])+0.7*kin['vbul']**2,kin['r'])
        gobs=vel2acc(kin['vobs']*np.abs(kin['vobs']),kin['r']); ok=(gbar>0)&(gobs>0)
        axs[2].plot(gbar[ok],gobs[ok],'.-',label=m['id'])
    axs[0].set(xlabel='R [kpc]',ylabel='disk M/L',title='Selected physical profiles')
    axs[1].set(xlabel='R [kpc]',ylabel='disk acceleration [m/s²]',title='vcdisk disk field')
    axs[2].set(xscale='log',yscale='log',xlabel='g_bar [m/s²]',ylabel='g_obs [m/s²]',title='RAR (positive points)')
    for ax in axs: ax.legend(fontsize=7); ax.grid(alpha=.2)
    fig.tight_layout(); image=OUT/'D564-8_force_only_contrast.png'; fig.savefig(image,dpi=160); plt.close(fig)
    meta['figure']=str(image); meta['figure_sha256']=sha(image); (OUT/'manifest.json').write_text(json.dumps(meta,indent=2)+'\n')
    print(f'Wrote {image} and {OUT/"manifest.json"}; no probability calculation or search performed.')

if __name__=='__main__': main()
