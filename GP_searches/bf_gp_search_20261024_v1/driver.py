#!/usr/bin/env python3
"""Isolated BF-aware GP search. Real execution requires exact approved compute batch."""
from pathlib import Path
import json,math,sys,time,signal,traceback,datetime,hashlib,os
from types import SimpleNamespace
import numpy as np
from scipy.linalg import solve_triangular
R=Path('/mnt/users/koe/RAR_mono_grid');sys.path[:0]=[str(R),str(R/'GP_searches'),str(R/'bayes_factors')]
from bf_gp_search_20261024_v1.core import BFRetention,JointEvaluator,Journal,StorageStop,positive_ml_interval,profile_hash,clean,native,a
import inverse_rar_diagnostic_corrected as inverse
import mono_SPARC_RAR as mono
from utils import load_sparc_table,lognormal_quadrature,vel2acc,_tie_groups
from gp_candidate_retention import limited_json
HERE=Path(__file__).resolve().parent

def save(path,value,limit=8*1024**2):
 text=limited_json(clean(value),limit,indent=2)
 with Path(path).open('x') as f:f.write(text)
def check(pins):
 for p,h in pins.items():a.require(a.file_sha(p)==h,'source/input pin drift '+p)
def prepare(m,g):
 check(m['pins']);check(m['galaxies'][g]['pins']);b=json.loads(Path(m['galaxies'][g]['baseline_receipt']).read_text());s=m['settings'];gate=inverse.cache_gate(g,Path(m['validation_root']));rad,sb,r,vobs,errv,gas,bul,zd,excluded=inverse.galaxy_data(g,Path(m['data_root']),Path(m['dens_root']));norm=json.loads(Path(m['normalization']).read_text())['galaxies'][g]
 a.require(np.array_equal(r,b['r_kpc']) and np.array_equal(vobs,b['vobs']) and np.array_equal(errv,b['sigma_v']),'baseline mask/observations differ');a.require(norm['kappa']==b['kappa_applied_once'],'normalization differs');a.require(math.isclose(norm['Ldisk_1e9']+norm['Lbulge_1e9'],norm['Ltotal_1e9'],rel_tol=1e-12),'luminosity identity')
 raw=np.atleast_2d(np.loadtxt(Path(m['data_root'])/(g+'_rotmod.dat')));from sparc_rotmod_columns import read_rotmod
 a.require(np.array_equal(raw,read_rotmod(Path(m['data_root'])/(g+'_rotmod.dat'))),'independent column schema disagrees');mask=np.isfinite(raw[:,0])&(raw[:,0]>=rad[0])&(raw[:,0]<=rad[-1]);a.require(np.array_equal(gas,raw[mask,3]) and np.array_equal(bul,raw[mask,5]),'gas/bulge column disagreement');row=load_sparc_table().loc[lambda t:t.Galaxy==g];a.require(len(row)==1,'galaxy join mismatch')
 return dict(rad=rad,sb=sb*norm['kappa'],r=r,vobs=vobs,errv=errv,gas=gas,bul=bul,zd=zd,row_ids=np.flatnonzero(mask).tolist(),rd=float(row.iloc[0]['Rdisk']),norm=norm,gate=gate,excluded=excluded)
def groups(g):return tuple(tuple(int(v) for v in z) for z in _tie_groups(np.asarray(g),0.))
def node_fields(acc,x,gn,bn):
 W=acc/(1e6/(x['r']*3.086e19));return [vel2acc(W+x['gas']*np.abs(x['gas'])*gas+x['bul']**2*bul,x['r']) for gas in gn for bul in bn]

def run_galaxy(m,g,stage,out):
 a.require(not out.exists(),'output exists: no overwrite/retry');out.mkdir(parents=True);settings=m['settings'];config=m['stages'][stage];start=time.perf_counter();deadline=time.monotonic()+config['galaxy_wall_cap_seconds'];force_calls=0;scalar=orders=None;retention=None;x=None;last_descriptor=None
 def timed(*_):raise TimeoutError('declared GP process/force timeout; no automatic retry')
 signal.signal(signal.SIGALRM,timed);signal.alarm(config['galaxy_wall_cap_seconds'])
 try:
  x=prepare(m,g);knots=np.asarray(x['rad']);mean=np.full(len(knots),.5);H=inverse.hat_matrix(x['rad'],knots);maxforces=len(knots)+6
  def force(surface):
   nonlocal force_calls
   a.require(force_calls<maxforces,'force-call cap exceeded');a.require(time.monotonic()<deadline,'galaxy deadline reached');force_calls+=1;signal.alarm(min(300,max(1,int(deadline-time.monotonic()))));W,warnings,error=mono.solve_disk_basis(x['rad'],surface,x['zd'],x['r']);signal.alarm(max(1,int(deadline-time.monotonic())));a.require(not error and not warnings and W is not None and np.isfinite(W).all(),'force failure/warnings '+str((error,warnings)));return vel2acc(W,x['r'])
  D=np.column_stack([force(x['sb']*H[:,j]) for j in range(len(knots))]);Lref=np.linalg.cholesky(inverse.kernel(knots,settings['length_fraction']));gn,gw=lognormal_quadrature(settings['gas_ml'],settings['gas_ml_dex'],settings['ml_samples']);bn,bw=lognormal_quadrature(settings['bulge_ml'],settings['bulge_ml_dex'],settings['ml_samples']);nuisance=np.array([vel2acc(x['gas']*np.abs(x['gas'])*gas+x['bul']**2*bul,x['r']) for gas in gn for bul in bn]);quad=[dict(gas_ml=float(gas),bulge_ml=float(bul),weight=float(wg)*float(wb)) for gas,wg in zip(gn,gw) for bul,wb in zip(bn,bw)]
  projections={};rngstates={}
  for i,ell in enumerate(settings['ell_over_rd']):
   projections[str(i)]=native.proposal_cholesky(knots,x['rd'],ell)[0];rngstates[str(ell)]=native.nested_direction_rng(settings['seed'],g,ell,settings['replicate']).bit_generator.state
  np.savez_compressed(out/'response_basis.npz',rad=x['rad'],brightness_normalized_once=x['sb'],r=x['r'],H=H,D_acceleration=D,L_ref=Lref,**{'L_proposal_'+k:v for k,v in projections.items()});a.require((out/'response_basis.npz').stat().st_size<=config['force_basis_bytes'],'force basis storage cap')
  context=dict(schema='bf-aware-gp-context-v1',galaxy=g,stage=stage,manifest_sha256=a.file_sha(HERE/'manifest.json'),input_pins=m['galaxies'][g]['pins'],source_pins=m['pins'],parser='sparc_rotmod_header_v1_gas3_disk4',normalization=dict(policy='disc_to_catalogue_total_minus_published_bulge_v1',applications=1,**x['norm']),radii_kpc=x['r'].tolist(),row_ids=x['row_ids'],vobs=x['vobs'].tolist(),sigma_v=x['errv'].tolist(),vgas=x['gas'].tolist(),vbul=x['bul'].tolist(),knots_kpc=knots.tolist(),rd_kpc=x['rd'],settings=settings,quadrature=quad,domain='signed_finite_gbar_full_nuisance_mass',positive_ML=True,proposal_bounds=[settings['min_ml'],settings['max_ml']],proposal_endpoint_roundoff_tolerance=64*np.finfo(float).eps*settings['max_ml'],J_ref='fixed reference kernel at length_fraction of full knot span; same metric for all proposal ell/Rd',tie_policy=a.TIE_POLICY,denominator_method=a.METHOD,RNG_states=rngstates,RNG_bit_generator='PCG64',profile_hash_encoding='SHA256 little-endian float64 knot vector',profile_replay='same pinned numerical libraries, persisted proposal matrices, nested RNG direction descriptor; verify profile_sha256 before reuse',basis_sha256=a.file_sha(out/'response_basis.npz'),cache_gate=x['gate'],scientific_validated=False,cut_accuracy='unverified')
  save(out/'context.json',context);scalar=Journal(out/'all_candidates.jsonl.gz',config['scalar_raw_bytes'],config['scalar_compressed_bytes']);orders=Journal(out/'paired_order_cache.jsonl.gz',config['order_raw_bytes'],config['order_compressed_bytes']);args=SimpleNamespace(gas_ml=settings['gas_ml'],bulge_ml=settings['bulge_ml'],marginalise_gas_bulge=True,prob_grid_size=settings['prob_grid_size'],zpadding=settings['zpadding'],tie_atol=0.,max_tie_permutations=720,boundary_log_tolerance=settings['boundary_log_tolerance'],boundary_max_refinements=settings['boundary_max_refinements']);evaluator=JointEvaluator(D,x['r'],x['vobs'],x['errv'],x['gas'],x['bul'],args,gn,gw,bn,bw,nuisance,-math.inf,order_journal=orders,max_orders=config['max_orders']);retention=BFRetention(scalar,max_retained=config['max_retained_profiles']);direct_qc={}
  def qc(xi):
   key=profile_hash(xi)
   if key in direct_qc:return direct_qc[key]
   cached=D@xi;direct=force(x['sb']*(H@xi));scale=max(float(np.max(np.abs(direct))),1e-30);delta=float(np.max(np.abs(direct-cached)));a.require(delta<=1e-8*scale,'direct/cache field mismatch')
   lhs=node_fields(direct,x,gn,bn);rhs=node_fields(cached,x,gn,bn);a.require(all(groups(u)==groups(v) for u,v in zip(lhs,rhs)),'direct/cache nuisance node order mismatch: stop, no silent probability repair');result=dict(relative_to_peak=delta/scale,maximum_absolute_acceleration_delta=delta,nuisance_order_signatures_match=True,node_count=len(lhs),force_tolerance=1e-8);direct_qc[key]=result;return result
  qc(mean)
  def add(xi,source,ell=None,direction=None,amplitude=None):
   nonlocal last_descriptor
   a.require(time.monotonic()<deadline,'galaxy deadline reached');xi=np.asarray(xi,float);z=solve_triangular(Lref,xi-mean,lower=True);last_descriptor=dict(source=source,ell_over_rd=ell,direction=direction,amplitude_ml=amplitude,J_ref=float(z@z));values=H@xi;tol=context['proposal_endpoint_roundoff_tolerance']
   if not np.isfinite(values).all() or np.any(values<=0) or np.min(xi)<settings['min_ml']-tol or np.max(xi)>settings['max_ml']+tol:score=dict(valid=False,reason='nonpositive_nonfinite_or_outside_proposal_ML')
   else:score=evaluator.evaluate(xi)
   a.require(time.monotonic()<deadline,'galaxy deadline reached after candidate evaluation');entry=dict(last_descriptor,**score);row=retention.add(entry,xi)
   if retention.total%1000==0:scalar.flush();orders.flush();limited_json(clean(retention.export()[3]),config['retained_json_bytes']);print('GP_PROGRESS',g,retention.total,len(retention.frontier),flush=True)
   return row
  add(mean,'SPARC_constant_0.5')
  for i,ell in enumerate(settings['ell_over_rd']):
   rng=native.nested_direction_rng(settings['seed'],g,ell,settings['replicate']);Lprop=projections[str(i)]
   for direction_id in range(config['directions']):
    direction=Lprop@rng.standard_normal(len(knots));direction/=np.max(np.abs(direction));interval=positive_ml_interval(mean,direction,settings['min_ml'],settings['max_ml']);a.require(interval is not None,'positive-ML direction interval missing')
    for amplitude in native.amplitude_grid(interval,settings['amplitude_steps']):
     if amplitude!=0:add(mean+amplitude*direction,'GP_direction',ell,direction_id,float(amplitude))
   records,vectors,frontier,summary=retention.export();save(out/('progress_ell_'+str(i)+'.json'),dict(status='partial_search_not_resumable',completed_ell_over_rd=ell,summary=summary,records=records,profiles=[v.tolist() for v in vectors]),config['retained_json_bytes']);scalar.flush();orders.flush()
  expected=1+len(settings['ell_over_rd'])*config['directions']*(settings['amplitude_steps']-1);a.require(retention.total==expected and scalar.rows==expected,'candidate scalar count/expected schedule mismatch');records,vectors,frontier,summary=retention.export();important={summary['baseline_record_index'],summary['maximum_probability_record_index'],summary['maximum_BF_record_index'],*[v['record_index'] for v in summary['threshold_winners'].values()]};important.discard(None)
  for index in sorted(important):qc(vectors[index])
  profiles=[dict(record_index=i,xi=v.tolist(),profile_on_photometric_grid=(H@v).tolist(),profile_sha256=profile_hash(v),direct_cache_qc=direct_qc.get(profile_hash(v),dict(status='not_direct_checked'))) for i,v in enumerate(vectors)];scalar.close();orders.close();check(m['pins']);check(m['galaxies'][g]['pins']);result=dict(schema='bf-aware-gp-result-v1',galaxy=g,status='finite_search_complete' if not summary['invalid'] else 'finite_search_with_candidate_failures',diagnostic_only=True,scientific_validated=False,cut_accuracy='unverified',settings=settings,search_budget=config['directions'],retention=summary,records=records,retained_profiles=profiles,frontier_indices=frontier,storage=dict(all_candidate_scalars_saved=True,all_candidate_profile_vectors_saved=False,reconstructible_descriptor_and_hash_for_discarded=True,all_candidates_path='all_candidates.jsonl.gz',paired_order_cache_path='paired_order_cache.jsonl.gz',context_path='context.json',force_basis_path='response_basis.npz',no_automatic_resume=True),force_solves=force_calls,force_solve_cap=maxforces,numerator_order_cache_entries=len(evaluator.memo),denominator_requests=evaluator.prior_requests,paired_order_entries=len(evaluator.ids),cache_hits=evaluator.cache_hits,seconds=time.perf_counter()-start,claim_limit='finite signed-domain BF-aware search; no global optimum, no theory evidence or calibrated sigma')
  a.require(time.monotonic()<deadline,'galaxy deadline reached before publication');save(out/'result.json',result,config['retained_json_bytes']);save(out/'completion.json',dict(status=result['status'],galaxy=g,evaluated=expected,invalid=summary['invalid'],force_solves=force_calls,seconds=result['seconds'],files={p.name:a.file_sha(p) for p in out.iterdir() if p.is_file()},manifest_sha256=a.file_sha(HERE/'manifest.json'),scientific_validated=False));print('GP_DONE',g,expected,summary['retained'],summary['invalid'],flush=True)
 except BaseException as e:
  signal.alarm(0)
  for journal in [scalar,orders]:
   if journal:
    try:journal.close()
    except Exception:pass
  save(out/'failure.json',dict(status='incomplete_requires_review_no_retry',type=type(e).__name__,message=str(e),traceback=traceback.format_exc(),evaluated=None if retention is None else retention.total,pending_record=None if retention is None else retention.pending_record,last_descriptor=last_descriptor,force_solves=force_calls,science_validated=False));raise
 finally:signal.alarm(0)

def main():
 m=json.loads((HERE/'manifest.json').read_text());check(m['pins'])
 if sys.argv[1:]==['--preflight']:
  report=[]
  for g in m['galaxies']:
   x=prepare(m,g);report.append(dict(galaxy=g,radii=len(x['r']),knots=len(x['rad']),force_cap=len(x['rad'])+6,once_kappa=x['norm']['kappa']))
  print('GP_INPUT_PREFLIGHT_OK',report,'forces=0 numerator_evaluations=0');return
 a.require(len(sys.argv)==3 and sys.argv[1] in m['stages'] and sys.argv[2] in m['batches'][sys.argv[1]],'use --preflight or <stage> <batch>');stage,batch=sys.argv[1:];a.require(os.environ.get('SLURM_JOB_ID'),'real GP execution requires approved compute job');approval=json.loads((HERE/('approved_'+stage+'.json')).read_text());a.require(approval.get('status')=='approved_exact_GP_batch' and approval.get('manifest_sha256')==a.file_sha(HERE/'manifest.json'),'exact GP batch approval receipt required')
 if stage=='production':
  review=json.loads((HERE/'reviewed_pilot.json').read_text());a.require(review.get('status')=='approved_for_production' and review.get('manifest_sha256')==a.file_sha(HERE/'manifest.json'),'reviewed real-galaxy GP pilot required')
 root=Path(m['output_root'])/stage/batch;a.require(not root.exists(),'batch output exists; no overwrite');root.mkdir(parents=True);done=[];failed=[];job_start=time.monotonic()
 for g in m['batches'][stage][batch]:
  if time.monotonic()-job_start>=m['stages'][stage]['job_wall_cap_seconds']:failed.append(dict(galaxy=g,reason='batch deadline before galaxy'));break
  try:run_galaxy(m,g,stage,root/g);done.append(g)
  except Exception as e:failed.append(dict(galaxy=g,reason=str(e)))
 save(root/'batch_receipt.json',dict(completed=done,failed=failed,unprocessed=[g for g in m['batches'][stage][batch] if g not in done and g not in [v['galaxy'] for v in failed]],no_auto_retry=True))
 if failed or len(done)!=len(m['batches'][stage][batch]):raise SystemExit(1)
if __name__=='__main__':main()
