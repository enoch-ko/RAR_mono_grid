"""Joint fixed-profile BF scoring and audited BF-aware retention; no search on import."""
from pathlib import Path
import collections,gzip,hashlib,json,math,os
import numpy as np
import gp_rar_search_retained_corrected as native
import mono_SPARC_RAR as mono
from utils import _tie_groups
from bf_joint_adapter_20261023_v1 import adapter as a
CUTS=(1.,1.5,2.)
class StorageStop(RuntimeError):pass

def profile_hash(x):return hashlib.sha256(np.asarray(x,dtype='<f8').tobytes()).hexdigest()
def clean(v):
 if isinstance(v,dict):return {str(k):clean(z) for k,z in v.items()}
 if isinstance(v,(list,tuple)):return [clean(z) for z in v]
 if isinstance(v,np.ndarray):return clean(v.tolist())
 if isinstance(v,(float,np.floating)):return float(v) if math.isfinite(v) else None
 if isinstance(v,np.integer):return int(v)
 return v
class Journal:
 def __init__(self,path,raw_limit,compressed_limit):
  self.path=Path(path);self.raw=self.path.open('xb');self.gz=gzip.GzipFile(filename='',mode='wb',fileobj=self.raw,mtime=0);self.raw_limit=raw_limit;self.compressed_limit=compressed_limit;self.bytes=0;self.rows=0;self.closed=False
 def append(self,row):
  data=(json.dumps(clean(row),separators=(',',':'),allow_nan=False)+chr(10)).encode();
  if self.bytes+len(data)>self.raw_limit:raise StorageStop('journal raw byte cap; pending record preserved in failure receipt')
  self.gz.write(data);self.bytes+=len(data);self.rows+=1
  if self.rows%1000==0:self.flush()
 def flush(self):
  self.gz.flush();self.raw.flush()
  if self.raw.tell()>self.compressed_limit:raise StorageStop('journal compressed byte cap; no truncation/retry')
 def close(self):
  if not self.closed:
   self.gz.close();self.raw.close();self.closed=True
  if self.path.stat().st_size>self.compressed_limit:raise StorageStop('closed journal exceeds compressed cap')

class BFRetention:
 def __init__(self,journal=None,max_retained=2048):
  self.journal=journal;self.max_retained=max_retained;self.total=self.valid=0;self.baseline=self.maximum_p=self.maximum_bf=None;self.winners={c:None for c in CUTS};self.frontier=[];self.reasons=collections.Counter();self.examples=[];self.pending_record=None
 def add(self,record,xi):
  xi=np.asarray(xi,float);a.require(xi.ndim==1 and np.isfinite(xi).all(),'nonfinite candidate profile');row=dict(record,original_candidate_index=self.total,profile_sha256=profile_hash(xi));self.pending_record=row
  if row.get('valid'):
   required=['J_ref','conditional_logp','marginal_logp','log_cbar','log_BF','log10_BF'];a.require(all(math.isfinite(row.get(k,math.nan)) for k in required),'valid BF record lacks finite scores');a.require(row['J_ref']>=0 and row['marginal_logp']<=1e-12 and row['conditional_logp']<=1e-12 and row['log_cbar']<=1e-12,'invalid probability/deviation');a.require(abs(row['log_BF']-(row['marginal_logp']-row['log_cbar']))<1e-10,'ratio-of-averages BF mismatch')
  if self.journal:self.journal.append(row)
  self.total+=1;self.pending_record=None;item=(row,xi.copy());index=row['original_candidate_index']
  if row.get('source')=='SPARC_constant_0.5':
   a.require(self.baseline is None,'duplicate baseline');self.baseline=item
  if not row.get('valid'):
   self.reasons[str(row.get('reason','unspecified'))]+=1
   if len(self.examples)<8:self.examples.append(row)
   return row
  self.valid+=1;j=row['J_ref'];bf=row['log_BF'];lp=row['marginal_logp']
  for c in CUTS:
   old=self.winners[c]
   if row['log10_BF']>c and (old is None or (j,index)<(old[0]['J_ref'],old[0]['original_candidate_index'])):self.winners[c]=item
  if self.maximum_p is None or (-lp,j,index)<(-self.maximum_p[0]['marginal_logp'],self.maximum_p[0]['J_ref'],self.maximum_p[0]['original_candidate_index']):self.maximum_p=item
  if self.maximum_bf is None or (-bf,j,index)<(-self.maximum_bf[0]['log_BF'],self.maximum_bf[0]['J_ref'],self.maximum_bf[0]['original_candidate_index']):self.maximum_bf=item
  if not any(r['J_ref']<=j and r['log_BF']>=bf for r,x in self.frontier):
   self.frontier=[(r,x) for r,x in self.frontier if not(j<=r['J_ref'] and bf>=r['log_BF'])];self.frontier.append(item);self.frontier.sort(key=lambda z:(z[0]['J_ref'],-z[0]['log_BF']))
  if len(self.frontier)>self.max_retained:raise StorageStop('BF frontier cap: do not truncate a finite-search frontier')
  return row
 def export(self):
  items=self.frontier+[self.baseline,self.maximum_p,self.maximum_bf]+list(self.winners.values());unique={r['original_candidate_index']:(r,x) for item in items if item is not None for r,x in [item]};ids=sorted(unique);dense={j:i for i,j in enumerate(ids)};records=[unique[i][0] for i in ids];vectors=[unique[i][1] for i in ids];frontier=[dense[r['original_candidate_index']] for r,x in self.frontier]
  def idx(item):return None if item is None else dense[item[0]['original_candidate_index']]
  summary=dict(evaluated=self.total,valid=self.valid,invalid=self.total-self.valid,retained=len(records),discarded_profile_vectors=self.total-len(records),all_candidate_scalars_saved=self.journal is not None,frontier_size=len(frontier),failure_counts=dict(self.reasons),failure_examples=self.examples,baseline_record_index=idx(self.baseline),maximum_probability_record_index=idx(self.maximum_p),maximum_BF_record_index=idx(self.maximum_bf),threshold_winners={str(c):dict(operator='strict >',log10_BF_cut=c,status='not_found_in_finite_search' if self.winners[c] is None else 'point_found_accuracy_unverified',record_index=idx(self.winners[c])) for c in CUTS})
  return records,vectors,frontier,summary

def positive_ml_interval(mean,direction,lower,upper):
 mean=np.asarray(mean);direction=np.asarray(direction);a.require(0<lower<.5<upper and np.isfinite([*mean,*direction]).all(),'invalid positive M/L proposal bounds');lo=0.;hi=math.inf
 for v,d in zip(mean,direction):
  if d>0:hi=min(hi,(upper-v)/d)
  elif d<0:hi=min(hi,(lower-v)/d)
  elif not lower<=v<=upper:return None
 return None if not math.isfinite(hi) or hi<0 else (lo,hi)

class JointEvaluator(native.ProbabilityEvaluator):
 """Capture exact numerator-call fields/orders; calculate c before retention."""
 def __init__(self,*args,order_journal=None,max_orders=50000,**kwargs):
  super().__init__(*args,**kwargs);self.floor=-math.inf;self.order_journal=order_journal;self.max_orders=max_orders;self.prior=a.DenominatorCache();self.ids={};self.events=[];self.fatal=None;self.weights=[float(g)*float(b) for g in self.gas_weights for b in self.bulge_weights];self.prior_requests=0
 def _ordered_logp(self,gbar,*values,**options):
  result=super()._ordered_logp(gbar,*values,**options);groups=tuple(tuple(int(j) for j in z) for z in _tie_groups(np.asarray(gbar),options.get('tie_atol',0.)));event=None
  if math.isfinite(float(result)) and groups in self.memo:
   try:
    lc=self.prior.evaluate(self.r,list(range(len(self.r))),gbar,tie_atol=options.get('tie_atol',0.),max_permutations=options.get('max_tie_permutations',720));self.prior_requests+=1
    if groups not in self.ids:
     if len(self.ids)>=self.max_orders:raise StorageStop('paired order-cache cap; stop search without discard')
     index=len(self.ids);entry=dict(order_id=index,tie_groups=[list(z) for z in groups],log_Q=float(result),log_c=lc,tie_policy=a.TIE_POLICY,denominator_method=a.METHOD)
     if self.order_journal:self.order_journal.append(entry)
     self.ids[groups]=index
    event=(self.ids[groups],float(result),lc,float(np.min(gbar)))
   except StorageStop as e:self.fatal=e;raise
  self.events.append(event);return result
 def evaluate(self,xi):
  self.events=[];self.fatal=None
  try:result=super().evaluate(xi)
  except StorageStop:raise
  except Exception as e:return dict(valid=False,reason=type(e).__name__+': '+str(e))
  if self.fatal:raise self.fatal
  if not result.get('valid'):return dict(result,node_order_ids=[None if z is None else z[0] for z in self.events])
  a.require(len(self.events)==len(self.weights)+1 and all(z is not None for z in self.events),'numerator call/quad/order audit incomplete')
  nodes=self.events[1:];logp=mono.log_weighted_average([z[1] for z in nodes],self.weights);a.require(logp==result['marginal_logp'],'native nuisance weights/call order mismatch');lc=mono.log_weighted_average([z[2] for z in nodes],self.weights);bf=logp-lc
  return dict(result,log_cbar=lc,log_BF=bf,log10_BF=bf/math.log(10),conditional_order_id=self.events[0][0],node_order_ids=[z[0] for z in nodes],minimum_nuisance_gbar_m_s2=min(z[3] for z in nodes),cut_accuracy='unverified',point_cuts={str(c):bf/math.log(10)>c for c in CUTS},scientific_validated=False)
