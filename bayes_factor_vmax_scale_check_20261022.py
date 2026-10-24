#!/usr/bin/env python3
"""Denominator scale-invariance diagnostic; numerator fixed at broad-prior limit."""
import csv, hashlib, json, time
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from bounded_uniform_bayes_factor_20261021 import _finite_support_probability

ROOT = Path(__file__).resolve().parent
SRC = ROOT / 'outputs/bayes_factor_full_sample_20261022'
OUT = ROOT / 'outputs/bayes_factor_vmax_scale_check_20261022'
VMAXES = (1.0, 10.0, 100.0, 1000.0)
NGRID = 32001

def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()

def evaluate(receipts, results):
    rows, timings = [], []
    for rp in receipts:
        rec = json.loads(rp.read_text())
        res = json.loads((results / rp.name).read_text())
        if res['status'] != 'complete_numerical_not_yet_claim_validated':
            raise RuntimeError(f"Unexpected status for {rec['galaxy']}")
        args = (np.asarray(rec['r_kpc']), np.asarray(rec['gbar']), np.asarray(rec['vobs']), np.asarray(rec['sigma_v']))
        q = float(rec['saved_conditional_logP'])
        for vmax in VMAXES:
            t0 = time.perf_counter()
            out = _finite_support_probability(*args, vmax, NGRID, float(rec['tie_atol']), False)
            elapsed = time.perf_counter() - t0
            timings.append(elapsed)
            rows.append({'galaxy': rec['galaxy'], 'vmax_km_s': vmax, 'saved_logQ_fixed': q, 'logc': out.log_probability, 'logBF': q-out.log_probability, 'seconds': elapsed, 'receipt': str(rp.relative_to(ROOT))})
    return rows, timings

def main():
    import sys
    pilot = '--pilot' in sys.argv
    if OUT.exists(): raise FileExistsError(f"Refusing to overwrite {OUT}")
    receipts = sorted((SRC/'receipts').glob('*.json')); results = SRC/'results'
    if len(receipts)!=175 or len(list(results.glob('*.json')))!=175: raise RuntimeError('Expected exactly 175 receipt/result JSON files')
    if pilot: receipts=receipts[:4]
    start=time.perf_counter(); rows,timings=evaluate(receipts,results); runtime=time.perf_counter()-start
    print(f"pilot={pilot} galaxies={len(receipts)} calls={len(rows)} runtime_s={runtime:.3f} max_call_s={max(timings):.3f}")
    if pilot: return
    if runtime>=120: raise TimeoutError(f"Full pass took {runtime:.1f}s; output not written (120s limit)")
    OUT.mkdir(parents=True,exist_ok=False)
    base={r['galaxy']:r['logBF'] for r in rows if r['vmax_km_s']==1.0}
    ranked={g:i+1 for i,(g,_) in enumerate(sorted(base.items(),key=lambda kv:kv[1]))}
    for r in rows:
        r['rank_by_vmax1_logBF']=ranked[r['galaxy']]; r['delta_logBF_vs_vmax1']=r['logBF']-base[r['galaxy']]; r['delta_logc_vs_vmax1']=-r['delta_logBF_vs_vmax1']
    with (OUT/'vmax_scale_check.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    fig,(ax,residual)=plt.subplots(2,1,figsize=(13,8),sharex=True,gridspec_kw={'height_ratios':[3,1]})
    colors=['black','#2878b5','#e07a1f','#4b9b58']
    for vmax,color in zip(VMAXES,colors):
        ss=sorted((r for r in rows if r['vmax_km_s']==vmax),key=lambda r:r['rank_by_vmax1_logBF']); x=[r['rank_by_vmax1_logBF'] for r in ss]
        ax.plot(x,[r['logBF']/np.log(10) for r in ss],color=color,lw=1.1,label=f'Vmax={vmax:g} km/s')
        residual.plot(x,[r['delta_logc_vs_vmax1'] for r in ss],color=color,lw=1.0,label=f'{vmax:g}')
    ax.axhline(0,color='0.4',ls='--',lw=1,label='BF = 1'); ax.set_ylabel('log10 BF (fixed saved numerator)'); ax.legend(ncol=3,fontsize=9)
    residual.axhline(0,color='0.4',lw=.8); residual.set_ylabel('Δ log c'); residual.set_xlabel('Galaxy rank (sorted once by Vmax=1 broad-limit log BF)')
    fig.suptitle('Denominator scale-invariance check; numerator fixed at broad-prior limit'); fig.tight_layout()
    fig.savefig(OUT/'vmax_scale_check.png',dpi=180); fig.savefig(OUT/'vmax_scale_check.pdf'); plt.close(fig)
    diffs=[abs(r['delta_logc_vs_vmax1']) for r in rows]
    manifest={'diagnostic':'denominator scale-invariance check; numerator fixed at broad-prior limit','not_a_claim':'not finite-prior BF sensitivity or full-prior robustness','galaxies':len(receipts),'vmax_km_s':VMAXES,'grid_points':NGRID,'runtime_seconds':runtime,'max_abs_delta_logc_from_vmax1':max(diffs),'mean_seconds_per_probability_call':float(np.mean(timings)),'input_results_dir':str(SRC),'helper':'bounded_uniform_bayes_factor_20261021._finite_support_probability','helper_sha256':sha(ROOT/'bounded_uniform_bayes_factor_20261021.py'),'script_sha256':sha(Path(__file__)),'files':{p.name:sha(p) for p in OUT.iterdir() if p.is_file()}}
    (OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n'); print(json.dumps(manifest,indent=2))

if __name__=='__main__': main()
