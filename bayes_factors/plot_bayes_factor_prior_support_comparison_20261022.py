#!/mnt/users/koe/RAR_mono_grid/.venv-osqp/bin/python
"""Restyle saved denominator-only scale comparison; no scientific evaluation."""
import collections,csv,datetime,hashlib,json,math,itertools
import numpy as np
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
SRC=ROOT/'outputs/bayes_factor_vmax_scale_check_20261022'
OUT=ROOT/'outputs/bayes_factor_prior_support_comparison_20261022_v1'
TITLE='Bayes factors (uniform priors)'
STYLES={1:dict(kind='line',color='black',linewidth=1.2,alpha=1.0),10:dict(kind='scatter',color='#2878b5',marker='o',size=28,alpha=0.4,hollow=True),100:dict(kind='scatter',color='#e07a1f',marker='^',size=30,alpha=0.4,hollow=True),1000:dict(kind='scatter',color='#4b9b58',marker='x',size=20,alpha=0.4,hollow=False)}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def identity_checks():
    from bounded_uniform_bayes_factor_20261021 import _finite_support_probability
    proof=('For independent v_i~Uniform(-V,V), abs(v_i)~Uniform(0,V), exactly the positive prior. '
           'The event depends only on v_i^2/R_i and fixed baryonic order/tie averaging. '
           'Hence c_positive=c_symmetric at every common V; for fixed saved Q, Q/c is identical. '
           'The Gaussian likelihood is not sign invariant, so this does not equate posterior numerators.')
    cases=[]
    u=(np.arange(400)+0.5)/400
    a,b=np.meshgrid(u,u,indexing='ij')
    for radii,expected in [([1.,4.],0.75),([4.,1.],0.25)]:
        r=np.array(radii)
        grid=float(np.mean(a*a/r[0]>b*b/r[1]))
        assert abs(grid-expected)<0.002
        helper=float(np.exp(_finite_support_probability(r,np.array([2.,1.]),np.zeros(2),np.ones(2),1.,32001,0.,False).log_probability))
        assert abs(helper-expected)<1e-4
        cases.append(dict(radii=radii,analytic_probability=expected,independent_midpoint_probability=grid,helper_probability=helper))
    u=(np.arange(24)+0.5)/24
    magnitudes=np.array(list(itertools.product(u,repeat=3)))
    r=np.array([1.,2.3,5.7])
    orders=[(0,1,2),(1,0,2)]
    def scores(v):
        g=v*v/r
        return sum(np.all(np.diff(g[:,order],axis=1)<0,axis=1).astype(float) for order in orders)/2
    base=scores(magnitudes)
    sign_checks=0
    for vmax in [1.,10.,100.,1000.]:
        positive=scores(vmax*magnitudes)
        assert np.array_equal(base,positive)
        for signs in itertools.product([-1.,1.],repeat=3):
            assert np.array_equal(positive,scores(vmax*magnitudes*np.array(signs)))
            sign_checks+=1
    return dict(analytic_identity=proof,independent_midpoint_cases=cases,midpoints_per_axis=24,magnitude_vectors=len(magnitudes),sign_patterns_per_scale=8,scales=[1,10,100,1000],tie_orders=orders,sign_and_scale_checks=sign_checks,max_event_score_difference=0.,production_numerator_evaluations=0,production_denominator_evaluations=0,synthetic_helper_calls=2,helper_event='strict descending v_i squared / R_i after descending g_bar ordering; arithmetic mean over tie permutations',positive_matching_posterior_numerator='not evaluated',all_checks_passed=True)

def main():
    if OUT.exists():raise FileExistsError('Preserve previous rendering: '+str(OUT))
    original=json.load(open(SRC/'manifest.json'))
    for f,h in original['files'].items():
        if sha(SRC/f)!=h:raise ValueError('Source changed: '+f)
    rows=list(csv.DictReader(open(SRC/'vmax_scale_check.csv')))
    groups=collections.defaultdict(list)
    for r in rows:groups[r['galaxy']].append(r)
    assert len(rows)==700 and len(groups)==175
    assert all(len(a)==4 and {float(r['vmax_km_s']) for r in a}=={1,10,100,1000} for a in groups.values())
    assert all(len({r['saved_logQ_fixed'] for r in a})==1 for a in groups.values())
    fig,(ax,residual,ratio)=plt.subplots(3,1,figsize=(13,10),sharex=True,gridspec_kw={'height_ratios':[3,1,1]})
    check=identity_checks()
    artists=[]
    for vmax,style in STYLES.items():
        selected=sorted((r for r in rows if float(r['vmax_km_s'])==vmax),key=lambda r:int(r['rank_by_vmax1_logBF']))
        x=[int(r['rank_by_vmax1_logBF']) for r in selected]
        values=([float(r['logBF'])/math.log(10) for r in selected],[float(r['delta_logc_vs_vmax1']) for r in selected])
        ratio.plot(x,[1.0]*len(x),color=style['color'],lw=1,alpha=0.7)
        for axis,y in zip((ax,residual),values):
            label=rf'$V_{{\max}}={vmax:g}$ km/s'
            if style['kind']=='line':artist,=axis.plot(x,y,color=style['color'],lw=style['linewidth'],alpha=style['alpha'],label=label,zorder=2)
            else:
                kw=dict(marker=style['marker'],s=style['size'],alpha=style['alpha'],linewidths=0.85,label=label,zorder=3)
                if style['hollow']:kw.update(facecolors='none',edgecolors=style['color'])
                else:kw.update(color=style['color'])
                artist=axis.scatter(x,y,**kw)
            artists.append(artist)
    ax.axhline(0,color='0.4',ls='--',lw=0.9,label='BF = 1',zorder=1)
    ax.set_ylabel('log10 BF (fixed saved numerator)')
    ax.legend(ncol=3,fontsize=9,scatterpoints=1)
    residual.axhline(0,color='0.4',lw=0.8,zorder=1)
    residual.set_ylabel('Δ log c')
    ratio.set_xlabel('Galaxy rank (sorted once by Vmax=1 broad-limit log BF)')
    ax.text(0.02,0.80,'Denominators: Uniform(0, Vmax) and Uniform(-Vmax, Vmax)\nExact equality for the squared-velocity ordering event',transform=ax.transAxes,fontsize=9,va='top')
    ratio.set_ylabel('c positive / c symmetric')
    ratio.set_ylim(0.98,1.02)
    ratio.set_yticks([0.99,1.0,1.01])
    ratio.text(0.02,0.78,'Analytic ratio = 1 for all 175 galaxies and all four Vmax values',transform=ratio.transAxes,fontsize=9)
    fig.suptitle(TITLE)
    fig.text(0.5,0.012,'Fixed saved Pmono numerator approximation: identical Pmono/c curves for both denominators.\nA positive-support posterior numerator was not evaluated; this is not a validated positive-prior Bayes factor.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,0.055,1,0.98))
    assert fig._suptitle.get_text()==TITLE
    assert len(ax.collections)==3 and len(residual.collections)==3
    assert all(a.get_alpha()==0.4 for a in ax.collections+residual.collections)
    OUT.mkdir(parents=True,exist_ok=False)
    png=OUT/'prior_support_comparison.png';pdf=OUT/'prior_support_comparison.pdf'
    fig.savefig(png,dpi=180);fig.savefig(pdf);plt.close(fig)
    with (OUT/'identity_check_receipt.json').open('x') as f:json.dump(check,f,indent=2);f.write(chr(10))
    manifest=dict(artifact='denominator prior-support identity comparison with fixed saved numerator',created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),title=TITLE,source_csv=str(SRC/'vmax_scale_check.csv'),source_csv_sha256=sha(SRC/'vmax_scale_check.csv'),source_manifest=str(SRC/'manifest.json'),source_manifest_sha256=sha(SRC/'manifest.json'),renderer=str(Path(__file__).resolve()),renderer_sha256=sha(__file__),style_settings=STYLES,rows=700,galaxies=175,numbers_changed=False,scientific_evaluations=0,synthetic_denominator_checks=2,positive_support_numerator_evaluated=False,positive_denominator_method="analytic equality with saved symmetric denominator, not a new numerical integration",identity_check_receipt_sha256=sha(OUT/"identity_check_receipt.json"),helper_sha256=sha(ROOT/"bounded_uniform_bayes_factor_20261021.py"),diagnostic_scope=original['diagnostic'],claim_limit=original['not_a_claim'],figure_inches=[13,10],png_dpi=180,residual_panel_retained=True,files={p.name:sha(p) for p in [png,pdf]})
    with (OUT/'rendering_manifest.json').open('x') as f:json.dump(manifest,f,indent=2);f.write(chr(10))
    print('RENDERED',png,pdf,'SOURCE_CSV_UNCHANGED',sha(SRC/'vmax_scale_check.csv')==original['files']['vmax_scale_check.csv'])
if __name__=='__main__':main()
