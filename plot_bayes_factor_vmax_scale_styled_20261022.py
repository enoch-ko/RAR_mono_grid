#!/mnt/users/koe/RAR_mono_grid/.venv-osqp/bin/python
"""Restyle saved denominator-only scale comparison; no scientific evaluation."""
import collections,csv,datetime,hashlib,json,math
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
SRC=ROOT/'outputs/bayes_factor_vmax_scale_check_20261022'
OUT=ROOT/'outputs/bayes_factor_vmax_scale_styled_20261022_v1'
TITLE='Bayes factors (uniform priors)'
STYLES={1:dict(kind='line',color='black',linewidth=1.2,alpha=1.0),10:dict(kind='scatter',color='#2878b5',marker='o',size=28,alpha=0.4,hollow=True),100:dict(kind='scatter',color='#e07a1f',marker='^',size=30,alpha=0.4,hollow=True),1000:dict(kind='scatter',color='#4b9b58',marker='x',size=20,alpha=0.4,hollow=False)}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
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
    fig,(ax,residual)=plt.subplots(2,1,figsize=(13,8),sharex=True,gridspec_kw={'height_ratios':[3,1]})
    artists=[]
    for vmax,style in STYLES.items():
        selected=sorted((r for r in rows if float(r['vmax_km_s'])==vmax),key=lambda r:int(r['rank_by_vmax1_logBF']))
        x=[int(r['rank_by_vmax1_logBF']) for r in selected]
        values=([float(r['logBF'])/math.log(10) for r in selected],[float(r['delta_logc_vs_vmax1']) for r in selected])
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
    residual.set_xlabel('Galaxy rank (sorted once by Vmax=1 broad-limit log BF)')
    fig.suptitle(TITLE)
    fig.tight_layout()
    assert fig._suptitle.get_text()==TITLE
    assert len(ax.collections)==3 and len(residual.collections)==3
    assert all(a.get_alpha()==0.4 for a in ax.collections+residual.collections)
    OUT.mkdir(parents=True,exist_ok=False)
    png=OUT/'vmax_scale_check_styled.png';pdf=OUT/'vmax_scale_check_styled.pdf'
    fig.savefig(png,dpi=180);fig.savefig(pdf);plt.close(fig)
    manifest=dict(artifact='rendering-only restyle',created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),title=TITLE,source_csv=str(SRC/'vmax_scale_check.csv'),source_csv_sha256=sha(SRC/'vmax_scale_check.csv'),source_manifest=str(SRC/'manifest.json'),source_manifest_sha256=sha(SRC/'manifest.json'),renderer=str(Path(__file__).resolve()),renderer_sha256=sha(__file__),style_settings=STYLES,rows=700,galaxies=175,numbers_changed=False,scientific_evaluations=0,diagnostic_scope=original['diagnostic'],claim_limit=original['not_a_claim'],figure_inches=[13,8],png_dpi=180,residual_panel_retained=True,files={p.name:sha(p) for p in [png,pdf]})
    with (OUT/'rendering_manifest.json').open('x') as f:json.dump(manifest,f,indent=2);f.write(chr(10))
    print('RENDERED',png,pdf,'SOURCE_CSV_UNCHANGED',sha(SRC/'vmax_scale_check.csv')==original['files']['vmax_scale_check.csv'])
if __name__=='__main__':main()
