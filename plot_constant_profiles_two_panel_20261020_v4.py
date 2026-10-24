
#!/usr/bin/env python3
"""Create a two-panel readability version of saved constant-profile selections; no new analysis."""
import csv, hashlib, json, math, sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
ROOT=Path('/mnt/users/koe/RAR_mono_grid')
BASE=ROOT/'outputs/saved_profile_comparison_20261019_exactcut_v1'
OUT=BASE/'convention_review_20261020_v4'
CONSTANT_CSV=BASE/'constant_rawP_accepted_grid_runs_source_scoped_v7.csv'
AUDIT_JSON=BASE/'grid_figure_manifest_20261020_source_scoped_v7.json'
sys.path.insert(0,str(ROOT))
from utils import load_sparc_table

def read_csv(path):
    with path.open(newline='') as f:
        return list(csv.DictReader(f))
def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()
if OUT.exists():
    raise FileExistsError(f'Refusing to overwrite existing output: {OUT}')
constant=read_csv(CONSTANT_CSV)
audit=json.loads(AUDIT_JSON.read_text())
per_galaxy=audit['constant_cell_audit_by_galaxy']
table=load_sparc_table().copy()
if len(table)!=175 or len(per_galaxy)!=175:
    raise ValueError('Expected exactly 175 catalogue/audit galaxies')
if any(int(v.get('invalid',0)) or int(v.get('nonfinite_log_probability',0)) or int(v.get('malformed_log_probability',0)) for v in per_galaxy.values()):
    raise ValueError('Constant audit contains invalid or malformed cells')
if set(table['Galaxy'].astype(str))!=set(per_galaxy):
    raise ValueError('Audit galaxy set differs from catalogue')
# L is in 10^9 Lsun; fiducial stellar M/L is 0.5 Msun/Lsun.
table['stellar_mass_1e9_msun']=0.5*table['L'].astype(float)
table=table.sort_values(['stellar_mass_1e9_msun','Galaxy'],ascending=[True,True]).reset_index(drop=True)
qualifies={g:int(v['qualifying'])>0 for g,v in per_galaxy.items()}
rows_by_galaxy={}
for row in constant:
    rows_by_galaxy.setdefault(row['galaxy'],[]).append(row)
if {g for g,ok in qualifies.items() if ok}!=set(rows_by_galaxy):
    raise ValueError('Accepted interval CSV and qualifying-cell audit disagree')
# Split contiguous ranks in the single global ascending-mass order.
subsets=[table.iloc[:88],table.iloc[88:]]
if [len(s) for s in subsets] not in ([88,87],[87,88]):
    raise ValueError('Unexpected two-panel split')
all_min=min(float(r['log10_ml_min_sample']) for r in constant)
all_max=max(float(r['log10_ml_max_sample']) for r in constant)
center=math.log10(0.5)
xlo=min(all_min,center-0.1); xhi=max(all_max,center+0.1)
pad=0.025*(xhi-xlo)
xlim=(xlo-pad,xhi+pad)
norm=Normalize(vmin=float(table.stellar_mass_1e9_msun.min()),vmax=float(table.stellar_mass_1e9_msun.max()))
cmap=plt.get_cmap('viridis')
OUT.mkdir(parents=True)
fig,axes=plt.subplots(1,2,figsize=(18,17),sharex=True,layout='constrained',gridspec_kw={'wspace':0.04})
for panel,subset in zip(axes,subsets):
    for y,row in enumerate(subset.itertuples()):
        galaxy=str(row.Galaxy)
        color=cmap(norm(float(row.stellar_mass_1e9_msun)))
        for interval in rows_by_galaxy.get(galaxy,[]):
            panel.hlines(y,float(interval['log10_ml_min_sample']),float(interval['log10_ml_max_sample']),color=color,linewidth=2.2)
    panel.axvspan(center-0.1,center+0.1,color='0.45',alpha=0.18,zorder=0)
    panel.axvline(center,color='0.35',linewidth=0.9,zorder=1)
    panel.set_xlim(*xlim)
    panel.set_ylim(len(subset)-0.5,-0.5)
    panel.set_yticks(range(len(subset)))
    panel.set_yticklabels([str(r.Galaxy) for r in subset.itertuples()],fontsize=8.5,rotation=45,ha='right',rotation_mode='anchor')
    for label,row in zip(panel.get_yticklabels(),subset.itertuples()):
        label.set_fontweight('bold' if not qualifies[str(row.Galaxy)] else 'normal')
    panel.grid(axis='x',alpha=0.18)
    panel.tick_params(axis='x',labelsize=10)
bs=chr(92)
fig.suptitle('Constant profiles: saved raw $P '+bs+'geq '+bs+'Phi(-3)$ selections',fontsize=22)
x_label='$'+bs+'log_{10} '+bs+'Upsilon_{'+bs+'rm disk} (M_{'+bs+'odot}/L_{'+bs+'odot})$'
fig.supxlabel(x_label,fontsize=13)
fig.supylabel('SPARC galaxies ordered by fiducial stellar mass (ascending)',fontsize=13)
sm=ScalarMappable(norm=norm,cmap=cmap); sm.set_array([])
cb=fig.colorbar(sm,ax=axes,location='right',pad=0.015,fraction=0.035)
cb_label='Stellar mass ($10^9'+bs+','+'M_{'+bs+'odot}$; SPARC fiducial $M/L$)'
cb.set_label(cb_label,fontsize=12)
cb.ax.tick_params(labelsize=10)
png=OUT/'constant_mass_ordered_intervals_two_panel_v4.png'
pdf=OUT/'constant_mass_ordered_intervals_two_panel_v4.pdf'
fig.savefig(png,dpi=240); fig.savefig(pdf); plt.close(fig)
ordered_names=[str(x) for s in subsets for x in s['Galaxy']]
if len(ordered_names)!=175 or len(set(ordered_names))!=175:
    raise ValueError('Plot split does not contain each of 175 galaxies exactly once')
outputs=[png,pdf]
manifest={'purpose':'format/readability revision of existing source-scoped v7 constant-profile selections only; no fit/probability recomputation','inputs':{str(p):{'sha256':sha256(p),'bytes':p.stat().st_size} for p in (CONSTANT_CSV,AUDIT_JSON)},'figure':{'panels':2,'galaxies_per_panel':[len(s) for s in subsets],'galaxies_total':175,'unique_galaxies':len(set(ordered_names)),'global_order':'ascending fiducial stellar mass = 0.5 x SPARC catalogue total luminosity','panel_1_range':[str(subsets[0].iloc[0].Galaxy),str(subsets[0].iloc[-1].Galaxy)],'panel_2_range':[str(subsets[1].iloc[0].Galaxy),str(subsets[1].iloc[-1].Galaxy)],'qualifying_galaxies':sum(qualifies.values()),'no_qualifier_galaxies':sum(not v for v in qualifies.values()),'bold':'galaxies with no qualifying cell','raw_probability_threshold':'raw P >= Phi(-3)','shared_x_limits':list(xlim),'shared_colour_normalization':'viridis, linear physical fiducial stellar mass in units of 10^9 Msun','colourbar_label':r'Stellar mass ($10^9\,M_{\odot}$; SPARC fiducial M/L)','shared_x_label':r'$\log_{10}\,\Upsilon_{\rm disk}\;(M_{\odot}/L_{\odot})$','shared_y_label':'SPARC galaxies ordered by fiducial stellar mass (ascending)','y_tick_rotation_degrees':45},'galaxies_in_plot_order':ordered_names,'outputs':[{'path':str(p),'bytes':p.stat().st_size,'sha256':sha256(p)} for p in outputs]}
origin_figure=BASE/'convention_review_20261020_v3/constant_mass_ordered_intervals_two_panel_v3.png'
origin_manifest=BASE/'convention_review_20261020_v3/manifest_two_panel_v3.json'
manifest['origin']={'figure':str(origin_figure),'figure_sha256':sha256(origin_figure),'manifest':str(origin_manifest),'manifest_sha256':sha256(origin_manifest)}
manifest['style_change']='Title font increased from 16 to 22 pt only; figure labels, data, order, axes, and layout unchanged.'
manifest['figure']['shared_title_fontsize_pt']=22
manifest_path=OUT/'manifest_two_panel_v4.json'
manifest_path.write_text(json.dumps(manifest,indent=2))
print(json.dumps({'output_directory':str(OUT),'files':[p.name for p in outputs]+[manifest_path.name],'panels':[len(s) for s in subsets],'galaxies':len(ordered_names),'qualifying':sum(qualifies.values()),'no_qualifier':sum(not v for v in qualifies.values()),'shared_x_limits':xlim,'output_hashes':[sha256(p) for p in outputs]},indent=2))