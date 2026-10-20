#!/usr/bin/env python3
"""Render saved SPARC M/L selections for review; performs no new fits or probabilities."""
import csv, hashlib, json, math, sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
ROOT=Path('/mnt/users/koe/RAR_mono_grid')
BASE=ROOT/'outputs/saved_profile_comparison_20261019_exactcut_v1'
OUT=BASE/'convention_review_20261020_v1'
CONSTANT_CSV=BASE/'constant_rawP_accepted_grid_runs_source_scoped_v7.csv'
LINEAR_CSV=BASE/'linear_marginalised_rawP_leastJrad_selected_profiles_v7.csv'
AUDIT_JSON=BASE/'grid_figure_manifest_20261020_source_scoped_v7.json'
sys.path.insert(0,str(ROOT))
from utils import load_sparc_table

def read_csv(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle))
def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

if OUT.exists():
    raise FileExistsError(f'Refusing to overwrite existing output: {OUT}')
constant=read_csv(CONSTANT_CSV)
linear=read_csv(LINEAR_CSV)
audit=json.loads(AUDIT_JSON.read_text())
per_galaxy=audit['constant_cell_audit_by_galaxy']
table=load_sparc_table().copy()
if len(table)!=175 or len(constant)!=67 or len(linear)!=77 or len(per_galaxy)!=175:
    raise ValueError('Saved-data coverage differs from expected v7 inputs')
if any(int(x.get('invalid',0)) or int(x.get('nonfinite_log_probability',0)) or int(x.get('malformed_log_probability',0)) for x in per_galaxy.values()):
    raise ValueError('Constant audit contains failed/malformed cells; inspect before plotting')
if set(table['Galaxy'].astype(str))!=set(per_galaxy):
    raise ValueError('Audit galaxy set does not match SPARC catalogue')
# The agreed ordering/colour proxy is 0.5 times catalogue total luminosity.
table['mass_proxy']=0.5*table['L'].astype(float)
table=table.sort_values(['mass_proxy','Galaxy'],ascending=[True,True]).reset_index(drop=True)
mass_by_galaxy={str(r.Galaxy):float(r.mass_proxy) for r in table.itertuples()}
rows_by_galaxy={}
for row in constant:
    rows_by_galaxy.setdefault(row['galaxy'],[]).append(row)
qualifies={g:int(v['qualifying'])>0 for g,v in per_galaxy.items()}
if {g for g,ok in qualifies.items() if ok}!=set(rows_by_galaxy):
    raise ValueError('Qualifying-cell audit and accepted-interval CSV disagree')
norm=Normalize(vmin=float(table.mass_proxy.min()),vmax=float(table.mass_proxy.max()))
cmap=plt.get_cmap('viridis')
fig,ax=plt.subplots(figsize=(13,39),layout='constrained')
for y,row in enumerate(table.itertuples()):
    galaxy=str(row.Galaxy)
    for interval in rows_by_galaxy.get(galaxy,[]):
        ax.hlines(y,float(interval['log10_ml_min_sample']),float(interval['log10_ml_max_sample']),color=cmap(norm(float(row.mass_proxy))),linewidth=2.2)
center=math.log10(0.5)
ax.axvspan(center-0.1,center+0.1,color='0.45',alpha=0.18,label='SPARC prior band: log10(0.5) ± 0.1 dex')
ax.axvline(center,color='0.35',linewidth=0.8)
ax.set_yticks(range(len(table)))
ax.set_yticklabels([str(r.Galaxy) for r in table.itertuples()],fontsize=6,rotation=45,ha='right',rotation_mode='anchor')
for label,row in zip(ax.get_yticklabels(),table.itertuples()):
    label.set_fontweight('bold' if not qualifies[str(row.Galaxy)] else 'normal')
ax.set_ylim(len(table)-0.5,-0.5)
ax.set_xlabel('log10 Υdisk (solar units)')
ax.set_ylabel('SPARC galaxies ordered by 0.5 × catalogue total L (ascending)')
ax.set_title('Constant profiles: raw P ≥ Φ(−3); 67 qualifying intervals, 108 no-qualifier galaxies')
sm=ScalarMappable(norm=norm,cmap=cmap); sm.set_array([])
fig.colorbar(sm,ax=ax,label='0.5 × catalogue total L (10⁹ L☉)')
ax.legend(loc='lower right',fontsize=8); ax.grid(axis='x',alpha=0.2)
OUT.mkdir(parents=True)
constant_png=OUT/'constant_mass_ordered_intervals_v1.png'
constant_pdf=OUT/'constant_mass_ordered_intervals_v1.pdf'
fig.savefig(constant_png,dpi=180); fig.savefig(constant_pdf); plt.close(fig)
# One saved minimum-J_rad point per galaxy; its mean is at that galaxy's own B/A pivot.
fig,ax=plt.subplots(figsize=(10,7),layout='constrained')
xs=[float(r['mean_ml_parameter']) for r in linear]
ys=[float(r['disk_slope']) for r in linear]
colors=[mass_by_galaxy[r['galaxy']] for r in linear]
sc=ax.scatter(xs,ys,c=colors,cmap=cmap,norm=norm,s=27,edgecolors='black',linewidths=0.3)
lo,hi=0.5*10**(-0.1),0.5*10**(0.1)
ax.axvspan(lo,hi,color='0.45',alpha=0.18,label='SPARC constant-prior band at pivot: ±0.1 dex')
ax.axvline(0.5,color='0.35',linewidth=0.8); ax.axhline(0,color='0.35',linewidth=0.8)
ax.set_xlabel('Mean Υdisk parameter (value at each galaxy’s saved pivot B/A)')
ax.set_ylabel('Linear slope s [Υdisk kpc⁻¹]')
ax.set_title('Linear profiles: saved least-J_rad selection (n=77); SPARC-only, pivots vary, external K_s overlay withheld')
ax.legend(fontsize=8); ax.grid(alpha=0.2)
fig.colorbar(sc,ax=ax,label='0.5 × catalogue total L (10⁹ L☉)')
linear_png=OUT/'linear_leastJrad_points_v1.png'
linear_pdf=OUT/'linear_leastJrad_points_v1.pdf'
fig.savefig(linear_png,dpi=180); fig.savefig(linear_pdf); plt.close(fig)
statuses=[{'galaxy':str(r.Galaxy),'status':'qualifying_interval' if qualifies[str(r.Galaxy)] else 'no_qualifying_cell','qualifying_cells':int(per_galaxy[str(r.Galaxy)]['qualifying']),'marginalised_cells':int(per_galaxy[str(r.Galaxy)]['marginalised_rows']),'invalid_cells':int(per_galaxy[str(r.Galaxy)]['invalid'])} for r in table.itertuples()]
outputs=[constant_png,constant_pdf,linear_png,linear_pdf]
manifest={'purpose':'visual review of existing saved selections only; no fits or probabilities recomputed','source_version':'source-scoped v7; preserved unchanged','inputs':{str(p):{'sha256':sha256(p),'bytes':p.stat().st_size} for p in (CONSTANT_CSV,LINEAR_CSV,AUDIT_JSON)},'constant':{'catalogue_rows':len(table),'qualifying_galaxies':sum(qualifies.values()),'no_qualifying_galaxies':sum(not x for x in qualifies.values()),'audit_invalid_nonfinite_malformed_cells':0,'ordering':'ascending 0.5 × SPARC catalogue total luminosity L','colour':'viridis by 0.5 × catalogue total L','threshold':'raw P >= exact Phi(-3); no sigma calibration implied','prior_band':'log10(0.5) +/- 0.1 dex','label_rotation_degrees':45,'bold_labels':'only galaxies with zero qualifying cells','status_by_galaxy':statuses},'linear':{'selected_points':len(linear),'selection':'existing saved minimum-J_rad qualifying profile per galaxy','axes':'saved mean parameter at each galaxy pivot B/A versus disk_slope in kpc^-1','external_overlay':'withheld; no verified shared coefficient/pivot mapping or band/component compatibility'},'external_radius_note':'No R_bar/R_d rescaling applied; R_bar is not assumed to equal R_d.','figure_files':[{'path':str(p),'bytes':p.stat().st_size,'sha256':sha256(p)} for p in outputs]}
manifest_path=OUT/'manifest_convention_review_v1.json'
manifest_path.write_text(json.dumps(manifest,indent=2))
print(json.dumps({'output_directory':str(OUT),'files':[p.name for p in outputs]+[manifest_path.name],'constant_galaxies':len(table),'constant_qualifying':sum(qualifies.values()),'constant_no_qualifier':sum(not x for x in qualifies.values()),'linear_points':len(linear)},indent=2))
