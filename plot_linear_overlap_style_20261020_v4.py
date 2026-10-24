#!/usr/bin/env python3
"""Readability-only re-render of saved linear selections; no fits or probabilities are recomputed."""
import csv, hashlib, json, sys
from collections import Counter
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
ROOT=Path('/mnt/users/koe/RAR_mono_grid')
BASE=ROOT/'outputs/saved_profile_comparison_20261019_exactcut_v1'
OUT=BASE/'convention_review_20261020_v4'
CSV=BASE/'linear_marginalised_rawP_leastJrad_selected_profiles_v7.csv'
AUDIT=BASE/'grid_figure_manifest_20261020_source_scoped_v7.json'
ORIGIN=BASE/'convention_review_20261020_v4'
PNG=OUT/'linear_leastJrad_points_overlap_v4.png'
PDF=OUT/'linear_leastJrad_points_overlap_v4.pdf'
MANIFEST=OUT/'manifest_linear_overlap_v4.json'
SOURCE=ROOT/'plot_linear_overlap_style_20261020_v3.py'
def sha256(path): return hashlib.sha256(path.read_bytes()).hexdigest()
for path in (PNG,PDF,MANIFEST):
    if path.exists(): raise FileExistsError(f'Refusing to overwrite existing output: {path}')
with CSV.open(newline='') as f: rows=list(csv.DictReader(f))
audit=json.loads(AUDIT.read_text())
if len(rows)!=77: raise ValueError(f'Expected 77 saved rows, found {len(rows)}')
sys.path.insert(0,str(ROOT))
from utils import load_sparc_table
table=load_sparc_table().copy()
if len(table)!=175: raise ValueError(f'Expected 175 SPARC galaxies, found {len(table)}')
table['mass_proxy']=0.5*table['L'].astype(float)
norm=Normalize(vmin=float(table.mass_proxy.min()),vmax=float(table.mass_proxy.max()))
cmap=plt.get_cmap('viridis')
xs=[float(r['mean_ml_parameter']) for r in rows]
ys=[float(r['disk_slope']) for r in rows]
colors=[float(table.loc[table.Galaxy.astype(str)==r['galaxy'],'mass_proxy'].iloc[0]) for r in rows]
counts=Counter(zip(xs,ys))
baseline=(0.5,0.0)
baseline_n=sum(abs(x-baseline[0])<=1e-12 and abs(y-baseline[1])<=1e-12 for x,y in zip(xs,ys))
if len(counts)!=18 or baseline_n!=58: raise ValueError(f'Coordinate audit changed: unique={len(counts)}, baseline multiplicity={baseline_n}')
fig,ax=plt.subplots(figsize=(12,8),layout='constrained')
sc=ax.scatter(xs,ys,c=colors,cmap=cmap,norm=norm,s=38,edgecolors='black',linewidths=0.35,alpha=0.75,zorder=3)
for (x,y),n in counts.items():
    if n>1:
        is_baseline=abs(x-baseline[0])<=1e-12 and abs(y)<=1e-12
        ring_size=80 if is_baseline else 50+18*n
        ax.scatter([x],[y],s=ring_size,facecolors='none',edgecolors='#b2182b',linewidths=1.5,zorder=4)
        dx,dy=(-55,-35) if is_baseline else (12,28)
        ax.annotate(f'{n} profiles overlap',xy=(x,y),xytext=(dx,dy),textcoords='offset points',fontsize=9,color='#8c1426',arrowprops={'arrowstyle':'-','color':'#8c1426','lw':0.8})
lo,hi=0.5*10**(-0.1),0.5*10**(0.1)
ax.axvspan(lo,hi,color='0.45',alpha=0.18,label='SPARC constant-prior band at pivot: ±0.1 dex')
ax.axvline(0.5,color='0.35',linewidth=1); ax.axhline(0,color='0.35',linewidth=1)
ax.set_xlabel('Mean Υdisk at each galaxy’s saved pivot B/A',fontsize=12)
ax.set_ylabel('Linear slope s [Υdisk kpc⁻¹]',fontsize=12)
ax.set_title('Saved least-Jrad profiles (n=77; 18 unique' \
' coordinates)',fontsize=14)
ax.text(0.02,0.98,'All 77 rows plotted at recorded coordinates; red rings mark duplicate coordinates; no jitter.',transform=ax.transAxes,va='top',fontsize=9)
ax.legend(fontsize=9); ax.grid(alpha=0.2)
fig.colorbar(sc,ax=ax,label='0.5 × catalogue total L (10⁹ L☉)')
fig.savefig(PNG,dpi=240); fig.savefig(PDF); plt.close(fig)
manifest={'purpose':'Final readability-only render. User-reviewed label wording omits exactly; the 58-profile callout is shallower and points down-left, closer to the angle of the other callouts. No profile values or science content changed.','origin':{'png':str(ORIGIN/'linear_leastJrad_points_overlap_v3.png'),'png_sha256':sha256(ORIGIN/'linear_leastJrad_points_overlap_v3.png'),'manifest':str(ORIGIN/'manifest_linear_overlap_v3.json'),'manifest_sha256':sha256(ORIGIN/'manifest_linear_overlap_v3.json'),'source_script':str(SOURCE),'source_script_sha256':sha256(SOURCE)},'data_inputs':{str(p):{'sha256':sha256(p),'bytes':p.stat().st_size} for p in (CSV,AUDIT)},'style_change':{'baseline_overlap_label_offset_points':[-55,-35],'other_duplicate_label_offset_points':[12,28],'label_text':'N profiles overlap','all_77_rows_at_original_coordinates':True,'coordinate_jitter':False},'audit':{'saved_rows':len(rows),'unique_exact_coordinates':len(counts),'baseline_coordinate':[0.5,0.0],'baseline_multiplicity':baseline_n,'duplicate_coordinate_groups':{f'{x:.16g},{y:.16g}':n for (x,y),n in counts.items() if n>1}},'outputs':[{'path':str(p),'bytes':p.stat().st_size,'sha256':sha256(p)} for p in (PNG,PDF)]}
MANIFEST.write_text(json.dumps(manifest,indent=2))
print(json.dumps({'output_directory':str(OUT),'png':str(PNG),'pdf':str(PDF),'manifest':str(MANIFEST),'rows':len(rows),'unique_coordinates':len(counts),'baseline_multiplicity':baseline_n,'png_sha256':sha256(PNG),'pdf_sha256':sha256(PDF)},indent=2))
