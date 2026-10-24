#!/usr/bin/env python3
"""Isolated four-normalization comparison; one raw Casertano solve, unit disc M/L.
No force-source edits, Vobs fit, empirical normalization or inclination adjustment.
Published bulge luminosity/distance are explicit supplied parameters with source.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import numpy as np
from utils import DEFAULT_DENS_ROOT, DEFAULT_TABLE_PATH, compute_A_B, load_sparc_table, vel2acc
from sparc_scale_proxies import sha256_file
ROOT = Path(__file__).resolve().parent


def factors(disk_l, bulge_l, total_l, published_bulge_l):
    if not all(math.isfinite(x) for x in (disk_l,bulge_l,total_l,published_bulge_l)):
        raise ValueError("nonfinite luminosity")
    if disk_l <= 0 or bulge_l < 0 or total_l <= 0 or not 0 <= published_bulge_l < total_l:
        raise ValueError("invalid component luminosities")
    return {"raw_absolute":1., "disk_to_catalogue_total":total_l/disk_l,
            "joint_disk_bulge_light":total_l/(disk_l+bulge_l),
            "disk_to_total_minus_published_bulge":(total_l-published_bulge_l)/disk_l}


def summarize(raw_w, r, reference_v, scales):
    raw_w,r,reference_v = (np.asarray(x,float) for x in (raw_w,r,reference_v))
    if raw_w.shape != r.shape or reference_v.shape != r.shape or np.any(r <= 0):
        raise ValueError("invalid comparison geometry")
    if not all(np.all(np.isfinite(x)) for x in (raw_w,r,reference_v)):
        raise ValueError("nonfinite comparison")
    target = vel2acc(reference_v*np.abs(reference_v),r)
    peak = float(np.max(np.abs(target)))
    if peak <= 0: raise ValueError("zero disc reference")
    nonzero = target != 0
    summaries,radial = [],[]
    for label,scale in scales.items():
        w = raw_w*scale; acceleration = vel2acc(w,r);delta = acceleration-target
        ratio = acceleration[nonzero]/target[nonzero]
        corr = (float(np.corrcoef(acceleration,target)[0,1])
                if np.std(acceleration)>0 and np.std(target)>0 else None)
        summaries.append({"case":label,"density_multiplier":scale,
                          "acceleration_rmse_m_s2":float(np.sqrt(np.mean(delta*delta))),
                          "acceleration_nrmse_to_peak":float(np.sqrt(np.mean(delta*delta))/peak),
                          "ratio_median":float(np.median(ratio)),
                          "ratio_q16":float(np.quantile(ratio,.16)),
                          "ratio_q84":float(np.quantile(ratio,.84)),
                          "acceleration_shape_correlation":corr,
                          "zero_target_ratio_exclusions":int(np.count_nonzero(~nonzero))})
        for i in range(len(r)):
            radial.append({"case":label,"row_index":i,"W_disc_kms2":float(w[i]),
                           "signed_acceleration_m_s2":float(acceleration[i]),
                           "residual_acceleration_m_s2":float(delta[i]),
                           "acceleration_ratio":float(acceleration[i]/target[i]) if nonzero[i] else None})
    return summaries,radial


def write_csv(path,rows):
    with path.open("x",newline="") as stream:
        w=csv.DictWriter(stream,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def self_test():
    scales=factors(2.,1.,6.,1.5)
    assert scales==dict(raw_absolute=1.,disk_to_catalogue_total=3.,joint_disk_bulge_light=2.,disk_to_total_minus_published_bulge=2.25)
    summaries,rows=summarize(np.array([4.,-9.,16.]),np.array([1.,2.,3.]),np.array([2.,-3.,4.]),scales)
    assert summaries[0]["acceleration_nrmse_to_peak"]==0
    assert math.isclose(summaries[1]["ratio_median"],3) and rows[1]["W_disc_kms2"]==-9
    try: factors(2.,1.,6.,7.)
    except ValueError: pass
    else: raise AssertionError("negative inferred disc luminosity accepted")
    print("Normalization/sign/scaling checks passed; no Casertano solve.")


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--self-test",action="store_true")
    p.add_argument("--galaxy")
    p.add_argument("--published-bulge-luminosity",type=float,help="1e9 Lsun from named source")
    p.add_argument("--published-distance",type=float,help="Mpc; must match catalogue and rotmod")
    p.add_argument("--bulge-source",help="Published table URL and/or version identifier")
    p.add_argument("--dens-root",type=Path,default=Path(DEFAULT_DENS_ROOT))
    p.add_argument("--data-root",type=Path,default=Path('/mnt/users/koe/data'))
    p.add_argument("--table",type=Path,default=Path(DEFAULT_TABLE_PATH))
    p.add_argument("--output-root",type=Path)
    p.add_argument("--dry-run",action="store_true")
    p.add_argument("--plot",action="store_true")
    args=p.parse_args()
    if args.self_test: self_test();return
    if not args.galaxy or args.published_bulge_luminosity is None or args.published_distance is None or not args.bulge_source:
        p.error("explicit galaxy, published luminosity/distance and source required")
    if not args.dry_run and args.output_root is None: p.error("fresh --output-root required")
    if args.output_root is not None and args.output_root.exists(): raise FileExistsError("no overwrite")
    table=load_sparc_table(args.table);cat=table.loc[table.Galaxy==args.galaxy]
    if len(cat)!=1: raise ValueError("galaxy absent/nonunique")
    cat=cat.iloc[0];dens=args.dens_root/(args.galaxy+'.dens');rotmod=args.data_root/(args.galaxy+'_rotmod.dat')
    density=np.loadtxt(dens,skiprows=1,ndmin=2);kin=np.loadtxt(rotmod,comments='#',ndmin=2)
    if density.shape[1]!=3 or kin.shape[1]<6 or len(density)<2:
        raise ValueError("invalid density/kinematic schema")
    if not np.all(np.isfinite(density)) or not np.all(np.isfinite(kin[:,:6])):
        raise ValueError("nonfinite inputs; no dropped bad rows")
    rad,sd,sb=density.T;sd=sd*1e6;sb=sb*1e6
    if np.any(rad<0) or np.any(np.diff(rad)<=0) or np.any(sd<0) or np.any(sb<0):
        raise ValueError("invalid photometric inputs")
    if np.any(kin[:,0]<=0) or np.any(np.diff(kin[:,0])<=0): raise ValueError("invalid kinematic radii")
    header=rotmod.read_text().splitlines()[0]
    distance=float(header.split('=')[1].split()[0])
    if not math.isclose(distance,float(cat.D),abs_tol=.005) or not math.isclose(args.published_distance,float(cat.D),abs_tol=.005):
        raise ValueError("distance mismatch; no automatic luminosity rescaling")
    positive=rad[rad>0]
    if len(positive)<2: raise ValueError("insufficient positive force support")
    # Restrict comparisons to solver support to prevent endpoint clamping.
    keep=(kin[:,0]>=positive[0])&(kin[:,0]<=positive[-1]);k=kin[keep]
    if len(k)<2: raise ValueError("insufficient common comparison support")
    ad=compute_A_B(rad,sd)[0]/1e9;ab=compute_A_B(rad,sb)[0]/1e9
    scales=factors(ad,ab,float(cat.L),args.published_bulge_luminosity)
    files={'dens':dens,'rotmod':rotmod,'table':args.table}
    sources={name:ROOT/name for name in ('mono_SPARC_RAR.py','utils.py','sparc_scale_proxies.py',Path(__file__).name)}
    hashes={key:sha256_file(path) for key,path in {**files,**sources}.items()}
    info={'galaxy':args.galaxy,'disc_ML':1.,'bulge_reference_ML':1.,'distance_Mpc':float(cat.D),
          'inclination_deg':float(cat.Inc),'zd_kpc':.196*float(cat.Rdisk)**.633,
          'raw_disk_luminosity_1e9':ad,'raw_bulge_light_integral_1e9':ab,
          'catalogue_total_luminosity_1e9':float(cat.L),'published_bulge_luminosity_1e9':args.published_bulge_luminosity,
          'published_source':args.bulge_source,'published_value_policy':'explicit parameter copied from cited table; full public table not archived by this script',
          'bulge_distance_assumption':'explicit supplied SPARC reference distance; bulge table itself has no distance column',
          'density_header':dens.read_text().splitlines()[0],'density_units':'Lsun/pc2 converted to Lsun/kpc2',
          'inclination_calibration':'UNRESOLVED; no cos(i) correction inferred or applied',
          'normalization_factors':scales,'force_support':[float(positive[0]),float(positive[-1])],
          'zero_radius_nodes_removed_by_common_solver':int(np.count_nonzero(rad==0)),
          'kinematic_rows_excluded':int(np.count_nonzero(~keep)),'compared_rows':len(k),
          'input_source_sha256':hashes,'numerical_policy':'one raw solve; signed W and signed-force interpolation; identical height/radii/support in every case',
          'SPARC_outer_extrapolation':'not recreated; finite .dens support is a separate limitation',
          'global_scale_comparison_not_physical_model_selection':True}
    if args.dry_run:
        print(json.dumps({'status':'prepared_no_force_solve','output':str(args.output_root),'configuration':info},indent=2,allow_nan=False));return
    import mono_SPARC_RAR as mono
    import inspect
    vcdisk_path=Path(inspect.getfile(mono.vcdisk));vcdisk_hash=sha256_file(vcdisk_path)
    info.update(vcdisk_file=str(vcdisk_path),vcdisk_sha256=vcdisk_hash,numpy_version=np.__version__)
    raw_w,warnings,error=mono.solve_disk_basis(rad,sd,info['zd_kpc'],k[:,0])
    if error or raw_w is None: raise RuntimeError('raw solver failed: '+error)
    if sha256_file(vcdisk_path)!=vcdisk_hash: raise RuntimeError('vcdisk source changed')
    summaries,rows=summarize(raw_w,k[:,0],k[:,4],scales)
    info.update(solver_warnings=warnings,raw_Casertano_solve_count=1)
    if hashes!={key:sha256_file(path) for key,path in {**files,**sources}.items()}:
        raise RuntimeError('source/input changed during comparison')
    base=[{'row_index':i,'radius_kpc':float(row[0]),'SPARC_Vdisk_unitML_kms':float(row[4]),
           'SPARC_Vbul_unitML_kms':float(row[5]),'SPARC_Vgas_kms':float(row[3]),
           'raw_Casertano_W_kms2':float(raw_w[i])} for i,row in enumerate(k)]
    args.output_root.mkdir(parents=True)
    write_csv(args.output_root/'reference_radial.csv',base)
    write_csv(args.output_root/'case_radial.csv',rows)
    write_csv(args.output_root/'summary.csv',summaries)
    with (args.output_root/'manifest.json').open('x') as stream:
        json.dump(info,stream,indent=2,allow_nan=False);stream.write('\n')
    if args.plot:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig,axes=plt.subplots(1,2,figsize=(12,4),layout='constrained')
        axes[0].plot(k[:,0],k[:,4],color='black',lw=2,label='SPARC disc, ML=1')
        axes[1].plot(k[:,0],vel2acc(k[:,4]*np.abs(k[:,4]),k[:,0]),color='black',lw=2,label='SPARC disc')
        for label,scale in scales.items():
            w=raw_w*scale
            axes[0].plot(k[:,0],np.sign(w)*np.sqrt(np.abs(w)),alpha=.5,label=label)
            axes[1].plot(k[:,0],vel2acc(w,k[:,0]),alpha=.5,label=label)
        axes[0].set_ylabel('Signed disc velocity [km/s]');axes[1].set_ylabel('Signed disc acceleration [m/s2]')
        for ax in axes: ax.set_xlabel('Radius [kpc]');ax.legend(fontsize=7);ax.grid(alpha=.15)
        fig.suptitle(args.galaxy+': normalization comparison; inclination calibration unresolved')
        fig.savefig(args.output_root/'disc_normalizations.png',dpi=150);plt.close(fig)
    print('ONE_RAW_SOLVE',args.galaxy,'cases',len(scales),'output',args.output_root)


if __name__=='__main__':
    main()
