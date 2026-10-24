#!/usr/bin/env python3
"""Plot saved controlled mocks by deterministic replay, never rerun probabilities."""
import argparse
import csv
import json
import math
from pathlib import Path
import numpy as np
import scipy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import mock_fixed_profile_N as producer
from plot_gp_rar_search import observed_interval, Q
from sparc_scale_proxies import sha256_file
from utils import vel2acc
ROOT = Path(__file__).resolve().parent


def reconstruct(manifest):
    settings=manifest['settings'];design=manifest['design']
    counts=sorted(set(settings['counts']));size=max(counts)
    arrays=producer.latent_design(size,settings['fractional_velocity_error'])
    keys=('radius_kpc','disk_g_m_s2','gas_unit_g_m_s2','bulge_unit_g_m_s2',
          'latent_velocity_kms','error_velocity_kms')
    for key,array in zip(keys,arrays):
        if not np.array_equal(array,np.asarray(design[key],float)):
            raise ValueError('replayed parent design differs: '+key)
    order=producer.nested_indices(size)
    for n in counts:
        if np.sort(order[:n]).tolist()!=design['subset_indices'][str(n)]:
            raise ValueError('saved nested subset differs')
    r,d,g,b,vtrue,errors=arrays
    gbar=d+g+.7*b;latent=vel2acc(vtrue*vtrue,r)
    if np.any(gbar<=0) or np.any(latent<=0) or not np.all(np.diff(latent[np.argsort(gbar)])>0):
        raise ValueError('latent RAR is not strictly increasing')
    replay=[]
    for rep in range(settings['replicates']):
        rng=np.random.default_rng(np.random.SeedSequence([settings['seed'],rep]))
        observed=vtrue+rng.normal(size=size)*errors
        interval=observed_interval(observed,errors,r)
        if interval.shape!=(2,size) or not np.all(np.isfinite(interval)) or np.any(interval<0):
            raise ValueError('invalid squared-Gaussian intervals')
        replay.append({'replicate':rep,'rng':rng.bit_generator.__class__.__name__,
                       'observed_velocity_kms':observed.tolist(),
                       'gobs_central_m_s2':vel2acc(observed*observed,r).tolist(),
                       'gobs_interval_m_s2':interval.tolist()})
    return gbar,latent,replay


def plot_saved(input_root,output_root):
    if output_root.exists(): raise FileExistsError('fresh output root required')
    paths={'manifest':input_root/'manifest.json','rows':input_root/'mock_rows.csv',
           'plot_script':Path(__file__).resolve(),'interval_helper':ROOT/'plot_gp_rar_search.py'}
    hashes={key:sha256_file(path) for key,path in paths.items()}
    manifest=json.loads(paths['manifest'].read_text())
    if manifest.get('scope')!='controlled_fixed_profile_N_dependence_only' or manifest.get('failures')!=0:
        raise ValueError('requires completed controlled null without failed rows')
    for name,expected in manifest['source_sha256'].items():
        if sha256_file(ROOT/name)!=expected: raise ValueError('producer dependency changed: '+name)
    rows=list(csv.DictReader(paths['rows'].open()))
    settings=manifest['settings'];counts=sorted(set(settings['counts']))
    if len(counts)!=4: raise ValueError('this layout requires four saved N values')
    expected={(rep,n) for rep in range(settings['replicates']) for n in counts}
    keyed={}
    for row in rows:
        key=(int(row['replicate']),int(row['bins']))
        if key in keyed or row['status']!='ok' or not math.isfinite(float(row['logp'])) or float(row['logp'])>0:
            raise ValueError('duplicate/failed/nonfinite saved row')
        keyed[key]=row
    if set(keyed)!=expected: raise ValueError('incomplete saved coverage')
    gbar,latent,replay=reconstruct(manifest)
    positives=[latent]
    for item in replay:
        values=np.r_[item['gobs_central_m_s2'],np.asarray(item['gobs_interval_m_s2']).ravel()]
        positives.append(values[values>0])
    values=np.concatenate(positives)
    xlim=(float(np.min(gbar)*.8),float(np.max(gbar)*1.25))
    ylim=(float(np.min(values)*.8),float(np.max(values)*1.25))
    output_root.mkdir(parents=True)
    for item in replay:
        rep=item['replicate'];central=np.asarray(item['gobs_central_m_s2'])
        low,high=np.asarray(item['gobs_interval_m_s2'])
        fig,axes=plt.subplots(2,2,figsize=(10,8),sharex=True,sharey=True)
        panels=[]
        for ax,n in zip(axes.ravel(),counts):
            ix=np.asarray(manifest['design']['subset_indices'][str(n)],int)
            row=keyed[(rep,n)];lp=float(row['logp'])
            ax.set_xscale('log');ax.set_yscale('log');ax.set_xlim(xlim);ax.set_ylim(ylim)
            order=np.argsort(gbar)
            ax.plot(gbar[order],latent[order],color='#0077bb',lw=1.7)
            ax.scatter(gbar[ix],latent[ix],s=12,color='#0077bb',alpha=.25)
            # Draw intervals independently: the noisy central square need not lie
            # inside an equal-tailed squared-Gaussian interval at low SNR.
            shown_low=np.maximum(low[ix],ylim[0]);shown_high=np.maximum(high[ix],ylim[0])
            ax.vlines(gbar[ix],shown_low,shown_high,color='#ee7733',alpha=.6,lw=.9)
            ax.scatter(gbar[ix],shown_low,marker='_',s=18,color='#ee7733',alpha=.6)
            ax.scatter(gbar[ix],shown_high,marker='_',s=18,color='#ee7733',alpha=.6)
            positive=central[ix]>0
            ax.plot(gbar[ix][positive],central[ix][positive],'o-',color='#ee7733',ms=3,lw=.6,alpha=.8)
            zero=~positive
            if np.any(zero): ax.scatter(gbar[ix][zero],np.full(np.count_nonzero(zero),ylim[0]),marker='v',color='#ee7733')
            zero_lower=low[ix]==0
            if np.any(zero_lower): ax.scatter(gbar[ix][zero_lower],np.full(np.count_nonzero(zero_lower),ylim[0]),marker='v',color='#777777')
            ax.set_title(f'N={n}: saved P={math.exp(lp):.6g}, lnP={lp:.5g}',fontsize=10)
            ax.grid(alpha=.18)
            panels.append({'bins':n,'saved_row':row,'indices':ix.tolist(),'gbar_reference_m_s2':gbar[ix].tolist(),
                           'latent_gobs_m_s2':latent[ix].tolist(),'noisy_gobs_central_m_s2':central[ix].tolist(),
                           'conditional_interval_lower_m_s2':low[ix].tolist(),
                           'conditional_interval_upper_m_s2':high[ix].tolist(),'zero_lower_indices':ix[zero_lower].tolist(),
                           'zero_central_indices':ix[zero].tolist()})
        for ax in axes[:,0]: ax.set_ylabel('g_obs [m/s²]')
        for ax in axes[-1,:]: ax.set_xlabel('g_bar,ref [m/s²]')
        handles=[Line2D([0],[0],color='#0077bb',label='Strictly monotonic latent RAR'),
                 Line2D([0],[0],color='#ee7733',marker='o',lw=.6,label='Noisy central values (connected for display)'),
                 Line2D([0],[0],color='#ee7733',marker='|',ls='none',label='Pointwise 68.27% squared-Gaussian intervals')]
        fig.suptitle(f'Controlled monotonic null — realization {rep}\nSame parent/noise, nested subsets; ordering probabilities are NOT calibrated significance',fontsize=12)
        fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.045),ncol=1,fontsize=9,frameon=False)
        fig.text(.5,.015,'x: synthetic nuisance-reference coordinate; saved P is GH-marginalized (shared component shapes preserve ranking).',ha='center',fontsize=8)
        fig.subplots_adjust(left=.10,right=.98,bottom=.22,top=.87,hspace=.23,wspace=.10)
        name=f'replicate_{rep:02d}'
        fig.savefig(output_root/(name+'.png'),dpi=160);fig.savefig(output_root/(name+'.pdf'));plt.close(fig)
        sidecar={'source_hashes':hashes,'producer_source_hashes':manifest['source_sha256'],'settings':settings,
                 'numpy_version':np.__version__,'scipy_version':scipy.__version__,'replicate':rep,'replayed_parent':item,
                 'axis_limits':{'x':xlim,'y':ylim},'interval_quantiles':Q.tolist(),
                 'interval_policy':'conditional pointwise folded/squared Gaussian velocity; both signs; not GP/profile or simultaneous uncertainty',
                 'probabilities_recomputed':False,'new_realizations_generated':False,
                 'reconstruction':'recorded SeedSequence([master,replicate]) and hash-verified producer; parent arrays/subsets match exactly; original noisy vectors were not saved for independent byte comparison',
                 'zero_log_policy':'true zeros retained in sidecar; shown as down-arrows at log boundary, not silently replaced in data',
                 'gbar_policy':'median gas1/bulge0.7 reference; no horizontal nuisance intervals; all synthetic components share shape',
                 'panels':panels}
        with (output_root/(name+'.plotmeta.json')).open('x') as stream:
            json.dump(sidecar,stream,indent=2,allow_nan=False);stream.write('\n')
    if hashes!={key:sha256_file(path) for key,path in paths.items()}:
        raise RuntimeError('input/plot source changed during rendering')
    for name,expected_hash in manifest['source_sha256'].items():
        if sha256_file(ROOT/name)!=expected_hash: raise RuntimeError('producer dependency changed')
    print('PLOTTED_SAVED_REALIZATIONS',len(replay),'panels',len(rows),'output',output_root)


def self_test():
    from scipy.special import ndtri
    r=np.array([1.,2.,3.]);v=np.array([0.,.3,-.8]);e=np.array([1.,1.2,.9])
    plus=observed_interval(v,e,r);minus=observed_interval(-v,e,r)
    assert plus.shape==(2,len(v)) and len(Q)==2
    assert np.array_equal(plus,minus) and np.all(np.diff(plus,axis=0)>0)
    expected=ndtri((1+Q)/2)**2*vel2acc(np.ones(1),np.ones(1))
    assert np.allclose(plus[:,0],expected,rtol=1e-10,atol=0)
    assert np.all(plus[:,0]>0)  # central square is zero, interval is not
    print('Conditional interval sign/zero/central-Gaussian checks passed.')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--self-test',action='store_true')
    p.add_argument('--input-root',type=Path)
    p.add_argument('--output-root',type=Path)
    args=p.parse_args()
    if args.self_test: self_test();return
    if args.input_root is None or args.output_root is None: p.error('explicit input and fresh output roots required')
    plot_saved(args.input_root.resolve(),args.output_root.resolve())


if __name__=='__main__':
    main()
