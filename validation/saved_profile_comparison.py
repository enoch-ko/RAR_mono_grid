#!/usr/bin/env python3
"""Isolated saved-product comparison tools; never solve forces or probabilities.

All outputs are new, exclusive-create artifacts. Historical inputs remain immutable.
GP covariance J and full-density uniform-dR grid J_rad are DISTINCT metrics.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parent
POLICY = 'disc_to_catalogue_total_minus_published_bulge_v1'
PARSER = 'sparc_rotmod_header_v1_gas3_disk4'
RAW_CUT = 1.35e-3


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def write(path, document):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as stream:
        json.dump(document, stream, indent=2, allow_nan=False)
        stream.write('\n')


def profile_grid(rad, sb, family, mean, slope=0.):
    rad, sb = np.asarray(rad), np.asarray(sb)
    weight = sb * rad
    area = np.trapz(weight, rad)
    if area <= 0 or not np.isfinite(area):
        raise ValueError('invalid full-density luminosity integral')
    if family == 'constant':
        values = np.full_like(rad, mean, dtype=float)
    elif family == 'linear':
        pivot = np.trapz(weight * rad, rad) / area
        values = mean + slope * (rad - pivot)
    elif family == 'log_linear':
        exponent = slope * rad
        shape = np.exp(exponent - exponent.max())
        values = mean * shape * area / np.trapz(weight * shape, rad)
    else:
        raise ValueError('unknown profile family')
    if not np.isfinite(values).all() or np.any(values <= 0):
        raise ValueError('nonphysical profile')
    return values


def metrics(rad, values, family='piecewise_linear', slope=0.):
    """Exact continuous integrals, not equal weighting of irregular density knots."""
    r, u = np.asarray(rad, float), np.asarray(values, float)
    if len(r) < 2 or len(r) != len(u) or not np.isfinite([*r, *u]).all():
        raise ValueError('invalid profile arrays')
    if np.any(np.diff(r) <= 0) or np.any(r < 0) or np.any(u <= 0):
        raise ValueError('invalid full-density support or M/L')
    span = r[-1] - r[0]
    if family == 'log_linear' and slope != 0:
        x = abs(slope) * span
        top = float(np.max(u))
        m1 = top * (-math.expm1(-x)) / x
        m2 = top * top * (-math.expm1(-2 * x)) / (2 * x)
    else:
        dr = np.diff(r)
        m1 = float(np.sum(dr * (u[:-1] + u[1:]) / 2) / span)
        m2 = float(np.sum(dr * (u[:-1]**2 + u[:-1]*u[1:] + u[1:]**2) / 3) / span)
    # Compute fractional-square integral directly to avoid cancellation near .5.
    if family != 'log_linear' or slope == 0:
        z = (u - .5) / .5
        j = float(np.sum(np.diff(r)*(z[:-1]**2+z[:-1]*z[1:]+z[1:]**2)/3)/span)
    elif abs(slope)*span<1e-3:
        nodes,weights=np.polynomial.legendre.leggauss(8)
        x=(nodes+1)*span/2
        z=2*u[0]*np.expm1(slope*x)+(2*u[0]-1)
        j=float(np.dot(weights,z*z)/2)
    else:
        j = max(0., 4*m2 - 4*m1 + 1)
    return dict(J_rad=j, uniform_radial_mean_ml=m1, minimum_ml=float(u.min()),
                maximum_ml=float(u.max()), log10_range=float(np.log10(u.max())-np.log10(u.min())),
                rms_departure_ml=.5*math.sqrt(j), support_kpc=[float(r[0]), float(r[-1])],
                support='FULL original supplied .dens, including R=0; uniform continuous dR')


def external(path):
    bins, groups, rejected = [], {}, []
    with Path(path).open() as stream:
        for i, row in enumerate(csv.DictReader(stream), 2):
            try:
                radius, ml = float(row['Radius_kpc']), float(row['ml_K'])
                if not np.isfinite([radius, ml]).all() or radius < 0 or ml <= 0:
                    raise ValueError('invalid radius or M/L')
                label, galaxy = row['Sample'], row['Galaxy']
                if not galaxy or not label:
                    raise ValueError('missing galaxy or Sample')
                bins.append(ml)
                key = label + ':' + galaxy
                rec = groups.setdefault(key, dict(galaxy=galaxy, sample=label, points=[]))
                try:
                    rd = float(row['R_bar'])
                    if not math.isfinite(rd) or rd <= 0:
                        rd = None
                except (ValueError, KeyError):
                    rd = None
                rec['points'].append(dict(radius_kpc=radius, ml=ml, R_bar_kpc=rd))
            except (ValueError, KeyError) as error:
                rejected.append(dict(csv_line=i, reason=str(error)))
    for rec in groups.values():
        scales={p['R_bar_kpc'] for p in rec['points'] if p['R_bar_kpc'] is not None}
        rec['rd_status']='provisional_consistent' if len(scales)==1 else 'missing_or_inconsistent'
        if len(scales)!=1:
            for point in rec['points']:
                point['R_bar_kpc']=None
    if not bins:
        raise ValueError('no valid raw bins')
    return dict(source=str(path), sha256=sha(path), bin_count=len(bins),
                galaxy_sample_groups=len(groups), population_band_95=np.percentile(bins,[2.5,97.5]).tolist(),
                weighting='each valid raw CSV bin once; NOT equal-galaxy or bootstrap uncertainty',
                bandpass='external Ks; SPARC 3.6um; no conversion applied',
                radius='kpc; R_bar provisional Rd; absent/nonpositive scale excludes only scaled point',
                profiles=list(groups.values()), rejected=rejected)


def source_guard(mapping):
    for path, expected in mapping.items():
        path = Path(path) if Path(path).is_absolute() else ROOT/path
        if not path.is_file() or sha(path) != expected:
            raise ValueError('frozen source/input mismatch: ' + str(path))


def independent_inputs(galaxy, normalization, data_root, dens_root, table):
    """Reconstruct evidence from raw arrays, not from the reference's invariant dict."""
    import utils
    import inspect
    import mono_SPARC_RAR as mono
    import common_isotonic_reference as common
    denspath = Path(dens_root)/(galaxy+'.dens')
    rotpath = Path(data_root)/(galaxy+'_rotmod.dat')
    dens = np.loadtxt(denspath, ndmin=2)
    raw = np.loadtxt(rotpath, ndmin=2)
    if dens.shape[1] != 3 or raw.shape[1] != 8:
        raise ValueError('independent schema mismatch')
    rad = dens[:,0]
    if not np.isfinite(dens).all() or np.any(np.diff(rad)<=0) or np.any(rad<0):
        raise ValueError('invalid original density support')
    mask = np.isfinite(raw[:,0]) & (raw[:,0]>0) & (raw[:,0]>=rad[0]) & (raw[:,0]<=rad[-1])
    arr = raw[mask]
    if len(arr)<2 or not np.isfinite(arr[:,:6]).all() or np.any(arr[:,2]<=0):
        raise ValueError('invalid retained kinematics')
    if arr[:,0].min() < rad[rad>0].min():
        raise ValueError('kinematics would require inner force clamping')
    # Header-aware loader provides an additional independent gas3/disc4 sentinel check.
    named, parse = common.read_rotmod(rotpath, (rad[0],rad[-1]))
    for k, col in [('Rad',0),('Vobs',1),('errV',2),('Vgas',3),('Vbul',5)]:
        if not np.array_equal(named[k],arr[:,col]):
            raise ValueError('independent raw-column/index mismatch: '+k)
    norm = read(normalization)
    if norm['policy'] != POLICY:
        raise ValueError('unsupported normalization')
    rec = norm['galaxies'][galaxy]
    cat = utils.load_sparc_table(table)
    rows = cat.loc[cat.Galaxy==galaxy]
    if len(rows)!=1:
        raise ValueError('catalogue join not unique')
    cat = rows.iloc[0]
    if not np.isclose(rec['Ldisk_1e9']+rec['Lbulge_1e9'],float(cat.L),rtol=1e-12):
        raise ValueError('component luminosity mismatch')
    if not np.isclose(rec['kappa'],rec['Ldisk_1e9']/float(cat.L),rtol=1e-12):
        raise ValueError('kappa mismatch')
    geom = dict(distance_Mpc=float(cat.D),inclination_deg=float(cat.Inc),
                zd_kpc=.196*float(cat.Rdisk)**.633,vertical_profile='exp')
    force = dict(geometry=geom,force_source_sha256={
        'mono_SPARC_RAR.py':sha(ROOT/'mono_SPARC_RAR.py'),
        'utils.py':sha(ROOT/'utils.py'),'vcdisk':sha(inspect.getfile(mono.vcdisk))},
        density_support='full supplied profile; remove R=0 before Casertano',
        kinematic_support='overlap only; refuse outside positive force support',
        interpolation='signed W/R force interpolation then multiply observed R',
        signed_W='V*abs(V)',units=dict(radius='kpc',velocity='km/s',acceleration='m/s^2'),
        kpc_metres=3.086e19)
    inv = dict(force_model=force,galaxy=galaxy,parser_id=PARSER,
        input_sha256=dict(rotmod=sha(rotpath),dens=sha(denspath),table=sha(table),normalization=sha(normalization)),
        normalization_id=POLICY,nuisance_model=None,retained_mask=mask.tolist(),
        kinematics_sha256=common.digest({k:arr[:,c].tolist() for k,c in [('Rad',0),('Vobs',1),('errV',2),('Vgas',3),('Vbul',5)]}))
    return inv, rad, dens[:,1], float(cat.Rdisk)


def gp_document(path, manifest):
    d, frozen = read(path), read(manifest)
    source_guard(frozen['source_sha256'])
    source_guard(frozen['input_sha256'])
    prov = d['provenance']
    source_guard(prov['source_sha256'])
    if prov.get('input_parser',{}).get('id') != PARSER:
        raise ValueError('wrong or missing native GP parser evidence')
    if sha(ROOT/prov['script_file']) != prov['script_sha256']:
        raise ValueError('GP entrypoint mismatch')
    recovery = prov['recovery_configuration']
    if recovery['recovery_manifest_sha256'] != sha(manifest):
        raise ValueError('GP recovery manifest mismatch')
    inv, rad, sb, rd = independent_inputs(d['galaxy'],frozen['normalization'],
        frozen['data_root'],frozen['dens_root'],frozen['table'])
    if not np.array_equal(rad,np.asarray(d['knots_kpc'])):
        raise ValueError('GP profile knots differ from full original density support')
    raw_inputs = {k:v for k,v in inv['input_sha256'].items() if k!='normalization'}
    if raw_inputs != prov['input_sha256']:
        raise ValueError('native GP input hashes disagree')
    if d['normalization']['normalization_sha256'] != sha(frozen['normalization']):
        raise ValueError('native GP normalization mismatch')
    s = d['settings']
    inv['nuisance_model'] = dict(gas_median=s['gas_ml'],gas_dex=s['gas_ml_dex'],
        bulge_median=s['bulge_ml'],bulge_dex=s['bulge_ml_dex'],samples=s['quadrature_nodes_per_term'])
    inv['probability_settings'] = dict(ngrid=s['probability_starting_grid'],zpadding=s['zpadding'],
        log_tolerance=s['boundary_log_tolerance'],max_refinements=s['boundary_max_refinements'],
        floor_relative=s['signed_gbar_floor_relative'],max_tie_permutations=720)
    inv['probability_sha256'] = prov['source_sha256']['probability_recursion.py']
    if not np.array_equal(np.asarray(d['retained_radii_kpc']), np.loadtxt(Path(frozen['data_root'])/(d['galaxy']+'_rotmod.dat'))[np.asarray(inv['retained_mask']),0]):
        raise ValueError('native GP retained radii mismatch')
    if d['physical_domain']['policy'] != 'signed net gbar at every gas/bulge quadrature node must meet the existing numerical floor':
        raise ValueError('unreviewed native GP physical policy')
    retained = {p['record_index']:p for p in d['retained_profiles']}
    records = []
    for i, row in enumerate(d['records']):
        if i not in retained or not row['valid']:
            continue
        u = retained[i]['profile_on_photometric_grid']
        m = metrics(rad,u)
        logp, j = float(row['marginal_logp']), float(row['J_ref'])
        if math.isnan(logp) or logp>0 or not math.isfinite(j) or j<0:
            raise ValueError('invalid saved score')
        records.append(dict(id=i,valid=True,warnings=[],probability_kind='marginalised_gas_bulge',
            log_probability=logp,J_ref=j,radii_kpc=rad.tolist(),profile_ml=u,
            physical_summary=m,rd_kpc=rd))
    if not records:
        raise ValueError('no retained valid GP profiles')
    return dict(schema='common-reference-candidates-v1',invariants=inv,
        parser_source_sha256=sha(ROOT/'common_isotonic_reference.py'),
        provenance=dict(native=str(path),native_sha256=sha(path),manifest_sha256=sha(manifest),
                        parser_evidence=prov['input_parser'],independent_raw_column_checks=True),
        deviation_metric=dict(id='GP_covariance_J_unchanged',settings={k:s[k] for k in ['reference_length_fraction_of_span','reference_nugget']}),
        diagnostic_only=bool(d.get('diagnostic_only',True)),records=records,
        selection_scope='saved unchanged-J frontier only; not global physical-excursion minimum')


def compatible_settings(args):
    keys = ['grid_size','prob_grid_size','zpadding','tie_atol','max_tie_permutations',
        'gas_ml','bulge_ml','gas_ml_dex','bulge_ml_dex','disk_ml','disk_ml_dex',
        'disk_grid_half_width_dex','slope_min','slope_max','ml_samples','probability_method',
        'boundary_log_tolerance','boundary_max_refinements','marginalise_gas_bulge']
    return {k:args[k] for k in keys}


def choose_sources(normalization, old_root, corrected_root, family):
    norm = read(normalization)
    if norm['policy'] != POLICY:
        raise ValueError('unsupported normalization')
    groups = {'original':{},'corrected':{}}
    settings = None
    for kind, root in [('original',old_root),('corrected',corrected_root)]:
        for meta in sorted((Path(root)/family).glob('run_metadata.shard_*.json')):
            d = read(meta)
            if d['profile']!=family:
                raise ValueError('family mismatch')
            args = d['arguments']
            current = compatible_settings(args)
            if settings is not None and current!=settings:
                raise ValueError('incompatible scientific grid settings: '+str(meta))
            settings = current
            numeric = json.loads(args['numerical_source_sha256'])
            required={'mono_SPARC_RAR.py','mono_SPARC_RAR_boundary.py','probability_recursion.py','utils.py'}
            if not required.issubset(numeric):
                raise ValueError('missing numerical source evidence')
            source_guard({k:v for k,v in numeric.items() if k.endswith('.py')})
            if kind=='corrected' and (d.get('normalization',{}).get('policy')!=POLICY or d.get('normalization',{}).get('normalization_sha256')!=sha(normalization)):
                raise ValueError('corrected grid normalization mismatch')
            if kind=='corrected' and d['normalization']['launcher_sha256']!=sha(ROOT/'run_bulge_reruns_20261015.py'):
                raise ValueError('normalization launcher mismatch')
            suffix = meta.stem.split('shard_')[1]
            grid = meta.with_name('profile_grid.shard_'+suffix+'.csv')
            full = meta.with_name('fully_marginalised.shard_'+suffix+'.csv')
            if not grid.is_file() or not full.is_file():
                raise ValueError('incomplete final shard set')
            names = set()
            with full.open() as stream:
                for row in csv.DictReader(stream):
                    names.add(row['galaxy'])
            if len(names)!=d['galaxy_count']:
                raise ValueError('galaxy count mismatch')
            for galaxy in names:
                if galaxy not in norm['galaxies'] or galaxy in groups[kind]:
                    raise ValueError('unknown/duplicate galaxy '+galaxy)
                groups[kind][galaxy] = dict(grid=str(grid),metadata=str(meta),
                    fully_marginalised=str(full),metadata_sha256=sha(meta),numeric_source_sha256=numeric,
                    historical_input_hashes='not present in ordinary shard metadata; do not infer from current files')
    ledger=[]
    for galaxy, rec in norm['galaxies'].items():
        bulged = float(rec['Lbulge_1e9'])>0
        if not bulged and float(rec['kappa'])!=1.:
            raise ValueError('control is not unchanged kappa=1')
        kind = 'corrected' if bulged else 'original'
        source = groups[kind].get(galaxy)
        ledger.append(dict(galaxy=galaxy,family=family,required_source=kind,
            status='saved_grid_available' if source else 'pending_required_shard',
            normalization=rec,source=source,settings=settings,
            common_reference_status='pending_independent_grid_policy_evidence'))
    return ledger


def grid_summaries(ledger, dens_root):
    by_path={}
    for item in ledger:
        if item['source']:
            by_path.setdefault(item['source']['grid'],{})[item['galaxy']]=item
    selected=[]
    for path, targets in by_path.items():
        before=sha(path)
        for galaxy,item in targets.items():
            denspath=Path(dens_root)/(galaxy+'.dens')
            meta=read(item['source']['metadata'])
            rotpath=Path(meta['arguments']['data_root'])/(galaxy+'_rotmod.dat')
            raw=np.loadtxt(rotpath,ndmin=2)
            rad=np.loadtxt(denspath,ndmin=2)[:,0]
            mask=np.isfinite(raw[:,0]) & (raw[:,0]>0) & (raw[:,0]>=rad[0]) & (raw[:,0]<=rad[-1])
            item['independent_current_input_evidence']=dict(dens_sha256=sha(denspath),rotmod_sha256=sha(rotpath),retained_mask=mask.tolist(),full_density_support_kpc=[float(rad[0]),float(rad[-1])],proof_scope='current inputs and source-inferred mask, NOT historical input identity')
        state={g:dict(seen=set(),errors=0,count=0,best=None) for g in targets}
        arrays={g:np.loadtxt(Path(dens_root)/(g+'.dens'),ndmin=2) for g in targets}
        with Path(path).open() as stream:
            for row in csv.DictReader(stream):
                galaxy=row['galaxy']
                if galaxy not in targets or row['probability_kind']!='marginalised_gas_bulge':
                    continue
                item, st = targets[galaxy],state[galaxy]
                key=(row['grid_i'],row['grid_j'])
                if key in st['seen']:
                    raise ValueError('duplicate native grid cell '+galaxy)
                st['seen'].add(key)
                if row['valid']!='True':
                    if row['invalid_reason']!='nonpositive_or_nonfinite_disk_ml':
                        st['errors']+=1
                    continue
                logp=float(row['log_probability'])
                if math.isnan(logp) or logp>0 or row.get('warnings'):
                    st['errors']+=1
                    continue
                st['count']+=1
                if logp < math.log(RAW_CUT):
                    continue
                dens=arrays[galaxy]
                slope=0. if item['family']=='constant' else float(row['disk_slope'])
                try:
                    u=profile_grid(dens[:,0],dens[:,1],item['family'],float(row['upsilon_disk_mean']),slope)
                    m=metrics(dens[:,0],u,item['family'],slope)
                except ValueError:
                    st['errors']+=1
                    continue
                if st['best'] is None or m['J_rad']<st['best']['physical_summary']['J_rad']:
                    st['best']=dict(galaxy=galaxy,family=item['family'],log_probability=logp,
                        grid_i=row['grid_i'],grid_j=row['grid_j'],radii_kpc=dens[:,0].tolist(),
                        profile_ml=u.tolist(),physical_summary=m,source_grid=path,
                        selection='minimum found grid J_rad at fixed raw-P threshold; not common-q selection')
        if before!=sha(path):
            raise ValueError('input shard changed during read')
        for galaxy,item in targets.items():
            st=state[galaxy]
            expected=int(item['settings']['grid_size'])**(1 if item['family']=='constant' else 2)
            item.update(grid_sha256=before,observed_cells=len(st['seen']),expected_cells=expected,
                valid_cells=st['count'],numerical_invalid_cells=st['errors'])
            if len(st['seen'])!=expected or st['errors']:
                item['status']='numerical_or_incomplete_grid'
            elif st['best']:
                item['status']='raw_threshold_selected'
                selected.append(st['best'])
            else:
                item['status']='no_qualifying_evaluated_cell'
    return selected


def plot_summary(document, output, coordinate='kpc'):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,3,figsize=(15,4),sharey=True)
    band=document['external']['population_band_95']
    labels=['Grids: raw P>=1.35e-3 (partial)','GP: common q>=0.05 (3 pilots)','External raw bins (Ks)']
    for ax,label in zip(axes,labels):
        ax.axhspan(*band,color='orange',alpha=.12,label='pooled external raw-bin95%')
        ax.axhspan(.5*10**(-.1),.5*10**.1,color='grey',alpha=.15,label='SPARC .1dex prior')
        ax.set(title=label,xlabel='Radius [kpc]' if coordinate=='kpc' else 'R/Rd',yscale='log',ylabel='disc M/L')
    colors={'constant':'tab:blue','linear':'tab:green','log_linear':'tab:red'}
    for row in document.get('grid_selected',[]):
        axes[0].plot(np.asarray(row['radii_kpc'])/(row.get('rd_kpc',1.) if coordinate=='rd' else 1.),row['profile_ml'],color=colors[row['family']],alpha=.12,lw=.7)
    for row in document.get('gp_selected',[]):
        axes[1].plot(np.asarray(row['radii_kpc'])/(row.get('rd_kpc',1.) if coordinate=='rd' else 1.),row['profile_ml'],color='purple',alpha=.4,lw=.8)
    palette={'XMM':'tab:blue','COSMOS_kin':'tab:orange','COSMOS_DR1':'tab:green','LADUMA':'tab:red'}
    for row in document['external']['profiles']:
        points=sorted([p for p in row['points'] if coordinate=='kpc' or p['R_bar_kpc'] is not None],key=lambda p:p['radius_kpc'])
        axes[2].plot([p['radius_kpc']/(p['R_bar_kpc'] if coordinate=='rd' else 1.) for p in points],[p['ml'] for p in points],
                     color=palette.get(row['sample'],'grey'),alpha=.2,lw=.6)
    from matplotlib.lines import Line2D
    axes[0].legend(handles=[Line2D([0],[0],color=c,label=f) for f,c in colors.items()],fontsize=7)
    axes[2].legend(handles=[Line2D([0],[0],color=c,label=f) for f,c in palette.items()],fontsize=7)
    from matplotlib.patches import Patch
    axes[1].legend(handles=[Patch(facecolor='orange',alpha=.12,label='External pooled-bin95% spread'),Patch(facecolor='grey',alpha=.15,label='SPARC .5 +/- .1dex prior')],fontsize=7)
    fig.suptitle('Saved-product diagnostic; partial coverage; thresholds/metrics differ; Ks vs 3.6um unresolved',fontsize=9)
    fig.tight_layout()
    if Path(output).exists():
        raise FileExistsError(output)
    fig.savefig(output,dpi=130)
    plt.close(fig)


def run(args):
    import time
    started=time.monotonic()
    out=Path(args.output_root)
    out.mkdir(parents=True,exist_ok=False)
    norm=ROOT/'outputs/bulge_normalized_20261015/inputs/normalization.json'
    manifest=ROOT/'outputs/corrected_recovery_prep_20261016/manifest.json'
    frozen=read(manifest)
    source_guard(frozen['source_sha256'])
    source_guard(frozen['input_sha256'])
    document=dict(schema='saved-profile-comparison-v1',script_sha256=sha(__file__),
        diagnostic_only=True,grid_selected=[],gp_selected=[],gp_raw_selected=[],ledger=[],gp_status=[],
        external=external(ROOT/'ML_Andreea.csv'),raw_probability_cut=RAW_CUT,
        caveats=['No sigma significance or population inference',
                 'Grid J_rad and GP covariance J are distinct; no frontier reranking',
                 'GP pilot B64 and tiny finite search are not convergence evidence',
                 'Centrally ordered disc QP is superseded and not a dependency',
                 'Current input hashes do not prove historical grid input immutability'])
    for family in ['constant','linear','log_linear']:
        ledger=choose_sources(norm,ROOT/'outputs/rerun_boundary_20261006',ROOT/'outputs/bulge_normalized_20261015',family)
        document['grid_selected']+=grid_summaries(ledger,frozen['dens_root'])
        document['ledger']+=ledger
    import common_reference_driver as driver
    import utils
    catalogue=utils.load_sparc_table(frozen['table']).set_index('Galaxy')
    for row in document['grid_selected']:
        row['rd_kpc']=float(catalogue.loc[row['galaxy'],'Rdisk'])
    for path in sorted((ROOT/'outputs/corrected_recovery_pilots_20261016/gp').glob('*/*_gp_search.json')):
        galaxy=path.parent.name
        try:
            native=gp_document(path,manifest)
            reference_path=ROOT/'outputs/corrected_recovery_pilots_20261016/reference'/(galaxy+'_B64.json')
            ref=read(reference_path)
            # Existing strict annotator enforces independently built invariant equality.
            annotation=driver.annotate_document(ref,native)
            raw_records=[r for r in annotation['records'] if r['log_probability']>=math.log(RAW_CUT)]
            if raw_records:
                document['gp_raw_selected'].append(dict(min(raw_records,key=lambda r:r['J_ref']),galaxy=galaxy))
            index=annotation['minimum_found_J_record_index']
            if index is not None:
                row=annotation['records'][index]
                document['gp_selected'].append(dict(row,galaxy=galaxy))
            write(out/(galaxy+'_gp_candidates.json'),native)
            write(out/(galaxy+'_gp_annotations.json'),annotation)
            document['gp_status'].append(dict(galaxy=galaxy,status='annotated_diagnostic_only',reference_sha256=sha(reference_path)))
        except (ValueError,KeyError,FileNotFoundError) as error:
            document['gp_status'].append(dict(galaxy=galaxy,status='incompatible_or_missing',reason=str(error)))
    for item in document['ledger']:
        if item['galaxy'] in [p.parent.name for p in (ROOT/'outputs/corrected_recovery_pilots_20261016/gp').glob('*/*_gp_search.json')]:
            try:
                grid_native_document(item,None)
            except ValueError as error:
                item['common_reference_blocker']=str(error)
    source_guard(frozen['source_sha256'])
    source_guard(frozen['input_sha256'])
    document['elapsed_saved_processing_seconds']=time.monotonic()-started
    write(out/'comparison.json',document)
    if args.plot:
        plot_summary(document,out/'paired_clouds_diagnostic.png')
        plot_summary(document,out/'paired_clouds_rd_diagnostic.png','rd')
    print('SAVED',out,'grid selected',len(document['grid_selected']),'GP selected',len(document['gp_selected']))
    print('GP STATUS',document['gp_status'])


def grid_native_document(item, reference):
    """Explicit fail-closed native grid adapter: do not fabricate policy equivalence.

    Ordinary CSVs do not retain per-node signed gbar/floor certificates or complete
    historical input hashes. Independent current masks alone cannot supply those.
    Compatible grid q-selection therefore needs a validated evidence export first.
    Raw-P/J_rad selections above remain available, with separate provenance labels.
    """
    if item.get('source') is None:
        raise ValueError('missing required component-normalized grid shard')
    if item.get('common_reference_status')!='independently_verified_grid_invariants':
        raise ValueError('native grid reference annotation blocked: historical input identity and all-node floor/support equivalence unproved; no rescore authorized')
    raise ValueError('unrecognized evidence version; no automatic reference-invariant copying')


class Tests(unittest.TestCase):
    def test_constant(self):
        self.assertEqual(metrics([0,1,9],[.5,.5,.5])['J_rad'],0.)
        self.assertAlmostEqual(metrics([0,1,9],[1,1,1])['J_rad'],1.)
    def test_irregular_linear(self):
        self.assertAlmostEqual(metrics([0,.1,2],[.5,.55,1.5])['J_rad'],4/3)
    def test_full_support(self):
        self.assertNotEqual(metrics([0,1,2],[.5,1,1])['J_rad'],metrics([1,2],[1,1])['J_rad'])
    def test_invalid(self):
        for r,u in [([0,0],[1,1]),([0,1],[0,1]),([0,1],[1,float('nan')])]:
            with self.assertRaises(ValueError): metrics(r,u)
    def test_loglinear(self):
        r=np.array([0,.1,2.]);u=.5*np.exp(.2*r)
        m=metrics(r,u,'log_linear',.2)
        x=np.linspace(0,2,100001)
        self.assertAlmostEqual(m['J_rad'],float(np.trapz((np.exp(.2*x)-1)**2,x)/2),places=9)
    def test_normalization(self):
        r=np.array([0,.1,2.]);sb=np.array([10,9,1.])
        for family in ['constant','linear','log_linear']:
            u=profile_grid(r,sb,family,.5,.05)
            self.assertAlmostEqual(float(np.trapz(sb*r*u,r)/np.trapz(sb*r,r)),.5)
    def test_reject_grid_evidence(self):
        with self.assertRaisesRegex(ValueError,'missing required'):
            grid_native_document(dict(source=None),{})
        with self.assertRaisesRegex(ValueError,'unproved'):
            grid_native_document(dict(source={'grid':'x'}),{'context':{'invariants':{}}})
    def test_small_loglinear(self):
        r=np.array([0,.1,2.]);u=.5*np.exp(1e-8*r)
        self.assertLess(abs(metrics(r,u,'log_linear',1e-8)['J_rad']),1e-14)
    def test_external_bin_weighting(self):
        import tempfile
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp)/'raw.csv'
            with p.open('w') as stream:
                stream.write('Galaxy,Sample,Radius_kpc,ml_K,R_bar\nA,XMM,1,1,2\nA,XMM,2,2,2\nB,LADUMA,1,9,3\nB,LADUMA,2,-1,3\n')
            d=external(p)
            self.assertEqual(d['bin_count'],3)
            self.assertEqual(len(d['rejected']),1)
            self.assertEqual(d['population_band_95'],np.percentile([1,2,9],[2.5,97.5]).tolist())
    def test_exclusive_output(self):
        import tempfile
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp)/'a.json'
            write(p,{'x':1})
            with self.assertRaises(FileExistsError):write(p,{'x':2})
    def test_source_guard(self):
        with self.assertRaisesRegex(ValueError,'mismatch'):
            source_guard({str(ROOT/'saved_profile_comparison.py'):'wrong'})
    def test_negative_slope(self):
        r=np.array([0,.1,2.]);u=np.exp(-.2*r)
        x=np.linspace(0,2,100001)
        self.assertAlmostEqual(metrics(r,u,'log_linear',-.2)['J_rad'],float(np.trapz(((np.exp(-.2*x)-.5)/.5)**2,x)/2),places=9)
    def fixture(self, directory):
        root=Path(directory)
        rec=dict(Lbulge_1e9=0.,kappa=1.)
        norm=root/'normalization.json'
        write(norm,dict(policy=POLICY,galaxies={'control':rec,'bulged':dict(Lbulge_1e9=1.,kappa=.5)}))
        a=dict(grid_size='2',prob_grid_size='400',zpadding='8',tie_atol='0',max_tie_permutations='720',gas_ml='1',bulge_ml='.7',gas_ml_dex='.04',bulge_ml_dex='.1',disk_ml='.5',disk_ml_dex='.1',disk_grid_half_width_dex='1',slope_min='-5',slope_max='5',ml_samples='9',probability_method='boundary_log_cumulative_v1',boundary_log_tolerance='2e-5',boundary_max_refinements='8',marginalise_gas_bulge='True')
        a['numerical_source_sha256']=json.dumps({k:sha(ROOT/k) for k in ['mono_SPARC_RAR.py','mono_SPARC_RAR_boundary.py','probability_recursion.py','utils.py']})
        p=root/'old/constant'
        p.mkdir(parents=True)
        for name,text in [('profile_grid.shard_01.csv','galaxy\ncontrol\nbulged\n'),('fully_marginalised.shard_01.csv','galaxy\ncontrol\nbulged\n')]:
            with (p/name).open('w') as stream:stream.write(text)
        write(p/'run_metadata.shard_01.json',dict(profile='constant',arguments=a,galaxy_count=2))
        return norm,root/'old',root/'corrected'
    def test_no_stale_bulge_fallback(self):
        import tempfile
        with tempfile.TemporaryDirectory() as temp:
            norm,old,new=self.fixture(temp)
            rows=choose_sources(norm,old,new,'constant')
            self.assertEqual(rows[0]['status'],'saved_grid_available')
            self.assertEqual(rows[1]['status'],'pending_required_shard')
            self.assertIsNone(rows[1]['source'])
    def test_duplicate_galaxy_refusal(self):
        import tempfile
        with tempfile.TemporaryDirectory() as temp:
            norm,old,new=self.fixture(temp)
            p=old/'constant'
            for name in ['run_metadata','profile_grid','fully_marginalised']:
                suffix='json' if name=='run_metadata' else 'csv'
                with (p/(name+'.shard_02.'+suffix)).open('w') as stream:
                    stream.write((p/(name+'.shard_01.'+suffix)).read_text())
            with self.assertRaisesRegex(ValueError,'duplicate'):
                choose_sources(norm,old,new,'constant')
    def test_settings_refusal(self):
        import tempfile
        with tempfile.TemporaryDirectory() as temp:
            norm,old,new=self.fixture(temp)
            p=old/'constant'
            d=read(p/'run_metadata.shard_01.json')
            d['arguments']['prob_grid_size']='800'
            write(p/'run_metadata.shard_02.json',d)
            with self.assertRaisesRegex(ValueError,'scientific grid settings'):
                choose_sources(norm,old,new,'constant')
    def test_reference_policy_refusal(self):
        import common_reference_driver as driver
        reference={'context':{'invariants':{'floor':1}}}
        document={'schema':'common-reference-candidates-v1','invariants':{'floor':0}}
        with self.assertRaisesRegex(ValueError,'invariant provenance mismatch'):
            driver.annotate_document(reference,document)
    def test_metric_not_gp_ranking(self):
        self.assertNotIn('J_rad',gp_document.__code__.co_consts)
    def test_no_reference_copy(self):
        self.assertNotIn('reference',gp_document.__code__.co_varnames)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-root',type=Path)
    parser.add_argument('--plot',action='store_true')
    parser.add_argument('--self-test',action='store_true')
    args=parser.parse_args()
    if args.self_test:
        suite=unittest.defaultTestLoader.loadTestsFromTestCase(Tests)
        result=unittest.TextTestRunner(verbosity=2).run(suite)
        if args.output_root is not None:
            write(args.output_root/'self_test_receipt.json',dict(tests=result.testsRun,failures=len(result.failures),errors=len(result.errors),successful=result.wasSuccessful(),script_sha256=sha(__file__)))
        sys.exit(0 if result.wasSuccessful() else 1)
    if args.output_root is None:
        parser.error('--output-root required (must not already exist)')
    run(args)
