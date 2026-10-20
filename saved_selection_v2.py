#!/usr/bin/env python3
"""Versioned evidence-led saved grid/GP export; no force, P, mocks or search calls.

The score, reference-availability and native-admissibility contracts are separate.
All outputs are exclusive-create. Unproved identity is WITHHELD, not disproved.
"""
import argparse
import bisect
from collections import Counter
import copy
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import sys
import time
import unittest
from unittest.mock import patch
import numpy as np
import saved_profile_comparison as legacy

ROOT = Path(__file__).resolve().parent
LEGACY_SHA = '2bc5fd807e662ce976c2626c151aa020e3e22ceb51e76cdbec356ad0a2331122'
RAW_P = .5 * math.erfc(3 / math.sqrt(2))
RAW_LOGP = math.log(RAW_P)
THRESHOLDS = {
    'raw_Phi_minus3': {'kind':'raw','probability':RAW_P,'log_probability':RAW_LOGP,
        'label':'raw P >= Phi(-3); Gaussian-tail numerical comparison, not calibrated significance'},
    'common_q05': {'kind':'rank','numerator':1,'denominator':20,
        'label':'main operational common-reference q >= 0.05 (inclusive)'},
    'common_q10': {'kind':'rank','numerator':1,'denominator':10,
        'label':'sensitivity operational common-reference q >= 0.10 (inclusive)'}}
VERIFIED = 'verified_score_compatible'
WITHHELD = 'evidence_incomplete_withheld'
MISMATCH = 'actual_scoring_mismatch'
INVALID = 'native_invalid'
WARNING = 'warning_tie_exception'


def validate_reference(reference):
    import common_isotonic_reference as common
    if reference.get('status')!='complete' or reference.get('failures'):
        raise ValueError('reference unavailable or failed mock')
    scores = reference.get('sorted_mock_logp')
    if not isinstance(scores,list) or not scores or len(scores)!=reference.get('B'):
        raise ValueError('reference B/score-count mismatch')
    if any(not math.isfinite(x) for x in scores) or scores!=sorted(scores):
        raise ValueError('reference mock scores must be finite and sorted')
    if common.digest(reference['context'])!=reference['context_sha256']:
        raise ValueError('reference context identity mismatch')
    return scores


def rank(logp, scores):
    if math.isnan(logp) or logp>0:
        raise ValueError('invalid saved logP')
    count = bisect.bisect_right(scores,logp)
    return {'count':count,'B':len(scores),'q':count/len(scores)}


def qualifies(row, threshold):
    if threshold['kind']=='raw':
        return row['log_probability']>=threshold['log_probability']
    r = row.get('common_reference_rank')
    return r is not None and r['count']*threshold['denominator']>=r['B']*threshold['numerator']


def cell_tier(evidence, native_valid=True, warnings=False):
    if not native_valid:
        return INVALID
    if warnings:
        return WARNING
    if evidence['mismatches']:
        return MISMATCH
    if evidence['missing']:
        return WITHHELD
    return VERIFIED


def select_all(records, metric, complete=True, tainted=False, allow_main=True):
    """Select from ALL supplied cells, never from old threshold winners alone."""
    results={}
    for name, threshold in THRESHOLDS.items():
        pool=[r for r in records if r.get('native_valid',True) and not r.get('warnings')
              and r.get(metric) is not None and qualifies(r,threshold)]
        winner=min(pool,key=lambda r:(r[metric],-r['log_probability'],str(r['id']))) if pool else None
        baseline=next((r for r in pool if r.get('is_baseline')),None)
        if tainted:
            status='failed_or_warning_tainted_galaxy'
        elif not complete:
            status='withheld_incomplete_galaxy'
        elif not allow_main:
            status='withheld_evidence_or_frontier'
        elif baseline is not None:
            status='baseline_already_acceptable'
        elif winner is not None:
            status='selected_minimum_found'
        else:
            status='no_qualifying_evaluated_candidate'
        proof_clean=all(r.get('tier')==VERIFIED for r in records if r.get('native_valid',True) and not r.get('warnings'))
        if allow_main and not proof_clean and not tainted and complete:
            status='withheld_mixed_candidate_evidence'
        main = complete and not tainted and allow_main and proof_clean
        results[name]={'status':status,'winner':winner if main else None,
            'diagnostic_finite_pool_winner':winner if not main else None,
            'selection_metric':metric,'global_optimum_claim':False,
            'baseline_qualifies_in_supplied_pool':baseline is not None}
    return results


def bound_manifest(meta, directory):
    """Use an output-bound historical digest, not today's arbitrary input manifest."""
    binding=meta.get('normalization',{}).get('manifest_sha256')
    path=Path(directory)/'batch_manifest.json'
    if not binding:
        return None,{'status':'missing','reason':'no historical manifest binding in native metadata'}
    if not path.is_file():
        return None,{'status':'missing','reason':'bound historical manifest unavailable'}
    if legacy.sha(path)!=binding:
        return None,{'status':'mismatch','reason':'bound manifest digest differs'}
    d=legacy.read(path)
    # Only the already-audited checksum-enforcing launcher is supported here.
    launcher=ROOT/'run_bulge_reruns_20261015.py'
    if meta['normalization'].get('launcher_sha256')!=d['source_sha256'].get(str(launcher)):
        return None,{'status':'mismatch','reason':'historical launcher binding differs'}
    if legacy.sha(launcher)!=meta['normalization']['launcher_sha256']:
        return None,{'status':'mismatch','reason':'checksum-enforcement launcher bytes differ'}
    return d,{'status':'verified_historical_binding','path':str(path),'sha256':binding,
        'enforcement':'worker load() checks pinned inputs/sources before each galaxy and at completion'}


def score_contract(inv):
    """Reference baseline floor is availability QC, NOT a candidate-score gate."""
    return {'galaxy':inv['galaxy'],'parser_id':inv['parser_id'],
        'input_sha256':inv['input_sha256'],'retained_mask':inv['retained_mask'],
        'kinematics_sha256':inv['kinematics_sha256'],'force_model':inv['force_model'],
        'normalization_id':inv['normalization_id'],'nuisance_model':inv['nuisance_model'],
        'probability_sha256':inv['probability_sha256'],
        'probability_settings':{k:v for k,v in inv['probability_settings'].items() if k!='floor_relative'}}


def grid_evidence(item, reference, normalization):
    e={'schema':'grid-score-evidence-v2','missing':[],'mismatches':[],
       'native_admissibility':{'disc_ml':'positive finite on full supplied density profile',
            'net_gbar':'finite signed values sorted as supplied; NO GP positive-floor screen',
            'nuisance_nodes':'all nodes; no failed-node removal'},
       'score_scope':'warning-free, exact-tie cases within 720 permutations; complete galaxy QC',
       'reference_availability':{'status':'missing_saved_reference'},'historical_proof':None}
    if reference is not None:
        try:
            validate_reference(reference)
            e['reference_availability']={'status':'complete_saved_CDF','B':reference['B'],
                'context_sha256':reference['context_sha256'],
                'baseline_floor_qc_only':reference['context']['probability_settings']['floor_relative'],
                'production_precision_established':False}
        except (ValueError,KeyError) as error:
            e['reference_availability']={'status':'invalid_saved_CDF','reason':str(error)}
    if item.get('source') is None:
        e['missing'].append('required native grid shard')
        return e
    source=item['source']; meta_path=Path(source['metadata']); meta=legacy.read(meta_path)
    e['native_metadata']={'path':str(meta_path),'sha256':legacy.sha(meta_path)}
    a=meta['arguments']; numeric=json.loads(a['numerical_source_sha256'])
    required={'mono_SPARC_RAR.py','mono_SPARC_RAR_boundary.py','probability_recursion.py','utils.py'}
    if not required.issubset(numeric):
        e['missing'].append('native numerical source evidence incomplete')
        return e
    try:
        legacy.source_guard({k:v for k,v in numeric.items() if k.endswith('.py')})
    except ValueError as error:
        e['missing'].append('archived numerical source bytes unavailable: '+str(error))
    historical,proof=bound_manifest(meta,meta_path.parent.parent)
    e['historical_proof']=proof
    if proof['status']=='mismatch':
        e['mismatches'].append(proof['reason'])
    elif historical is None:
        e['missing'].append('historical raw inputs for this run (later snapshots do not bind earlier jobs)')
    if reference is None:
        e['missing'].append('saved fixed baseline reference for galaxy')
        return e
    try:
        validate_reference(reference)
    except ValueError as error:
        e['missing'].append(str(error))
        return e
    e['reference_availability']={'status':'complete_saved_CDF','B':reference['B'],
        'context_sha256':reference['context_sha256'],
        'baseline_floor_qc_only':reference['context']['probability_settings']['floor_relative'],
        'production_precision_established':False}
    try:
        inv,rad,sb,rd=legacy.independent_inputs(item['galaxy'],normalization,
            Path(a['data_root']),Path(a['dens_root']),Path(a['table_path']))
    except (ValueError,KeyError,FileNotFoundError) as error:
        e['missing'].append('independent raw/schema/support evidence unavailable: '+str(error))
        return e
    # These are separately reconstructed CURRENT arrays, never historical proof alone.
    e['current_reconstruction_diagnostic']={'retained_mask':inv['retained_mask'],
        'kinematics_sha256':inv['kinematics_sha256'],'rd_kpc':rd}
    paths={'rotmod':Path(a['data_root'])/(item['galaxy']+'_rotmod.dat'),
        'dens':Path(a['dens_root'])/(item['galaxy']+'.dens'),'table':Path(a['table_path']),
        'normalization':Path(normalization)}
    recorded={}
    for key,path in paths.items():
        digest=None if historical is None else historical.get('input_sha256',{}).get(str(path))
        recorded[key]=digest
        if digest is None:
            e['missing'].append('historical '+key+' digest')
        elif digest!=reference['context']['input_sha256'][key]:
            e['mismatches'].append('historical/reference '+key+' identity differs')
        elif legacy.sha(path)!=digest:
            e['missing'].append('bound historical '+key+' version unavailable at declared path')
    solver=[] if historical is None else [v for k,v in historical.get('source_sha256',{}).items() if 'vcdisk' in k.lower()]
    if len(set(solver))!=1:
        e['missing'].append('historical installed vcdisk identity bound to this grid run')
        historical_solver=None
    else:
        historical_solver=solver[0]
    # Build a score contract from native settings and historical hashes, not reference copying.
    inv['input_sha256']=recorded
    inv['force_model']['force_source_sha256']['vcdisk']=historical_solver
    inv['force_model']['force_source_sha256']['mono_SPARC_RAR.py']=numeric['mono_SPARC_RAR.py']
    inv['force_model']['force_source_sha256']['utils.py']=numeric['utils.py']
    inv['nuisance_model']={'gas_median':float(a['gas_ml']),'gas_dex':float(a['gas_ml_dex']),
        'bulge_median':float(a['bulge_ml']),'bulge_dex':float(a['bulge_ml_dex']),'samples':int(a['ml_samples'])}
    inv['probability_settings']={'ngrid':int(a['prob_grid_size']),'zpadding':float(a['zpadding']),
        'log_tolerance':float(a['boundary_log_tolerance']),'max_refinements':int(a['boundary_max_refinements']),
        'max_tie_permutations':int(a['max_tie_permutations'])}
    inv['probability_sha256']=numeric['probability_recursion.py']
    e['native_score_contract']=score_contract(inv)
    if float(a['tie_atol'])!=0:
        e['mismatches'].append('native tie tolerance differs from exact reference ties')
    target=score_contract(reference['context']['invariants'])
    for key in ['galaxy','parser_id','retained_mask','kinematics_sha256','normalization_id',
                'nuisance_model','probability_settings','probability_sha256']:
        if e['native_score_contract'][key]!=target[key]:
            e['mismatches'].append('native/reference score contract differs: '+key)
    for key in ['geometry','interpolation','signed_W','units','kpc_metres']:
        if inv['force_model'][key]!=target['force_model'][key]:
            e['mismatches'].append('native/reference force contract differs: '+key)
    for key,digest in inv['force_model']['force_source_sha256'].items():
        if digest is not None and digest!=target['force_model']['force_source_sha256'][key]:
            e['mismatches'].append('historical/reference force source differs: '+key)
    if not e['missing'] and not e['mismatches'] and e['native_score_contract']!=target:
        e['mismatches'].append('unhandled score-contract difference')
    e['compatibility_tier']=cell_tier(e)
    return e


def frontier_sufficiency(native,records,removed_ids=(),metric_unchanged=True):
    problems=[]
    storage=native.get('checkpoint_storage',{})
    summary=storage.get('summary',{})
    ids={r['id'] for r in records}
    native_profiles={p['record_index']:p for p in native.get('retained_profiles',[])}
    by_record={r['id']:r for r in records}
    for point in native.get('frontier',[]):
        index=point['record_index']
        expected=native_profiles.get(index)
        rec=by_record.get(index)
        if expected is None or rec is None or point.get('profile_on_photometric_grid')!=expected.get('profile_on_photometric_grid') or expected.get('profile_on_photometric_grid')!=rec.get('profile_ml'):
            problems.append('frontier/profile identity missing or inconsistent')
    if len(ids)!=len(records):problems.append('duplicate retained record IDs')
    if summary.get('evaluated')!=native.get('physical_domain',{}).get('expected_evaluated_candidates'):
        problems.append('completed-search expected/evaluated count mismatch')
    allowed_native_rejections={'signed_gbar_below_floor','nonpositive_or_nonfinite_disk_ml'}
    if any(k not in allowed_native_rejections and v for k,v in summary.get('failure_counts',{}).items()):
        problems.append('unresolved numerical/unknown failed evaluations; no finite-search minimum certificate')
    fronts=[x['record_index'] for x in native.get('frontier',[])]
    base,maximum=summary.get('baseline_record_index'),summary.get('maximum_probability_record_index')
    if native.get('status')!='finite_search_completed':problems.append('search not complete')
    if storage.get('mode')!='retained-only' or storage.get('discard_policy')!='exact dominance; baseline, threshold minima and maximum P pinned':
        problems.append('exact full-frontier retention contract absent; old pinned winners insufficient')
    if not metric_unchanged:problems.append('metric changed; retained J frontier cannot be reranked')
    if len(fronts)!=summary.get('frontier_size') or len(set(fronts))!=len(fronts) or not fronts:
        problems.append('missing/duplicate full frontier')
    if not set(fronts).issubset(ids):problems.append('frontier profile/record unavailable')
    if set(fronts)&set(removed_ids):problems.append('eligibility removed a dominating frontier candidate')
    if len(native.get('records',[]))!=summary.get('retained') or summary.get('evaluated',-1)!=summary.get('valid',0)+summary.get('invalid',0):
        problems.append('retention/evaluation denominators inconsistent')
    if summary.get('discarded')!=summary.get('evaluated',0)-summary.get('retained',0):
        problems.append('discard count inconsistent')
    if base not in ids:problems.append('baseline missing or newly ineligible')
    if maximum not in ids:problems.append('maximum-P profile missing or newly ineligible')
    by_id={r['id']:r for r in records}
    frontier=[by_id[i] for i in fronts if i in by_id]
    for i,r in enumerate(frontier):
        for j,s in enumerate(frontier):
            if i!=j and s['J_ref']<=r['J_ref'] and s['log_probability']>=r['log_probability']:
                problems.append('saved frontier contains dominated/duplicate point')
    for r in records:
        if not any(f['J_ref']<=r['J_ref'] and f['log_probability']>=r['log_probability'] for f in frontier):
            problems.append('retained point not covered by full frontier')
    if maximum in by_id and any(r['log_probability']>by_id[maximum]['log_probability'] for r in records):
        problems.append('maximum-P pin inconsistent')
    return {'sufficient':not problems,'problems':sorted(set(problems)),
        'scope':'minimum FOUND in unchanged native eligible search, not global optimum',
        'baseline_id':base,'maximum_probability_id':maximum,'frontier_ids':fronts,
        'evaluated':summary.get('evaluated'),'valid':summary.get('valid'),'invalid':summary.get('invalid')}


def attach_metrics(record,rad,values,family,slope=0.):
    metrics=legacy.metrics(rad,values,family,slope)
    metrics['absolute_rms_departure_ml']=metrics.pop('rms_departure_ml')
    metrics['fractional_rms_departure']=math.sqrt(metrics['J_rad'])
    return dict(record,radii_kpc=np.asarray(rad).tolist(),profile_ml=np.asarray(values).tolist(),
                J_rad=metrics['J_rad'],physical_summary=metrics)


def grid_export(items,references,out,normalization):
    groups={}; summaries=[]
    for item in items:
        g,f=item['galaxy'],item['family']
        ref=references.get(g)
        evidence=grid_evidence(item,ref,normalization)
        entry={'galaxy':g,'family':f,'evidence':evidence,'source':item.get('source'),
            'reference_B':None if ref is None else ref.get('B'),
            'expected_cells':None if item['settings'] is None else int(item['settings']['grid_size'])**(1 if f=='constant' else 2),
            'observed_cells':0,'native_valid':0,'native_invalid':0,'numerical_failures':0,'warning_tie_exceptions':0,
            'records':[],'complete':False,'reconstruction_failures':0,'duplicate_cells':0}
        if item.get('source'):
            groups.setdefault(item['source']['grid'],{})[g]=(item,entry)
        else:
            entry['selections']=select_all([], 'J_rad',complete=False,allow_main=False)
            summaries.append(entry)
    fields=['galaxy','family','grid_i','grid_j','log_probability','raw_P_display','J_rad',
            'q_common_ref','q_count_inclusive','tier','native_valid','export_valid','reason']
    with (out/'grid_cell_annotations.csv').open('x',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader()
        for path,targets in groups.items():
            digest=legacy.sha(path)
            states={g:{'seen':set(),'pools':[],'best':{},'baseline':None} for g in targets}
            densities={g:np.loadtxt(Path(legacy.read(item['source']['metadata'])['arguments']['dens_root'])/(g+'.dens'),ndmin=2)
                       for g,(item,entry) in targets.items()}
            with Path(path).open() as native:
                for row in csv.DictReader(native):
                    g=row['galaxy']
                    if g not in targets or row['probability_kind']!='marginalised_gas_bulge':continue
                    item,entry=targets[g];state=states[g]
                    key=(row['grid_i'],row['grid_j'])
                    if key in state['seen']:
                        entry['duplicate_cells']+=1
                        continue
                    state['seen'].add(key);entry['observed_cells']+=1
                    valid=row['valid']=='True';warning=bool(row.get('warnings'))
                    native_invalid=not valid and row.get('invalid_reason')=='nonpositive_or_nonfinite_disk_ml'
                    if not valid:
                        entry['native_invalid' if native_invalid else 'numerical_failures']+=1
                        tier=INVALID
                        record=None
                    else:
                        lp=float(row['log_probability'])
                        if not math.isfinite(lp) or lp>0:
                            entry['numerical_failures']+=1;valid=False;tier=INVALID;record=None
                        else:
                            entry['native_valid']+=1
                            if warning:entry['warning_tie_exceptions']+=1
                            tier=cell_tier(entry['evidence'],warnings=warning)
                            ref=references.get(g)
                            ranks=None
                            if ref is not None and not warning and tier!=MISMATCH:
                                try:ranks=rank(lp,validate_reference(ref))
                                except ValueError:pass
                            slope=0. if item['family']=='constant' else float(row['disk_slope'])
                            baseline=math.isclose(float(row['upsilon_disk_mean']),.5,rel_tol=0,abs_tol=1e-14) and slope==0
                            record={'id':':'.join(key),'log_probability':lp,'native_valid':True,
                                'warnings':row.get('warnings',''),'tier':tier,'is_baseline':baseline,
                                'common_reference_rank':ranks,'native_admissibility':'positive disc ML; finite SIGNED net gbar; no GP gate'}
                            # Read all relevant cells, not the old raw-P winners.
                            needs_profile=ref is not None or baseline or any(qualifies(record,t) for t in THRESHOLDS.values())
                            if needs_profile:
                                try:
                                    dens=densities[g]
                                    u=legacy.profile_grid(dens[:,0],dens[:,1],item['family'],float(row['upsilon_disk_mean']),slope)
                                    record=attach_metrics(record,dens[:,0],u,item['family'],slope)
                                except ValueError:
                                    entry['reconstruction_failures']+=1;valid=False;tier=WITHHELD;record=None
                            if record is not None and not warning and 'J_rad' in record:
                                for name,cut in THRESHOLDS.items():
                                    if qualifies(record,cut):
                                        previous=state['best'].get(name)
                                        if previous is None or (record['J_rad'],-lp,record['id'])<(previous['J_rad'],-previous['log_probability'],previous['id']):
                                            state['best'][name]=record
                                if baseline:state['baseline']=record
                    # Full cell annotations for reference-covered galaxies; other galaxies remain denominator-ledger/raw-diagnostic entries.
                    if g in references or tier==VERIFIED:
                        rp=None if record is None else record.get('common_reference_rank')
                        lp=row['log_probability']
                        writer.writerow(dict(galaxy=g,family=item['family'],grid_i=key[0],grid_j=key[1],
                            log_probability=lp,raw_P_display='' if record is None else math.exp(record['log_probability']),
                            J_rad='' if record is None else record.get('J_rad',''),q_common_ref='' if rp is None else rp['q'],
                            q_count_inclusive='' if rp is None else rp['count'],tier=tier,native_valid=row['valid']=='True',export_valid=valid,
                            reason=row.get('invalid_reason') or row.get('warnings') or ''))
            if legacy.sha(path)!=digest:raise ValueError('native CSV changed during export')
            for g,(item,entry) in targets.items():
                state=states[g]
                pool=list({r['id']:r for r in list(state['best'].values())+([state['baseline']] if state['baseline'] else [])}.values())
                # Every cell has been examined; only the union of new-cut minima/baseline needs profile retention.
                entry['complete']=entry['observed_cells']==entry['expected_cells']
                entry['grid_sha256']=digest
                entry['selections']=select_all(pool,'J_rad',complete=entry['complete'],
                    tainted=bool(entry['numerical_failures'] or entry['warning_tie_exceptions'] or entry['reconstruction_failures'] or entry['duplicate_cells']),
                    allow_main=cell_tier(entry['evidence'])==VERIFIED)
                entry['reference_B']=references[g]['B'] if g in references else None
                entry['scan_scope']='all relevant saved marginalised cells scanned; no old-threshold prefilter'
                entry.pop('records')
                summaries.append(entry)
    return summaries


def gp_export(path,manifest,reference,out):
    native=legacy.read(path)
    if native.get('basis_warnings'):
        raise ValueError('native GP basis warnings; compatibility withheld')
    checkpoint_path=Path(native['checkpoint_directory'])/'manifest.json'
    checkpoint=legacy.read(checkpoint_path)
    for key in ['script_sha256','source_sha256','input_sha256']:
        if checkpoint.get(key)!=native['provenance'].get(key):
            raise ValueError('checkpoint/native provenance mismatch: '+key)
    # Existing independent adapter verifies source/input/parser/normalization/masks.
    document=legacy.gp_document(path,manifest)
    scores=validate_reference(reference)
    records=[]
    for row in document['records']:
        rec=dict(row,log_probability=float(row['log_probability']),native_valid=True,
                 common_reference_rank=rank(float(row['log_probability']),scores),
                 is_baseline=row['id']==native['checkpoint_storage']['summary']['baseline_record_index'],tier=VERIFIED)
        physical=rec['physical_summary']
        physical['absolute_rms_departure_ml']=physical.pop('rms_departure_ml')
        physical['fractional_rms_departure']=math.sqrt(physical['J_rad'])
        records.append(rec)
    import common_reference_driver as driver
    # Strict old GP proof is still checked; its minimum-q05 selector is NOT reused.
    driver.annotate_document(reference,document)
    suff=frontier_sufficiency(native,records)
    summary=native['checkpoint_storage']['summary']
    selection=select_all(records,'J_ref',complete=native['status']=='finite_search_completed',allow_main=suff['sufficient'])
    result={'schema':'retained-selection-view-v2','galaxy':native['galaxy'],
        'thresholds':THRESHOLDS,'sufficiency':suff,'selections':selection,'reference_B':reference['B'],
        'unchanged_metric':document['deviation_metric'],'retained_records_with_profiles':records,
        'full_frontier_ids':suff['frontier_ids'],'baseline_id':suff['baseline_id'],
        'maximum_probability_id':suff['maximum_probability_id'],'original_denominators':summary,
        'native_admissibility':native['physical_domain'],
        'reference_availability':{'B':reference['B'],'context_sha256':reference['context_sha256'],
            'baseline_floor_is_availability_QC':True},
        'provenance':dict(document['provenance'],checkpoint_manifest_sha256=legacy.sha(checkpoint_path),reference_sha256=legacy.sha(ROOT/'outputs/corrected_recovery_pilots_20261016/reference'/(native['galaxy']+'_B64.json'))),'diagnostic_only':True,
        'interpretation':'B64 and finite pilot; selected minimum FOUND, not a global optimum or calibrated significance'}
    legacy.write(out/(native['galaxy']+'_gp_retained_selection_v2.json'),result)
    return result


def output_summary(document,out):
    fields=['strand','galaxy','family','threshold','status','winner_id','raw_logP','raw_P_display',
            'q_common_ref','metric_id','J','reference_B','diagnostic_only']
    with (out/'selection_summary.csv').open('x',newline='') as stream:
        w=csv.DictWriter(stream,fieldnames=fields);w.writeheader()
        for strand,items in [('grid',document['grids']),('gp',document['gp'])]:
            for item in items:
                for name,selection in item.get('selections',{}).items():
                    winner=selection['winner'];rp=None if winner is None else winner.get('common_reference_rank')
                    metric=selection['selection_metric']
                    w.writerow(dict(strand=strand,galaxy=item['galaxy'],family=item.get('family','GP'),threshold=name,
                        status=selection['status'],winner_id='' if winner is None else winner['id'],
                        raw_logP='' if winner is None else winner['log_probability'],
                        raw_P_display='' if winner is None else math.exp(winner['log_probability']),
                        q_common_ref='' if rp is None else rp['q'],metric_id=metric,
                        J='' if winner is None else winner[metric],reference_B=item.get('reference_B',''),diagnostic_only=True))


def run(args):
    started=time.monotonic()
    if legacy.sha(ROOT/'saved_profile_comparison.py')!=LEGACY_SHA:
        raise ValueError('legacy helper changed; review rather than importing unpinned code')
    manifest=ROOT/'outputs/corrected_recovery_prep_20261016/manifest.json'
    frozen=legacy.read(manifest)
    legacy.source_guard(frozen['source_sha256']);legacy.source_guard(frozen['input_sha256'])
    out=args.output_root;out.mkdir(parents=True,exist_ok=False)
    refroot=ROOT/'outputs/corrected_recovery_pilots_20261016/reference'
    # Name extraction is explicit; never infer galaxy from a title or CSV row number.
    references={p.name.removesuffix('_B64.json'):legacy.read(p) for p in refroot.glob('*_B64.json')}
    items=[]
    for family in ['constant','linear','log_linear']:
        items+=legacy.choose_sources(frozen['normalization'],ROOT/'outputs/rerun_boundary_20261006',
                                     ROOT/'outputs/bulge_normalized_20261015',family)
    document={'schema':'saved-score-selection-v2','created_utc':datetime.now(timezone.utc).isoformat(),
        'thresholds':THRESHOLDS,'source_sha256':{Path(__file__).name:legacy.sha(__file__),
            'saved_profile_comparison.py':LEGACY_SHA},'grids':grid_export(items,references,out,frozen['normalization']),
        'gp':[],'gp_failures':[],'diagnostic_only':True,
        'saved_references':{g:{'path':str(refroot/(g+'_B64.json')),'sha256':legacy.sha(refroot/(g+'_B64.json')),'B':ref['B'],'context_sha256':ref['context_sha256']} for g,ref in references.items()},
        'integral_convention':'grid J_rad: full original density span uniform continuous dR; R0 included',
        'rms_convention':'absolute ML RMS=.5sqrt(J_rad); fractional RMS=sqrt(J_rad)',
        'small_slope_integral':'8-point Gauss-Legendre numerical approximation, not analytic exactness',
        'main_grid_rule':'only independently evidenced score-compatible, warning-free complete galaxies',
        'reference_rule':'one unchanged inclusive empirical CDF; no interpolated quantiles',
        'native_domains':'signed finite grid gbar and GP positive-floor domain are distinct and unchanged',
        'selection_scope':'minimum found under unchanged family metric; no GP frontier reranking by J_rad',
        'production_precision':'not established; all real saved CDFs here B64',
        'artifact_retention':'all-cell scan ledger; only new-threshold minima/baseline/maxP/full same-J GP frontier profiles retained'}
    for p in sorted((ROOT/'outputs/corrected_recovery_pilots_20261016/gp').glob('*/*_gp_search.json')):
        g=p.parent.name
        try:
            document['gp'].append(gp_export(p,manifest,references[g],out))
        except (ValueError,KeyError,FileNotFoundError) as error:
            document['gp_failures'].append({'galaxy':g,'status':WITHHELD,'reason':str(error)})
    legacy.source_guard(frozen['source_sha256']);legacy.source_guard(frozen['input_sha256'])
    document['elapsed_seconds']=time.monotonic()-started
    legacy.write(out/'selection_evidence_v2.json',document)
    output_summary(document,out)
    print('SAVED',out,'seconds',document['elapsed_seconds'])
    print('GRID TIERS',dict(Counter(cell_tier(x['evidence']) for x in document['grids'])))
    print('GP',len(document['gp']),'WITHHELD',document['gp_failures'])


class Tests(unittest.TestCase):
    def evidence(self,missing=(),mismatch=()):return {'missing':list(missing),'mismatches':list(mismatch)}
    def row(self,i,j,lp,baseline=False,tier=VERIFIED):
        return {'id':i,'J_rad':j,'J_ref':j,'log_probability':lp,'native_valid':True,'warnings':[],
                'is_baseline':baseline,'tier':tier,'common_reference_rank':rank(lp,[-9,-8,-7,-6])}
    def test_exact_raw_threshold(self):
        self.assertAlmostEqual(RAW_P,.0013498980316301,places=17)
        self.assertTrue(qualifies(self.row(0,0,RAW_LOGP),THRESHOLDS['raw_Phi_minus3']))
        self.assertFalse(qualifies(self.row(0,0,math.nextafter(RAW_LOGP,-math.inf)),THRESHOLDS['raw_Phi_minus3']))
    def test_inclusive_rank(self):
        self.assertEqual(rank(-3,[-4,-3,-3,-1]),{'count':3,'B':4,'q':.75})
    def test_integer_rank_thresholds(self):
        r=self.row(0,0,-9)
        for name,num in [('common_q05',4),('common_q10',7)]:
            r['common_reference_rank']={'count':num,'B':64,'q':num/64}
            self.assertTrue(qualifies(r,THRESHOLDS[name]))
            r['common_reference_rank']['count']-=1
            self.assertFalse(qualifies(r,THRESHOLDS[name]))
    def test_tiers(self):
        self.assertEqual(cell_tier(self.evidence()),VERIFIED)
        self.assertEqual(cell_tier(self.evidence(['solver'])),WITHHELD)
        self.assertEqual(cell_tier(self.evidence(mismatch=['P algorithm'])),MISMATCH)
        self.assertEqual(cell_tier(self.evidence(),False),INVALID)
        self.assertEqual(cell_tier(self.evidence(),warnings=True),WARNING)
    def test_native_signed_is_not_gp_gate(self):
        e=dict(self.evidence(),native_admissibility={'net_gbar':'signed finite; no GP floor'})
        self.assertEqual(cell_tier(e),VERIFIED)
    def test_all_cell_reselection(self):
        old_raw_winner=self.row(1,10,-2)
        previously_not_raw_winner=self.row(2,1,-7)
        result=select_all([old_raw_winner,previously_not_raw_winner],'J_rad')
        self.assertEqual(result['common_q05']['winner']['id'],2)
        self.assertEqual(result['raw_Phi_minus3']['winner']['id'],1)
    def test_incomplete_galaxy(self):
        r=select_all([self.row(1,1,-1)],'J_rad',complete=False)
        self.assertIsNone(r['common_q05']['winner'])
        self.assertIsNotNone(r['common_q05']['diagnostic_finite_pool_winner'])
    def test_tainted_galaxy(self):
        self.assertEqual(select_all([self.row(1,1,-1)],'J_rad',tainted=True)['common_q05']['status'],'failed_or_warning_tainted_galaxy')
    def test_baseline_label(self):
        self.assertEqual(select_all([self.row(0,0,-1,True)],'J_rad')['common_q05']['status'],'baseline_already_acceptable')
    def test_no_found(self):
        self.assertEqual(select_all([],'J_rad')['raw_Phi_minus3']['status'],'no_qualifying_evaluated_candidate')
    def test_withheld_evidence(self):
        self.assertIsNone(select_all([self.row(1,0,-1)],'J_rad',allow_main=False)['common_q05']['winner'])
    def test_floor_not_score_invariant(self):
        inv={'galaxy':'x','parser_id':'p','input_sha256':{},'retained_mask':[],
            'kinematics_sha256':'k','force_model':{},'normalization_id':'n','nuisance_model':{},
            'probability_sha256':'s','probability_settings':{'floor_relative':1e-8,'ngrid':400}}
        other=copy.deepcopy(inv);other['probability_settings']['floor_relative']=1e-4
        self.assertEqual(score_contract(inv),score_contract(other))
    def test_rms_units(self):
        r=attach_metrics({},[0,1],[1,1],'constant')
        self.assertEqual(r['physical_summary']['absolute_rms_departure_ml'],.5)
        self.assertEqual(r['physical_summary']['fractional_rms_departure'],1.)
    def native_fixture(self):
        records=[self.row(0,0,-10,True),self.row(1,1,-5)]
        for r in records:r['profile_ml']=[.5+r['J_ref'],.5+r['J_ref']]
        native={'status':'finite_search_completed','records':[{},{}],
            'physical_domain':{'expected_evaluated_candidates':10},
            'retained_profiles':[{'record_index':r['id'],'profile_on_photometric_grid':r['profile_ml']} for r in records],
            'frontier':[{'record_index':r['id'],'profile_on_photometric_grid':r['profile_ml']} for r in records],
            'checkpoint_storage':{'mode':'retained-only',
                'discard_policy':'exact dominance; baseline, threshold minima and maximum P pinned',
                'summary':{'evaluated':10,'valid':10,'invalid':0,'retained':2,'discarded':8,
                    'frontier_size':2,'baseline_record_index':0,'maximum_probability_record_index':1}}}
        return native,records
    def test_full_frontier_sufficient(self):
        n,r=self.native_fixture();self.assertTrue(frontier_sufficiency(n,r)['sufficient'])
    def test_old_winners_insufficient(self):
        n,r=self.native_fixture();n['checkpoint_storage'].pop('discard_policy')
        self.assertFalse(frontier_sufficiency(n,r)['sufficient'])
    def test_missing_profile(self):
        n,r=self.native_fixture();self.assertFalse(frontier_sufficiency(n,r[:1])['sufficient'])
    def test_removed_dominator(self):
        n,r=self.native_fixture();self.assertFalse(frontier_sufficiency(n,r,removed_ids=[1])['sufficient'])
    def test_metric_change(self):
        n,r=self.native_fixture();self.assertFalse(frontier_sufficiency(n,r,metric_unchanged=False)['sufficient'])
    def test_numerical_failed_search(self):
        n,r=self.native_fixture()
        n['checkpoint_storage']['summary']['failure_counts']={'probability_failure':1}
        self.assertFalse(frontier_sufficiency(n,r)['sufficient'])
    def test_mixed_evidence_withheld(self):
        a=self.row(1,2,-1);b=self.row(2,1,-1,tier=WITHHELD)
        self.assertIsNone(select_all([a,b],'J_rad')['common_q05']['winner'])
    def test_frontier_profile_mismatch(self):
        n,r=self.native_fixture();n['frontier'][1]['profile_on_photometric_grid']=[9,9]
        self.assertFalse(frontier_sufficiency(n,r)['sufficient'])
    def test_missing_inputs_withheld_not_invalid(self):
        self.assertEqual(cell_tier(self.evidence(['historical rotmod'])),WITHHELD)
        self.assertEqual(cell_tier(self.evidence(['historical vcdisk'])),WITHHELD)
    def collector_fixture(self, expected=2, extra=None):
        import tempfile
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);dens=root/'X.dens'
            with dens.open('w') as stream:stream.write('0 10 0\n1 1 0\n')
            meta=root/'meta.json';legacy.write(meta,{'arguments':{'dens_root':str(root)}})
            path=root/'grid.csv'
            rows=[dict(galaxy='X',profile='constant',probability_kind='marginalised_gas_bulge',grid_i='0',grid_j=str(i),upsilon_disk_mean=str(ml),disk_slope='nan',log_probability=str(lp),valid='True',invalid_reason='',warnings='') for i,ml,lp in [(0,.5,-7),(1,1.,-2)]]
            if extra is not None:rows.append(dict(rows[1],grid_j='2',**extra))
            with path.open('w',newline='') as stream:
                writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
            item={'galaxy':'X','family':'constant','settings':{'grid_size':str(expected)},'source':{'grid':str(path),'metadata':str(meta)}}
            with patch.dict(grid_export.__globals__,grid_evidence=lambda *a:self.evidence(),validate_reference=lambda *a:[-9,-8,-7,-6]):
                result=grid_export([item],{'X':{'B':4}},root,root/'norm.json')
            with (root/'grid_cell_annotations.csv').open() as stream:annotations=list(csv.DictReader(stream))
            return result[0],annotations
    def test_collector_scans_nonraw_winner(self):
        result,rows=self.collector_fixture()
        self.assertEqual(len(rows),2)
        self.assertEqual(result['selections']['common_q05']['winner']['id'],'0:0')
        self.assertEqual(result['selections']['raw_Phi_minus3']['winner']['id'],'0:1')
    def test_collector_incomplete_withheld(self):
        result,rows=self.collector_fixture(expected=3)
        self.assertEqual(result['observed_cells'],2)
        self.assertIsNone(result['selections']['common_q05']['winner'])
    def test_collector_native_invalid_denominator(self):
        result,rows=self.collector_fixture(expected=3,extra={'valid':'False','invalid_reason':'nonpositive_or_nonfinite_disk_ml'})
        self.assertEqual(result['observed_cells'],3)
        self.assertEqual(result['native_invalid'],1)
        self.assertEqual(result['numerical_failures'],0)
        self.assertEqual(rows[-1]['tier'],INVALID)
        self.assertIsNotNone(result['selections']['common_q05']['winner'])
    def test_collector_numerical_failure_taints(self):
        result,rows=self.collector_fixture(expected=3,extra={'valid':'False','invalid_reason':'nonfinite_probability'})
        self.assertEqual(result['numerical_failures'],1)
        self.assertIsNone(result['selections']['common_q05']['winner'])
    def test_collector_tie_warning_taints(self):
        result,rows=self.collector_fixture(expected=3,extra={'warnings':'tie cap exceeded; stable descending order'})
        self.assertEqual(result['warning_tie_exceptions'],1)
        self.assertEqual(rows[-1]['tier'],WARNING)
        self.assertIsNone(result['selections']['common_q05']['winner'])
    def test_bound_manifest_verified_chain(self):
        import tempfile
        with tempfile.TemporaryDirectory() as t:
            launcher=ROOT/'run_bulge_reruns_20261015.py'
            document={'source_sha256':{str(launcher):legacy.sha(launcher)},'input_sha256':{}}
            path=Path(t)/'batch_manifest.json';legacy.write(path,document)
            metadata={'normalization':{'manifest_sha256':legacy.sha(path),'launcher_sha256':legacy.sha(launcher)}}
            saved,proof=bound_manifest(metadata,Path(t))
            self.assertEqual(proof['status'],'verified_historical_binding')
            self.assertEqual(saved,document)
    def test_grid_evidence_missing_numerical(self):
        import tempfile
        with tempfile.TemporaryDirectory() as t:
            path=Path(t)/'meta.json'
            legacy.write(path,{'arguments':{'numerical_source_sha256':'{}'}})
            e=grid_evidence({'source':{'metadata':str(path)}},None,Path(t)/'norm')
            self.assertEqual(cell_tier(e),WITHHELD)
    def test_wrong_loader_cannot_be_rescued(self):
        bad={'provenance':{'source_sha256':{},'input_parser':{'id':'wrong_gas4'}}}
        fake_manifest={'source_sha256':{},'input_sha256':{}}
        with patch.object(legacy,'read',side_effect=[bad,fake_manifest]):
            with self.assertRaisesRegex(ValueError,'parser evidence'):
                legacy.gp_document('native','manifest')
    def evidence_fixture(self, mode='complete'):
        import tempfile
        import common_isotonic_reference as common
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);family=root/'constant';family.mkdir()
            files={'rotmod':root/'X_rotmod.dat','dens':root/'X.dens','table':root/'table.txt','normalization':root/'normalization.json'}
            for key,path in files.items():
                with path.open('w') as stream:stream.write('synthetic '+key)
            solver=root/'archived_vcdisk.py'
            with solver.open('w') as stream:stream.write('# synthetic identity fixture, never imported')
            launcher=ROOT/'run_bulge_reruns_20261015.py'
            rawhash={key:legacy.sha(path) for key,path in files.items()}
            numeric={key:legacy.sha(ROOT/key) for key in ['mono_SPARC_RAR.py','mono_SPARC_RAR_boundary.py','probability_recursion.py','utils.py']}
            args={'numerical_source_sha256':json.dumps(numeric),'data_root':str(root),'dens_root':str(root),'table_path':str(files['table']),
                'gas_ml':'1','gas_ml_dex':'.04','bulge_ml':'.7','bulge_ml_dex':'.1','ml_samples':'9',
                'prob_grid_size':'400','zpadding':'8','boundary_log_tolerance':'2e-5','boundary_max_refinements':'8','max_tie_permutations':'720','tie_atol':'0'}
            inv={'galaxy':'X','parser_id':legacy.PARSER,'input_sha256':dict(rawhash),'retained_mask':[True,True],'kinematics_sha256':'synthetic arrays',
                'normalization_id':legacy.POLICY,'nuisance_model':{'gas_median':1.,'gas_dex':.04,'bulge_median':.7,'bulge_dex':.1,'samples':9},
                'probability_sha256':numeric['probability_recursion.py'],
                'probability_settings':{'ngrid':400,'zpadding':8.,'log_tolerance':2e-5,'max_refinements':8,'max_tie_permutations':720,'floor_relative':1e-8},
                'force_model':{'force_source_sha256':{'mono_SPARC_RAR.py':numeric['mono_SPARC_RAR.py'],'utils.py':numeric['utils.py'],'vcdisk':legacy.sha(solver)},
                    'geometry':{},'interpolation':'signed W/R','signed_W':'V*abs(V)','units':{},'kpc_metres':3.086e19,'density_support':'full','kinematic_support':'used points inside force support'}}
            manifest={'source_sha256':{str(launcher):legacy.sha(launcher),str(solver):legacy.sha(solver)},
                'input_sha256':{str(files[k]):v for k,v in rawhash.items()}}
            if mode=='missing_input':manifest['input_sha256'].pop(str(files['rotmod']))
            if mode=='missing_solver':manifest['source_sha256'].pop(str(solver))
            mpath=root/'batch_manifest.json';legacy.write(mpath,manifest)
            metadata={'arguments':args,'normalization':{'manifest_sha256':legacy.sha(mpath),'launcher_sha256':legacy.sha(launcher)}}
            meta=family/'run_metadata.json';legacy.write(meta,metadata)
            target=copy.deepcopy(inv)
            if mode=='score_mismatch':target['probability_settings']['ngrid']=800
            context={'invariants':target,'input_sha256':target['input_sha256'],'probability_settings':target['probability_settings']}
            reference={'status':'complete','failures':[],'B':4,'sorted_mock_logp':[-9,-8,-7,-6],
                'context':context,'context_sha256':common.digest(context)}
            with patch.object(legacy,'independent_inputs',return_value=(inv,np.array([0,1]),np.array([1,1]),1.)):
                return grid_evidence({'galaxy':'X','family':'constant','source':{'metadata':str(meta)}},reference,files['normalization'])
    def test_full_evidence_graph_accepted(self):
        e=self.evidence_fixture()
        self.assertEqual(e['missing'],[])
        self.assertEqual(e['mismatches'],[])
        self.assertEqual(cell_tier(e),VERIFIED)
    def test_graph_missing_historical_input(self):
        e=self.evidence_fixture('missing_input')
        self.assertEqual(cell_tier(e),WITHHELD)
        self.assertEqual(e['mismatches'],[])
    def test_graph_missing_historical_solver(self):
        e=self.evidence_fixture('missing_solver')
        self.assertEqual(cell_tier(e),WITHHELD)
        self.assertEqual(e['mismatches'],[])
    def test_graph_actual_score_mismatch(self):
        self.assertEqual(cell_tier(self.evidence_fixture('score_mismatch')),MISMATCH)
    def test_reference_failures_not_discarded(self):
        with self.assertRaisesRegex(ValueError,'failed mock'):
            validate_reference({'status':'complete','failures':['mock failed']})
    def test_reference_available_when_grid_missing(self):
        import common_isotonic_reference as common
        context={'probability_settings':{'floor_relative':1e-8}}
        ref={'status':'complete','failures':[],'B':4,'sorted_mock_logp':[-9,-8,-7,-6],
             'context':context,'context_sha256':common.digest(context)}
        evidence=grid_evidence({'source':None},ref,Path('unused'))
        self.assertEqual(evidence['reference_availability']['status'],'complete_saved_CDF')
        self.assertEqual(cell_tier(evidence),WITHHELD)
    def test_summary_reference_B_not_hardcoded(self):
        import tempfile
        with tempfile.TemporaryDirectory() as t:
            item={'galaxy':'X','reference_B':4,'selections':select_all([self.row(0,0,-1,True)],'J_ref')}
            output_summary({'grids':[],'gp':[item]},Path(t))
            with (Path(t)/'selection_summary.csv').open() as stream:
                rows=list(csv.DictReader(stream))
            self.assertTrue(all(row['reference_B']=='4' for row in rows))
    def test_missing_manifest_binding(self):
        d,p=bound_manifest({},Path('/no_such_historical_manifest'))
        self.assertIsNone(d);self.assertEqual(p['status'],'missing')
    def test_tampered_manifest(self):
        import tempfile
        with tempfile.TemporaryDirectory() as t:
            p=Path(t)/'batch_manifest.json';legacy.write(p,{})
            d,proof=bound_manifest({'normalization':{'manifest_sha256':'wrong'}},Path(t))
            self.assertIsNone(d);self.assertEqual(proof['status'],'mismatch')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-root',type=Path,required=True)
    p.add_argument('--self-test',action='store_true')
    args=p.parse_args()
    if args.self_test:
        suite=unittest.TestSuite([unittest.defaultTestLoader.loadTestsFromTestCase(legacy.Tests),
            unittest.defaultTestLoader.loadTestsFromTestCase(Tests)])
        result=unittest.TextTestRunner(verbosity=2).run(suite)
        legacy.write(args.output_root/'test_receipt.json',{'schema':'saved-selection-tests-v2',
            'tests':result.testsRun,'failures':len(result.failures),'errors':len(result.errors),
            'successful':result.wasSuccessful(),'source_sha256':legacy.sha(__file__),'legacy_helper_sha256':legacy.sha(ROOT/'saved_profile_comparison.py')})
        sys.exit(0 if result.wasSuccessful() else 1)
    run(args)
