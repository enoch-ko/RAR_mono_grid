from pathlib import Path
import csv,json,math,time,datetime,collections
R=Path('/mnt/users/koe/RAR_mono_grid');P=R/'bf_accuracy_scoped_20261024_v1';assert not P.exists();P.mkdir();OUT=R/'outputs/bf_grid_rollout_20261023_v1';CUTS=(1.,1.5,2.);BANDS=(1e-5,1e-4,1e-3,1e-2);counts={str(c):{str(b):0 for b in BANDS} for c in CUTS};problems={str(c):collections.Counter() for c in CUTS};nearest={};selection={};smallgap=[];minimum=None;stats=collections.Counter();start=time.perf_counter()
def candidate(p,row,family,galaxy,n,reason,cut):
 k='|'.join([family,galaxy,row['grid_i'],row['grid_j']]);v=dict(candidate_id=k,family=family,galaxy=galaxy,grid_i=row['grid_i'],grid_j=row['grid_j'],context_path=str(p/'context.json'),force_bases_path=str(p/'force_bases.npz'),cells_path=str(p/'cells.csv'),row=row,retained_bins=n,distance_to_cut=abs(row['log10_BF']-cut),selection_reason=reason,cut=cut);return v
for cp in sorted(OUT.glob('*/*/cells.csv')):
 p=cp.parent;summary=json.loads((p/'summary.json').read_text());context=json.loads((p/'context.json').read_text());family=summary['family'];galaxy=summary['galaxy'];n=len(context['radii_kpc']);rows=[];stats['galaxy_family_results']+=1
 with cp.open() as f:
  for row in csv.DictReader(f):
   if row['status'].startswith('excluded'):stats['physical_exclusions']+=1;continue
   assert row['status']=='accuracy_unverified';row['log10_BF']=float(row['log10_BF']);row['J_rad']=float(row['J_rad']);row['log_cbar']=float(row['log_cbar']);row['minimum_gap_over_roundoff_screen']=float(row['minimum_gap_over_roundoff_screen']);rows.append(row);stats['valid_cells']+=1;gap=row['minimum_gap_over_roundoff_screen']
   v=candidate(p,row,family,galaxy,n,'minimum_order_gap',0.)
   if minimum is None or gap<minimum['row']['minimum_gap_over_roundoff_screen']:minimum=v
   for level in (1.,10.,100.,1000.,1e6):stats['order_gap_ratio_le_'+str(level)]+=gap<=level
   for cut in CUTS:
    delta=abs(row['log10_BF']-cut)
    for b in BANDS:counts[str(cut)][str(b)]+=delta<=b
    key=(family,cut)
    if key not in nearest or delta<nearest[key]['distance_to_cut']:nearest[key]=candidate(p,row,family,galaxy,n,'nearest_cell_for_family_cut',cut)
 for cut in CUTS:
  key=str(cut);winner=summary['point_cut_minimum_J_found'].get(key)
  # JSON keys retain 1,1.5,2 rather than float strings 1.0/2.0.
  if winner is None:winner=summary['point_cut_minimum_J_found'].get(str(int(cut)) if cut.is_integer() else key)
  if winner:
   problems[key]['point_winner_exists']+=1;wr=next(x for x in rows if x['grid_i']==winner['grid_i'] and x['grid_j']==winner['grid_j']);cand=[('minimum_J_point_winner',wr)]
   lower=[x for x in rows if x['J_rad']<wr['J_rad'] and x['log10_BF']<=cut]
   if lower:cand.append(('closest_cut_strictly_lower_J_failing_competitor',max(lower,key=lambda z:z['log10_BF'])))
   higher=[x for x in rows if x['J_rad']>wr['J_rad'] and x['log10_BF']>cut]
   if higher:cand.append(('next_higher_J_passing_alternative',min(higher,key=lambda z:z['J_rad'])))
  else:
   problems[key]['no_point_passing_grid_cell']+=1;cand=[('maximum_BF_when_no_pass',max(rows,key=lambda z:z['log10_BF']))]
  for reason,row in cand:
   v=candidate(p,row,family,galaxy,n,reason,cut);k=v['candidate_id'];entry=selection.setdefault(k,dict(**{k:z for k,z in v.items() if k not in ['selection_reason','cut','distance_to_cut']},reasons=[]));entry['reasons'].append(dict(reason=reason,cut=cut,distance_to_cut=v['distance_to_cut']))
   if reason in ['minimum_J_point_winner','closest_cut_strictly_lower_J_failing_competitor']:
    for band in BANDS:problems[key][reason+'_within_'+str(band)]+=abs(row['log10_BF']-cut)<=band
primary=[];used=set()
def pick(v):
 if v['candidate_id'] not in used:primary.append(v);used.add(v['candidate_id'])
# Cover each family/cut's nearest saved cell, then sensitive winner/competitor cases.
for v in sorted(nearest.values(),key=lambda z:(z['distance_to_cut'],z['retained_bins'])):pick(v)
critical=[]
for v in selection.values():
 for reason in v['reasons']:
  if reason['reason'] in ['minimum_J_point_winner','closest_cut_strictly_lower_J_failing_competitor']:critical.append(dict(**v,selection_reason=reason['reason'],cut=reason['cut'],distance_to_cut=reason['distance_to_cut']))
for v in sorted(critical,key=lambda z:(z['distance_to_cut'],z['retained_bins'])):
 if len(primary)>=12:break
 pick(v)
report=dict(status='saved_results_screen_only_no_probability_force_evaluation',created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),statistics=dict(stats),cell_distance_counts=counts,cut_problem_counts={k:dict(v) for k,v in problems.items()},nearest_by_family_cut=list(nearest.values()),selected_winner_competitor_candidate_count=len(selection),primary_precision_candidates=primary,minimum_saved_order_gap=minimum,selection_bands_are_not_error_bounds=True,screen_seconds=time.perf_counter()-start,force_solves=0,numerator_evaluations=0,denominator_evaluations=0)
(P/'screen.json').write_text(json.dumps(report,indent=2)+chr(10));(P/'selected_winners_competitors.json').write_text(json.dumps(list(selection.values()),indent=2)+chr(10));print('SCREEN_DONE',dict(stats),'SELECTED',len(selection),'PRIMARY',len(primary),'SECONDS',report['screen_seconds']);print('DISTANCE_COUNTS',counts);print('PRIMARY',[(v['candidate_id'],v['cut'],v['distance_to_cut'],v['retained_bins']) for v in primary]);print('MIN_ORDER_GAP',minimum['candidate_id'],minimum['row']['minimum_gap_over_roundoff_screen'])
