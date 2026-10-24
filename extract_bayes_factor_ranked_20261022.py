from pathlib import Path
import csv,json,math,hashlib,datetime,decimal
ROOT=Path('/mnt/users/koe/RAR_mono_grid')
SRC=ROOT/'outputs/bayes_factor_vmax_scale_check_20261022'
CANON=ROOT/'outputs/bayes_factor_full_sample_20261022/summary.csv'
OUT=ROOT/'outputs/bayes_factor_ranked_20261022_v1'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    assert not OUT.exists()
    original=json.loads((SRC/'manifest.json').read_text())
    assert sha(SRC/'vmax_scale_check.csv')==original['files']['vmax_scale_check.csv']
    assert sha(ROOT/'bayes_factor_vmax_scale_check_20261022.py')==original['script_sha256']
    assert sha(ROOT/'bounded_uniform_bayes_factor_20261021.py')==original['helper_sha256']
    raw=list(csv.DictReader((SRC/'vmax_scale_check.csv').open()))
    base=sorted([r for r in raw if float(r['vmax_km_s'])==1.],key=lambda r:float(r['logBF']))
    canon=list(csv.DictReader(CANON.open()));cm={r['galaxy']:r for r in canon}
    assert len(raw)==700 and len(base)==175 and len({r['galaxy'] for r in base})==175
    assert len(canon)==len(cm)==175 and set(cm)=={r['galaxy'] for r in base}
    result=[]
    decimal.getcontext().prec=17
    for i,r in enumerate(base,1):
        q,c,bf=map(float,[r['saved_logQ_fixed'],r['logc'],r['logBF']])
        assert all(math.isfinite(x) for x in [q,c,bf]) and abs(q-c-bf)<1e-10
        assert int(r['rank_by_vmax1_logBF'])==i
        cc=cm[r['galaxy']]
        assert abs(float(cc['log_BF'])-bf)<1e-10 and abs(float(cc['saved_conditional_logP'])-q)<1e-10 and abs(float(cc['log_c_32001'])-c)<1e-10
        result.append(dict(rank=i,galaxy=r['galaxy'],log10BF=bf/math.log(10),lnBF=bf,Q=str(decimal.Decimal(r['saved_logQ_fixed']).exp()),c=str(decimal.Decimal(r['logc']).exp()),lnQ=q,lnc=c,vmax_km_s=1.0))
    below=max((r for r in result if r['lnBF']<0),key=lambda r:r['lnBF'])
    above=min((r for r in result if r['lnBF']>0),key=lambda r:r['lnBF'])
    selections=dict(bottom5=result[:5],median5=result[85:90],top5=result[-5:],nearest_BF1_below=below,nearest_BF1_above=above)
    OUT.mkdir()
    with (OUT/'ranked_galaxies.csv').open('x',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(result[0]));w.writeheader();w.writerows(result)
    lines=['# Galaxies ranked by increasing log10 BF','', '175 unique galaxies; saved Vmax = 1 km/s denominator and unchanged saved Pmono numerator approximation.', '', 'Natural-log source fields: logBF = saved_logQ_fixed - logc; log10 BF = ln BF / ln(10). Q and c below are exponentiated saved natural logs, not new probability calculations.', '', 'The ratio is a numerical diagnostic; no new posterior numerator, scientific validation, force calculation, or search was performed. The positive/symmetric denominator identity does not establish matching posterior numerators.', '', '| Rank | Galaxy | log10 BF | ln BF | Q (saved numerator) | c (denominator) |','|---:|:---|---:|---:|---:|---:|']
    for r in result:lines.append('| {rank} | {galaxy} | {log10BF:.10g} | {lnBF:.10g} | {Q} | {c} |'.format(**r))
    lines+=['','Median: rank 88. Candidate examples are numerical rank selections only.','', 'Source: '+str(SRC/'vmax_scale_check.csv'),'Source SHA256: '+sha(SRC/'vmax_scale_check.csv'),'Cross-checked all 175 rows against: '+str(CANON)]
    (OUT/'ranked_galaxies.md').write_text(chr(10).join(lines)+chr(10))
    (OUT/'example_selections.json').write_text(json.dumps(selections,indent=2)+chr(10))
    check=list(csv.DictReader((OUT/'ranked_galaxies.csv').open()))
    assert len(check)==175 and [r['galaxy'] for r in check]==[r['galaxy'] for r in result]
    manifest=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),source_csv=str(SRC/'vmax_scale_check.csv'),source_csv_sha256=sha(SRC/'vmax_scale_check.csv'),source_manifest_sha256=sha(SRC/'manifest.json'),canonical_summary=str(CANON),canonical_summary_sha256=sha(CANON),rows=175,sort='ascending saved logBF at Vmax=1; verified existing ranks 1 through 175',source_log_base='natural',reported_log10_conversion='lnBF / ln(10)',canonical_all175_match=True,source_original_hashes_match=True,probability_evaluations=0,numerator='fixed saved conditional Pmono approximation',claim_limit='not a new positive-support posterior BF or scientific validation',files={p.name:sha(p) for p in OUT.iterdir() if p.is_file()})
    (OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+chr(10))
    print('VERIFIED_175_RANKS_AND_CANONICAL_MATCH',OUT)
    for label,rs in selections.items():
        print(label)
        for r in (rs if isinstance(rs,list) else [rs]):print(r['rank'],r['galaxy'],format(r['log10BF'],'.10g'))
if __name__=='__main__':main()
