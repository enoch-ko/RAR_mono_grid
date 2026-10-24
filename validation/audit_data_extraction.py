#!/usr/bin/env python3
"""Read-only input extraction audit; writes NEW audit artifacts only, no fits/forces."""
import argparse
import ast
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import numpy as np
import pandas as pd
import utils

ROOT=Path(__file__).resolve().parent

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def inventory():
    result=[]; notebooks=[]; skips=[]
    excluded={"outputs",".venv-osqp",".git","__pycache__","slurm_logs","plots"}
    for parent,dirs,files in os.walk(ROOT):
        dirs[:]=[d for d in dirs if d not in excluded and not d.startswith(".")]
        for name in sorted(files):
            path=Path(parent)/name
            if path.suffix not in {".py",".ipynb",".sh",".yaml",".yml",".toml",".ini"}:continue
            if path.stat().st_size>2*1024*1024:
                skips.append(str(path));continue
            text=path.read_text(errors="replace")
            hits=[{"line":i,"text":line.strip()} for i,line in enumerate(text.splitlines(),1)
                  if any(token in line for token in ("loadtxt","read_csv","read_fwf","genfromtxt","DictReader","csv.reader","np.load(","json.loads","json.load(",
                    "read_table","galaxy_data(","get_SBdisk(","MHI","helium","raw[","data[keep","kin[:,"))]
            imports=[]
            if path.suffix==".py":
                try:
                    tree=ast.parse(text)
                    imports=[ast.unparse(n) for n in ast.walk(tree) if isinstance(n,(ast.Import,ast.ImportFrom))]
                except SyntaxError as exc:
                    skips.append(str(path)+":syntax:"+str(exc))
            if path.suffix==".ipynb":notebooks.append(str(path.relative_to(ROOT)))
            result.append({"path":str(path.relative_to(ROOT)),"sha256":sha(path),
                           "hits":hits,"imports":imports,"test_or_archive":name.startswith("test_") or "ARCHIVED" in name,
                           "shell_payload":text if path.suffix==".sh" else None})
    return {"files":result,"notebooks":notebooks,"skipped":skips,"excluded_directories":sorted(excluded)}

def master_check(table):
    path=Path(utils.DEFAULT_TABLE_PATH);lines=path.read_text().splitlines();spec=[]
    pattern=re.compile(r"^\s*(\d+)(?:-\s*(\d+))?\s+(\S+)\s+(\S+)\s+(\S+)")
    for line in lines[:98]:
        match=pattern.match(line)
        if match and match[3][0] in "AIF" and match[3][1:2].isdigit():spec.append((int(match[1])-1,int(match[2] or match[1]),match[3],match[4],match[5]))
    names=[{"L[3.6]":"L","e_L[3.6]":"e_L"}.get(item[-1],item[-1]) for item in spec]
    if len(names)!=19 or set(names)!=set(table.columns):raise ValueError("declared schema field names mismatch")
    records=[];mismatches=[]
    for index,line in enumerate(lines[98:],99):
        if not line.strip():continue
        values=line.split(maxsplit=18)
        if len(values)!=len(names):raise ValueError("catalogue token count mismatch")
        if values[0] not in set(table.Galaxy):
            mismatches.append({"line":index,"reason":"unknown galaxy","value":values[0]});continue
        row=table.loc[table.Galaxy==values[0]]
        if len(row)!=1:raise ValueError("duplicate catalogue galaxy")
        row=row.iloc[0]
        for column,value,(_,_,fmt,_,_) in zip(names,values,spec):
            if fmt.startswith("A"):
                equal=str(row[column]).strip()==value
            else:
                equal=np.isclose(float(row[column]),float(value),rtol=0,atol=1e-12)
            if not equal:mismatches.append({"galaxy":values[0],"column":column,"byte_value":value,"reader_value":str(row[column])})
        records.append(values[0])
    return {"path":str(path),"sha256":sha(path),"declared_byte_schema":spec,"comparison_parser":"independent ordered whitespace fields; actual row padding differs from declared bytes", "row_count":len(records),
            "unique_galaxies":len(set(records)),"mismatches":mismatches}

def galaxy_check(galaxy,row):
    rot=Path('/mnt/users/koe/data')/(galaxy+'_rotmod.dat')
    dens=Path(utils.DEFAULT_DENS_ROOT)/(galaxy+'.dens')
    header=next(l.lstrip('# ').split() for l in rot.read_text().splitlines() if l.startswith('#') and 'Vobs' in l)
    expected=['Rad','Vobs','errV','Vgas','Vdisk','Vbul','SBdisk','SBbul']
    if header!=expected:raise ValueError(galaxy+':unexpected rotmod header')
    a=np.loadtxt(rot,ndmin=2);b=np.loadtxt(dens,ndmin=2)
    if a.shape[1]!=8 or b.shape[1]!=3:raise ValueError(galaxy+':column count')
    rad,brightness=utils.get_SBdisk(galaxy)
    expected_brightness=b[:,1]*1e6*float(row.L)*1e9/utils.compute_A_B(b[:,0],b[:,1]*1e6)[0]
    rows_equal=np.array_equal(rad,b[:,0]) and np.allclose(brightness,expected_brightness,rtol=1e-13,atol=0)
    positive=b[b[:,0]>0,0]
    keep=np.isfinite(a[:,0])&(a[:,0]>=b[0,0])&(a[:,0]<=b[-1,0])
    distance=float(rot.read_text().splitlines()[0].split('=')[1].split()[0])
    return {"galaxy":galaxy,"rotmod_sha256":sha(rot),"dens_sha256":sha(dens),
        "density_header":dens.read_text().splitlines()[0],"n_rotmod":len(a),"n_dens":len(b),
        "finite_rotmod":bool(np.isfinite(a).all()),"finite_dens":bool(np.isfinite(b).all()),
        "strictly_increasing_rotmod":bool(np.all(np.diff(a[:,0])>0)),
        "strictly_increasing_dens":bool(np.all(np.diff(b[:,0])>0)),
        "nonpositive_error_rows":int(np.count_nonzero(a[:,2]<=0)),
        "nonpositive_kinematic_radius_rows":int(np.count_nonzero(a[:,0]<=0)),
        "negative_gas_rows":int(np.count_nonzero(a[:,3]<0)),
        "negative_disk_rows":int(np.count_nonzero(a[:,4]<0)),
        "negative_bulge_rows":int(np.count_nonzero(a[:,5]<0)),
        "negative_density_rows":int(np.count_nonzero(np.any(b[:,1:]<0,axis=1))),
        "zero_density_radius_rows":int(np.count_nonzero(b[:,0]==0)),
        "joint_overlap_rows":int(np.count_nonzero(keep)),
        "retained_below_positive_solver_support":int(np.count_nonzero(keep&(a[:,0]<positive[0]))),
        "excluded_rows":np.flatnonzero(~keep).tolist(),
        "density_reader_independently_matches":bool(rows_equal),
        "distance_difference_Mpc":distance-float(row.D),
        "disk_luminosity_ratio_to_catalogue_total":float(utils.compute_A_B(rad,brightness)[0]/(float(row.L)*1e9)),
        "gas_disk_distinct_rows":int(np.count_nonzero(a[:,3]!=a[:,4]))}

def run(output):
    if output.exists():raise FileExistsError('new audit directory required')
    table=utils.load_sparc_table()
    if len(table)!=175 or table.Galaxy.nunique()!=175:raise ValueError('catalogue identity/coverage')
    master=master_check(table)
    galaxies=[galaxy_check(str(row.Galaxy),row) for _,row in table.iterrows()]
    external=ROOT/'ML_Andreea.csv';frame=pd.read_csv(external)
    grouped=frame.groupby(['Sample','Galaxy'],dropna=False)
    ext={"sha256":sha(external),"columns":list(frame.columns),"rows":len(frame),
         "unique_galaxies":int(frame.Galaxy.nunique()),"sample_galaxy_groups":len(grouped),
         "samples":{str(k):int(v) for k,v in frame.Sample.value_counts(dropna=False).items()},
         "missing_identifiers":int(frame[['Sample','Galaxy']].isna().any(axis=1).sum()),
         "nonpositive_radius":int((frame.Radius_kpc<=0).sum()),
         "nonfinite_radius_or_ml":int((~np.isfinite(frame[['Radius_kpc','ml_K']]).all(axis=1)).sum()),
         "nonpositive_ml":int((frame.ml_K<=0).sum()),
         "nonpositive_Rbar":int((frame.R_bar<=0).sum()),
         "inconsistent_Rbar_groups":int(grouped.R_bar.nunique(dropna=False).gt(1).sum()),
         "same_radius_only_groups":[str(k) for k,g in grouped if g.Radius_kpc.nunique()<2],
         "Sample_vs_sample_mismatches":int((frame.Sample!=frame['sample']).sum())}
    totals={key:sum(int(g[key]) for g in galaxies) for key in ['nonpositive_error_rows',
        'nonpositive_kinematic_radius_rows','negative_gas_rows','negative_disk_rows',
        'negative_bulge_rows','negative_density_rows','zero_density_radius_rows',
        'retained_below_positive_solver_support']}
    result={"created_utc":datetime.now(timezone.utc).isoformat(),"script_sha256":sha(__file__),
            "master":master,"galaxies":galaxies,"totals":totals,"external":ext,
            "inventory":inventory(),"scope":"extraction/schema audit only; no fits, force solves or production modifications"}
    output.mkdir(parents=True)
    with (output/'data_extraction_audit.json').open('x') as stream:json.dump(result,stream,indent=2,allow_nan=False)
    print(json.dumps({"master_mismatches":len(master['mismatches']),"galaxies":len(galaxies),
        "totals":totals,"density_matches":sum(g['density_reader_independently_matches'] for g in galaxies),
        "distance_mismatches":sum(abs(g['distance_difference_Mpc'])>.005 for g in galaxies),
        "external":ext,"inventory_files":len(result['inventory']['files']),"skipped":result['inventory']['skipped']},indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    run(parser.parse_args().output)
