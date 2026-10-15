#!/usr/bin/env python3
"""Controlled fixed-profile monotonic nulls; NOT GP-selection calibration.
Independent Gaussian velocity noise; fixed geometry; nested radial subsets.
A phenomenological increasing latent RAR is used, not a fitted MOND law.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import numpy as np
from scipy.special import logsumexp
from utils import vel2acc, lognormal_quadrature, _tie_groups
from probability_recursion import boundary_log_probability_for_gbar
from sparc_scale_proxies import sha256_file
ROOT = Path(__file__).resolve().parent


def nested_indices(size):
    """Endpoints first, then farthest index: every N is a prefix of one design."""
    if size < 2: raise ValueError("parent size must be >=2")
    chosen = [0, size-1]
    while len(chosen) < size:
        remaining = [i for i in range(size) if i not in chosen]
        chosen.append(max(remaining, key=lambda i: (min(abs(i-j) for j in chosen), -i)))
    return chosen


def latent_design(size, fractional_error, a0=1.2e-10):
    if not math.isfinite(fractional_error) or fractional_error <= 0:
        raise ValueError("positive fractional velocity error required")
    r = np.geomspace(1., 20., size)
    gbar = np.geomspace(3e-10, 3e-12, size)
    gtrue = gbar+np.sqrt(a0*gbar)  # strictly increasing function of gbar
    vtrue = np.sqrt(gtrue/vel2acc(np.ones(size), r))
    return r, .7*gbar, .25*gbar, .05*gbar/.7, vtrue, fractional_error*vtrue


def fixed_probability(disk, gas_unit, bulge_unit, observed, errors, radius,
                      samples=9, ngrid=400, log_tolerance=2e-5, max_refinements=8):
    """Same gas/bulge lognormal nodes and ordered velocity-error integration.
    Cache identical tie/order groups, since gbar enters this statistic via order.
    No nuisance-node failure is dropped or reweighted.
    """
    arrays = [np.asarray(x, float) for x in (disk, gas_unit, bulge_unit, observed, errors, radius)]
    if any(a.ndim != 1 for a in arrays) or len({a.shape for a in arrays}) != 1:
        raise ValueError("fixed components and kinematics must share a 1D geometry")
    if not all(np.all(np.isfinite(a)) for a in arrays) or not len(arrays[0]):
        raise ValueError("nonempty finite fixed arrays required")
    disk, gas_unit, bulge_unit, observed, errors, radius = arrays
    gas, wg = lognormal_quadrature(1., .04, samples)
    bulge, wb = lognormal_quadrature(.7, .1, samples)
    values, weights, cache = [], [], {}
    for g, gw in zip(gas, wg):
        for b, bw in zip(bulge, wb):
            total = disk+g*gas_unit+b*bulge_unit
            if np.any(total <= 0):
                raise ValueError("nonpositive net acceleration at nuisance node")
            key = tuple(tuple(int(i) for i in group) for group in _tie_groups(total, 0.))
            if key not in cache:
                cache[key] = boundary_log_probability_for_gbar(
                    total, observed, errors, radius, ngrid=ngrid, zpadding=8.,
                    log_tolerance=log_tolerance, max_refinements=max_refinements)
            values.append(cache[key]);weights.append(gw*bw)
    answer = float(logsumexp(np.asarray(values)+np.log(weights)))
    if not math.isfinite(answer) or answer > 1e-10:
        raise ValueError("invalid marginal ordering probability")
    return answer, len(cache), len(weights)


def run_null(counts, replicates, seed, fractional_error, samples, ngrid):
    counts = sorted(set(counts))
    if not counts or counts[0] < 2 or replicates < 1 or seed < 0:
        raise ValueError("counts>=2, replicates>=1 and seed>=0 required")
    parent = max(counts);order = nested_indices(parent)
    r, disk, gas, bulge, vtrue, errors = latent_design(parent, fractional_error)
    rows = []
    for rep in range(replicates):
        rng = np.random.default_rng(np.random.SeedSequence([seed, rep]))
        observed = vtrue+rng.normal(size=parent)*errors
        for n in counts:
            ix = np.sort(order[:n])
            row = {"replicate":rep, "bins":n, "status":"ok", "logp":None,
                   "unique_order_integrations":None, "quadrature_nodes":None, "error":None}
            try:
                lp, unique, nodes = fixed_probability(disk[ix],gas[ix],bulge[ix],observed[ix],
                                                     errors[ix],r[ix],samples=samples,ngrid=ngrid)
                row.update(logp=lp, unique_order_integrations=unique, quadrature_nodes=nodes)
            except (ValueError, ArithmeticError, RuntimeError) as error:
                row.update(status="numerical_or_domain_failure",error=str(error))
            rows.append(row)
    return rows, {"radius_kpc":r.tolist(), "disk_g_m_s2":disk.tolist(),
                  "gas_unit_g_m_s2":gas.tolist(), "bulge_unit_g_m_s2":bulge.tolist(),
                  "latent_velocity_kms":vtrue.tolist(), "error_velocity_kms":errors.tolist(),
                  "subset_indices":{str(n):np.sort(order[:n]).tolist() for n in counts}}


def self_test():
    order = nested_indices(9)
    assert len(set(order)) == 9 and order[:2] == [0,8]
    assert set(order[:3]) <= set(order[:5]) <= set(order[:9])
    r,d,g,b,v,e = latent_design(3,.05)
    assert np.all(np.diff(d+g+.7*b)<0) and np.all(np.diff(vel2acc(v*v,r))<0)
    a = run_null([2,3],1,20261014,.05,1,64)
    c = run_null([2,3],1,20261014,.05,1,64)
    assert a == c and all(x["status"] == "ok" for x in a[0])
    lp, unique, nodes = fixed_probability(np.array([2.,1.]),np.zeros(2),np.zeros(2),
                                         np.full(2,10.),np.ones(2),np.ones(2),ngrid=64)
    assert abs(lp-math.log(.5)) < 1e-4 and unique == 1 and nodes == 81
    try:
        fixed_probability(-np.ones(2),np.zeros(2),np.zeros(2),np.ones(2),np.ones(2),np.ones(2))
    except ValueError:
        pass
    else:
        raise AssertionError("nonpositive acceleration accepted")
    print("Fixed-profile mock synthetic checks passed; no GP fitting or galaxy run.")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--self-test",action="store_true")
    p.add_argument("--counts",type=int,nargs="+",default=[5,10,20,40])
    p.add_argument("--replicates",type=int)
    p.add_argument("--seed",type=int,default=20261014)
    p.add_argument("--fractional-velocity-error",type=float,default=.05)
    p.add_argument("--ml-samples",type=int,default=9)
    p.add_argument("--prob-grid-size",type=int,default=400)
    p.add_argument("--output-root",type=Path)
    args = p.parse_args()
    if args.self_test: self_test();return
    if args.replicates is None or args.output_root is None:
        p.error("explicit --replicates and fresh --output-root required")
    if args.output_root.exists(): raise FileExistsError("no overwrite")
    paths = {name:ROOT/name for name in (Path(__file__).name,"probability_recursion.py","utils.py","sparc_scale_proxies.py")}
    hashes = {k:sha256_file(v) for k,v in paths.items()}
    rows, design = run_null(args.counts,args.replicates,args.seed,args.fractional_velocity_error,args.ml_samples,args.prob_grid_size)
    if hashes != {k:sha256_file(v) for k,v in paths.items()}: raise RuntimeError("source changed")
    manifest = {"schema_version":1, "scope":"controlled_fixed_profile_N_dependence_only",
                "settings":{k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()},
                "design":design,"source_sha256":hashes,
                "noise":"independent Gaussian velocity noise; no rejected/resampled draws",
                "seed_recipe":"SeedSequence([master,replicate]); same parent draw for all N",
                "nesting":"within declared parent size; changing maxN changes design",
                "gas_prior":"median1.0,0.04dex; lognormal GH", "bulge_prior":"median0.7,0.1dex; lognormal GH",
                "probability_controls":{"starting_grid":args.prob_grid_size,"zpadding":8.,"log_tolerance":2e-5,"max_refinements":8},
                "failures":sum(x["status"]!="ok" for x in rows),
                "calibrated_GP_thresholds":False,"selection_repeated":False,
                "galaxy_matched":False,"distance_inclination_fixed":True,
                "limitation":"component fractions share shape, so nuisance rescaling normally preserves ordering; not full galaxy-shaped nuisance sensitivity",
                "interpretation":"descriptive ordering probabilities, not p-values or physical MG classification"}
    args.output_root.mkdir(parents=True)
    with (args.output_root/"mock_rows.csv").open("x",newline="") as stream:
        w=csv.DictWriter(stream,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    with (args.output_root/"manifest.json").open("x") as stream:
        json.dump(manifest,stream,indent=2,allow_nan=False);stream.write("\n")
    print("Fixed null rows",len(rows),"failures",manifest["failures"],"output",args.output_root)


if __name__ == "__main__":
    main()
