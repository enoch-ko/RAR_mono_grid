#!/usr/bin/env python3
"""Independent folded-normal pair/triple quadrature and ordering bounds.
Uses saved fixed-nuisance UGC11914 fit; no recursion for reference integrals.
"""
import json
import hashlib
import signal
import time
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.special import ndtr
import probability_recursion as boundary
import utils

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "outputs/inverse_rar_diagnostic/comparison_cholesky/checks_1_4_20261006/UGC11914_fit.json"
OUTPUT = SOURCE.with_name("UGC11914_independent.json")

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def normal_mass(a, b, mu, sigma):
    a, b = (a-mu)/sigma, (b-mu)/sigma
    return float(ndtr(-a)-ndtr(-b) if a >= 0 else ndtr(b)-ndtr(a))

def main():
    if OUTPUT.exists():
        raise FileExistsError(OUTPUT)
    start = time.monotonic()
    def timeout(signum, frame):
        raise TimeoutError("40-second login-node limit")
    signal.signal(signal.SIGALRM, timeout)
    signal.alarm(40)
    d = json.loads(SOURCE.read_text())
    g = np.array(d["fit"]["g_bar"])
    order = np.argsort(-g, kind="stable")
    r, v, e = [np.array(d[k])[order] for k in ("r", "vobs", "errv")]
    mu, sigma = v/np.sqrt(r), e/np.sqrt(r)
    def cdf(i, x):
        return normal_mass(-x, x, mu[i], sigma[i])
    def sf(i, x):
        return float(ndtr((mu[i]-x)/sigma[i])+ndtr((-mu[i]-x)/sigma[i]))
    def reference(ids):
        middle = ids[0] if len(ids) == 2 else ids[1]
        def integrand(z):
            x = abs(mu[middle]+sigma[middle]*z)
            factor = cdf(ids[-1], x)
            if len(ids) == 3:
                factor *= sf(ids[0], x)
            return np.exp(-z*z/2)/np.sqrt(2*np.pi)*factor
        # Integrate in standardized Gaussian coordinates, independent of either mesh.
        points = [-12., 12., -mu[middle]/sigma[middle]]
        for i in ids:
            for sign in (-1, 1):
                points.append((sign*mu[i]-mu[middle])/sigma[middle])
        points = sorted(set(max(-12., min(12., x)) for x in points))
        values = [quad(integrand, a, b, epsabs=1e-12, epsrel=1e-10, limit=150)
                  for a, b in zip(points[:-1], points[1:]) if b > a]
        return float(sum(x[0] for x in values)), float(sum(x[1] for x in values))
    pairs = []
    triples = []
    for size, records in [(2, pairs), (3, triples)]:
        for j in range(len(mu)-size+1):
            ids = list(range(j,j+size))
            p, error = reference(ids)
            records.append({"start":j,"radii_kpc":r[ids].tolist(),"probability":p,
                            "quad_error_estimate":error})
    # Full ordering implies each independent disjoint block is ordered.
    upper = min(sum(np.log(x["probability"]) for x in records[offset::size])
                for size, records in [(2,pairs),(3,triples)] for offset in range(size))
    # Conversely, all Xi in these nonoverlapping intervals guarantees full ordering.
    cuts = (mu[:-1]+mu[1:])/2
    masses = []
    for i in range(len(mu)):
        lo = cuts[i] if i < len(cuts) else 0.
        hi = cuts[i-1] if i > 0 else np.inf
        masses.append(normal_mass(lo,hi,mu[i],sigma[i]) + normal_mass(-hi,-lo,mu[i],sigma[i]))
    lower = float(np.sum(np.log(masses)))
    selected = []
    for size, records in [(2,pairs),(3,triples)]:
        picks = sorted(set([int(np.argmin([x["probability"] for x in records])),0,len(records)//2,len(records)-1]))
        for j in picks:
            ids = np.arange(j,j+size)
            item = dict(records[j],size=size)
            for name, fun in [("legacy",utils.log_monotonic_probability_for_gbar),
                              ("boundary",boundary.boundary_log_probability_for_gbar)]:
                lp = float(fun(g[order][ids],v[ids],e[ids],r[ids],ngrid=400,zpadding=8))
                item[name+"_logp"] = lp
                item[name+"_probability"] = float(np.exp(lp))
            selected.append(item)
    result = {"galaxy":"UGC11914","diagnostic_only":True,"points":len(mu),
              "exact_ties":len(g)-len(np.unique(g)),"central_reversals":int(np.sum(np.diff(mu)>0)),
              "pairs":pairs,"triples":triples,"selected_comparisons":selected,
              "full_logp_lower_bound":lower,"full_logp_upper_bound":float(upper),
              "pair_probability_range":[min(x["probability"] for x in pairs),max(x["probability"] for x in pairs)],
              "triple_probability_range":[min(x["probability"] for x in triples),max(x["probability"] for x in triples)],
              "maximum_quad_error_estimate":max(x["quad_error_estimate"] for x in pairs+triples),
              "normal_tail_omission_bound":float(2*ndtr(-12)),"source_sha256":sha(SOURCE),
              "script_sha256":sha(__file__),"utils_sha256":sha(utils.__file__),
              "boundary_sha256":sha(boundary.__file__),"seconds":time.monotonic()-start}
    OUTPUT.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    signal.alarm(0)
    print(json.dumps({k:v for k,v in result.items() if k not in ("pairs","triples")}))

if __name__ == "__main__":
    main()
