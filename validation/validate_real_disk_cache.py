#!/usr/bin/env python3
"""Light real-SPARC direct-vs-cached Casertano validation; no production outputs."""
import argparse
import hashlib
import json
import time
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import mono_SPARC_RAR as mono

GALAXIES = ("IC4202", "NGC0891", "UGC07399", "D631-7")
MEANS = (0.5, 0.8)
SLOPES = (-0.01, 0.0, 0.01)

def solve(rad, density, zd, r):
    value, caught, error = mono.solve_disk_basis(rad, density, zd, r)
    if error or not np.all(np.isfinite(value)):
        raise RuntimeError(f"vcdisk failure: {error}; warnings={caught}")
    return value, caught

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--galaxy", action="append")
    p.add_argument("--range", nargs=2, type=int, metavar=("START", "END"))
    p.add_argument("--prob-grid-size", type=int, default=101)
    p.add_argument("--output", type=Path, default=Path("outputs/real_disk_cache_validation_20261005.json"))
    args = p.parse_args()
    table = mono.load_sparc_table()
    if args.range is not None and (args.galaxy or not (1 <= args.range[0] <= args.range[1] <= len(table))):
        p.error("choose --galaxy or a valid inclusive --range within the SPARC catalogue")
    prob_args = SimpleNamespace(prob_grid_size=args.prob_grid_size, zpadding=8.0,
                                tie_atol=0.0, max_tie_permutations=720)
    records = []
    qc = []
    t0 = time.monotonic()
    for galaxy in (args.galaxy or (list(table.Galaxy.iloc[args.range[0] - 1:args.range[1]]) if args.range else GALAXIES)):
        rad, sb = mono.get_SBdisk(galaxy, table=table)
        data = np.loadtxt(mono.DEFAULT_DATA_ROOT / (galaxy + "_rotmod.dat"))
        keep = np.isfinite(data[:, 0]) & (data[:, 0] >= rad[0]) & (data[:, 0] <= rad[-1])
        r, vobs, errv, vgas, vbul = (data[keep, j] for j in (0, 1, 2, 3, 5))
        if len(r) < 2 or not np.all(r > 0):
            raise RuntimeError(f"{galaxy}: invalid retained radius sample")
        rdisk = float(table.loc[table.Galaxy == galaxy, "Rdisk"].iloc[0])
        zd = 0.196 * rdisk ** 0.633
        A, B = mono.compute_A_B(rad, sb)
        center = B / A
        base, base_warnings = solve(rad, sb, zd, r)
        slope_basis, slope_warnings = solve(rad, sb * (rad - center), zd, r)
        log_basis = {}
        log_shape = {}
        for slope in SLOPES:
            _, beta = mono.compute_A_beta(rad, sb, slope)
            shape = A / beta * np.exp(slope * rad)
            log_shape[slope] = shape
            log_basis[slope], _ = solve(rad, sb * shape, zd, r)
        q = dict(galaxy=galaxy, photometric_points=int(len(rad)), photometric_zero_start=bool(rad[0] == 0),
                 photometric_range_kpc=[float(rad[0]), float(rad[-1])], kinematic_input=int(len(data)),
                 kinematic_retained=int(len(r)), kinematic_excluded=int(len(data) - len(r)),
                 retained_range_kpc=[float(r[0]), float(r[-1])], center_kpc=float(center),
                 zd_kpc=float(zd), basis_warning_count=len(base_warnings) + len(slope_warnings))
        qc.append(q)
        for profile in ("constant", "linear", "log_linear"):
            for mu in MEANS:
                for slope in ((0.0,) if profile == "constant" else SLOPES):
                    if profile == "constant":
                        ups = np.full_like(rad, mu)
                        cached = mu * base
                    elif profile == "linear":
                        ups = mu + slope * (rad - center)
                        cached = mu * base + slope * slope_basis
                    else:
                        ups = mu * log_shape[slope]
                        cached = mu * log_basis[slope]
                    if not np.all(np.isfinite(ups)) or np.min(ups) <= 0:
                        q.setdefault("invalid_cells", []).append(dict(profile=profile, mean_ml=mu, slope=slope, reason="nonphysical_on_full_photometric_grid"))
                        continue
                    direct, direct_warnings = solve(rad, sb * ups, zd, r)
                    g_direct = mono.vel2acc(direct, r)
                    g_cached = mono.vel2acc(cached, r)
                    delta = np.abs(g_direct - g_cached)
                    scale = float(np.max(np.abs(g_direct)))
                    bound = 1e-10 * scale + 1e-8 * np.abs(g_direct)
                    gas = vgas * np.abs(vgas) + 0.7 * vbul ** 2
                    gbar_direct = mono.vel2acc(direct + gas, r)
                    gbar_cached = mono.vel2acc(cached + gas, r)
                    rank_agrees = bool(np.array_equal(np.argsort(gbar_direct), np.argsort(gbar_cached)))
                    order = np.argsort(gbar_direct)
                    gaps = np.diff(gbar_direct[order])
                    min_gap = float(np.min(gaps)) if gaps.size else None
                    p_direct = mono.evaluate_probability(gbar_direct, vobs, errv, r, prob_args)
                    p_cached = mono.evaluate_probability(gbar_cached, vobs, errv, r, prob_args)
                    rec = dict(galaxy=galaxy, profile=profile, mean_ml=mu, slope=slope,
                               slope_units="ML_per_kpc" if profile == "linear" else "inverse_kpc" if profile == "log_linear" else None,
                               min_profile_ml=float(np.min(ups)), max_profile_ml=float(np.max(ups)),
                               force_pass=bool(np.all(delta <= bound)), max_abs_accel_delta=float(np.max(delta)),
                               max_delta_over_scale=float(np.max(delta) / scale),
                               max_delta_over_bound=float(np.max(delta / bound)),
                               sign_agrees=bool(np.array_equal(np.sign(direct), np.sign(cached))),
                               rank_agrees=rank_agrees, min_gbar_gap=min_gap,
                               max_gbar_delta=float(np.max(np.abs(gbar_cached - gbar_direct))),
                               logp_direct=p_direct[0] if p_direct[1] else None,
                               logp_cached=p_cached[0] if p_cached[1] else None,
                               prob_valid_agrees=bool(p_direct[1] == p_cached[1]),
                               prob_reason_direct=p_direct[2], prob_reason_cached=p_cached[2],
                               logp_abs_delta=abs(p_direct[0] - p_cached[0]) if p_direct[1] and p_cached[1] else None,
                               direct_warning_count=len(direct_warnings))
                    records.append(rec)
        local = [x for x in records if x["galaxy"] == galaxy]
        print(galaxy, "cases", len(local), "force_pass", sum(x["force_pass"] for x in local),
              "max_scaled_delta", max(x["max_delta_over_scale"] for x in local),
              "rank_agree", sum(x["rank_agrees"] for x in local),
              "max_logp_delta", max((x["logp_abs_delta"] or 0) for x in local), flush=True)
    payload = dict(status="diagnostic_only", galaxies=qc, records=records,
                   criterion="abs(g_direct-g_cached) <= 1e-10*g_scale + 1e-8*abs(g_direct)",
                   probability=dict(grid_size=args.prob_grid_size, zpadding=8.0, tie_atol=0.0,
                                    gas_ml=1.0, bulge_ml=0.7),
                   code_sha256=hashlib.sha256(Path(mono.__file__).read_bytes()).hexdigest(),
                   diagnostic_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                   runtime_seconds=time.monotonic() - t0)
    if args.output.exists():
        raise FileExistsError(args.output)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + chr(10))
    print("Saved", args.output, "seconds", round(payload["runtime_seconds"], 2), flush=True)

if __name__ == "__main__":
    main()
