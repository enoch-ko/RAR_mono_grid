"""Exploratory fixed-baryonic Bayes-factor prototype for symmetric bounded velocities.

For supplied radii and g_bar, calculate the monotonic-order probability under
(i) independent Uniform(-Vmax,Vmax) velocities and (ii) the corresponding
Gaussian likelihood posterior truncated to the SAME bounds. Probabilities use
a finite-support ordering recursion; the prior also has a seeded Monte Carlo
cross-check. This is not the production RAR solver and has no real-data loader.
"""
from dataclasses import dataclass
from itertools import permutations, product

import math
import numpy as np
from scipy.special import logsumexp
from scipy.stats import norm, truncnorm






MAX_TIE_ORDERINGS = 5040


@dataclass(frozen=True)
class ProbabilityResult:
    log_probability: float
    grid_points: int
    orderings_averaged: int

    @property
    def probability(self):
        return float(np.exp(self.log_probability)) if self.log_probability > -745.0 else 0.0


@dataclass(frozen=True)
class PriorMonteCarloCheck:
    draws: int
    successes: int
    probability: float
    standard_error: float
    zero_success_upper_95: float | None
    z_score_vs_recursion: float | None


@dataclass(frozen=True)
class BayesFactorResult:
    prior: ProbabilityResult
    posterior: ProbabilityResult
    log_bayes_factor: float
    bayes_factor: float | None
    prior_monte_carlo: PriorMonteCarloCheck


def _validate(radii, g_bar, v_obs, sigma_v, vmax, ngrid, mc_draws, seed, tie_atol):
    arrays = [np.asarray(x, dtype=float) for x in (radii, g_bar, v_obs, sigma_v)]
    if any(a.ndim != 1 for a in arrays) or len({a.size for a in arrays}) != 1:
        raise ValueError("radii, g_bar, v_obs, sigma_v must be equal-length 1D arrays")
    r, gb, vo, ev = arrays
    if r.size < 2 or not all(np.all(np.isfinite(a)) for a in arrays):
        raise ValueError("need at least two rows and all inputs must be finite")
    if np.any(r <= 0.0) or np.any(ev <= 0.0):
        raise ValueError("radii and sigma_v must be strictly positive")
    if not np.isfinite(vmax) or vmax <= 0.0:
        raise ValueError("an explicit finite positive Vmax is required")
    if isinstance(ngrid, bool) or int(ngrid) != ngrid or ngrid < 32:
        raise ValueError("ngrid must be an integer >= 32")
    if isinstance(mc_draws, bool) or int(mc_draws) != mc_draws or mc_draws < 1:
        raise ValueError("mc_draws must be a positive integer")
    if seed is None:
        raise ValueError("an explicit RNG seed is required")
    if not np.isfinite(tie_atol) or tie_atol < 0.0:
        raise ValueError("tie_atol must be finite and nonnegative")
    return r, gb, vo, ev, float(vmax), int(ngrid), int(mc_draws), int(seed), float(tie_atol)


def _orderings(g_bar, tie_atol):
    order = np.argsort(-g_bar, kind="stable")
    values = g_bar[order]
    groups = []
    start = 0
    while start < len(order):
        stop = start + 1
        while stop < len(order) and abs(values[stop] - values[start]) <= tie_atol:
            stop += 1
        groups.append(order[start:stop])
        start = stop
    count = 1
    for group in groups:
        if len(group) > 1:
            count *= int(math.factorial(len(group)))
    if count > MAX_TIE_ORDERINGS:
        raise ValueError(f"tied g_bar values imply {count} orderings; cap is {MAX_TIE_ORDERINGS}")
    choices = [list(permutations(group)) if len(group) > 1 else [tuple(group)] for group in groups]
    return [np.asarray(sum((list(x) for x in choice), []), dtype=int) for choice in product(*choices)]


def _cumtrapz(y, x):
    out = np.zeros_like(y, dtype=float)
    out[1:] = np.cumsum(0.5 * (y[:-1] + y[1:]) * np.diff(x))
    return out


def _one_order_log_probability(r, vo, ev, vmax, ordering, ngrid, posterior):
    r = r[ordering]
    vo = vo[ordering]
    ev = ev[ordering]
    bounds = vmax / np.sqrt(r)
    xmax = float(np.max(bounds))
    x = np.unique(np.concatenate((np.linspace(0.0, xmax, ngrid), bounds)))
    next_h = np.ones_like(x)
    log_scale = 0.0

    for i in range(len(r) - 1, -1, -1):
        supported = x <= bounds[i]
        if posterior:
            a = (-vmax - vo[i]) / ev[i]
            b = (vmax - vo[i]) / ev[i]
            xp = np.sqrt(r[i]) * x[supported]
            lp = truncnorm.logpdf(xp, a, b, loc=vo[i], scale=ev[i])
            lm = truncnorm.logpdf(-xp, a, b, loc=vo[i], scale=ev[i])
            logf = 0.5 * np.log(r[i]) + logsumexp(np.vstack((lp, lm)), axis=0)
            density = np.zeros_like(x)
            density[supported] = np.exp(logf)
        else:
            density = np.where(supported, np.sqrt(r[i]) / vmax, 0.0)
        h = _cumtrapz(density * next_h, x)
        scale = float(h[-1])
        if not np.isfinite(scale) or scale <= 0.0:
            raise ArithmeticError("finite-support recursion failed; increase ngrid")
        next_h = h / scale
        log_scale += np.log(scale)
    return float(log_scale)


def _finite_support_probability(r, gb, vo, ev, vmax, ngrid, tie_atol, posterior):
    orders = _orderings(gb, tie_atol)
    logs = [_one_order_log_probability(r, vo, ev, vmax, order, ngrid, posterior) for order in orders]
    logp = float(logsumexp(logs) - np.log(len(logs)))
    return ProbabilityResult(logp, ngrid, len(logs))


def _unbounded_gaussian_cutoff(r, vo, ev, tail_tolerance):
    if not np.isfinite(tail_tolerance) or not 0.0 < tail_tolerance < 1.0:
        raise ValueError("tail_tolerance must lie strictly between 0 and 1")
    z = float(norm.isf(tail_tolerance / (2.0 * len(r))))
    velocity_cutoffs = np.abs(vo) + ev * z
    xmax = float(np.max(velocity_cutoffs / np.sqrt(r)))
    if not np.isfinite(xmax) or xmax <= 0.0:
        raise ArithmeticError("could not construct a finite Gaussian-tail integration bound")
    return xmax


def _one_order_log_probability_unbounded_gaussian(r, vo, ev, ordering, ngrid, xmax):
    r, vo, ev = r[ordering], vo[ordering], ev[ordering]
    x = np.linspace(0.0, xmax, ngrid)
    next_h = np.ones_like(x)
    log_scale = 0.0
    for i in range(len(r) - 1, -1, -1):
        xp = np.sqrt(r[i]) * x
        lp = norm.logpdf(xp, loc=vo[i], scale=ev[i])
        lm = norm.logpdf(-xp, loc=vo[i], scale=ev[i])
        logf = 0.5 * np.log(r[i]) + logsumexp(np.vstack((lp, lm)), axis=0)
        density = np.exp(logf)
        h = _cumtrapz(density * next_h, x)
        scale = float(h[-1])
        if not np.isfinite(scale) or scale <= 0.0:
            raise ArithmeticError("Gaussian ordering recursion failed; increase ngrid")
        next_h = h / scale
        log_scale += np.log(scale)
    return float(log_scale)


def estimate_broad_limit_bayes_factor(
    radii, g_bar, v_obs, sigma_v, *, seed, ngrid=16001, mc_draws=100_000,
    tie_atol=0.0, tail_tolerance=1e-10,
):
    """Estimate the broad symmetric-prior conditional BF at fixed baryonic inputs.

    The denominator uses iid Uniform(-1,1) velocities; its event probability
    equals that for any common finite Vmax. The numerator uses the unbounded
    signed Gaussian posterior for Vobs. Gaussian integration is truncated only
    numerically at a union-bound tail mass <= tail_tolerance, which is returned
    as an additive upper bound on the omitted ordering probability.
    """
    r, gb, vo, ev, _, ngrid, mc_draws, seed, tie_atol = _validate(
        radii, g_bar, v_obs, sigma_v, 1.0, ngrid, mc_draws, seed, tie_atol
    )
    prior = _finite_support_probability(r, gb, vo, ev, 1.0, ngrid, tie_atol, posterior=False)
    prior_mc = _prior_mc(r, gb, 1.0, mc_draws, tie_atol, seed, prior.log_probability)
    xmax = _unbounded_gaussian_cutoff(r, vo, ev, tail_tolerance)
    orders = _orderings(gb, tie_atol)
    logqs = [
        _one_order_log_probability_unbounded_gaussian(r, vo, ev, order, ngrid, xmax)
        for order in orders
    ]
    logq_lower = float(logsumexp(logqs) - np.log(len(logqs)))
    q_lower = float(np.exp(logq_lower)) if logq_lower > -745.0 else 0.0
    q_upper = min(1.0, q_lower + tail_tolerance)
    logq_upper = min(0.0, float(np.logaddexp(logq_lower, np.log(tail_tolerance))))
    return {
        "prior": prior,
        "posterior_log_probability_lower": logq_lower,
        "posterior_probability_lower": q_lower,
        "posterior_probability_upper": q_upper,
        "posterior_probability_tail_upper_bound": float(tail_tolerance),
        "log_bayes_factor_lower": logq_lower - prior.log_probability,
        "log_bayes_factor_upper": logq_upper - prior.log_probability,
        "prior_monte_carlo": prior_mc,
        "tail_tolerance": float(tail_tolerance),
        "integration_xmax": xmax,
    }


def _prior_mc(r, gb, vmax, draws, tie_atol, seed, exact_logp):
    rng = np.random.default_rng(seed)
    orders = _orderings(gb, tie_atol)
    order = orders[0]
    sorted_g = gb[order]
    groups = []
    start = 0
    while start < len(order):
        stop = start + 1
        while stop < len(order) and abs(sorted_g[stop] - sorted_g[start]) <= tie_atol:
            stop += 1
        groups.append((start, stop))
        start = stop
    hits = 0
    batch = min(20_000, draws)
    completed = 0
    while completed < draws:
        n = min(batch, draws - completed)
        v = rng.uniform(-vmax, vmax, size=(n, len(r)))
        gobs = (v * v / r[None, :])[:, order]
        if len(orders) == 1:
            hits += int(np.count_nonzero(np.all(np.diff(gobs, axis=1) < 0.0, axis=1)))
        else:
            indicators = np.empty(n, dtype=float)
            for row in range(n):
                row_g = gobs[row].copy()
                for start, stop in groups:
                    if stop - start > 1:
                        row_g[start:stop] = row_g[start:stop][rng.permutation(stop - start)]
                indicators[row] = np.all(np.diff(row_g) < 0.0)
            hits += int(np.count_nonzero(indicators))
        completed += n
    p = hits / draws
    se = float(np.sqrt(p * (1.0 - p) / draws))
    upper = float(1.0 - 0.05 ** (1.0 / draws)) if hits == 0 else None
    exact = float(np.exp(exact_logp)) if exact_logp > -745.0 else 0.0
    z = (p - exact) / se if se > 0.0 else None
    return PriorMonteCarloCheck(draws, hits, p, se, upper, z)


def estimate_bayes_factor(
    radii, g_bar, v_obs, sigma_v, *, vmax, seed, ngrid=4001, mc_draws=100_000, tie_atol=0.0
):
    """Return fixed-input conditional BF and a prior Monte Carlo cross-check.

    Prior: independent Uniform(-Vmax,Vmax) velocities. Posterior: independent
    N(v_obs,sigma_v) likelihood factors truncated and renormalised on precisely
    [-Vmax,Vmax]. The folded posterior acceleration density has both velocity
    branches: sqrt(R) * (p(+sqrt(R)*x|D) + p(-sqrt(R)*x|D)).

    Vmax and seed are mandatory; no catalogue-derived default is provided.
    The recursion is deterministic quadrature on finite support. Monte Carlo is
    used only as an independent prior-probability cross-check. This conditional
    BF is not a nuisance-integrated model comparison or a scientific result.
    """
    r, gb, vo, ev, vmax, ngrid, mc_draws, seed, tie_atol = _validate(
        radii, g_bar, v_obs, sigma_v, vmax, ngrid, mc_draws, seed, tie_atol
    )
    prior = _finite_support_probability(r, gb, vo, ev, vmax, ngrid, tie_atol, posterior=False)
    posterior = _finite_support_probability(r, gb, vo, ev, vmax, ngrid, tie_atol, posterior=True)
    logbf = posterior.log_probability - prior.log_probability
    bf = float(np.exp(logbf)) if -745.0 < logbf < 709.0 else None
    mc = _prior_mc(r, gb, vmax, mc_draws, tie_atol, seed, prior.log_probability)
    return BayesFactorResult(prior, posterior, float(logbf), bf, mc)


def _main():
    import argparse
    import json

    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--plan", action="store_true", help="print the input-only proposed real-data pilot receipt")
    modes.add_argument("--mocktest", action="store_true", help="run the broad-limit API on a synthetic fixture only")
    modes.add_argument("--synthetic-demo", action="store_true", help="run the legacy finite-bound synthetic example")
    parser.add_argument("--vmax", type=float, help="explicit common velocity bound for --synthetic-demo only")
    parser.add_argument("--seed", type=int, help="explicit random seed for a synthetic run")
    parser.add_argument("--draws", type=int, help="explicit prior Monte Carlo draws for synthetic runs")
    parser.add_argument("--ngrid", type=int, default=16001, help="quadrature grid size")
    parser.add_argument("--tail-tolerance", type=float, default=1e-10, help="upper bound on omitted Gaussian posterior tail mass")
    args = parser.parse_args()

    if args.plan:
        print(json.dumps({
            "scope": "input-only proposed pilot; not executed",
            "preferred_target": "NGC4068",
            "profile_conditioning": {"disk_ML": 0.5, "gas_ML": 1.0, "bulge": "absent / bulgeless"},
            "target_retained_radii": 6,
            "model": "fixed-theta conditional broad-prior-limit BF = Q(theta)/c(theta)",
            "denominator": "c(theta) under iid Uniform(-1,1) velocities; common finite velocity-bound scale cancels from the ordering event",
            "numerator": "Q(theta) from the matching unbounded signed-Gaussian Vobs likelihood, folded into V^2/R",
            "must_pin_before_real_data": [
                "exact source files and hashes for R, Vobs, sigmaV and fixed baryonic inputs",
                "retained-radius mask, ordering, and tie policy",
                "seed, prior-MC draw budget, quadrature grid and refinement schedule",
                "rare-event/zero-hit reporting, tail tolerance, and runtime cap",
            ],
            "status": "settings are proposals only; no real-data calculation, queue job, or manuscript change",
            "warning": "Do not substitute nuisance-marginalised saved GP probabilities for this fixed-theta conditional numerator without proving equivalence.",
        }, indent=2))
        return

    if args.seed is None or args.draws is None:
        parser.error("--seed and --draws are required for every synthetic execution mode")
    if args.draws < 1:
        parser.error("--draws must be positive")

    if args.mocktest:
        result = estimate_broad_limit_bayes_factor(
            [1.0, 2.0, 3.0], [3.0, 2.0, 1.0], [70.0, 55.0, 35.0], [5.0, 5.0, 5.0],
            seed=args.seed, ngrid=args.ngrid, mc_draws=args.draws,
            tail_tolerance=args.tail_tolerance,
        )
        output = {
            "scope": "synthetic fixture only; no real-galaxy data",
            "seed": args.seed, "draws": args.draws, "ngrid": args.ngrid,
            "prior_probability": result["prior"].probability,
            "prior_log_probability": result["prior"].log_probability,
            "prior_monte_carlo": result["prior_monte_carlo"].__dict__,
            "posterior_probability_lower": result["posterior_probability_lower"],
            "posterior_probability_upper": result["posterior_probability_upper"],
            "posterior_probability_tail_upper_bound": result["posterior_probability_tail_upper_bound"],
            "log_bayes_factor_lower": result["log_bayes_factor_lower"],
            "log_bayes_factor_upper": result["log_bayes_factor_upper"],
            "integration_xmax": result["integration_xmax"],
            "warning": "synthetic smoke test only; not a real-data result or convergence claim",
        }
        print(json.dumps(output, indent=2))
        return

    if args.vmax is None or not np.isfinite(args.vmax) or args.vmax <= 0.0:
        parser.error("--synthetic-demo requires an explicit positive --vmax")
    result = estimate_bayes_factor(
        [1.0, 2.0, 3.0], [3.0, 2.0, 1.0], [70.0, 55.0, 35.0], [5.0, 5.0, 5.0],
        vmax=args.vmax, seed=args.seed, ngrid=args.ngrid, mc_draws=args.draws,
    )
    print(json.dumps({
        "scope": "synthetic only; fixed inputs; no real-galaxy data",
        "vmax": args.vmax, "seed": args.seed, "draws": args.draws, "ngrid": args.ngrid,
        "prior_log_probability": result.prior.log_probability,
        "prior_probability": result.prior.probability,
        "posterior_log_probability": result.posterior.log_probability,
        "posterior_probability": result.posterior.probability,
        "log_bayes_factor": result.log_bayes_factor,
        "bayes_factor": result.bayes_factor,
        "prior_mc": result.prior_monte_carlo.__dict__,
        "warning": "conditional fixed-input diagnostic, not nuisance-marginalised evidence or a science result",
    }, indent=2))


if __name__ == "__main__":
    _main()
