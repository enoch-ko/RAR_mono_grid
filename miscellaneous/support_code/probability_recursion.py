"""Boundary-resolving quadrature for the exact Gaussian ordering recursion.

No probability floor, indicator smoothing, failed-node omission or renormalisation
is used. This is the manuscript's continuous recursion evaluated cumulatively.
The positive and negative Z branches are summed exactly before integration.
"""
import numpy as np
from scipy.special import log_ndtr

METHOD = "boundary_log_cumulative_v1"
LOG_SQRT_2PI = 0.5 * np.log(2.0 * np.pi)


class ProbabilityConvergenceError(RuntimeError):
    """Numerical refinement failed; callers must retain the invalid status."""


def _log_exprel(d):
    """log((exp(d)-1)/d), including d=0, without overflow/cancellation."""
    d = np.asarray(d, dtype=float)
    out = np.empty_like(d)
    tiny = np.abs(d) < 1e-4
    high = d > 50.0
    low = d < -50.0
    other = ~(tiny | high | low)
    out[tiny] = d[tiny] / 2.0 + d[tiny]**2 / 24.0
    out[high] = d[high] + np.log1p(-np.exp(-d[high])) - np.log(d[high])
    out[low] = np.log(-np.expm1(d[low])) - np.log(-d[low])
    out[other] = np.log(np.expm1(d[other]) / d[other])
    return out


def _grid(mu, logscale, size, padding):
    # x=|V|/sqrt(r); y=log(x). Scale x by the largest error scale, which
    # cancels from all orderings. This is only an integration coordinate.
    low = np.min(logscale - np.log1p(np.abs(mu))) - 30.0
    high = np.max(logscale + np.log(np.abs(mu) + padding))
    pieces = [np.linspace(low, high, size)]
    offsets = np.linspace(-padding, padding, size)
    for centre, scale in zip(np.abs(mu), logscale):
        u = centre + offsets
        u = u[u > 0]
        pieces.append(scale + np.log(u))
        # Resolve rare-event overlap between separated Gaussian peaks too.
        full = np.linspace(0.0, centre + padding, size)[1:]
        pieces.append(scale + np.log(full))
    y = np.unique(np.concatenate(pieces))
    return y[(y >= low) & (y <= high)]


def _at_resolution(mu, logscale, size, padding):
    y = _grid(mu, logscale, size, padding)
    dy = np.diff(y)
    log_dy = np.log(dy)
    logH = np.zeros_like(y)  # empty continuation has probability one
    log_origin = 0.0
    for count, i in enumerate(range(len(mu)-1, -1, -1), start=1):
        u = np.exp(y - logscale[i])
        logf = np.logaddexp(-0.5*(u-mu[i])**2,
                           -0.5*(u+mu[i])**2) - LOG_SQRT_2PI - logscale[i]
        logh = y + logf + logH
        # Integrate the exponential interpolant in y over EVERY cell.
        # Narrow ordering intervals never become empty sets of sample points.
        cells = log_dy + logh[:-1] + _log_exprel(np.diff(logh))
        # At x->0, H_i(x)=prod_j f_j(0) * x^count/count! + O(x^(count+2)).
        # The lower grid bound puts every |mu_j|*x/scale_j below exp(-30).
        # This analytic origin cell is not a floor or a fictitious point mass.
        logf0 = np.log(2.0) - 0.5*mu[i]**2 - LOG_SQRT_2PI - logscale[i]
        log_origin += logf0 + y[0] - np.log(count)
        logH = np.logaddexp.accumulate(np.r_[log_origin, cells])
    upper_u = np.exp(y[-1] - logscale)
    log_tail_bound = float(np.logaddexp.reduce(np.logaddexp(
        log_ndtr(mu-upper_u), log_ndtr(-mu-upper_u))))
    return float(logH[-1]), len(y), log_tail_bound


def boundary_log_probability(Vobs, errV, r, ngrid=400, zpadding=8.0,
                             log_tolerance=2e-5, max_refinements=8,
                             return_diagnostics=False):
    """Compute log ordering probability with mandatory resolution refinement.

    ngrid is the STARTING resolution per Gaussian neighbourhood and logarithmic
    background grid, not the old uniform-grid point count. Refinement doubles it.
    log_tolerance bounds successive logP changes, NOT an exact error guarantee.
    zpadding controls the finite upper domain; check it independently for science.
    A convergence failure raises instead of substituting zero probability.
    """
    v, e, r = (np.asarray(a, dtype=float) for a in (Vobs, errV, r))
    if any(a.ndim != 1 for a in (v, e, r)) or not (len(v) == len(e) == len(r)):
        raise ValueError("velocity, uncertainty and radius arrays must match")
    if not len(v) or not all(np.all(np.isfinite(a)) for a in (v, e, r)):
        raise ValueError("finite, nonempty arrays required")
    if np.any(e <= 0) or np.any(r <= 0):
        raise ValueError("positive uncertainties and radii required")
    if int(ngrid) < 2 or not np.isfinite(zpadding) or zpadding <= 0:
        raise ValueError("invalid resolution or tail padding")
    if not np.isfinite(log_tolerance) or log_tolerance <= 0 or max_refinements < 1:
        raise ValueError("positive tolerance and at least one refinement required")
    if len(v) == 1:
        info = {"method": METHOD, "log_probability": 0.0, "levels": [],
                "log_tolerance": log_tolerance, "zpadding": zpadding}
        return info if return_diagnostics else 0.0
    mu = v/e
    logscale = np.log(e) - 0.5*np.log(r)
    logscale -= np.max(logscale)
    levels = []
    previous = None
    for level in range(int(max_refinements)+1):
        size = max(32, int(ngrid)) * 2**level
        value, points, log_tail = _at_resolution(mu, logscale, size, zpadding)
        if not np.isfinite(value):
            raise ProbabilityConvergenceError("nonfinite boundary integral")
        delta = None if previous is None else abs(value-previous)
        levels.append({"resolution": size, "points": points,
                       "log_probability": value, "log_delta": delta,
                       "log_upper_tail_union_bound": log_tail})
        if delta is not None and delta <= log_tolerance:
            if value > log_tolerance:
                raise ProbabilityConvergenceError("quadrature probability exceeds one")
            value = min(value, 0.0)  # only roundoff-sized overshoot after refinement
            info = {"method": METHOD, "log_probability": value, "levels": levels,
                    "log_tolerance": log_tolerance, "zpadding": zpadding,
                    "tail_note": "absolute union bound; not a relative-error certificate"}
            return info if return_diagnostics else value
        previous = value
    raise ProbabilityConvergenceError(
        f"logP refinement did not converge: delta={delta:.6g}, tolerance={log_tolerance}")


def boundary_log_probability_for_gbar(g_bar, Vobs, errV, r, ngrid=400,
                                      zpadding=8.0, tie_atol=0.0,
                                      max_tie_permutations=720, **kwargs):
    """Ordered-gbar entry point; retain arithmetic probability averaging over ties."""
    import math
    import warnings
    from itertools import permutations, product
    from scipy.special import logsumexp
    from utils import _tie_groups, _validate_probability_inputs
    g_bar = np.asarray(g_bar, dtype=float)
    v, e, r = _validate_probability_inputs(Vobs, errV, r, ngrid, zpadding)
    if g_bar.ndim != 1 or len(g_bar) != len(v) or not np.all(np.isfinite(g_bar)):
        raise ValueError("g_bar must match the finite kinematic arrays")
    if not np.isfinite(tie_atol) or tie_atol < 0 or int(max_tie_permutations) < 1:
        raise ValueError("invalid tie tolerance or permutation cap")
    if kwargs.get("return_diagnostics", False):
        raise ValueError("request per-order diagnostics from boundary_log_probability")
    groups = _tie_groups(g_bar, tie_atol)
    # Count actual permutations before allocation; sizes are not factorials.
    total = math.prod(math.factorial(len(g)) for g in groups)
    if total > int(max_tie_permutations):
        warnings.warn(f"{total} tied-order permutations exceed max_tie_permutations="
                      f"{max_tie_permutations}; using stable descending order", RuntimeWarning)
        orders = [np.concatenate(groups)]
    else:
        choices = [list(permutations(g)) for g in groups]
        orders = (np.concatenate(parts) for parts in product(*choices))
    values = [boundary_log_probability(v[o], e[o], r[o], ngrid, zpadding, **kwargs)
              for o in orders]
    return float(logsumexp(values)-np.log(len(values)))
