#!/usr/bin/env python
"""Shared numerical and SPARC data utilities for the RAR monotonicity analysis."""

import warnings
from itertools import permutations, product

import numpy as np
import pandas as pd
from scipy.integrate import trapezoid
from scipy.special import erf, logsumexp

DEFAULT_PROJECT_ROOT = "/mnt/users/koe/RAR_mono_grid"
DEFAULT_TABLE_PATH = f"{DEFAULT_PROJECT_ROOT}/SPARC_table.txt"
DEFAULT_DENS_ROOT = "/mnt/users/koe/SPARC_RAR/BulgeDiskDec_LTG"


def load_sparc_table(file=None):
    """Load the SPARC catalogue table used by the active workspace."""
    # These names follow the fixed-width column layout distributed with SPARC.
    columns = ["Galaxy", "T", "D", "e_D", "f_D", "Inc", "e_Inc", "L",
               "e_L", "Reff", "SBeff", "Rdisk", "SBdisk", "MHI", "RHI",
               "Vflat", "e_Vflat", "Q", "Ref."]
    file = DEFAULT_TABLE_PATH if file is None else file
    return pd.read_fwf(file, skiprows=98, names=columns)


def vel2acc(vel_sq, rad):
    """Convert velocity-squared in (km/s)^2 and radius in kpc to m/s^2."""
    # 1 (km/s)^2 = 10^6 (m/s)^2 and 1 kpc = 3.086e19 m.
    return (vel_sq * 1e6) / (rad * 3.086e19)


def get_SBdisk(galaxy: str, dens_root=DEFAULT_DENS_ROOT, table=None):
    """Read a galaxy's disk profile and normalise it to its SPARC luminosity."""
    file = f"{dens_root}/{galaxy}.dens"
    columns = ["Rad", "SBdisk", "SBbul"]
    data = pd.read_csv(file, names=columns, skiprows=1, sep="\t")
    rad = data["Rad"].to_numpy()
    SBdisk_raw = data["SBdisk"].to_numpy() * 1e6
    Ltot, _ = compute_A_B(rad, SBdisk_raw)
    table = load_sparc_table() if table is None else table
    matches = np.flatnonzero(table["Galaxy"].to_numpy() == galaxy)
    if len(matches) != 1:
        raise ValueError(f"Expected one SPARC table row for {galaxy!r}, found {len(matches)}")
    L = table.iloc[matches[0]]["L"] * 1e9
    if not np.isfinite(Ltot) or Ltot <= 0:
        raise ValueError(f"Invalid disk luminosity integral for {galaxy!r}: {Ltot}")
    return rad, SBdisk_raw * L / Ltot


def compute_A_B(R, I_R):
    """Return the two radial integrals needed to normalise a disk profile."""
    # A is the luminosity integral; B/A is the luminosity-weighted mean radius.
    fac = 2 * np.pi
    return fac * trapezoid(I_R * R, R), fac * trapezoid(I_R * R**2, R)


def compute_A_beta(R, I_R, log_slope):
    """Return the normalisation integrals for an exponential disk M/L profile."""
    fac = 2 * np.pi
    A = fac * trapezoid(I_R * R, R)
    beta = fac * trapezoid(I_R * np.exp(log_slope * R) * R, R)
    return A, beta


def lognormal_quadrature(median, sigma_dex, samples=9):
    """Return Gauss-Hermite nodes and weights for lognormal M/L uncertainty.

    The quoted SPARC uncertainty is interpreted as a Gaussian standard deviation
    in log10(M/L).  The returned positive nodes therefore preserve the SPARC
    convention while allowing deterministic marginalisation over that uncertainty.
    """
    median = float(median)
    sigma_dex = float(sigma_dex)
    samples = int(samples)
    if not np.isfinite(median) or median <= 0:
        raise ValueError("the lognormal median must be finite and positive")
    if not np.isfinite(sigma_dex) or sigma_dex < 0:
        raise ValueError("the lognormal scatter must be finite and non-negative")
    if samples < 1:
        raise ValueError("the number of quadrature samples must be positive")
    if sigma_dex == 0 or samples == 1:
        return np.array([median]), np.array([1.0])
    # Hermite quadrature integrates exp(-x^2) f(x); sqrt(2) rescales a standard
    # normal in log10(M/L), and 1/sqrt(pi) normalises the weights.
    nodes, weights = np.polynomial.hermite.hermgauss(samples)
    log10_nodes = np.log10(median) + np.sqrt(2.0) * sigma_dex * nodes
    return 10.0 ** log10_nodes, weights / np.sqrt(np.pi)


def log_weighted_average(log_values, weights):
    """Compute log(sum_i weights_i exp(log_values_i)) stably."""
    log_values = np.asarray(log_values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    if log_values.shape != weights.shape or log_values.ndim != 1:
        raise ValueError("log_values and weights must be one-dimensional and equally sized")
    if len(log_values) == 0 or not np.all(np.isfinite(weights)) or np.any(weights < 0):
        raise ValueError("weights must be finite, non-negative, and non-empty")
    if not np.isclose(np.sum(weights), 1.0):
        weights = weights / np.sum(weights)
    if np.any(weights == 0):
        keep = weights > 0
        log_values, weights = log_values[keep], weights[keep]
    return float(logsumexp(log_values + np.log(weights)))


LOG_SQRT_2PI = 0.5 * np.log(2.0 * np.pi)


def log_normal_pdf(x, mu):
    """Evaluate the standard-normal log density at x-mu."""
    return -0.5 * (x - mu)**2 - LOG_SQRT_2PI


def normal_cdf(x):
    """Evaluate the standard-normal cumulative distribution function."""
    return 0.5 * (1.0 + erf(x / np.sqrt(2.0)))


def _validate_probability_inputs(Vobs, errV, r, ngrid, zpadding):
    """Validate and convert inputs shared by the probability routines."""
    Vobs = np.asarray(Vobs, dtype=float)
    errV = np.asarray(errV, dtype=float)
    r = np.asarray(r, dtype=float)
    if Vobs.ndim != 1 or errV.ndim != 1 or r.ndim != 1:
        raise ValueError("Vobs, errV, and r must be one-dimensional arrays")
    if not (len(Vobs) == len(errV) == len(r)):
        raise ValueError("Vobs, errV, and r must have the same length")
    if len(Vobs) == 0:
        raise ValueError("At least one data point is required")
    if not (np.all(np.isfinite(Vobs)) and np.all(np.isfinite(errV)) and np.all(np.isfinite(r))):
        raise ValueError("Vobs, errV, and r must contain only finite values")
    if np.any(errV <= 0) or np.any(r <= 0):
        raise ValueError("errV and r must be strictly positive")
    if int(ngrid) < 2:
        raise ValueError("ngrid must be at least 2")
    if not np.isfinite(zpadding) or zpadding <= 0:
        raise ValueError("zpadding must be finite and positive")
    return Vobs, errV, r


def log_monotonic_probability(Vobs, errV, r, ngrid=400, zpadding=8.0):
    """Compute log P(g_1 > g_2 > ... > g_N), where g_i=V_i^2/r_i."""
    Vobs, errV, r = _validate_probability_inputs(Vobs, errV, r, ngrid, zpadding)
    N = len(Vobs)
    if N < 2:
        return 0.0

    # Transform each velocity to a standard normal and integrate the ordered
    # region recursively, avoiding a high-dimensional Monte Carlo integral.
    mu = Vobs / errV
    a = errV**2 / r
    k = np.sqrt(a[:-1] / a[1:])
    zmax = np.max(np.abs(mu)) + zpadding
    z = np.linspace(-zmax, zmax, int(ngrid))
    dz = z[1] - z[0]
    log_dz = np.log(dz)

    # Start with the innermost point and work outward through the recursion.
    upper = k[-1] * np.abs(z)
    TN = normal_cdf(upper - mu[-1]) - normal_cdf(-upper - mu[-1])
    TN = np.maximum(TN, 1e-300)
    logT = np.log(TN)

    # At each step, integrate only over values that preserve the required order.
    for i in range(N - 2, 0, -1):
        logphi = log_normal_pdf(z, mu[i])
        logT_new = np.empty_like(z)
        for j, zprev in enumerate(z):
            limit = k[i - 1] * abs(zprev)
            mask = np.abs(z) <= limit
            if not np.any(mask):
                logT_new[j] = -np.inf
                continue
            logT_new[j] = logsumexp(logphi[mask] + logT[mask]) + log_dz
        logT = logT_new

    logP = logsumexp(log_normal_pdf(z, mu[0]) + logT) + log_dz
    return float(min(logP, 0.0))


def _tie_groups(g_bar, tie_atol):
    """Group consecutive descending baryonic accelerations treated as ties."""
    order = np.argsort(-g_bar, kind="mergesort")
    groups = []
    current = [int(order[0])]
    for index in order[1:]:
        index = int(index)
        reference = g_bar[current[-1]]
        tied = (g_bar[index] == reference) if tie_atol == 0 else abs(g_bar[index] - reference) <= tie_atol
        if tied:
            current.append(index)
        else:
            groups.append(tuple(current))
            current = [index]
    groups.append(tuple(current))
    return groups


def log_monotonic_probability_for_gbar(g_bar, Vobs, errV, r, ngrid=400,
                                       zpadding=8.0, tie_atol=0.0,
                                       max_tie_permutations=720):
    """Return log monotonicity probability after ordering points by decreasing g_bar.

    Exact or tolerance-defined ties are handled by averaging the probabilities for
    all allowed within-tie orderings.  The averaging is done in log space, but is
    the arithmetic mean in probability space requested for the edge case.
    """
    g_bar = np.asarray(g_bar, dtype=float)
    if g_bar.ndim != 1 or not np.all(np.isfinite(g_bar)):
        raise ValueError("g_bar must be a one-dimensional finite array")
    Vobs, errV, r = _validate_probability_inputs(Vobs, errV, r, ngrid, zpadding)
    if len(g_bar) != len(Vobs):
        raise ValueError("g_bar, Vobs, errV, and r must have the same length")
    if not np.isfinite(tie_atol) or tie_atol < 0:
        raise ValueError("tie_atol must be finite and non-negative")
    if int(max_tie_permutations) < 1:
        raise ValueError("max_tie_permutations must be positive")

    groups = _tie_groups(g_bar, tie_atol)
    total = int(np.prod([len(group) for group in groups], dtype=np.int64))
    if total > int(max_tie_permutations):
        warnings.warn(
            f"{total} tied-order permutations exceed max_tie_permutations="
            f"{max_tie_permutations}; using stable descending order",
            RuntimeWarning,
        )
        order = np.concatenate([np.asarray(group, dtype=int) for group in groups])
        return log_monotonic_probability(Vobs[order], errV[order], r[order], ngrid, zpadding)

    # Enumerate only permutations within tied groups; all strictly ordered groups
    # retain their descending g_bar order, making the treatment deterministic.
    ordering_choices = [list(permutations(group)) for group in groups]
    log_probabilities = []
    for choice in product(*ordering_choices):
        order = np.concatenate([np.asarray(part, dtype=int) for part in choice])
        log_probabilities.append(log_monotonic_probability(
            Vobs[order], errV[order], r[order], ngrid, zpadding))
    return float(logsumexp(log_probabilities) - np.log(len(log_probabilities)))
