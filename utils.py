#!/usr/bin/env python
# from jax import numpy as jnp
import pandas as pd
import numpy as np
from scipy.special import erf, logsumexp
from scipy.integrate import trapezoid

def load_sparc_table():
    """
    Load the SPARC table from the .mrt file.
    Returns a pandas DataFrame with the relevant columns.
    """
    SPARC_c = [ "Galaxy", "T", "D", "e_D", "f_D", "Inc",
                "e_Inc", "L", "e_L", "Reff", "SBeff", "Rdisk",
                "SBdisk", "MHI", "RHI", "Vflat", "e_Vflat", "Q", "Ref."]
    file = "/mnt/users/koe/mono_SPARC_RAR/SPARC_table.txt"
    table = pd.read_fwf(file, skiprows=98, names=SPARC_c)
    return table


def vel2acc(vel_sq, rad):
    """
    Calculate the gravitational acceleration (g_obs or g_bar).
    In the data, vel is in km/s, and rad is in kpc.
    Returns g_x in m/s².
    """
    return (vel_sq * 1e6) / (rad * 3.086e19)


def get_SBdisk(galaxy:str):
    """
    Extract surface brightness profile of disk component
    from SPARC .dens file for vcdisk calculation.
    
    Parameters
    ----------
    galaxy : str
        Name of the galaxy.
    inc : float
        Inclination in degrees.
    
    Returns
    -------
    rad (kpc), SBdisk (Lsun/kpc^2) : 1D np arrays
    """

    file =  f"/mnt/users/koe/SPARC_RAR/BulgeDiskDec_LTG/{galaxy}.dens"
    columns = [ "Rad", "SBdisk", "SBbul" ]
    data = pd.read_csv(file, names=columns, skiprows=1, sep="\t")

    rad = data["Rad"].to_numpy()    # in kpc

    # Normalise SBdisk such that total disk luminosity matches L[3.6] in SPARC data (see email from Federico)
    SBdisk_raw = data["SBdisk"].to_numpy() * 1e6
    Ltot, _ = compute_A_B(rad, SBdisk_raw)

    table = load_sparc_table()
    i_table = np.where(table["Galaxy"] == galaxy)[0][0]
    L = table["L"][i_table] * 1e9   # in Lsun

    SBdisk = SBdisk_raw * L / Ltot

    return rad, SBdisk


def compute_A_B(R, I_R):
    """
    Compute A = 2π ∫ I(R) R dR,
            B = 2π ∫ I(R) R^2 dR
    for arbitrary nonuniform SPARC radii R.
    Used in linear M/L profile.
    
    Parameters
    ----------
    R   : radii (1D array) in kpc
    I_R : surface brightness (Lsun/kpc^2) at each radius
    
    Returns
    -------
    A, B : floats
    """

    # geometric factor
    fac = 2 * np.pi

    # trapezoidal integrals
    A = fac * trapezoid(I_R * R, R)
    B = fac * trapezoid(I_R * R**2, R)

    return A, B


def compute_A_beta(R, I_R, log_slope):
    """
    Compute A = 2π ∫ I(R) R dR,
            beta = 2π ∫ I(R) e^(bR) R dR
    for arbitrary nonuniform SPARC radii R.
    Used in linear M/L profile.
    
    Parameters
    ----------
    R   : radii (1D array) in kpc
    I_R : surface brightness (Lsun/kpc^2) at each radius
    log_slope : float
        Logarithmic slope of the M/L profile.
    
    Returns
    -------
    A, beta : floats
    """
    
    # geometric factor
    fac = 2 * np.pi

    # trapezoidal integrals
    A = fac * trapezoid(I_R * R, R)
    beta = fac * trapezoid(I_R * np.exp(log_slope * R) * R, R)

    return A, beta


# Monotonic likelihood functions for SPARC RAR data
LOG_SQRT_2PI = 0.5 *np.log(2.0 *np.pi)

def log_normal_pdf(x, mu):
    """
    log N(x | mu, 1)
    """
    return -0.5 * (x - mu)**2 - LOG_SQRT_2PI

def normal_cdf(x):
    return 0.5 * (1.0 + erf(x /np.sqrt(2.0)))

def log_monotonic_probability( Vobs, errV, r, ngrid=400, zpadding=8.0 ):
    """
    Computes
        log P(g1 > g2 > ... > gN)
    where
        g_i = V_i^2 / r_i
    and
        V_i ~ N(Vobs_i, errV_i^2)

    Parameters
    ----------
    Vobs : array_like
    errV : array_like
    r : array_like

    ngrid : int
        Number of grid points

    zpadding : float
        Number of sigma beyond largest |mu|
        included in integration grid

    Returns
    -------
    logP : float
    """

    Vobs = np.asarray(Vobs)
    errV = np.asarray(errV)
    r = np.asarray(r)

    N = len(Vobs)

    if N < 2:
        return 0.0

    mu = Vobs / errV
    a = errV**2 / r
    k = np.sqrt(a[:-1] / a[1:])

    # ------------------------------------------
    # Global z-grid
    # ------------------------------------------

    zmax = np.max(np.abs(mu)) + zpadding
    z = np.linspace(-zmax, zmax, ngrid)
    dz = z[1] - z[0]
    log_dz = np.log(dz)

    # ------------------------------------------
    # Base case: T_N(z_{N-1})
    #
    # T_N(z) =
    # P(|Z_N| < k |z|)
    # ------------------------------------------

    upper = k[-1] * np.abs(z)

    TN = (
        normal_cdf(upper - mu[-1])
        - normal_cdf(-upper - mu[-1])
    )

    TN = np.maximum(TN, 1e-300)

    logT = np.log(TN)

    # ------------------------------------------
    # Backward recursion
    # ------------------------------------------

    for i in range(N - 2, 0, -1):
        logphi = log_normal_pdf(z, mu[i])
        logT_new = np.empty_like(z)

        for j, zprev in enumerate(z):
            limit = k[i - 1] * abs(zprev)
            mask = np.abs(z) <= limit

            if not  np.any(mask):
                logT_new[j] = -np.inf
                continue

            log_integrand = (
                logphi[mask]
                + logT[mask]
            )

            logT_new[j] = (
                logsumexp(log_integrand)
                + log_dz
            )

        logT = logT_new

    # ------------------------------------------
    # Final integral
    #
    # P = ∫ phi(z1) T2(z1) dz1
    # ------------------------------------------

    logphi1 = log_normal_pdf(z, mu[0])
    log_integrand = logphi1 + logT
    logP = logsumexp(log_integrand) + log_dz

    return logP
