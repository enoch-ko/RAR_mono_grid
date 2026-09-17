#!/usr/bin/env python
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from vcdisk import vcdisk
from utils import load_sparc_table, get_SBdisk, compute_A_B, compute_A_beta, vel2acc, log_monotonic_probability

testing = False
profile = "constant"  # Available profiles are "log_linear", "linear", and "constant"

valid_profiles = ["log_linear", "linear", "constant"]
if profile not in valid_profiles:
    raise ValueError(f"Invalid profile: {profile}. Must be one of {valid_profiles}")

# grid_size = 21 if (linear_profile and not testing) else 51
grid_size = 51

table = load_sparc_table()
galaxy_count = len(table["Galaxy"])
columns = [ "Rad", "Vobs", "errV", "Vgas",
            "Vdisk", "Vbul", "SBdisk", "SBbul" ]

Upsilon_gas = 1.0
Upsilon_bulge = 0.7

log_disk_SPARC = np.log10(0.5)
log_disk = np.linspace(log_disk_SPARC-1.0, log_disk_SPARC+1.0, grid_size, endpoint=True)
ups_disk = 10**log_disk

if testing: galaxy_count = 1

# for i in tqdm(range(galaxy_count), desc="SPARC galaxies"):
for i in range(galaxy_count):
    g = table["Galaxy"][i]
    print(f"Processing galaxy {i+1}/{galaxy_count}: {g}")

    # if not testing and i == 0: continue

    data_path = "/mnt/users/koe/data/"+g+"_rotmod.dat"
    rawdata = np.loadtxt(data_path)
    data = pd.DataFrame(rawdata, columns=columns)

    r = data["Rad"].to_numpy()
    Vobs, errV = data["Vobs"].to_numpy(), data["errV"].to_numpy()
    Vgas, Vbul = data["Vgas"].to_numpy(), data["Vbul"].to_numpy()

    if profile == "constant":
        log_prob_array = np.zeros(grid_size)
        Vdisk = data["Vdisk"].to_numpy()

        for p, pdisk in enumerate(ups_disk):
            Vbar_squared = Vdisk**2 * pdisk + Vgas*np.abs(Vgas) * Upsilon_gas + Vbul**2 * Upsilon_bulge
            g_bar = vel2acc(Vbar_squared, r)

            idx = np.flip(np.argsort(g_bar))
            Vobs_sorted, errV_sorted, r_sorted = Vobs[idx], errV[idx], r[idx]

            log_prob = log_monotonic_probability(Vobs_sorted, errV_sorted, r_sorted)
            log_prob_array[p] = log_prob

        # Plot log_prob_array as coloured 1D plot
        log_prob_plot = np.maximum(np.minimum(log_prob_array, 0.0), -50.0)
        fig, ax = plt.subplots(figsize=(7, 5))
        ax.plot(log_disk, log_prob_plot, linestyle='-', color='blue')
        ax.set_xlabel('log10(Upsilon_disk)')
        ax.set_ylabel('log(probability)')
        ax.set_title(f'Log Probability - {g}')
    
        plt.tight_layout()
        # plt.show()
        plt.savefig(f'/mnt/users/koe/mono_SPARC_RAR/plots/constant_ML_plots/{g}.png', dpi=150, bbox_inches='tight')
        plt.close()

    else:
        i_tab = np.where(table["Galaxy"] == g)[0][0]
        Rdisk = table["Rdisk"][i_tab]
        zd = 0.196 * Rdisk**0.633   # SPARC scale height relation (see Lelli+2016)

        rad, SBdisk = get_SBdisk(g)
        log_prob_matrix = np.zeros((grid_size, grid_size))
        negative_upsilon_mask = np.zeros((grid_size, grid_size), dtype=bool)

        disk_slope = np.linspace(-5.0, 5.0, grid_size, endpoint=True)
        A, B = compute_A_B(rad, SBdisk)

        for p, pdisk in enumerate(ups_disk):
            for m, slope in enumerate(disk_slope):

                if profile == "linear":
                    y_intercept = pdisk - slope * (B / A)
                    Upsilon_disk = slope * rad + y_intercept
                elif profile == "log_linear":
                    _, beta = compute_A_beta(rad, SBdisk, slope)
                    norm_scale = pdisk / (beta / A)
                    Upsilon_disk = norm_scale * np.exp(slope * rad)

                negative_upsilon_mask[m, p] = np.any(Upsilon_disk < 0)
                surface_density = SBdisk * Upsilon_disk

                Vdisk_dens = vcdisk(rad, surface_density, z0=zd, rhoz='exp')
                Vdisk = np.interp(r, rad, Vdisk_dens)
                Vbar_squared = Vdisk**2 + Vgas*np.abs(Vgas) * Upsilon_gas + Vbul**2 * Upsilon_bulge
                g_bar = vel2acc(Vbar_squared, r)
                
                idx = np.flip(np.argsort(g_bar))
                Vobs_sorted, errV_sorted, r_sorted = Vobs[idx], errV[idx], r[idx]

                log_prob = log_monotonic_probability(Vobs_sorted, errV_sorted, r_sorted)
                log_prob_matrix[m, p] = log_prob
    
        # Plot log_prob_matrix as coloured 2D histogram with contours
        log_prob_plot = np.maximum(np.minimum(log_prob_matrix, 0.0), -50.0)
        log_levels = np.linspace(np.min(log_prob_plot), 0.0, 11)
        fig, ax = plt.subplots(figsize=(10, 8))
        im = ax.contourf(log_disk, disk_slope, log_prob_plot, levels=log_levels, cmap='viridis', vmin=log_levels[0], vmax=0.0)

        # Shade region where the model predicts a negative mass-to-light ratio for the stellar disk (Upsilon_disk <= 0).
        neg_mask_plot = np.ma.masked_where(~negative_upsilon_mask, np.zeros_like(negative_upsilon_mask, dtype=float))
        ax.imshow(neg_mask_plot, origin='lower', extent=[log_disk[0], log_disk[-1], disk_slope[0], disk_slope[-1]],
                  cmap='gray', alpha=0.5, aspect='auto', interpolation='nearest', zorder=2)
        
        # Add contours
        contours = ax.contour(log_disk, disk_slope, log_prob_plot, levels=log_levels, colors='black', alpha=0.3, linewidths=0.5)
        ax.clabel(contours, inline=True, fontsize=8)

        # Overlay red star at (log_disk_SPARC, 0.0)
        ax.plot(log_disk_SPARC, 0.0, marker='*', color='red', markersize=14, markeredgecolor='black')
        
        ax.set_xlabel('log10(Upsilon_disk)')
        if profile == "linear": ax.set_ylabel('disk_slope')
        elif profile == "log_linear": ax.set_ylabel('log_disk_slope')
        ax.set_title(f'Log Probability Matrix - {g}')
        
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('log(probability)')

        plt.tight_layout()
        if profile == "linear":
            plt.savefig(f'/mnt/users/koe/mono_SPARC_RAR/plots/linear_2D_grids/{g}.png', dpi=150, bbox_inches='tight')
        elif profile == "log_linear":
            plt.savefig(f'/mnt/users/koe/mono_SPARC_RAR/plots/log_linear_2D_grids/{g}.png', dpi=150, bbox_inches='tight')
        plt.close()
