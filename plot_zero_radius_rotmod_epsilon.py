
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import mono_SPARC_RAR as rar

GALAXIES = ("IC4202", "NGC0891", "UGC09133", "UGC11914")
EPSILON_FRACTIONS = (1e-3, 1e-6, 1e-9)
ROT_COLUMNS = ["Rad", "Vobs", "errV", "Vgas", "Vdisk", "Vbul", "SBdisk", "SBbul"]
OUT_DIR = Path(__file__).resolve().parent / "outputs" / "zero_radius_rotmod_ml1_plots"


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    table = rar.load_sparc_table()
    print("galaxy case n MSE_km2_s-2 RMSE_km_s")
    for galaxy in GALAXIES:
        info = table.loc[table["Galaxy"] == galaxy].iloc[0]
        radius, sb = rar.get_SBdisk(galaxy, table=table)
        radius = np.asarray(radius, dtype=float)
        sb = np.asarray(sb, dtype=float)
        if radius[0] != 0 or len(radius) < 2:
            raise ValueError(f"{galaxy}: expected a .dens grid beginning at R=0")
        rdisk = float(info["Rdisk"])
        z0 = 0.196 * rdisk**0.633
        rot_path = Path(rar.DEFAULT_DATA_ROOT) / f"{galaxy}_rotmod.dat"
        rot = pd.read_csv(rot_path, sep=r"\s+", comment="#", header=None, names=ROT_COLUMNS)
        r_obs = rot["Rad"].to_numpy(dtype=float)
        v_sparc = rot["Vdisk"].to_numpy(dtype=float)
        fig, ax = plt.subplots(figsize=(8.2, 5.4))
        plotted = False
        finite_reference = np.isfinite(r_obs) & np.isfinite(v_sparc)
        ax.scatter(r_obs[finite_reference], v_sparc[finite_reference], s=20, marker="s", zorder=5,
                   label=r"raw SPARC $V_{disk}$ ($\Upsilon_d=1$)")

        # Reference case: remove the R=0 density sample and solve on the untouched
        # positive-radius grid, without synthesizing an inner radius.
        r_model = radius[1:]
        sb_model = sb[1:]
        velocity = np.asarray(rar.vcdisk(r_model, sb_model, z0=z0, rhoz="exp"), dtype=float)
        force = np.asarray(rar.signed_disk_term(velocity), dtype=float)
        good = np.isfinite(r_model) & np.isfinite(force)
        covered = (finite_reference & (r_obs >= r_model[good].min())
                   & (r_obs <= r_model[good].max())) if good.any() else np.zeros_like(finite_reference)
        if good.sum() >= 2 and covered.any():
            r_common = r_obs[covered]
            force_common = rar.interpolate_disk_force(r_model[good], force[good], r_common)
            v_model = np.sign(force_common) * np.sqrt(np.abs(force_common))
            v_ref = v_sparc[covered]
            valid = np.isfinite(v_model) & np.isfinite(v_ref)
            r_common, v_model, v_ref = r_common[valid], v_model[valid], v_ref[valid]
            if not len(r_common):
                continue
            mse = float(np.mean((v_model - v_ref) ** 2))
            ax.plot(r_common, v_model, marker="o", markersize=3, linewidth=1.6, alpha=0.9,
                    label=fr"$R=0$ removed; MSE={mse:.2f} (km/s)$^2$")
            plotted = True
            print(f"{galaxy} R0_removed {len(r_common)} {mse:.6g} {np.sqrt(mse):.6g}")
        else:
            print(f"SKIP {galaxy} R0_removed: fewer than two finite radii or no shared coverage")

        for frac in EPSILON_FRACTIONS:
            r_model = np.r_[frac * radius[1], radius[1:]]
            sb_model = np.r_[sb[0], sb[1:]]
            velocity = np.asarray(rar.vcdisk(r_model, 1.0 * sb_model, z0=z0, rhoz="exp"), dtype=float)
            force = np.asarray(rar.signed_disk_term(velocity), dtype=float)
            good = np.isfinite(r_model) & (r_model > 0) & np.isfinite(force)
            if good.sum() < 2:
                print(f"SKIP {galaxy} epsilon={frac:g}: fewer than two finite model radii")
                continue
            covered = (np.isfinite(r_obs) & np.isfinite(v_sparc)
                       & (r_obs >= r_model[good].min()) & (r_obs <= r_model[good].max()))
            if not covered.any():
                print(f"SKIP {galaxy} epsilon={frac:g}: no shared radial coverage")
                continue
            r_common = r_obs[covered]
            force_common = rar.interpolate_disk_force(r_model[good], force[good], r_common)
            v_model = np.sign(force_common) * np.sqrt(np.abs(force_common))
            valid = np.isfinite(v_model) & np.isfinite(v_sparc[covered])
            r_common, v_model, v_ref = r_common[valid], v_model[valid], v_sparc[covered][valid]
            if not len(r_common):
                continue
            mse = float(np.mean((v_model - v_ref) ** 2))
            ax.plot(r_common, v_model, marker="o", markersize=3, linewidth=1.4, alpha=0.5,
                    label=fr"$\epsilon={frac:g}R_1$; MSE={mse:.2f} (km/s)$^2$")
            plotted = True
            print(f"{galaxy} {frac:g} {len(r_common)} {mse:.6g} {np.sqrt(mse):.6g}")
        if not plotted:
            plt.close(fig)
            raise RuntimeError(f"{galaxy}: no finite comparison could be plotted")
        ax.set(title=f"{galaxy}: SPARC disk speed vs vcdisk (disk M/L=1.0)", xlabel="Radius (kpc)", ylabel=r"$V_{disk}$ (km s$^{-1}$)")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(OUT_DIR / f"{galaxy}_vcdisk_vs_sparc_ml1.png", dpi=180)
        plt.close(fig)


if __name__ == "__main__":
    main()
