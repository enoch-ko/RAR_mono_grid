#!/usr/bin/env python3
"""Three-panel inverse-RAR diagnostic plot; exploratory, not a manuscript figure.

For a fitted profile, both RAR tracks use direct Casertano disc solves. The
SPARC tabulated disc is a separate reference. If no bounded fit is identified,
only SPARC-based curves are shown; no optimized M/L curve is invented.
"""
import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from utils import vel2acc


def signed_velocity(signed_v2):
    """Display signed W=V|V| without silently losing outward-force signs."""
    return np.sign(signed_v2) * np.sqrt(np.abs(signed_v2))


def direct_disc(diagnostic, rad, brightness, zd, r, ml_profile):
    response, warnings, error = diagnostic.mono.solve_disk_basis(
        rad, brightness * ml_profile, zd, r)
    if error or response is None or warnings or not np.all(np.isfinite(response)):
        raise ValueError(f"direct disc solve failed: {error}; warnings={warnings}")
    return np.asarray(response, float)


def render_plot(result, output, data_root, dens_root, diagnostic,
                replace_existing=False):
    """Save one PNG and return auditable plot QC for the diagnostic JSON."""
    output = Path(output)
    if output.exists() and not replace_existing:
        raise FileExistsError(output)
    galaxy = result["galaxy"]
    rad, brightness, r, vobs, errv, vgas, vbul, zd, excluded = diagnostic.galaxy_data(
        galaxy, data_root, dens_root)
    np.testing.assert_allclose(r, result["retained_radii_kpc"], rtol=0, atol=1e-12)
    raw = np.loadtxt(Path(data_root) / f"{galaxy}_rotmod.dat", ndmin=2)
    mask = np.isfinite(raw[:, 0]) & (raw[:, 0] >= rad[0]) & (raw[:, 0] <= rad[-1])
    np.testing.assert_allclose(raw[mask, 0], r, rtol=0, atol=1e-12)

    reference_ml = float(result["reference_ml"])
    gas_ml = float(result["gas_ml_qp"])
    bulge_ml = float(result["bulge_ml_qp"])
    gas_w = gas_ml * vgas * np.abs(vgas)
    bulge_w = bulge_ml * vbul**2
    # SPARC rotmod Vdisk is tabulated at M/L=1; rescale its signed V|V| term.
    sparc_disc = raw[mask, 3] * np.sqrt(reference_ml)
    sparc_vbar2 = sparc_disc * np.abs(sparc_disc) + gas_w + bulge_w
    sparc_valid = np.isfinite(sparc_vbar2) & (sparc_vbar2 > 0)
    gobs = vel2acc(vobs**2, r)
    has_fit = "xi" in result.get("optimization", {})
    qc = {"status": "created", "path": str(output),
          "kind": "optimized" if has_fit else "sparc_reference_only",
          "retained_points": len(r), "excluded_kinematic_rows": excluded,
          "sparc_rar_valid_points": int(np.count_nonzero(sparc_valid)),
          "sparc_disc_table_ml": 1.0, "sparc_disc_display_ml": reference_ml,
          "rar_reference": "Casertano constant M/L" if has_fit else "SPARC rotmod constant M/L"}

    fig = plt.figure(figsize=(16, 8.4))
    grid = fig.add_gridspec(2, 2, width_ratios=[1.12, 1],
                            height_ratios=[2.25, 1], wspace=0.27, hspace=0.06)
    axv = fig.add_subplot(grid[0, 0])
    axml = fig.add_subplot(grid[1, 0], sharex=axv)
    axrar = fig.add_subplot(grid[:, 1])
    blue, violet, green = "#2874A6", "#8E44AD", "#148F77"
    axv.errorbar(r, vobs, yerr=errv, fmt="ko", ms=4, capsize=2,
                 label=r"Observed $V_{obs}$", zorder=7)
    axv.plot(r, vgas, color="#16A085", lw=1.6, label=r"Gas $V_{gas}$")
    axv.plot(r, signed_velocity(bulge_w), color="#D68910", lw=1.6,
             label=r"Bulge $V_{bulge}$")
    axv.plot(r, sparc_disc, color="0.45", lw=1.8, ls=":",
             label=f"SPARC disc, M/L={reference_ml:g}")
    axml.axhline(reference_ml, color=blue if has_fit else "0.45", ls="--",
                 lw=1.8, label=f"Reference M/L={reference_ml:g}")
    axml.set(xlabel="Radius [kpc]", ylabel="Disc M/L")
    axml.set_xlim(float(rad[0]), float(rad[-1]))
    axv.set(ylabel="Circular speed [km s$^{-1}$]",
            title="Rotation-curve decomposition")
    axv.tick_params(labelbottom=False)
    axv.grid(alpha=0.2)
    axml.grid(alpha=0.2)
    axrar.set(xscale="log", yscale="log",
              xlabel=r"Baryonic acceleration $g_{bar}$ [m s$^{-2}$]",
              ylabel=r"Observed acceleration $g_{obs}$ [m s$^{-2}$]",
              title="RAR tracks, connected in radius order")
    axrar.grid(alpha=0.25, which="both")

    if has_fit:
        knots = np.asarray(result["knots_kpc"], float)
        xi = np.asarray(result["optimization"]["xi"], float)
        profile = diagnostic.hat_matrix(rad, knots) @ xi
        before_w = direct_disc(diagnostic, rad, brightness, zd, r,
                               np.full_like(rad, reference_ml))
        after_w = direct_disc(diagnostic, rad, brightness, zd, r, profile)
        before_vbar2 = before_w + gas_w + bulge_w
        after_vbar2 = after_w + gas_w + bulge_w
        if np.any(before_vbar2 <= 0) or np.any(after_vbar2 <= 0):
            raise ValueError("non-positive baryonic velocity square in fitted RAR")
        before_gbar = vel2acc(before_vbar2, r)
        after_gbar = vel2acc(after_vbar2, r)
        expected = np.asarray(result["optimization"]["gbar_m_s2"], float)
        delta = float(np.max(np.abs(after_gbar - expected)) /
                      np.max(np.abs(expected)))
        if delta > 1e-8:
            raise ValueError(f"optimized plot differs from saved gbar: {delta}")
        qc["direct_vs_saved_gbar_relative_to_peak"] = delta
        qc["optimized_ml_range"] = [float(np.min(profile)), float(np.max(profile))]
        axv.plot(r, signed_velocity(before_w), color=blue, lw=2, ls="--",
                 label=f"Casertano disc, M/L={reference_ml:g}")
        axv.plot(r, signed_velocity(after_w), color=violet, lw=2.2,
                 label="Casertano disc, optimized M/L(R)")
        axv.plot(r, np.sqrt(before_vbar2), color=green, lw=2, ls="--",
                 label=r"Before $V_{bar}$")
        axv.plot(r, np.sqrt(after_vbar2), color="#5B2C6F", lw=2.3,
                 label=r"After $V_{bar}$")
        axml.plot(rad, profile, color=violet, lw=2, label="Optimized M/L(R)")
        axrar.plot(before_gbar, gobs, "o--", color=green, ms=4, lw=1.7,
                   label=f"Before: Casertano M/L={reference_ml:g}")
        axrar.plot(after_gbar, gobs, "o-", color="#5B2C6F", ms=4, lw=2,
                   label="After: optimized M/L(R)")
        ties = len(result["optimization"].get("active_order_pairs", []))
        note = f"Weak-order fit; active ties: {ties}"
        subtitle = "Casertano before/after tracks; SPARC disc is a separate tabulated reference."
    else:
        axv.plot(r, np.where(sparc_valid, np.sqrt(np.clip(sparc_vbar2, 0, None)), np.nan),
                 color="#5B2C6F", lw=2.3, label=r"SPARC-based $V_{bar}$")
        axml.text(0.5, 0.78, "No feasible M/L profile identified",
                  transform=axml.transAxes, ha="center", va="center",
                  fontsize=11, color="#A93226")
        if np.any(sparc_valid):
            axrar.plot(vel2acc(sparc_vbar2[sparc_valid], r[sparc_valid]),
                       gobs[sparc_valid], "o-", color="#5B2C6F", ms=4, lw=2,
                       label=f"SPARC reference, M/L={reference_ml:g}")
        if result["status"] == "infeasible":
            note = ("No feasible M/L profile identified "
                    f"(upper bound {result['bounds_ml'][1]:g})")
        else:
            note = f"No optimized M/L profile: {result['status']}"
        subtitle = "SPARC-based reference only; no Casertano fit or before/after comparison."
    axrar.legend(loc="best", fontsize=9, frameon=False)
    axrar.text(0.97, 0.03, note, transform=axrar.transAxes,
               ha="right", va="bottom", fontsize=9,
               bbox={"facecolor": "white", "alpha": 0.88, "edgecolor": "none"})
    fig.suptitle(f"{galaxy}: inverse disc-M/L diagnostic", fontsize=15)
    handles, labels = axv.get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, 0.08),
               ncol=4, fontsize=8, frameon=False)
    fig.text(0.5, 0.025, subtitle + " Diagnostic only; not a converged science result.",
             ha="center", fontsize=8, color="0.35")
    fig.subplots_adjust(left=0.07, right=0.98, top=0.90, bottom=0.24)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    plt.close(fig)
    return qc


def main():
    import inverse_rar_diagnostic as diagnostic
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", type=Path, default=Path(
        "outputs/inverse_rar_diagnostic/pilots/inverse_UGC07399_densknots_ml10_final_20261005.json"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--replace-existing", action="store_true")
    args = parser.parse_args()
    if args.output.with_suffix(".plotmeta.json").exists() and not args.replace_existing:
        raise FileExistsError(args.output.with_suffix(".plotmeta.json"))
    result = json.loads(args.result.read_text())
    qc = render_plot(result, args.output, diagnostic.mono.DEFAULT_DATA_ROOT,
                     diagnostic.DEFAULT_DENS_ROOT, diagnostic, args.replace_existing)
    provenance = {"created_utc": datetime.now(timezone.utc).isoformat(),
                  "galaxy": result["galaxy"], "result": str(args.result),
                  "result_sha256": hashlib.sha256(args.result.read_bytes()).hexdigest(),
                  "plot_script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  **qc}
    args.output.with_suffix(".plotmeta.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(args.output)


if __name__ == "__main__":
    main()
