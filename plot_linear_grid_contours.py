#!/usr/bin/env python3
"""Plot sigma contours on existing linear-profile probability grids.

This only reads profile_grid.csv; it does not rerun the RAR analysis. Sigma levels
are one-sided Gaussian lower-tail probabilities of monotonicity, so contouring
log(P_mono) at log(Phi(-sigma)).
"""
import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

SIGMAS = (2, 5)
KINDS = ("conditional", "marginalised_gas_bulge")
REQUIRED = ["galaxy", "profile", "probability_kind", "grid_i", "grid_j",
            "upsilon_disk_mean", "disk_slope", "log_probability", "valid"]

def log_tail(sigma):
    return math.log(0.5 * math.erfc(float(sigma) / math.sqrt(2.0)))

def valid_mask(series):
    return series.astype(str).str.strip().str.lower().isin(("true", "1", "yes"))

def make_grid(group):
    group = group.copy()
    group["grid_i"] = pd.to_numeric(group["grid_i"], errors="raise").astype(int)
    group["grid_j"] = pd.to_numeric(group["grid_j"], errors="raise").astype(int)
    group["upsilon_disk_mean"] = pd.to_numeric(group["upsilon_disk_mean"], errors="coerce")
    group["disk_slope"] = pd.to_numeric(group["disk_slope"], errors="coerce")
    group["log_probability"] = pd.to_numeric(group["log_probability"], errors="coerce")
    group["_valid"] = valid_mask(group["valid"]) & np.isfinite(group["log_probability"])
    if (group["upsilon_disk_mean"] <= 0).any():
        raise ValueError(f"non-positive mean M/L for {group.galaxy.iloc[0]}")
    dup = group.duplicated(["grid_i", "grid_j"], keep=False)
    if dup.any():
        raise ValueError(f"duplicate grid cells for {group.galaxy.iloc[0]} / {group.probability_kind.iloc[0]}")
    x = np.log10(group.groupby("grid_j")["upsilon_disk_mean"].first().sort_index().to_numpy(float))
    y = group.groupby("grid_i")["disk_slope"].first().sort_index().to_numpy(float)
    zdf = group.pivot(index="grid_i", columns="grid_j", values="log_probability").sort_index().sort_index(axis=1)
    vdf = group.pivot(index="grid_i", columns="grid_j", values="_valid").sort_index().sort_index(axis=1)
    z = zdf.to_numpy(float)
    valid = vdf.to_numpy(bool) & np.isfinite(z)
    if z.shape != (len(y), len(x)) or min(z.shape) < 2:
        raise ValueError(f"incomplete/degenerate grid for {group.galaxy.iloc[0]}")
    return x, y, z, valid

def contour_segments(ax, x, y, z, valid, level, color, linestyle, linewidth=1.5):
    finite = z[valid]
    if finite.size == 0 or level < finite.min() or level > finite.max():
        return [], "not_bracketed"
    masked = np.ma.masked_where(~valid, z)
    try:
        cs = ax.contour(x, y, masked, levels=[level], colors=[color],
                        linestyles=[linestyle], linewidths=linewidth)
        segments = [np.asarray(s, dtype=float) for s in cs.allsegs[0] if len(s) >= 2]
    except (ValueError, TypeError):
        segments = []
    if not segments:
        return [], "no_segment"
    xtol = max(abs(x[-1] - x[0]), 1.0) * 1e-6
    ytol = max(abs(y[-1] - y[0]), 1.0) * 1e-6
    touches = any(np.any(np.isclose(s[:, 0], x[0], atol=xtol, rtol=0) |
                          np.isclose(s[:, 0], x[-1], atol=xtol, rtol=0) |
                          np.isclose(s[:, 1], y[0], atol=ytol, rtol=0) |
                          np.isclose(s[:, 1], y[-1], atol=ytol, rtol=0)) for s in segments)
    return segments, ("touches_grid_boundary" if touches else "contained")

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path("outputs/production_linear_51/profile_grid.csv"))
    parser.add_argument("--output-dir", type=Path, default=Path("outputs/production_linear_51/contour_plots"))
    parser.add_argument("--vmin", type=float, default=-10.0, help="lower colour limit for log probability; contours use unclipped values")
    parser.add_argument("--probability-kind", choices=(*KINDS, "all"), default="all")
    parser.add_argument("--galaxy", action="append", help="optional galaxy selection; may be repeated")
    args = parser.parse_args()
    if not args.vmin < 0:
        parser.error("--vmin must be below zero")
    if not args.input.is_file():
        parser.error(f"input does not exist: {args.input}")
    data = pd.read_csv(args.input, usecols=REQUIRED, low_memory=False)
    missing = set(REQUIRED) - set(data.columns)
    if missing:
        parser.error(f"missing input columns: {sorted(missing)}")
    data = data[data["profile"].astype(str).str.lower().eq("linear")].copy()
    selected_kinds = KINDS if args.probability_kind == "all" else (args.probability_kind,)
    data = data[data["probability_kind"].isin(selected_kinds)]
    if args.galaxy:
        data = data[data["galaxy"].isin(args.galaxy)]
        missing_galaxies = sorted(set(args.galaxy) - set(data["galaxy"]))
        if missing_galaxies:
            parser.error(f"galaxies not found in selected probability grids: {missing_galaxies}")
    if data.empty:
        parser.error("no matching linear grid rows")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    status_rows, vertex_rows = [], []
    kind_groups = {}
    for (kind, galaxy), group in data.groupby(["probability_kind", "galaxy"], sort=True):
        x, y, z, valid = make_grid(group)
        kind_groups.setdefault(kind, []).append((str(galaxy), x, y, z, valid))
        subdir = args.output_dir / kind
        heat_dir = subdir / "heatmaps"
        heat_dir.mkdir(parents=True, exist_ok=True)
        fig, ax = plt.subplots(figsize=(8, 6))
        cmap = plt.get_cmap("viridis").copy()
        cmap.set_bad("#d9d9d9")
        heat = np.ma.masked_where(~valid, np.clip(z, args.vmin, 0.0))
        mesh = ax.pcolormesh(x, y, heat, shading="nearest", cmap=cmap, vmin=args.vmin, vmax=0.0)
        for sigma, color, style in ((2, "deepskyblue", "solid"), (5, "orangered", "dashed")):
            level = log_tail(sigma)
            segments, status = contour_segments(ax, x, y, z, valid, level, color, style)
            outside = status in ("touches_grid_boundary", "not_bracketed")
            status_rows.append({"galaxy": galaxy, "probability_kind": kind, "sigma": sigma,
                                "log_probability_level": level, "probability_level": math.exp(level),
                                "status": status, "outside_plot_margin": outside,
                                "valid_min_log_probability": float(z[valid].min()) if valid.any() else math.nan,
                                "valid_max_log_probability": float(z[valid].max()) if valid.any() else math.nan})
            for seg_i, segment in enumerate(segments):
                for vertex_i, (xx, yy) in enumerate(segment):
                    vertex_rows.append({"galaxy": galaxy, "probability_kind": kind, "sigma": sigma,
                                        "segment": seg_i, "vertex": vertex_i,
                                        "log10_mean_upsilon_disk": xx, "disk_slope": yy})
        ax.set(xlabel=r"$\log_{10}(\Upsilon_{\rm disk})$", ylabel="Linear disk slope",
               title=f"{galaxy} — {kind.replace('_', ' ')}")
        ax.legend(handles=[Line2D([0], [0], color="deepskyblue", ls="-", label="2σ (one-sided)"),
                           Line2D([0], [0], color="orangered", ls="--", label="5σ (one-sided)")], loc="best")
        fig.colorbar(mesh, ax=ax, label="log P(monotonic), colour clipped below vmin")
        fig.tight_layout()
        fig.savefig(heat_dir / f"{galaxy}.png", dpi=140)
        plt.close(fig)
    status_df = pd.DataFrame(status_rows)
    vertices_df = pd.DataFrame(vertex_rows, columns=["galaxy", "probability_kind", "sigma", "segment",
                                                     "vertex", "log10_mean_upsilon_disk", "disk_slope"])
    status_df.to_csv(args.output_dir / "contour_status.csv", index=False)
    vertices_df.to_csv(args.output_dir / "contour_vertices.csv", index=False)
    outside_df = status_df[status_df["outside_plot_margin"]].copy()
    outside_df.to_csv(args.output_dir / "outside_plot_margins.csv", index=False)
    for kind, grids in kind_groups.items():
        fig, ax = plt.subplots(figsize=(10, 8))
        for galaxy, x, y, z, valid in grids:
            for sigma, color, style in ((2, "#0077bb", "solid"), (5, "#cc3311", "dashed")):
                segs = vertices_df[(vertices_df.galaxy == galaxy) & (vertices_df.probability_kind == kind) & (vertices_df.sigma == sigma)]
                for _, seg in segs.groupby("segment", sort=True):
                    ax.plot(seg.log10_mean_upsilon_disk, seg.disk_slope, color=color, ls=style, lw=0.8, alpha=0.42)
        ax.set_xlim(min(g[1][0] for g in grids), max(g[1][-1] for g in grids))
        ax.set_ylim(min(g[2][0] for g in grids), max(g[2][-1] for g in grids))
        kind_status = status_df[status_df.probability_kind == kind]
        summaries = []
        for sigma in SIGMAS:
            n_out = int(kind_status[(kind_status.sigma == sigma) & kind_status.outside_plot_margin].shape[0])
            summaries.append(f"{sigma}σ outside/edge: {n_out}/{len(grids)}")
        ax.text(0.02, 0.98, "\n".join(summaries) + "\nNames/status: outside_plot_margins.csv",
                transform=ax.transAxes, va="top", ha="left", fontsize=9,
                bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none"})
        ax.set(xlabel=r"$\log_{10}(\Upsilon_{\rm disk})$", ylabel="Linear disk slope",
               title=f"Linear-profile {kind.replace('_', ' ')}: 2σ and 5σ contours ({len(grids)} galaxies)")
        ax.legend(handles=[Line2D([0], [0], color="#0077bb", ls="-", label="2σ"),
                           Line2D([0], [0], color="#cc3311", ls="--", label="5σ")], loc="best")
        fig.tight_layout()
        fig.savefig(args.output_dir / f"all_galaxies_{kind}_contours.png", dpi=180)
        plt.close(fig)
    for kind, grids in kind_groups.items():
        kind_status = status_df[status_df.probability_kind == kind]
        for sigma, color, style in ((2, "#0077bb", "solid"), (5, "#cc3311", "dashed")):
            fig, ax = plt.subplots(figsize=(10, 8))
            for galaxy, x, y, z, valid in grids:
                segs = vertices_df[(vertices_df.galaxy == galaxy) &
                                   (vertices_df.probability_kind == kind) &
                                   (vertices_df.sigma == sigma)]
                for _, seg in segs.groupby("segment", sort=True):
                    ax.plot(seg.log10_mean_upsilon_disk, seg.disk_slope,
                            color=color, ls=style, lw=0.8, alpha=0.42)
            ax.set_xlim(min(g[1][0] for g in grids), max(g[1][-1] for g in grids))
            ax.set_ylim(min(g[2][0] for g in grids), max(g[2][-1] for g in grids))
            n_out = int(kind_status[(kind_status.sigma == sigma) &
                                    kind_status.outside_plot_margin].shape[0])
            ax.text(0.02, 0.98,
                    f"{sigma}σ outside/edge: {n_out}/{len(grids)}\nNames/status: outside_plot_margins.csv",
                    transform=ax.transAxes, va="top", ha="left", fontsize=9,
                    bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none"})
            ax.set(xlabel="log10(mean disk M/L)", ylabel="Linear disk slope",
                   title=f"Linear-profile {kind.replace('_', ' ')}: {sigma}σ contours ({len(grids)} galaxies)")
            ax.legend(handles=[Line2D([0], [0], color=color, ls=style, label=f"{sigma}σ")],
                      loc="best")
            fig.tight_layout()
            fig.savefig(args.output_dir / f"all_galaxies_{kind}_{sigma}sigma_contours.png", dpi=180)
            plt.close(fig)
    meta = {"input": str(args.input), "profile": "linear", "galaxy_count": int(data.galaxy.nunique()),
            "probability_kinds": sorted(kind_groups), "heatmap_log_probability_vmin": args.vmin,
            "contours_use_unclipped_log_probability": True, "sigma_convention": "one-sided Gaussian lower tail P(monotonic)=Phi(-sigma)",
            "sigma_log_probability_levels": {str(s): log_tail(s) for s in SIGMAS},
            "outside_statuses": ["touches_grid_boundary", "not_bracketed"],
            "created_utc": datetime.now(timezone.utc).isoformat()}
    (args.output_dir / "run_metadata.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(f"Plotted {data.galaxy.nunique()} galaxies; kinds: {', '.join(sorted(kind_groups))}")
    print(f"Heat maps and contour products: {args.output_dir}")
    print(f"Outside/edge records: {len(outside_df)} (see outside_plot_margins.csv)")

if __name__ == "__main__":
    main()
