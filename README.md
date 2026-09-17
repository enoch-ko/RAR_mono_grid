# Grid-based Monotonicity Probability Mapping of SPARC RARs

Purpose: Instead of MCMC, we perform exhaustive grid searches over M/L parameter space to
compute the log-probability of monotonicity in the RAR, i.e., P(g_obs increasing with 
g_bar), for all 175 SPARC galaxies.

Key Components:
- mono_SPARC_RAR.py — Main grid search script. Three M/L profile modes:
  - constant: 1D grid over Upsilon_disk (51 points) → 1D log-probability curve per galaxy.
  - linear: 2D grid over (Upsilon_disk, disk_slope) → 2D log-probability contour map.
  - log_linear: 2D grid over (Upsilon_disk, log_disk_slope) → 2D contour map with exponential
M/L profile.
  - Uses vcdisk (from separate module) for extracting disk rotation curves from surface brightness.
  - log_monotonic_probability — Computes P(g_1 > g_2 > ... > g_N) via numerical integration
over velocity error distributions (backward recursion on a z-grid).
- utils.py — SPARC table loading, vel2acc, get_SBdisk, compute_A_B, compute_A_beta, and the
core log_monotonic_probability implementation (exact probabilistic monotonicity test).
- RAR_plot.py — Interactive matplotlib demo for visualizing monotonicity ranking
(Spearman/Pearson) — appears to be a standalone exploratory tool.
- SPARC_table.txt — Local copy of the SPARC catalog.
- ML_Andreea.csv — MIGHTEE-HI + LADUMA M/L reference data from Andreea Vărăşteanu (paper: https://arxiv.org/abs/2608.03576)

Output: Per-galaxy PNG plots in plots/constant_ML_plots/, plots/linear_2D_grids/,
plots/log_linear_2D_grids/ showing log-probability landscapes with the SPARC fiducial point
(Upsilon_disk = 0.5, slope=0) marked as a red star.
