# Grid-search analysis

## Purpose and scope

The grid analysis evaluates monotonic radial-acceleration ordering for a finite set of stellar-disc mass-to-light (M/L) profiles in SPARC galaxies. It compares constant, linear, and log-linear radial profiles. A profile candidate specifies the stellar-disc M/L as a function of radius; it changes the disc surface density and hence the Newtonian baryonic field.

The shared calculation engine is in `../miscellaneous/support_code/mono_SPARC_RAR.py`; supporting disk, catalogue, and probability routines are in `../miscellaneous/support_code/`. The optional boundary-resolving probability entry point is `../miscellaneous/support_code/mono_SPARC_RAR_boundary.py`, with its integrator in `../miscellaneous/support_code/probability_recursion.py`. The folders here are not three complete standalone solvers: `constant/` and `linear/` currently contain plotters, while `log-linear/` contains dated recovery and checkpoint workflows.

## Status — 24 October 2026

At the latest recorded 24 October coverage snapshot, grid-family Bayes-factor processing covered 175/175 galaxies for constant profiles, 175/175 for linear profiles, and 130/175 for log-linear profiles. The remaining 45 log-linear cases are assigned to existing jobs; their final valid outputs still need to be published and processed. These are coverage counts over evaluated profiles, not evidence that the finite grids locate a global optimum or validate a scientific conclusion.

## Method

For each candidate disc profile, the stellar field is calculated from the SPARC surface-brightness profile and combined with the signed gas contribution and the bulge contribution at the kinematic radii. The observed acceleration is inferred from the measured rotation velocities and radii. The monotonicity statistic, `P_mono`, is the probability of the required ordering under the stated velocity-error model and nuisance assumptions; it is not a conventional rotation-curve likelihood or a probability assigned to a gravity theory.

The grid driver can write conditional probabilities, probabilities marginalised over gas and bulge M/L, and a separately prior-marginalised grid summary. Gas and bulge M/L are integrated with the configured lognormal quadrature. For the prior-marginalised summary, the disc prior is Gaussian in `log10(M/L)` around the specified SPARC median, and radial slopes use the configured uniform interval (default `[-5, 5]`). The priors are approximated on the finite configured grid; the final average is normalized over the valid prior mass represented there, and the retained mass is recorded. The result is therefore limited to the sampled profile support. Non-positive or non-finite disc profiles are excluded from the physical prior domain. Calculation failures remain marked invalid rather than receiving an artificially low probability.

The solver records run settings and source provenance alongside the profile-grid and marginalised tables. Use those records to identify the exact probability implementation, mask, priors, grids, and tie policy for any result. A finite grid only establishes what was found among the evaluated candidates. Downstream candidate selection and Bayes-factor cuts are described in `../bayes_factors/README.md`; they are separate from the raw grid probability and the supplementary prior-marginalised grid statistic.

## Files

### `constant/`

- `plot_constant_profiles_two_panel_20261020_v3.py` — renders saved constant-profile selections as a two-panel readability figure; it does not recompute forces or probabilities.
- `plot_constant_profiles_two_panel_20261020_v4.py` — later rendering variant of the same saved selections; visual presentation only.

### `linear/`

- `plot_linear_grid_contours.py` — draws sigma contours from an existing linear-profile `profile_grid.csv`. The contours map `P_mono` to one-sided Gaussian lower-tail levels; they are probability visualisations, not Bayes-factor contours.
- `plot_linear_overlap_style_20261020_v3.py` — re-renders the saved set of 77 linear-profile selections for readability; no fit or probability is recomputed.
- `plot_linear_overlap_style_20261020_v4.py` — later visual variant of those saved selections; no new analysis.

### `log-linear/`

These dated utilities target specific recovery cohorts; they are not a general log-linear production interface. Their plans, manifests, and existing output receipts determine the intended cohort and shard.

- `prepare_loglinear_recovery_20261019.py` — audits and prepares metadata for the original 59-galaxy recovery scope; no galaxy solves or queue changes.
- `run_loglinear_recovery_20261019.py` — runs that isolated recovery scope using positional shard identifiers; `--plan` is read-only, and the script does not submit jobs or automatically retry them.
- `prepare_loglinear_checkpoint37_20261019.py` — metadata-only preparation and overlap audit for a 37-galaxy checkpointed replacement set.
- `run_loglinear_checkpoint37_20261019.py` — checkpointed worker for that fixed replacement set; it does not implement automatic retry or resume.
- `prepare_loglinear_bulgeless47_20261019.py` — prepares the manifest and overlap audit for a 47-galaxy bulgeless recovery set; no galaxy solves.
- `run_loglinear_bulgeless47_20261019.py` — isolated runner for that set, with one galaxy per exclusive shard claim; `--plan` is read-only.
- `run_loglinear_cohort10_checkpoint_20261022.py` — checkpointed runner for a fixed cohort of 10. Its source docstring still says 37; use the file’s cohort schema and roster for its scope.
- `run_loglinear_waiting8_checkpoint_20261022.py` — checkpointed runner for the fixed waiting cohort of eight. Its source docstring also contains stale 37-galaxy wording; use its schema and roster.
- `run_loglinear_checkpoint37_20261019.py.sh` — dated Slurm launcher for the checkpoint37 worker. It contains site-specific paths and is not a portable local launcher.

## Reuse and interpretation

Read the selected script, its manifest, and output metadata before reusing a dated recovery tool. Some scripts contain machine-specific paths and are tied to a particular output layout. The `--plan` modes described above are for inspection; they do not validate science or authorize queue actions. A plot or a completed shard does not by itself establish convergence, complete coverage, or a validated scientific conclusion. Saved family outputs are grouped under [`outputs/grid_searches/constant/`](../outputs/grid_searches/constant/), [`outputs/grid_searches/linear/`](../outputs/grid_searches/linear/), and [`outputs/grid_searches/log-linear/`](../outputs/grid_searches/log-linear/); inspect shard metadata and receipts before treating a set as complete.
