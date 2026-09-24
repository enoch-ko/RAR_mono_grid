# Grid-based Monotonicity Probability Mapping of SPARC RARs

## Project goal

This project tests the compatibility of the radial acceleration relation (RAR) with
monotonicity across SPARC galaxies under explicit baryonic, observational-uncertainty,
and stellar-disc mass-to-light (M/L) assumptions. The results will be used to assess
implications for the SPARC error model and for Modified Inertia (MI) versus Modified
Gravity (MG). The MI interpretation is conditional on isolated, approximately circular
systems and a one-to-one monotonic relation between observed and baryonic acceleration.
Numerical results remain provisional until convergence, quality-control, and provenance
checks are documented.

## Analysis design

### SPARC grid statistics

- Evaluate constant, linear, and log-linear radial stellar-disc M/L profiles on declared
  parameter grids for each galaxy.
- At each grid point, marginalise over gas and bulge M/L using the SPARC uncertainty
  model; retain the bin-by-bin results as well as per-galaxy summaries.
- Compute the probability of the specified monotonic ordering of observed acceleration
  with baryonic acceleration. These are monotonicity probabilities under the stated
  uncertainty and nuisance-parameter model, not conventional rotation-curve likelihoods.
- Use the signed SPARC gas contribution, retain only positive finite stellar-disc
  profiles, and average probabilities over tied baryonic-acceleration orderings.

### Fully marginalised statistics

Keep grid-based probabilities distinct from fully marginalised probabilities. The latter
integrate over disc M/L using the SPARC-informed lognormal prior (default: M/L centred
on 0.5 with width 0.1 dex in log10(M/L)). For radial profiles, validate both a finite uniform
slope prior over [-5, 5] as a diagnostic and the SPARC-informed Gaussian prior
N(0, 0.5) as the intended headline case. Keep distance and inclination fixed to the
SPARC catalogue values for this phase; defer population aggregation until per-galaxy
results are validated.

### MIGHTEE-HI/LADUMA comparison

Use the supplied resolved K-band M/L data to fit constant, linear, and log-linear radial
profiles, and estimate per-galaxy slope distributions and within-galaxy variation across
radial bins. Compare these results with SPARC; they are descriptive comparisons, not
calibrated SPARC priors. Resolve sample membership and validate units, radii, uncertainties,
and other data-quality issues before drawing conclusions.

## Key components

- `mono_SPARC_RAR.py` — SPARC grid evaluation, nuisance integration, radial-prior cases,
  and output writing.
- `merge_run_shards.py` — validates compatible shard outputs and combines them with
  provenance metadata.
- `utils.py` — data loading, radial integrals, quadrature, monotonicity calculation,
  and tied-ordering handling.
- `test_marginalisation.py` — lightweight tests that do not require SPARC data or a
  galaxy analysis run.
- `RAR_plot.py` — exploratory RAR plotting tool.
- `SPARC_table.txt` — local SPARC catalogue copy.
- `ML_Andreea.csv` — MIGHTEE-HI/LADUMA M/L reference data supplied by Andreea
  Vărăşteanu; see [Vărăşteanu et al. (2026)](https://arxiv.org/abs/2608.03576).

Analysis products include per-galaxy profile-grid and fully marginalised tables,
run metadata, and diagnostic plots. Do not treat a run as complete until its expected
shards have been validated and merged; do not make scientific claims from results that
have not passed the planned convergence and quality-control checks.
