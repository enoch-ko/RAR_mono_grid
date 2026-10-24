# Bayesian Assessment of RAR Monotonicity in SPARC Galaxies

## Project overview

This project assesses whether SPARC galaxy rotation-curve data are compatible with a monotonic relation between observed and baryonic acceleration under alternative radial stellar-disc mass-to-light (M/L) profiles. Grid-based and Gaussian-process (GP) searches examine constant, linear and log-linear profiles. The analysis records the monotonicity probability, P_mono, and Bayes factors comparing the monotonic model with an encompassing model that allows unrestricted orderings.

The physical motivation includes the increasing circular-orbit mappings predicted by some modified-inertia (MI) models. Implications for modified-gravity (MG) models depend on the specific theory, geometry and environment. This project does not provide a full comparison of MI, MG and ΛCDM theories.

A separate comparison workflow fits resolved K-band M/L profiles from MIGHTEE-HI/LADUMA data. These comparisons are descriptive and are not used as calibrated priors for SPARC.

## Analysis and interpretation

The grid calculations evaluate the stated profile families while accounting for gas and bulge M/L uncertainty. GP searches explore additional candidate profiles and retain probability and Bayes-factor information. These quantities describe compatibility under the declared data, error model and priors; they are not conventional rotation-curve likelihoods.

Numerical results remain provisional until the relevant jobs are complete and numerical convergence, quality-control and provenance checks are documented. Processing coverage alone is not scientific validation.

## Processing snapshot — 24 October 2026

Bayes-factor processing covers 175/175 galaxies for constant profiles, 175/175 for linear profiles, and 130/175 for log-linear profiles. The remaining 45 log-linear cases are assigned to existing jobs; their final valid outputs still need to be published and processed.

All three BF-aware GP pilot jobs completed on Glamdring with scheduler exit code 0 and complete batch receipts: pilot_01/UGC00634 (job 991244), pilot_02/NGC4068 (job 991246), and pilot_03/NGC7331 (job 991247). Each reports `finite_search_complete`, 7,681 candidate evaluations, and zero invalid candidates. The local copy matches pilots 1 and 2 but is incomplete for pilot 3; [GP_searches/README.md](GP_searches/README.md) lists the eight artifacts to synchronize.

Pilot completion and processing coverage describe workflow state, not scientific validation. The GP outputs retain `scientific_validated=false` and `cut_accuracy=unverified`; numerical convergence and inference checks remain separate.

## Directory guide

- archive/ — retired or superseded files retained for reference.
- grid_searches/constant/, grid_searches/linear/ and grid_searches/log-linear/ — grid-search code and material organized by profile family.
- GP_searches/ — Gaussian-process search code, inputs and run definitions.
- bayes_factors/ — Bayes-factor calculations, processing tools and comparisons.
- outputs/ — results grouped under outputs/archive/, outputs/grid_searches/{constant,linear,log-linear}/, outputs/GP_searches/, outputs/bayes_factors/, outputs/tests/, outputs/validation/ and outputs/miscellaneous/; there is no nested outputs/ directory.
- tests/ — automated tests for analysis and support code.
- validation/ — numerical checks, quality-control tools and validation records.
- miscellaneous/ — documentation, reference data and supporting utilities.
