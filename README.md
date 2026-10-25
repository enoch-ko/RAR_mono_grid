# Bayesian Assessment of RAR Monotonicity in SPARC Galaxies

## Project overview

This repository contains research code, inputs, tests and selected outputs for studying the radial acceleration relation (RAR) in SPARC galaxy rotation curves. The analysis asks how compatible the observed and baryonic accelerations are with a monotonic ordering when the stellar-disc mass-to-light (M/L) profile is allowed to vary with radius.

The grid searches cover constant, linear and log-linear disc M/L profiles. Gaussian-process (GP) searches explore more flexible radial profiles. The workflows account for gas and bulge M/L uncertainty and report the monotonicity probability, `P_mono`. Bayes-factor analyses compare the posterior probability of the monotonic ordering with its matched prior ordering probability; the denominator is candidate-specific and is averaged over the same nuisance settings as the numerator.

A separate workflow compares resolved K-band M/L profiles from MIGHTEE-HI and LADUMA data. This comparison is descriptive context and is not used as a calibrated prior for SPARC. The physical motivation includes increasing circular-orbit mappings in some modified-inertia models; implications for modified-gravity theories depend on the specific theory, geometry and environment. This repository does not provide a full comparison of modified inertia, modified gravity and ΛCDM.

## Search results and interpretation

Finite grid and GP searches report the best-found candidates and frontiers within the evaluated domain; they do not establish a global optimum. GP proposal directions are search devices, not posterior samples. `P_mono` and the Bayes factors are conditional on the declared data, error model, nuisance assumptions and priors. They are not conventional rotation-curve likelihoods or probabilities that a physical theory is true.

As of 25 October 2026, the three BF-aware GP pilots (UGC00634, NGC4068 and NGC7331) are complete, but remain operational checks: their outputs mark `scientific_validated=false` and `cut_accuracy=unverified`, and they used an earlier driver/manifest. The user-approved full-sample batch has now been submitted as 20 GP jobs (992046–992065); all were running at the latest check. Ten multi-galaxy shards are ordered by increasing radial sample count, and the ten largest-sample galaxies run separately. The batch covers 175 galaxies and brings the active account total to 41 tasks. No production results are complete yet; output reconciliation, numerical accuracy near the provisional BF cuts and convergence checks remain open. See [the GP search guide](GP_searches/README.md) for the live batch record and [the Bayes-factor guide](bayes_factors/README.md) for the statistic and its limits.

## Directory structure

```text
archive/                         retired or superseded material

grid_searches/
  constant/
  linear/
  log-linear/

GP_searches/                     Gaussian-process search code and run definitions
bayes_factors/                   Bayes-factor methods and workflow tools

outputs/
  archive/
  grid_searches/
    constant/
    linear/
    log-linear/
  GP_searches/
  bayes_factors/
  tests/
  validation/
  miscellaneous/

tests/                           automated tests
validation/                      validation tools and records
miscellaneous/                   documentation, reference data and support code
```

The folder-specific README files describe each workflow. Some dated scripts use site-specific paths or dependencies, so inspect their instructions and provenance before running them.
