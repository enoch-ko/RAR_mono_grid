# Bayesian Assessment of RAR Monotonicity in SPARC Galaxies

## Project overview

This repository contains research code, inputs, tests and selected outputs for studying the radial acceleration relation (RAR) in SPARC galaxy rotation curves. The analysis asks how compatible the observed and baryonic accelerations are with a monotonic ordering when the stellar-disc mass-to-light (M/L) profile is allowed to vary with radius.

The grid searches cover constant, linear and log-linear disc M/L profiles. Gaussian-process (GP) searches explore more flexible radial profiles. The workflows account for gas and bulge M/L uncertainty and report the monotonicity probability, `P_mono`. Bayes-factor analyses compare the posterior probability of the monotonic ordering with its matched prior ordering probability; the denominator is candidate-specific and is averaged over the same nuisance settings as the numerator.

A separate workflow compares resolved K-band M/L profiles from MIGHTEE-HI and LADUMA data. This comparison is descriptive context and is not used as a calibrated prior for SPARC. The physical motivation includes increasing circular-orbit mappings in some modified-inertia models; implications for modified-gravity theories depend on the specific theory, geometry and environment. This repository does not provide a full comparison of modified inertia, modified gravity and ΛCDM.

## Search results and interpretation

Finite grid and GP searches report the best-found candidates and frontiers within the evaluated domain; they do not establish a global optimum. GP proposal directions are search devices, not posterior samples. `P_mono` and the Bayes factors are conditional on the declared data, error model, nuisance assumptions and priors. They are not conventional rotation-curve likelihoods or probabilities that a physical theory is true.

At 22:58:20 UTC on 27 October, Charlie's read-only snapshot showed 28 jobs running overall. Of GP production jobs 992046–992065, 17 were running and 992059, 992060 and 992065 had completed successfully; none of the 20 was pending, failed or missing. This is a dated scheduler snapshot, not live state or scientific validation. The active manifest and receipts remain unchanged and match the submitted batches. The user reviewed a separate review-only candidate manifest and 20 receipt drafts; they are not active and are not needed to authorize the existing jobs. Optional promotion is not underway and would require separate explicit direction with matching production-valid receipts.

The three BF-aware GP pilots (UGC00634, NGC4068 and NGC7331) completed as operational checks, but their outputs mark `scientific_validated=false` and `cut_accuracy=unverified`, and they used an earlier driver/manifest. On 25 October, matched diagnostic grid atlases were rendered for 175 constant, 175 linear and 130 log-linear profiles (59, 59 and 44 pages), with per-panel colour scales, log10(P_mono), log10(cbar), log10(BF), and provisional guides at 1, 1.5 and 2. Numerical and scientific validation remain open; direct tick verification of the linear atlas is incomplete. The three atlas PDFs are now in the local repository. Figure 2 still needs verification and refresh from those files, and Figure 1 remains to be reviewed. At 23:45 UTC on 27 October, 11 log-linear jobs were still running (nine recovery shards and two cohort-10 checkpoints); Bayes-factor processing and plots for those runs were deferred until completion. A 28 October roster check could not proceed because the existing VS Code connection could not reconnect while the Mac was locked; no server or scheduler command was run. The expected-shard completeness and contribution to the 45-profile deferred gap remain unconfirmed. Charlie completed and tested the runtime-cap and checkpoint/resumability implementation before GP production; empirical runtime scaling and GP output/search validation remain open. See [the GP search guide](GP_searches/README.md) and [the Bayes-factor guide](bayes_factors/README.md) for details.

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
