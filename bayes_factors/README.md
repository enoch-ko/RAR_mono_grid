# Bayes-factor analysis

## Question and interpretation

For a fixed stellar-disc profile, the Bayes factor asks whether the observed velocity data favor the monotonic-ordering restriction over an encompassing model that permits unrestricted orderings. It is a comparison of two probability models under stated velocity and nuisance assumptions. It is not the probability that a galaxy, profile family, MI/MG theory, or the RAR is true.

The raw monotonicity probability, `P_mono`, and the Bayes factor answer different questions. A large BF can coexist with a small `P_mono` when the same ordering is rarer under the unrestricted prior. The provisional thresholds `log10(BF) > 1, 1.5, 2` are strict cutoffs under review; they are not validated physical boundaries.

## Status — 24 October 2026

At the latest recorded 24 October processing snapshot, grid Bayes-factor processing covered 175/175 constant-profile galaxies, 175/175 linear-profile galaxies, and 130/175 log-linear-profile galaxies. The remaining 45 log-linear cases are assigned to existing jobs; final valid outputs remain pending publication and processing.

A targeted independent denominator check covered 19 selected saved candidates. All agreed with the higher-precision calculation, with a maximum change below `1.6e-15` in `log10(BF)` and no threshold flips. This verifies those selected denominator calculations only; it is not a global error certificate. Numerator, force, nuisance-quadrature, and tail-error guarantees remain separate validation questions.

## Adopted broad-limit construction

Let `θ` denote a fixed disc profile and `η` the gas/bulge nuisance settings. At one nuisance node, `Q(θ,η)` is the posterior probability of the monotonic ordering under the unrestricted velocity model, while `c(θ,η)` is the prior probability of that same ordering under independent signed velocities with a common uniform bound. In the specified broad-bound limit, the matched quantities are

- `P_mono(θ) = Σ_k w_k Q(θ,η_k)`
- `cbar(θ) = Σ_k w_k c(θ,η_k)`
- `BF_inf(θ) = P_mono(θ) / cbar(θ)`

where the same normalized nuisance weights `w_k` are used in numerator and denominator. This is the ratio of the averaged ordering probabilities, not an average of node-wise `Q/c` ratios. It corresponds to conditioning the joint velocity/nuisance prior on the monotonic event; under that conditioning the nuisance distribution is reweighted in proportion to `c(θ,η)`.

The denominator is the prior probability of the ordering event, not separate evidence for an “unrestricted” theory. Since the ordering and hence `c` can depend on the disc candidate and nuisance node, denominators must be recomputed for the matched candidate. A previously saved probability frontier may therefore omit a BF-optimal candidate.

The numerator and denominator must use the same galaxy, radii, mask, baryonic profile, nuisance weights, physical domain, and tie convention. Ties are handled by averaging over the declared auxiliary strict orderings; a weak-order event or a different tie tolerance changes the model. For the current BF-aware GP implementation, exact adjacent ties are grouped at zero tolerance and excessive permutation counts fail rather than silently selecting an ordering.

## Finite priors, scores, and limits

`BF_inf` is the broad-limit statistic for the stated signed-uniform velocity family. A finite `V_max` calculation is a different proper-prior analysis and needs the matching posterior numerator under the same truncation. Equality of prior ordering denominators for positive and symmetric velocity support does not imply equality of their posterior numerators.

Candidate scores are also distinct from evidence. `J_rad` for grid profiles and `J_ref` for GP proposals use different metrics. In a finite grid or random-direction search, report a minimum found only within the evaluated candidates; do not treat it as a certified global optimum. Failure to find a profile passing a threshold is not proof that no such profile exists. These candidate-level BFs are not evidence for a whole profile family or for MI versus MG.

## Files

- `bayes_factor_denominator_pilot_20261021.py` — prepares a hash-checked ordering receipt for a selected canonical conditional row, then runs a denominator-only broad-limit pilot with grid refinement and a Monte Carlo cross-check. It does not supply a nuisance-integrated numerator.
- `bayes_factor_full_sample_20261022.py` — checkpointed 175-galaxy calculation for saved constant-profile conditional rows. It reuses saved conditional `log P_mono` and computes a unit-scale signed-uniform denominator; its output is explicitly a conditional ratio, not nuisance-integrated evidence.
- `bayes_factor_vmax_scale_check_20261022.py` — recomputes the prior denominator at several common velocity scales while holding the saved numerator fixed; a denominator scale-invariance diagnostic, not a matched finite-prior sensitivity analysis.
- `bounded_uniform_bayes_factor_20261021.py` — fixed-input finite-support ordering integrator and prior Monte Carlo helper with synthetic/demo use; it is a prototype, not the production SPARC data runner.
- `extract_bayes_factor_ranked_20261022.py` — verifies and ranks previously saved galaxy ratios, then writes tabular and example-selection summaries; it performs no new force or probability evaluation.
- `plot_bayes_factor_prior_support_comparison_20261022.py` — rendering-only restyle of saved denominator-scale results; it does not recalculate scientific values.
- `plot_bayes_factor_vmax_scale_styled_20261022.py` — alternative rendering of the same saved scale diagnostic; no new evaluation.
- `prepare_bf_grid_extension_20261024_v1.py` — prepares frozen inputs, a manifest, and a proposed batch for five deferred log-linear galaxies. This dated planner has module-level filesystem and queue-query side effects; do not import or run it casually.
- `screen_bf_accuracy_20261024_v1.py` — screens existing grid BF cells near the provisional cutoffs and selects cases for later precision checks. It performs no force/probability evaluation and is not an accuracy certificate; it also has module-level output-writing side effects.

The current BF-aware GP calculation is implemented by [`JointEvaluator` in `core.py`](../GP_searches/bf_gp_search_20261024_v1/core.py) and orchestrated by [`driver.py`](../GP_searches/bf_gp_search_20261024_v1/driver.py). In this local checkout, the adapter imported by `core.py` is absent, so the package is not runnable as staged; see the GP guide for details. Saved dated outputs are under [`outputs/bayes_factors/`](../outputs/bayes_factors/), and recorded BF-aware GP pilot artifacts are under [`outputs/GP_searches/bf_gp_search_20261024_v1/pilot/`](../outputs/GP_searches/bf_gp_search_20261024_v1/pilot/). The files in this directory are dated pilots, conditional calculations, plots, extraction, and grid-extension/accuracy-screen helpers; they are not one end-to-end current runner for every profile family. Check each script’s paths and receipts before reuse.
