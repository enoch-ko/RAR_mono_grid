# RAR monotonicity project guidance

## Scope and scientific goal

This repository supports a publishable paper on monotonicity of the RAR in SPARC galaxies. It quantifies each galaxy's compatibility with a monotonic RAR under explicit baryonic, uncertainty, and stellar-disc M/L assumptions, then assesses implications for the SPARC error model and for Modified Inertia (MI) versus Modified Gravity (MG). The MI argument is restricted to isolated, approximately circular systems whose predictions imply a one-to-one monotonic relation between observed and baryonic acceleration.

Do not turn numerical outputs into scientific claims until convergence checks, quality-control rules, and data provenance are recorded.

## Shared coordination rules

- This file is the durable coordination record for the repository on Glamdring. The local repository copy should be synchronised from it when the user is ready.
- All four project agents are local: planning/tracking, manuscript editing, code/analysis/validation, and literature/theory.
- Record durable methodological decisions, interfaces, validation outcomes, and handoffs in AGENTS.md. Use Notion only for concise status, priorities, decisions, blockers, and daily records.
- The Research Logbook contains only Next priorities and Daily record. Broader goal, dated status, analysis design, and open decisions belong in RAR - Current Status.
- Treat existing uncommitted changes as user-owned. Inspect and preserve them before editing.
- Before each Glamdring submission batch, obtain the user’s explicit approval of the exact job count, commands, and purpose. After that approval, the code, analysis, and validation agent may submit only those approved jobs. Without approval, prepare commands, paths, expected products, and resource estimates for the user to queue manually; any material change to the approved batch requires renewed approval.

## Planning and tracking agent

- Maintain the roadmap, priorities, decisions, blockers, and handoffs.
- Keep the Research Logbook concise, with each daily entry as a date toggle containing bullet points.
- Keep broader project context in RAR - Current Status.
- Update Notion routinely when appropriate, but do not modify manuscript or code without explicit permission.
- Check this file and the status page before assigning work or resolving conflicts.

## Manuscript editing agent

- Edit the Overleaf manuscript only after explicit user approval, stating intended changes and justification briefly.
- Use \gptdel{...} and \gptnew{...} for substantive tracked changes.
- Preserve the acknowledgement section unless explicitly asked otherwise.
- Maintain coherent section, equation, figure, table, and citation cross-references.
- Keep claims provisional until relevant code, convergence, QC, and provenance checks are complete.
- Flag unresolved physics, interpretation, or literature issues rather than silently deciding them.

## Code, analysis, and validation agent

- Work from the local repository and keep it synchronised with this Glamdring copy through the established workflow.
- Focus on physical disc-field calculations, grid statistics, SPARC-informed marginalisation, convergence, unit tests, and QC.
- Validate constant, linear, and log-linear profiles, including UGC07399 pilots, grid/probability-grid/quadrature convergence, tied-g_bar handling, invalid cells, and retained physical prior mass.
- Use Gauss-Hermite quadrature for gas and bulge M/L in log10(M/L), trapezoidal weights for disc/slope grids, the signed SPARC gas convention, and positive finite stellar-disc profiles only.
- Run both radial-prior cases during validation: finite uniform slope prior [-5,5] as a diagnostic sensitivity case, and Gaussian slope prior N(0,0.5) as the intended headline fully marginalised case.
- Do not promote numerical outputs to manuscript claims until validation criteria and provenance are documented.

## Literature and theory agent

- Maintain the citation and conceptual basis for RAR monotonicity, the SPARC error model, and MI versus MG.
- Review MIGHTEE-HI/LADUMA literature and data provenance, distinguishing published results, private communication, and project inference.
- Provide citation-ready recommendations and identify claims requiring further verification.
- Keep the MI argument explicitly conditional on isolated, approximately circular systems and a one-to-one monotonic MI relation between observed and baryonic acceleration.
- Do not edit the manuscript directly unless explicitly assigned and approved; hand recommendations to the manuscript agent through AGENTS.md.

## Statistical design

### Grid-based statistics

For each SPARC galaxy and constant, linear, or log-linear disc-M/L profile, scan a declared disc-M/L grid. At each point, marginalise gas and bulge M/L with the SPARC lognormal uncertainty model. Retain bin-by-bin tables; histograms are grid summaries, not posteriors unless explicit prior weights are applied.

### Fully marginalised statistics

Use a Gaussian prior in log10(M/L), equivalently lognormal in M/L, centred on 0.5 with width 0.1 dex by default. For radial profiles, integrate over slope using either the finite uniform prior [-5,5] or the SPARC-informed Gaussian prior N(0,0.5). Keep distance and inclination fixed to SPARC catalogue values for this phase, and defer population aggregation until per-galaxy results are validated.

Keep distinct: grid probabilities with gas/bulge marginalisation, and fully marginalised probabilities with disc M/L and slope integrated under named priors. These are monotonicity probabilities under the stated observational-error and nuisance-prior model, not conventional rotation-curve likelihoods.

## MIGHTEE-HI/LADUMA workstream

The CSV data were supplied by Andreea Vărăşteanu via private communication. Cite Vărăşteanu et al. (2026) for the joint resolved dynamics sample and Ponomareva et al. (2026) for the resolved stellar-M/L and mass-modelling methodology when used.

Validate columns, units, radii, uncertainties, missing values, covariance, and upper limits. Fit positive constant and linear profiles, preferably in R/R_d; preserve likelihoods or posterior samples; compare with SPARC only after checking bandpass, population, selection, and radial-definition compatibility. Do not silently turn this comparative work into a SPARC prior.

## Validation and production sequence

1. Run test_marginalisation.py on a login node.
2. Queue UGC07399 pilots for constant, linear, and log-linear profiles.
3. Repeat at increasing grid, probability-grid, and M/L quadrature resolutions.
4. Check both prior-labelled probabilities, retained physical prior mass, invalid cells, and hook-removing regions.
5. Queue the full SPARC sample only after validation.
6. Fit MIGHTEE-HI/LADUMA profiles independently.
7. Add figures and quantitative claims only after validation.

## Current repository components

- mono_SPARC_RAR.py: SPARC grid evaluation, nuisance integration, radial-prior cases, and output writing.
- merge_run_shards.py: validates compatible shard outputs and combines them into one CSV pair and provenance file per profile batch.
- utils.py: loading, radial integrals, quadrature, monotonicity recursion, and tie handling.
- test_marginalisation.py: lightweight tests without SPARC files or galaxy analysis.
- mono_SPARC_RAR.py.sh: launcher; verify argument forwarding before production.
- profile_grid.csv, fully_marginalised.csv, and run_metadata.json: analysis products and provenance outputs.
- ML_Andreea.csv: MIGHTEE-HI/LADUMA stellar-M/L data supplied by Andreea Vărăşteanu via private communication.

## Parallel batch outputs

- Use one shared output root per profile/settings batch and a unique `--shard-id` for every concurrent `mono_SPARC_RAR.py` job. Shards write flat, uniquely named CSV and metadata files into that root; galaxy diagnostic plots share the profile's plots directory.
- Do not have concurrent jobs write the common `profile_grid.csv`, `fully_marginalised.csv`, or `run_metadata.json` directly. After all shard jobs finish, run `merge_run_shards.py` once for that profile and pass all expected shard IDs. The merger checks settings, source hashes, schemas, galaxy counts, and duplicate coverage before writing the common outputs.
- Merged outputs are not complete until this merger has run. Preserve the shard artifacts as provenance; use `--replace-existing` only when intentionally rebuilding an existing aggregate.
- Run the merger whenever completed shard jobs are found awaiting reorganisation, and remind the user that this step remains outstanding until it has been run.

## Editing and provenance conventions

- Manuscript revisions use \gptdel{...} and \gptnew{...}; preserve acknowledgements.
- Comment non-trivial code and preserve provenance.
- Average probabilities over tied g_bar orderings.
- Record uncertainties, quadrature, grids, priors, exclusions, fixed-geometry policy, warnings, and source hash.
- Ask for clarification on unresolved physics or interpretation rather than silently deciding.

## 2026-10-03 signed-force cache validation (not production validation)
- mono_SPARC_RAR.py now caches signed W=V|V|=-R F_R: two Casertano basis solves per linear-profile galaxy, and one normalised-shape solve per physical log-linear slope. Constant profiles retain tabulated SPARC disk velocities with the same signed convention. Gas/bulge quadrature, priors, physical M/L mask, invalid-cell accounting, and shard outputs are unchanged.
- The disk response is interpolated as signed force W/R, then multiplied by observed radius. Out-of-grid observed radii retain endpoint clamping. A zero-radius photometric node is omitted only from force interpolation; it is still passed to Casertano, so pre-existing solver failures are not silently repaired.
- python3 -m unittest discover -p 'test_*.py': 24 tests passed. UGC07399 at mean M/L 0.5 and linear slope 0.01: cached versus direct signed disk term differed by at most 2.9e-13 (km/s)^2; changing old velocity interpolation to force interpolation changed it by at most 0.38% relative; the tested conditional log probability was unchanged. These are diagnostics, not a full-grid convergence result.
- Radius audit: 121/175 galaxies have at least one observed radius beyond the photometric-grid bounds; 32/175 photometric grids start at zero. For IC4202, the installed vcdisk returns non-finite velocities at all radii when given its zero-start grid; investigate separately before treating all-galaxy varying-profile outputs as validated.
- The previous driver is preserved as mono_SPARC_RAR_pre_speedup_20261003.py (SHA-256 4daee1b2b5e2e45200ec52da0098fe883625239778ecf6aca3c1b3fe66f0e468). No production output or queue job was changed/submitted in this validation.

## Open decisions

- Whether the finite uniform-slope case appears in the final paper as a sensitivity result alongside the SPARC-informed headline case.
- Sensitivity to slope bounds and Gaussian width.
- Population-histogram grid measure.
- MIGHTEE-HI/LADUMA uncertainty/covariance and radius-normalisation conventions.
- Whether MIGHTEE-HI/LADUMA remains comparative or later informs an empirical prior.
- Glamdring queue syntax and resource limits.

## Historical coordination note

Previous work established the physical disc-field treatment, SPARC gas/bulge marginalisation, positive stellar-disc M/L filtering, exact tie handling, local UGC07399 pilot tests, tracked manuscript edits, and the MIGHTEE-HI/LADUMA comparative workstream. No Glamdring production jobs have been submitted through this project guidance.

## 2026-10-05 first-order inverse-RAR diagnostic

- `inverse_rar_diagnostic.py` is a separate diagnostic: piecewise-linear disc M/L knots spanning the full photometric grid, signed Casertano basis accelerations, physical knot bounds, and a normalized squared-exponential-kernel Mahalanobis minimum from M/L = 0.5 under weak central-value RAR ordering. Fixed gas/bulge M/L enter the QP; the existing full velocity-error probability with gas/bulge quadrature is evaluated only after optimization, without re-optimization. An active baryonic tie is not a strict one-to-one RAR.
- The script requires the completed five-shard, 175-galaxy direct-versus-cached diagnostic for the same driver hash and no target-galaxy cache warnings. It also checks the baseline and any optimized profile against fresh direct force solves. JSON records feasibility status, active constraints, knot values, probabilities where available, input/source hashes, tolerances, and exclusions. It is not an uncertainty-aware global optimum or a GP-sampling analysis.
- `python3 -m unittest -v test_disk_force_cache test_marginalisation test_inverse_rar_diagnostic`: 20 tests passed on the login node. The 4-knot/101-probability-grid UGC07399 and NGC0289 smoke runs returned `infeasible` within that declared profile class/bounds; neither is an astrophysical exclusion. No full inverse run, knot/kernel/bound sensitivity, or probability-grid/quadrature convergence has been completed; no queue jobs were submitted for this diagnostic.
## 2026-10-05 inverse bound/knot pilot and UGC07399 visualisation

- `inverse_rar_diagnostic.py` now accepts `--knots-from-dens` to use every photometric radius as an M/L knot, and rescales only the SLSQP objective numerically while reporting the original Mahalanobis metric. `--max-ml 10` was an explicit exploratory override; the default upper bound was 2 at the time.
- With exact .dens knots and M/L in [0.01,10], UGC07399 gave a feasible weak-order fit with 59 knots and one active tie; NGC0289 remained LP-infeasible with 80 knots. UGC07399 used probability grid 101 and three gas/bulge quadrature nodes, so these are diagnostic results, not strict monotonicity or converged scientific claims.
- An earlier version of `plot_inverse_rar_comparison.py` created a two-panel rotation/RAR PNG and companion JSON provenance from the saved UGC07399 result. Both before/after RAR tracks use fresh Casertano solves; the SPARC tabulated disc, scaled from its M/L=1 convention to 0.5, is only a separate reference. The plotted optimized track is checked against the saved acceleration.

## 2026-10-05 automated inverse plots and infeasible-bound retry

- Current `inverse_rar_diagnostic.py` CLI defaults to exact .dens-radius knots, an exploratory first upper M/L bound of 10, and one retry at 20 only if the first linear feasibility test is infeasible. Optimizer failures are not relabelled as infeasibility or retried. JSON records requested/effective bounds and every attempt. These wide bounds are diagnostic, not a physical prior.
- By default each galaxy now writes a JSON/PNG pair under `outputs/inverse_rar_diagnostic/`. The three-panel PNG has rotation curves and disc M/L versus a shared radius axis, plus RAR tracks. Fitted before/after tracks use direct Casertano solves; if no fit is identified, only SPARC rotmod-based reference curves and the constant 0.5 M/L line appear, with an explicit warning. Plot QC and script hashes are in the result JSON. `--replace-existing` is required to regenerate an existing pair.
- Nine earlier inverse pilot JSON files were moved from the top-level `outputs/` to `outputs/inverse_rar_diagnostic/pilots/`; the five cache-validation shards remain in `outputs/` because the gate reads them there. The original UGC07399 two-panel figure remains under `plots/`, and its provenance points to the relocated source JSON.
- The 35-test lightweight suite passed. End-to-end 101-probability-grid/3-node pilots created both plot types: UGC07399 remained a feasible weak-order fit with one tie at the first bound 10; NGC0289 remained LP-infeasible at 10 and 20. These are exploratory diagnostics, not converged or manuscript-ready results.

## 2026-10-05 full inverse-diagnostic batch prepared (not submitted)

- `run_inverse_rar_batch.py` deterministically allocates the 175 unique SPARC galaxies to ten disjoint shards by squared .dens knot count. `--dry-run` is read-only; `--shard N` runs one shard; `--audit` checks all 175 current JSON/PNG pairs and matching code/settings. The output root is `outputs/inverse_rar_diagnostic/full175_densknots_p400_q9_ml10_retry20/`, with `galaxies/<name>_inverse.json/.png` and `shard_NN.json` reports. Existing mismatched pairs are not overwritten. Each galaxy failure is recorded and does not prevent later galaxies in the shard from running.
- Defaults for this exploratory batch are exact .dens knots, M/L [0.01,10] with one LP-infeasible retry to 20, probability grid 400, and gas/bulge quadrature 9. They are not physical-prior or convergence endorsements. The two earlier 101/3 pilots are outside this root.
- Remote read-only checks found 175 unique catalogue galaxies, 15-109 knots, all 175 rotmod inputs, and ten shards of 16-18 galaxies with balanced knot-squared weights (61018-61165). `python3 -m py_compile run_inverse_rar_batch.py` and the 35-test suite passed. Queue submission still requires the user's explicit approval; no batch was submitted. After results, audit coverage/statuses, source hashes, direct/cache QC, active ties and failed plots before population figures or claims.

## 2026-10-06 boundary-resolving probability recursion (experimental)

- Enoch approved implementation and requested Bob update the numerical methods. Continuous Gaussian-error ordering equations and nuisance priors are unchanged. New probability_recursion.py sums both signed Gaussian branches and cumulatively integrates in log(|V|/sqrt(r)) on a mesh combining logarithmic background, Gaussian neighbourhood and intervening points. Exponential cell interpolants and an analytic small-radius origin cell replace hard uniform-grid masks; no smoothing, probability floor or omitted nuisance weight is introduced.
- Starting resolution doubles until successive logP changes are <=2e-5, up to eight refinements; failure raises ProbabilityConvergenceError and remains numerically invalid. This is a refinement criterion, not a rigorous total-error certificate. zpadding remains an independent finite-tail convergence control; diagnostics include an absolute union tail bound. ngrid now denotes STARTING points per neighbourhood/background component, not the old uniform point count.
- The ordered-gbar interface retains tie averaging and the declared stable-order warning fallback when the cap is exceeded. It counts actual factorial permutations BEFORE allocation; the old utilities incorrectly multiply group sizes. No existing shared solver has been altered.
- python3 -m unittest discover -p 'test_*.py': 49 tests passed in 3.103s, including 14 new analytic/reference, narrow-interval, extreme-tail, invariance, refinement-failure, baryonic-ordering and tie-policy tests. The synthetic narrow-interval reference is P=2.341993260969728e-7; old400/800 grids returned -inf, while the new refined calculation agrees within the tested2e-5 logP tolerance.
- UGC02916's formerly failing gas/bulge node (0.6599182903397267,0.24764151581870858) returned logP=-8499.257884397108, final refinement delta1.668176e-5 at starting resolution400 refined to6400 (540541 common mesh points). The bounded login-node pilot took25.31s including rebuilding the missing fit. All81 nuisance nodes and all175 galaxies have NOT been validated with the new solver. Broader checks should be queued only after exact user approval.
- inverse_rar_boundary.py is a separate entry point: python3 inverse_rar_boundary.py --galaxy UGC02916 --dry-run validates arguments without calculation. It injects the new probability solver only in its own process, retains disk/cache physics checks and existing fit/plot logic, records numerical refinement and source hashes, and writes by default to outputs/inverse_rar_diagnostic/boundary_v1/<galaxy>_inverse.json/.png. Invalid numbers become JSON null while flags/reasons remain. Existing outputs are not overwritten without --replace-existing.
- Shared utils.py, mono_SPARC_RAR.py, inverse_rar_diagnostic.py and the old production default remain untouched because batch935179 was still launching galaxy subprocesses during this edit. Activating a shared default mid-shard would mix numerical algorithms. A deliberate switch and probability-only validation/reprocessing remain pending once the old batch is finished. No new queue jobs or manuscript edits were made by Charlie. Bob was sent the request and exact implementation/validation limits; his completed manuscript edit is not yet verified.
