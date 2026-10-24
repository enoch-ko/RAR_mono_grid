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

## Open decisions

- Whether the finite uniform-slope case appears in the final paper as a sensitivity result alongside the SPARC-informed headline case.
- Sensitivity to slope bounds and Gaussian width.
- Population-histogram grid measure.
- MIGHTEE-HI/LADUMA uncertainty/covariance and radius-normalisation conventions.
- Whether MIGHTEE-HI/LADUMA remains comparative or later informs an empirical prior.
- Glamdring queue syntax and resource limits.

## Historical coordination note

Previous work established the physical disc-field treatment, SPARC gas/bulge marginalisation, positive stellar-disc M/L filtering, exact tie handling, local UGC07399 pilot tests, tracked manuscript edits, and the MIGHTEE-HI/LADUMA comparative workstream. No Glamdring production jobs have been submitted through this project guidance.
