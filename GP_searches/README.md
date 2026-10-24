# Gaussian-process profile searches

## Purpose and scope

The GP searches explore more flexible, radially varying stellar-disc M/L profiles than the constant, linear, and log-linear grid families. The profiles are represented at photometric-radius knots. A Gaussian-process covariance supplies correlated search directions over those knots, and a finite set of amplitudes turns those directions into candidate profiles.

The proposal GP is a search device. Its candidates are neither posterior samples nor draws from an inferred population of stellar profiles. A finite run reports a best-found candidate or frontier within its sampled directions; it cannot certify a global optimum.

This folder contains two related generations of code: earlier probability-only searches and the BF-aware implementation in `bf_gp_search_20261024_v1/`. They should not be treated as interchangeable. The shared grid probability and force routines live under `../miscellaneous/support_code/`; the BF-aware package also depends on corrected profile/force support code and frozen inputs listed in its manifest.

## Status — 24 October 2026

All three BF-aware pilots completed on Glamdring. Scheduler accounting reports exit code 0, and each batch receipt lists its assigned galaxy as completed with no failed or unprocessed entries.

| Pilot | Job | Galaxy | Candidates | Invalid | Force solves |
|---|---:|---|---:|---:|---:|
| `pilot_01` | 991244 | UGC00634 | 7,681 | 0 | 61 |
| `pilot_02` | 991246 | NGC4068 | 7,681 | 0 | 60 |
| `pilot_03` | 991247 | NGC7331 | 7,681 | 0 | 75 |

Each result reports `finite_search_complete`, but the outputs retain `scientific_validated=false` and `cut_accuracy=unverified`; scheduler success does not validate the numerical thresholds or scientific interpretation.

The local mirror matches the remote pilot_01 and pilot_02 artifacts by filename, size, and modification time. Pilot_03 is not fully synchronized. From `/mnt/users/koe/RAR_mono_grid/outputs/GP_searches/bf_gp_search_20261024_v1/pilot/pilot_03/`, copy these six files for a complete local review set:

- `batch_receipt.json`
- `NGC7331/completion.json`, `NGC7331/result.json`
- `NGC7331/progress_ell_1.json`, `NGC7331/progress_ell_2.json`, `NGC7331/progress_ell_3.json`

Refresh these two local partial journals from their completed remote versions:

- `NGC7331/all_candidates.jsonl.gz` (local 166,465 bytes; remote 631,521 bytes)
- `NGC7331/paired_order_cache.jsonl.gz` (local 14,174 bytes; remote 82,480 bytes)

The local `progress_ell_0.json`, `context.json`, and `response_basis.npz` match; no replacement is needed. These eight files have not been copied by this update.

## Method

A candidate is a positive, finite disc M/L vector on the galaxy’s photometric knots. The search constructs correlated directions using squared-exponential covariance matrices at declared physical length scales, then tests amplitudes over the feasible interval. The new BF-aware implementation uses a common fixed reference covariance metric, `J_ref`, across proposal length scales. The grid analysis uses the distinct physical radial-departure metric `J_rad`. Neither score alone measures astrophysical plausibility.

For each valid candidate, the stellar response is combined with the signed gas field and the bulge contribution over the declared nuisance quadrature. The earlier GP programs calculate the monotonicity probability and retain a finite frontier in probability and `J_ref`. The BF-aware package additionally evaluates the posterior ordering mass and the matching prior ordering mass at the same nuisance nodes, then scores their ratio of averages. It retains the exact score frontier, baseline, maximum-probability and maximum-BF candidates, and minima passing the configured strict BF cutoffs.

The cutoffs `log10(BF) > 1, 1.5, 2` are provisional. The finite random-direction search does not make these cutoffs numerically certified, and “minimum found” is not a guaranteed global minimum. Check each output’s validity, failure, force-check, storage, source-hash, and accuracy fields before interpreting it.

## Files

### BF-aware implementation: `bf_gp_search_20261024_v1/`

- `__init__.py` — marks the directory as an importable Python package.
- `core.py` — fixed-profile joint Bayes-factor evaluation, tie/order journals, storage safeguards, and BF-aware candidate retention; importing it does not start a search.
- `driver.py` — isolated driver that validates pinned inputs, builds the disk-response basis, evaluates the finite candidate search, and writes per-galaxy results and provenance. It hard-codes `/mnt/users/koe/RAR_mono_grid`, so this local checkout cannot run it unchanged.
- `worker.sh` — cluster environment and Python-path wrapper for the driver. It contains an absolute site path and is not a portable local launcher.
- `manifest.json` — frozen galaxy roster, numerical settings, source/input hashes, stage definitions, and resource caps. Its saved status text conflicts with the separate approval receipt; treat it as a configuration record, not live execution status.
- `submission_plan.json` — proposed job commands, rosters, and resources. Its saved status is pre-submission and is not a live job record.
- `approved_pilot.json` — records approval of the exact pilot batch and manifest hash; approval does not prove that a job ran or completed.
- `README.md` — detailed dated settings, storage limits, validation gates, and proposal history. Its statements that no approval receipt exists are stale in this local copy; use the approval receipt for authorization and output/job receipts for execution status.

### Earlier probability-only GP workflow

- `gp_rar_search.py` — exploratory correlated GP-direction search scored by monotonicity probability and a common reference metric.
- `gp_rar_search.py.sh` — dated site-specific launcher for one earlier GP search.
- `gp_rar_search_retained.py` — retained-profile variant of the probability-only search.
- `gp_rar_search_retained_corrected.py` — corrected retained probability-only implementation and a dependency of the BF-aware package.
- `gp_candidate_retention.py` — exact finite-search retention by probability and `J_ref`; helper only, with no astrophysical calculation.
- `gp_frontier_reader.py` — bounded-memory reader for saved legacy frontier records; it inspects metadata and retained records but does not recompute them.
- `run_gp_rar_batch.py` — prepares or executes the earlier frozen 25-shard, seven-galaxy probability-only workflow; it does not submit jobs.
- `run_gp_rar_batch.py.sh` — dated site-specific launcher for that earlier batch workflow.
- `run_gp_rar_retained_batch.py` — batch counterpart for the retained probability-only search; it also does not submit jobs.
- `run_gp_rar_retained_batch.py.sh` — dated site-specific launcher for the retained batch.
- `plot_gp_rar_search.py` — plots saved probability-only frontiers and profile/RAR panels; no search or probability recalculation.
- `plot_gp_rar_search_corrected.py` — corrected-format plotting counterpart; it also reads saved results only.
- `summarize_gp_thresholds.py` — audits old pilot outputs and reports minimum-found `J_ref` at specified raw-probability levels; these are not the current BF cutoffs.

## Local-copy cautions

The local package is incomplete as staged: `core.py` imports `bf_joint_adapter_20261023_v1.adapter`, which is pinned in `manifest.json` but absent from this checkout, and the nested README references `tests/test_bf_gp_20261024_v1.py`, which is also absent. Restore those dependencies before treating the package as runnable. The BF-aware `driver.py`, `worker.sh`, and several older scripts contain absolute site paths; inspect and adapt them before any local execution. The nested pilot documentation and submission plan are dated records, not a source of live queue status. `approved_pilot.json` records authorization only. Saved pilot artifacts are under [`outputs/GP_searches/bf_gp_search_20261024_v1/pilot/`](../outputs/GP_searches/bf_gp_search_20261024_v1/pilot/); use matching batch and per-galaxy receipts to establish what those artifacts contain. This guide describes the local checkout and does not imply any current remote execution state.
