# GP performance review — 15 October 2026

Scope: read-only review and a bounded synthetic profiling test. No performance change has been deployed to the ongoing searches or to the prepared rerun sampler.

## Code to inspect

- Sampling loop and fixed-reference deviation: /mnt/users/koe/RAR_mono_grid/gp_rar_search_retained.py, galaxy_run (line266 onward).
- Order-probability cache: the same file, ProbabilityEvaluator (line107 onward).
- Repeated nuisance integration: /mnt/users/koe/RAR_mono_grid/mono_SPARC_RAR.py, evaluate_disk_term (line149 onward).
- Tie grouping: /mnt/users/koe/RAR_mono_grid/utils.py, _tie_groups (line182 onward).
- Numerical integrator: /mnt/users/koe/RAR_mono_grid/probability_recursion.py.
- Batch allocation, not sampling: /mnt/users/koe/RAR_mono_grid/run_gp_rar_retained_batch.py.

## Measurements and limitations

Synthetic fixture:50 retained radii,60 knots,9x9 nuisance nodes,constant xi=.5; one warm evaluation then200 identical candidate evaluations, single-threaded,cProfile,30second cap. Used the existing ProbabilityEvaluator, not a galaxy fit or a calibration run.

-200 evaluations triggered16400 conditional/nuisance-node calls (82 per candidate), about5.54million function calls and approximately3.1seconds cumulative evaluate time under instrumentation.
-About2.45seconds cumulative was in _ordered_logp, including0.86seconds tie grouping and substantial tuple/factorial-key construction. Profiling overhead matters; this is NOT a full-search speed prediction.
-Real completed retained results: NGC0024 took7.32h with67533 distinct cached orders and24530253 cache hits; UGC04278 took4.80h with53753orders; NGC0247 took4.79h with28592orders. UGC00634 took.21h with4orders and24600160hits. These are different galaxies, not a controlled timing comparison.

## Prioritized implementation candidates — NOT implemented

1. A fast no-tie ordering path: stable descending sort, adjacent-equality check and compact order keys. Preserve the existing grouping/permutation treatment whenever ties occur; never replace it with an arbitrary stable ordering.
2. Precompute nuisance terms/weights and observation-only quantities; batch nuisance-node ordering and cached lookups rather than rebuilding them in82 Python calls per candidate.
3. Batch direction/response/whitened-deviation matrix operations, preserving the exact per-galaxy/per-lengthscale RNG stream and amplitude schedule. Floating-point reassociation can change near-tied orderings and must be checked.
4. Profile the cache-miss integrator separately before proposing a safe numerical refactor. Tens of thousands of genuinely distinct orders cannot be removed merely by optimizing cache hits. No reduction in direction count, quadrature or integration accuracy is assumed.
5. Add small periodic progress counters/timestamps without restoring candidate JSON saves. Existing per-lengthscale checkpoints are too coarse for reliable short-term ETAs.

The disc response is already cached within a galaxy. A persistent force-basis cache would help restarts, but it is not demonstrated to remove the dominant sampling/integration cost. The saved historical cache-validation JSON files are QC summaries, not reusable response matrices.

## Checks required before adopting an optimization

Compare cached and uncached probabilities, exact/near ties, all nuisance weights, signed-force positivity floors, warning/invalid-cell handling, random streams and amplitudes, J scores, Pareto frontiers and threshold winners. Keep numerical failures visible. Benchmark representative hard and easy galaxies with bounded pilots; queue any heavy comparison only after approval.

## Prepared rerun boundary

The36-command review uses the unchanged sampler through a new isolated normalization launcher. Original source files remain untouched. All corrected results are separate from the old searches; see submission_review.md. QP uses the existing isolated OSQP1.1.3 environment. Its100000iteration/10second numerical budget is explicit and can still fail on hard cases; no automatic escalation or old-normalization warm-fit import is hidden.
