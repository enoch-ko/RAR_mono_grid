# BF-aware GP proposal — awaiting exact approval

Manifest SHA256: 577cb8ce344869e3b18163d818fcebc76d6dbce06a32b1cd7a6c12c9a9af1981

Implementation and preparation only. No real GP force/probability evaluation or submission has occurred. Existing scientific kernels and jobs are preserved. Sixteen synthetic tests pass, including a fully mocked end-to-end driver. All 175 input-only checks pass.

## Immediate approval request: three pilot jobs

One job each for UGC00634 (3 radii, 59 knots, at most 65 forces), NGC4068 (6 radii, 58 knots, at most 64 forces), NGC7331 (34 radii, 73 knots, at most 79 forces; kappa 0.953724000622429 applied once).

Each job: queue cmb, one CPU, 8 GiB, two-hour process deadline, per-force cap 300 seconds, no retry or overwrite. Combined upper budget: 6 CPU hours and 208 force solves. Runtime is unmeasured until this pilot. The addqueue -c field is a comment, not a walltime request.

Seed 20261024, replicate 0; ell/Rdisk = 0.5, 1, 2, 4; 128 nested directions per ell; 16 amplitude nodes, baseline once, giving 7,681 candidate evaluations per galaxy and 23,043 overall. At most 629,842 numerator callback/denominator requests per galaxy (82 per candidate, including conditional); actual Gaussian integrations are reduced by warning-free order caching. Paired order-cache cap 50,000 per galaxy.

Fixed reference J uses length_fraction=0.3 for every proposal ell. Positive M/L proposal box [0.01,10]; signed finite gbar accepted; no nuisance positivity floor. Full matched GH9x9 gas/bulge quadrature: gas 1, 0.04 dex; bulge 0.7, 0.1 dex. Probability grid 400, zpadding 8, adaptive log tolerance 2e-5, 8 refinements, exact adjacent tie grouping at atol 0, auxiliary-uniform averaging with cap 720; a larger tie multiplicity fails rather than selecting a stable order.

Exact commands below are proposals only, one invocation each, from the corresponding isolated submission cwd in submission_plan.json:

    addqueue -q cmb -n 1 -m 8 -c BF_GP_pilot_01 /mnt/users/koe/RAR_mono_grid/GP_searches/bf_gp_search_20261024_v1/worker.sh pilot pilot_01
    addqueue -q cmb -n 1 -m 8 -c BF_GP_pilot_02 /mnt/users/koe/RAR_mono_grid/GP_searches/bf_gp_search_20261024_v1/worker.sh pilot pilot_02
    addqueue -q cmb -n 1 -m 8 -c BF_GP_pilot_03 /mnt/users/koe/RAR_mono_grid/GP_searches/bf_gp_search_20261024_v1/worker.sh pilot pilot_03

## Storage and scientific limits

All candidate scalar rows are written before any profile-vector discard, including validity/failure, descriptor/hash, J_ref, ln P, ln cbar, ln BF, log10 BF and matched node-order IDs when valid. Keep the exact BF–J frontier, all three strict cut minimum-J records (log10 BF > 1, 1.5, 2), baseline, maximum P, maximum BF and full retained profiles. Discarded vectors can be replayed from persisted proposal matrices, nested RNG descriptors and pinned environment; verify SHA256 before reuse. This is not storage of every full profile vector. Paired order journals preserve Q and c with tie groups. Context includes full quadrature, normalization, mask, observed arrays, response matrices, RNG states, hashes and cache gate evidence.

Pilot disk reservation 112 MiB per galaxy, 336 MiB total, with explicit per-file caps (scalar gzip 16 MiB, raw 64 MiB; paired orders gzip 8 MiB, raw 16 MiB; basis 8 MiB; retained/progress/context/failure JSON each 8 MiB). Caps stop the search and preserve an incomplete receipt; no frontier truncation or automatic retry. Max retained frontier 2,048. No cache spilling or silent data loss. Shared filesystem has 56 TiB free, but the personal quota command is unavailable; this does not establish personal quota. Recheck storage before submission; production remains gated.

No real pilot accuracy validation has occurred. Point cuts remain accuracy_unverified; a finite random-direction search does not certify a global optimum or theory evidence. Current tests check implementation contracts, not scientific accuracy. Real pilot review must assess force/order QC, failures, per-candidate cost, memory/storage and numerical adequacy before production. Baseline and unique baseline/maxP/maxBF/cut winner profiles receive direct-force and all-node order comparison (at most six direct checks); other retained vectors are marked not_direct_checked.

## Gated production proposal — not ready for approval until pilot review

Full 175-galaxy roster is frozen in manifest.json: 17 batches (five of 11 galaxies and twelve of 10), 5,000 directions per ell, 300,001 candidates per galaxy, 52,500,175 overall. Same seed/settings, with pilot directions a per-ell prefix. One CPU and 8 GiB per job, six-hour per-galaxy deadline, five-day batch ceiling (largest serial cap 66 hours), no retry. Total upper CPU budget 1,050 hours, subject to pilot-driven revision before approval.

Reservation 176 MiB per galaxy, 30,800 MiB (~30.1 GiB) overall; scalar gzip 64 MiB/raw 512 MiB; order gzip 16 MiB/raw 64 MiB; 200,000 paired order cap; other caps as pilot. A measured quota/storage assessment and reviewed_pilot.json tied to this manifest are mandatory before production. Exact commands/rosters are in submission_plan.json. At the saved queue snapshot 21 existing allocations + 17 production jobs = 38; refresh queue and obtain exact production approval after pilot review. Pilot and production are not proposed to overlap.

## Evidence

Current post-migration verification: /mnt/users/koe/RAR_mono_grid/outputs/validation/bf_gp_search_20261024_v1/post_migration/receipt.json. The corresponding GP_synthetic.log, GP_input_preflight.log and worker_shell_syntax.log are in the same directory. Historical receipts/logs are preserved in the adjacent pre_migration directory. Exact replacement commands and rosters are in submission_plan.json. No approved_pilot.json exists.

## Reorganisation on 24 October 2026

This unsubmitted proposal now uses GP_searches/ for source, tests/test_bf_gp_20261024_v1.py for tests, outputs/GP_searches/ for real results and outputs/tests/ for synthetic fixtures. Original proposal bytes are preserved in archive/20261024_reorganisation_v1/gp_proposal_before_path_update. Historical test/preflight receipts are in outputs/validation/bf_gp_search_20261024_v1/pre_migration; new checks are recorded separately in post_migration. The previous manifest and commands are superseded and remain unapproved. Scientific settings and core scoring code are unchanged.
