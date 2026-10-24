# Log-linear recovery: ready for queue review

Prepared only; no jobs submitted, cancelled or calculated. Enoch approved preparation in Alice's 19 October exchange (turn 01a154b4-90c1-7e33-9eaa-21f68113ab00).

22 jobs cover 59 unique galaxies: all32 positive-bulge galaxies and27 bulgeless galaxies from failed original shards03/06. The eight bulged galaxies in those failed shards are already among32. Allocation:15 shards of3 and7 shards of2. Heuristic weights use positive density rows times retained kinematic rows squared; they are not measured runtimes.

One CPU and8GB per job;22new +7existing=29 concurrent jobs, below30. Comments say "2 days" based on observed log-linear expense; this is not a hard limit or a completion guarantee. Original completed17/18-galaxy shards took4–12.4days and four remain active beyond13days. New2/3-galaxy shards should reduce serialization, but individual galaxies can dominate.

Settings: log_linear51/400/9; disc median.5/.1dex and log-grid half-width1dex; slope[-5,5]; gas1/.04dex, bulge.7/.1dex; boundary tolerance2e-5/refinements8,padding8,tie_atol0,tie cap720. Both gas/bulge and disc marginalisation are enabled. Correct named gas column3, separate bulge template5, disc from vcdisk. The existing normalization adapter applies Ldisk/Ltotal once in each new worker; no frozen science source edits.

Science output root (not yet created): /mnt/users/koe/RAR_mono_grid/outputs/grid_recovery_20261019/log_linear_v1
Each shard writes profile_grid.shard_NN.csv, fully_marginalised.shard_NN.csv, run_metadata.shard_NN.json, plots, exclusive claim and a short progress journal. Queue output goes under preparation_v1/submission_logs. No automatic retries or merging.

## Verification

53 synthetic/functional tests passed in0.576s:7 new tests plus46 existing policy/retention/dispatcher tests. New tests exercise the real core argument parser without calculation, mocked once-only normalization (bulge and bulgeless), allocation and changed-source rejection, scheduler guard, exclusive output and failure journals. No real-galaxy force, probability or grid evaluation was executed.

All22 plans validated; real CLI --plan01 passed and selected NGC6789,UGC00634,UGC02953.21 source pins and362 input pins verified. Output root creation=false; submissions=0; cancellations=0.

Worker SHA256:7b6cb9f043ee983eee9aa5cb2f250b7e6e8033b47496d3a917a6c2a48dece048
Manifest SHA256:f23528f761118f6dd44ccb6063b29216cfff07317f060f587541b19d8d321651
Test source SHA256:bc657ec506422e84c63e13312a932ad7e7d08fe4df8c0c5e93074753cbcb223c

Machine-readable evidence: manifest.json,launch_pin.json,resolved_plans.json,preparation_receipt.json,overlap_progress_audit.json in this directory.

## Four old jobs: retain; no cancellation criterion established

Two resource observations at15:10:39 and15:20:00UTC show each worker gained561–562CPU seconds. All four Python workers were RN/running at99.9%CPU; RSS remained below1GB with8GB allocations. This rules out idle waiting over the interval, but does not prove mathematical progress/convergence. Logs contain no recorded failure or quota message.

|Job|Shard/current galaxy|New-batch overlap/all assigned|Useful bulgeless/all assigned|Useful current+remaining|
|---|---|---:|---:|---:|
|935526|05,NGC5371 (11/18)|8|10|5|
|935528|07,UGC02953 (3/17)|7|10|10|
|935529|08,UGC06787 (12/17)|4|13|4|
|935531|10,UGC09133 (5/17)|3|14|11|

Their27 preceding galaxies have27 matching plots, but the original code writes CSVs only after the whole shard finishes.17 of those preceding galaxies are useful bulgeless objects, whose numerical rows currently exist only in worker memory. Altogether47 bulgeless galaxies would need replacement if all four jobs were cancelled; preserving coverage would expand the new59-galaxy scope to106 and needs separate approval. No complete job is redundant. Old bulged outputs will be superseded; all existing logs/plots/dependencies are preserved.

Log last-update times before the audit:93552619Oct04:28:37UTC;93552816Oct11:29:29UTC;93552917Oct19:52:31UTC;93553111Oct07:20:38UTC. These sparse per-galaxy markers are not per-cell progress. The new wrapper adds named start/finish journals for future diagnosis.

Current science coverage remains69 saved bulgeless log-linear galaxies +47 assigned active bulgeless +27 failed bulgeless +32 new normalized bulges =175. Constant/linear replacement inventories are complete; GP/reference proposals remain on hold. No new manuscript claims follow from these packaging checks.

## Exact22 commands: explicit submission approval still required

```bash
cd /mnt/users/koe/RAR_mono_grid/outputs/grid_recovery_20261019/preparation_v1/submission_logs
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_recovery_20261019.py 01
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_recovery_20261019.py 02
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_recovery_20261019.py 03
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_recovery_20261019.py 04
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_recovery_20261019.py 05
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_recovery_20261019.py 06
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_recovery_20261019.py 07
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_recovery_20261019.py 08
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_recovery_20261019.py 09
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_recovery_20261019.py 10
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_recovery_20261019.py 11
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_recovery_20261019.py 12
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_recovery_20261019.py 13
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_recovery_20261019.py 14
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_recovery_20261019.py 15
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_recovery_20261019.py 16
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_recovery_20261019.py 17
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_recovery_20261019.py 18
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_recovery_20261019.py 19
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_recovery_20261019.py 20
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_recovery_20261019.py 21
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_recovery_20261019.py 22
```
