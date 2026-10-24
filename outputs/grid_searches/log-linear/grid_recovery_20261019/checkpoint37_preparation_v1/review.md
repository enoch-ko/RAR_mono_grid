# 37-galaxy log-linear scheduling replacement: PREPARED ONLY

No science/cancellation/submission/merge performed. Old dependencies and active jobs untouched.
37 zero-bulge kappa=1 controls; source column correction does not require their rerun.
Repeat 12 useful calculations held only in RAM; all37 are covered, not just remaining25.
One galaxy per shard plus checkpoints/GALAXY: conditional.partial.csv, marginalised_gas_bulge.partial.csv and profile_grid.partial.csv stream completed cells; final CSVs/grids_complete.json then fully_marginalised.csv/galaxy_complete.json are atomically published with hashes and exact context. Partial files never imply complete marginalisation. Invalid cells/reasons remain explicit, no clipping.
Settings identical to approved51/400/9, full nuisance/disc marginalisation, signed gas, same boundary/tie/physical-prior policies.
Reuse frozen science with process-local I/O hooks; no existing source file modified. Exact conditional and nuisance rows flush per completed cell, original interleaved evaluation order unchanged. Complete grids published before disc marginalisation; complete galaxy numeric CSVs/metadata before plot/shard completion. No automatic resume.

28 active now; KEEP935526, cancel only3 would leave25, allowing15 replacements under40 cap.
Recheck queue and old completed products immediately before any approved action; hold if changed.
935526 and all22 corrected-grid jobs973851,973864..973884 and GP969420/969421 remain untouched.
No automated submission/refill. Later22 require explicit exact-command approval in capacity-bounded tranches.
1CPU/8GiB each; -q cmb; -c 2 days is a comment, not a hard limit. Initial15 request15CPU/120GiB total.
Old completed mixed17/18-galaxy shards ran days; elapsed13+days on active jobs is not a pergalaxy ETA.
Smaller shards do not accelerate one difficult galaxy, and proxy ordering is not measured runtime.
CSV-only base estimate is in preparation_audit.json; checkpoints add approximately two further copies of numeric row volume plus provenance. Reserve150MiB provisionally, not a hard cap. Hard-link staging aliases do not duplicate physical bytes. No per-cell JSON files. Failure preserves prior flushed cells where storage permits; SIGKILL/in-progress writes or quota cannot be guaranteed.
Merging is deferred: remind user to run merge_run_shards.py after all37 have verified complete compatible products.
Do not combine with other roots blindly: manifest/source signatures intentionally differ.

## Exact allocations

shard | galaxy | oldjob | repeated RAM-only
--- | --- | --- | ---
01 | UGC11455 | 935531 | False
02 | UGCA444 | 935531 | False
03 | UGC04278 | 935528 | False
04 | UGC05721 | 935528 | False
05 | UGC09037 | 935531 | True
06 | UGC06446 | 935529 | True
07 | UGC12732 | 935531 | False
08 | UGC12506 | 935531 | False
09 | UGC12632 | 935531 | False
10 | UGC04305 | 935528 | False
11 | UGC11557 | 935531 | False
12 | UGC05829 | 935529 | True
13 | UGC06917 | 935529 | False
14 | UGC08550 | 935531 | True
15 | UGC05986 | 935529 | True
16 | UGC05764 | 935529 | True
17 | UGC11820 | 935531 | False
18 | UGC06930 | 935529 | False
19 | UGC04499 | 935528 | False
20 | UGC08837 | 935531 | True
21 | UGC04325 | 935528 | False
22 | UGC05005 | 935528 | False
23 | UGC06399 | 935529 | True
24 | UGC06818 | 935529 | False
25 | UGC05716 | 935528 | False
26 | UGC06628 | 935529 | True
27 | UGCA281 | 935531 | False
28 | UGC04483 | 935528 | False
29 | UGC05414 | 935528 | False
30 | UGC06923 | 935529 | False
31 | UGC10310 | 935531 | False
32 | UGC06667 | 935529 | True
33 | UGC05750 | 935528 | False
34 | UGC09992 | 935531 | False
35 | UGCA442 | 935531 | False
36 | UGC05918 | 935529 | True
37 | UGC05999 | 935529 | True

## Proposed initial cancellation and fifteen-job tranche (approval required)

~~~sh
scancel 935528 935529 935531
cd /mnt/users/koe/RAR_mono_grid/outputs/grid_recovery_20261019/checkpoint37_preparation_v1/submission_logs
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_checkpoint37_20261019.py 01
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_checkpoint37_20261019.py 02
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_checkpoint37_20261019.py 03
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_checkpoint37_20261019.py 04
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_checkpoint37_20261019.py 05
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_checkpoint37_20261019.py 06
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_checkpoint37_20261019.py 07
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_checkpoint37_20261019.py 08
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_checkpoint37_20261019.py 09
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_checkpoint37_20261019.py 10
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_checkpoint37_20261019.py 11
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_checkpoint37_20261019.py 12
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_checkpoint37_20261019.py 13
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_checkpoint37_20261019.py 14
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_checkpoint37_20261019.py 15
~~~

## Later22 (HOLD; exact staged approval required)

~~~sh
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_checkpoint37_20261019.py 16
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_checkpoint37_20261019.py 17
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_checkpoint37_20261019.py 18
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_checkpoint37_20261019.py 19
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_checkpoint37_20261019.py 20
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_checkpoint37_20261019.py 21
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_checkpoint37_20261019.py 22
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_checkpoint37_20261019.py 23
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_checkpoint37_20261019.py 24
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_checkpoint37_20261019.py 25
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_checkpoint37_20261019.py 26
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_checkpoint37_20261019.py 27
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_checkpoint37_20261019.py 28
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_checkpoint37_20261019.py 29
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_checkpoint37_20261019.py 30
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_checkpoint37_20261019.py 31
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_checkpoint37_20261019.py 32
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_checkpoint37_20261019.py 33
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_checkpoint37_20261019.py 34
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_checkpoint37_20261019.py 35
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_checkpoint37_20261019.py 36
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_checkpoint37_20261019.py 37
~~~
