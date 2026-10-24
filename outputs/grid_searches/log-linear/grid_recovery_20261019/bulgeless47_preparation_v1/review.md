# 47-galaxy log-linear scheduling replacement: PREPARED ONLY

No science/cancellation/submission/merge performed. Old dependencies and active jobs untouched.
47 zero-bulge kappa=1 controls; source column correction does not require their rerun.
Repeat 17 useful calculations held only in RAM; all47 are covered, not just remaining30.
One galaxy per shard gives durable CSV/fully-marginalised/metadata after each galaxy, plus exclusive claim and events.
Settings identical to approved51/400/9, full nuisance/disc marginalisation, signed gas, same boundary/tie/physical-prior policies.
Reuse private instance of frozen worker implementation; no existing source file modified.

28 active now; cancelling these4 would leave24, allowing6 replacements under30 cap.
Recheck queue and old completed products immediately before any approved action; hold if changed.
All22 corrected-grid jobs973851,973864..973884 and GP969420/969421 remain untouched.
No automated submission/refill. Later41 require explicit exact-command approval in capacity-bounded tranches.
1CPU/8GiB each; -q cmb; -c 2 days is a comment, not a hard limit. Initial6 request6CPU/48GiB total.
Old completed mixed17/18-galaxy shards ran days; elapsed13+days on active jobs is not a pergalaxy ETA.
Smaller shards do not accelerate one difficult galaxy, and proxy ordering is not measured runtime.
CSV-only storage estimate is in preparation_audit.json; no full per-candidate JSON output.
Merging is deferred: remind user to run merge_run_shards.py after all47 have verified complete compatible products.
Do not combine with other roots blindly: manifest/source signatures intentionally differ.

## Exact allocations

shard | galaxy | oldjob | repeated RAM-only
--- | --- | --- | ---
01 | NGC6015 | 935526 | False
02 | NGC6503 | 935526 | False
03 | UGC11455 | 935531 | False
04 | NGC5055 | 935526 | True
05 | NGC4559 | 935526 | True
06 | UGCA444 | 935531 | False
07 | UGC04278 | 935528 | False
08 | UGC05721 | 935528 | False
09 | NGC4183 | 935526 | True
10 | UGC09037 | 935531 | True
11 | NGC5907 | 935526 | False
12 | NGC5585 | 935526 | False
13 | UGC06446 | 935529 | True
14 | UGC12732 | 935531 | False
15 | UGC12506 | 935531 | False
16 | NGC5371 | 935526 | False
17 | UGC12632 | 935531 | False
18 | UGC04305 | 935528 | False
19 | NGC4214 | 935526 | True
20 | UGC11557 | 935531 | False
21 | UGC05829 | 935529 | True
22 | UGC06917 | 935529 | False
23 | UGC08550 | 935531 | True
24 | UGC05986 | 935529 | True
25 | UGC05764 | 935529 | True
26 | UGC11820 | 935531 | False
27 | UGC06930 | 935529 | False
28 | UGC04499 | 935528 | False
29 | UGC08837 | 935531 | True
30 | UGC04325 | 935528 | False
31 | UGC05005 | 935528 | False
32 | UGC06399 | 935529 | True
33 | UGC06818 | 935529 | False
34 | UGC05716 | 935528 | False
35 | UGC06628 | 935529 | True
36 | UGCA281 | 935531 | False
37 | UGC04483 | 935528 | False
38 | UGC05414 | 935528 | False
39 | UGC06923 | 935529 | False
40 | UGC10310 | 935531 | False
41 | UGC06667 | 935529 | True
42 | UGC05750 | 935528 | False
43 | NGC4389 | 935526 | True
44 | UGC09992 | 935531 | False
45 | UGCA442 | 935531 | False
46 | UGC05918 | 935529 | True
47 | UGC05999 | 935529 | True

## Proposed initial cancellation and six-job tranche (approval required)

~~~sh
scancel 935526 935528 935529 935531
cd /mnt/users/koe/RAR_mono_grid/outputs/grid_recovery_20261019/bulgeless47_preparation_v1/submission_logs
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_bulgeless47_20261019.py 01
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_bulgeless47_20261019.py 02
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_bulgeless47_20261019.py 03
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_bulgeless47_20261019.py 04
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_bulgeless47_20261019.py 05
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_bulgeless47_20261019.py 06
~~~

## Later41 (HOLD; exact staged approval required)

~~~sh
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_bulgeless47_20261019.py 07
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_bulgeless47_20261019.py 08
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_bulgeless47_20261019.py 09
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_bulgeless47_20261019.py 10
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_bulgeless47_20261019.py 11
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_bulgeless47_20261019.py 12
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_bulgeless47_20261019.py 13
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_bulgeless47_20261019.py 14
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_bulgeless47_20261019.py 15
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_bulgeless47_20261019.py 16
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_bulgeless47_20261019.py 17
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_bulgeless47_20261019.py 18
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_bulgeless47_20261019.py 19
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_bulgeless47_20261019.py 20
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_bulgeless47_20261019.py 21
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_bulgeless47_20261019.py 22
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_bulgeless47_20261019.py 23
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_bulgeless47_20261019.py 24
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_bulgeless47_20261019.py 25
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_bulgeless47_20261019.py 26
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_bulgeless47_20261019.py 27
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_bulgeless47_20261019.py 28
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_bulgeless47_20261019.py 29
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_bulgeless47_20261019.py 30
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_bulgeless47_20261019.py 31
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_bulgeless47_20261019.py 32
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_bulgeless47_20261019.py 33
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_bulgeless47_20261019.py 34
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_bulgeless47_20261019.py 35
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_bulgeless47_20261019.py 36
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_bulgeless47_20261019.py 37
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_bulgeless47_20261019.py 38
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_bulgeless47_20261019.py 39
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_bulgeless47_20261019.py 40
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_bulgeless47_20261019.py 41
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_bulgeless47_20261019.py 42
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_bulgeless47_20261019.py 43
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_bulgeless47_20261019.py 44
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_bulgeless47_20261019.py 45
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_bulgeless47_20261019.py 46
addqueue -q cmb -n 1 -m 8 -c "2 days" /mnt/users/koe/RAR_mono_grid/run_loglinear_bulgeless47_20261019.py 47
~~~
