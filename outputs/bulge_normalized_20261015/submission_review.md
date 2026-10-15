# Bulge-normalized rerun submission review

PREPARED ONLY; NO JOBS SUBMITTED. 32 galaxies per family; 36 jobs total.
1 CPU and 8 GB each. Max36 CPU/288 GB if all concurrent; intended staged20+16.
Comments are expectations, not time limits; GP uses 2 days instead of optimistic5hrs.
Constant/linear:51/400/9, logmean .5+/-1dex, prior.5/.1dex, slope[-5,5], nuisance1/.04dex and.7/.1dex.
Boundary tolerance2e-5,maxrefinements8,padding8; full marginalisation in same grid jobs.
QP uses existing .venv-osqp; OSQP100000iterations/10s,eps1e-9; .01..10 with proven-infeasible-only retry20.
QP numerical failures remain failures; no automatic unapproved escalation or old warm-profile import.
GP5000directions/ell,Rd scales.5/1/2/4,16amplitudes,seed20261013,retained-only8MiBcap.
GP conditional on corrected QP artifacts; only accepted audited fits within GP bounds imported; others explicitly omitted.
Shared output root: /mnt/users/koe/RAR_mono_grid/outputs/bulge_normalized_20261015
Separate GP root: /mnt/users/koe/RAR_mono_grid/outputs/bulge_normalized_20261015/gp_d5000_a16
Old searches and source dependencies preserved; no automatic merge/replacement.
Production is blocked until normalization_validation.json matches current manifest.
No speed optimization is deployed in these commands.

## Stage A: 20 grid/QP jobs
```bash
addqueue -q cmb -n 1 -m 8 -c '5 hrs' /mnt/users/koe/RAR_mono_grid/run_bulge_reruns_20261015.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/bulge_normalized_20261015 --kind constant --shard 1
addqueue -q cmb -n 1 -m 8 -c '5 hrs' /mnt/users/koe/RAR_mono_grid/run_bulge_reruns_20261015.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/bulge_normalized_20261015 --kind constant --shard 2
addqueue -q cmb -n 1 -m 8 -c '5 hrs' /mnt/users/koe/RAR_mono_grid/run_bulge_reruns_20261015.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/bulge_normalized_20261015 --kind constant --shard 3
addqueue -q cmb -n 1 -m 8 -c '5 hrs' /mnt/users/koe/RAR_mono_grid/run_bulge_reruns_20261015.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/bulge_normalized_20261015 --kind constant --shard 4
addqueue -q cmb -n 1 -m 8 -c '5 hrs' /mnt/users/koe/RAR_mono_grid/run_bulge_reruns_20261015.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/bulge_normalized_20261015 --kind linear --shard 1
addqueue -q cmb -n 1 -m 8 -c '5 hrs' /mnt/users/koe/RAR_mono_grid/run_bulge_reruns_20261015.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/bulge_normalized_20261015 --kind linear --shard 2
addqueue -q cmb -n 1 -m 8 -c '5 hrs' /mnt/users/koe/RAR_mono_grid/run_bulge_reruns_20261015.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/bulge_normalized_20261015 --kind linear --shard 3
addqueue -q cmb -n 1 -m 8 -c '5 hrs' /mnt/users/koe/RAR_mono_grid/run_bulge_reruns_20261015.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/bulge_normalized_20261015 --kind linear --shard 4
addqueue -q cmb -n 1 -m 8 -c '5 hrs' /mnt/users/koe/RAR_mono_grid/run_bulge_reruns_20261015.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/bulge_normalized_20261015 --kind linear --shard 5
addqueue -q cmb -n 1 -m 8 -c '5 hrs' /mnt/users/koe/RAR_mono_grid/run_bulge_reruns_20261015.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/bulge_normalized_20261015 --kind linear --shard 6
addqueue -q cmb -n 1 -m 8 -c '5 hrs' /mnt/users/koe/RAR_mono_grid/run_bulge_reruns_20261015.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/bulge_normalized_20261015 --kind linear --shard 7
addqueue -q cmb -n 1 -m 8 -c '5 hrs' /mnt/users/koe/RAR_mono_grid/run_bulge_reruns_20261015.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/bulge_normalized_20261015 --kind linear --shard 8
addqueue -q cmb -n 1 -m 8 -c '5 hrs' /mnt/users/koe/RAR_mono_grid/run_bulge_reruns_20261015.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/bulge_normalized_20261015 --kind qp --shard 1
addqueue -q cmb -n 1 -m 8 -c '5 hrs' /mnt/users/koe/RAR_mono_grid/run_bulge_reruns_20261015.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/bulge_normalized_20261015 --kind qp --shard 2
addqueue -q cmb -n 1 -m 8 -c '5 hrs' /mnt/users/koe/RAR_mono_grid/run_bulge_reruns_20261015.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/bulge_normalized_20261015 --kind qp --shard 3
addqueue -q cmb -n 1 -m 8 -c '5 hrs' /mnt/users/koe/RAR_mono_grid/run_bulge_reruns_20261015.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/bulge_normalized_20261015 --kind qp --shard 4
addqueue -q cmb -n 1 -m 8 -c '5 hrs' /mnt/users/koe/RAR_mono_grid/run_bulge_reruns_20261015.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/bulge_normalized_20261015 --kind qp --shard 5
addqueue -q cmb -n 1 -m 8 -c '5 hrs' /mnt/users/koe/RAR_mono_grid/run_bulge_reruns_20261015.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/bulge_normalized_20261015 --kind qp --shard 6
addqueue -q cmb -n 1 -m 8 -c '5 hrs' /mnt/users/koe/RAR_mono_grid/run_bulge_reruns_20261015.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/bulge_normalized_20261015 --kind qp --shard 7
addqueue -q cmb -n 1 -m 8 -c '5 hrs' /mnt/users/koe/RAR_mono_grid/run_bulge_reruns_20261015.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/bulge_normalized_20261015 --kind qp --shard 8
```

## Stage B: 16 GP jobs
Only after all32 corrected QP JSONs exist. Seek explicit submission approval for this stage.
```bash
addqueue -q cmb -n 1 -m 8 -c '2 days' /mnt/users/koe/RAR_mono_grid/run_bulge_reruns_20261015.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/bulge_normalized_20261015 --kind gp --shard 1
addqueue -q cmb -n 1 -m 8 -c '2 days' /mnt/users/koe/RAR_mono_grid/run_bulge_reruns_20261015.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/bulge_normalized_20261015 --kind gp --shard 2
addqueue -q cmb -n 1 -m 8 -c '2 days' /mnt/users/koe/RAR_mono_grid/run_bulge_reruns_20261015.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/bulge_normalized_20261015 --kind gp --shard 3
addqueue -q cmb -n 1 -m 8 -c '2 days' /mnt/users/koe/RAR_mono_grid/run_bulge_reruns_20261015.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/bulge_normalized_20261015 --kind gp --shard 4
addqueue -q cmb -n 1 -m 8 -c '2 days' /mnt/users/koe/RAR_mono_grid/run_bulge_reruns_20261015.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/bulge_normalized_20261015 --kind gp --shard 5
addqueue -q cmb -n 1 -m 8 -c '2 days' /mnt/users/koe/RAR_mono_grid/run_bulge_reruns_20261015.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/bulge_normalized_20261015 --kind gp --shard 6
addqueue -q cmb -n 1 -m 8 -c '2 days' /mnt/users/koe/RAR_mono_grid/run_bulge_reruns_20261015.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/bulge_normalized_20261015 --kind gp --shard 7
addqueue -q cmb -n 1 -m 8 -c '2 days' /mnt/users/koe/RAR_mono_grid/run_bulge_reruns_20261015.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/bulge_normalized_20261015 --kind gp --shard 8
addqueue -q cmb -n 1 -m 8 -c '2 days' /mnt/users/koe/RAR_mono_grid/run_bulge_reruns_20261015.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/bulge_normalized_20261015 --kind gp --shard 9
addqueue -q cmb -n 1 -m 8 -c '2 days' /mnt/users/koe/RAR_mono_grid/run_bulge_reruns_20261015.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/bulge_normalized_20261015 --kind gp --shard 10
addqueue -q cmb -n 1 -m 8 -c '2 days' /mnt/users/koe/RAR_mono_grid/run_bulge_reruns_20261015.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/bulge_normalized_20261015 --kind gp --shard 11
addqueue -q cmb -n 1 -m 8 -c '2 days' /mnt/users/koe/RAR_mono_grid/run_bulge_reruns_20261015.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/bulge_normalized_20261015 --kind gp --shard 12
addqueue -q cmb -n 1 -m 8 -c '2 days' /mnt/users/koe/RAR_mono_grid/run_bulge_reruns_20261015.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/bulge_normalized_20261015 --kind gp --shard 13
addqueue -q cmb -n 1 -m 8 -c '2 days' /mnt/users/koe/RAR_mono_grid/run_bulge_reruns_20261015.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/bulge_normalized_20261015 --kind gp --shard 14
addqueue -q cmb -n 1 -m 8 -c '2 days' /mnt/users/koe/RAR_mono_grid/run_bulge_reruns_20261015.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/bulge_normalized_20261015 --kind gp --shard 15
addqueue -q cmb -n 1 -m 8 -c '2 days' /mnt/users/koe/RAR_mono_grid/run_bulge_reruns_20261015.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/bulge_normalized_20261015 --kind gp --shard 16
```

## Allocation
constant 1: NGC4013, NGC4138, NGC5985, NGC7814, UGC02487, UGC02953, UGC03546, UGC06614
constant 2: IC4202, NGC4157, NGC4217, NGC7331, UGC02916, UGC06786, UGC06973, UGC09133
constant 3: NGC2683, NGC2841, NGC2955, NGC5033, NGC6674, UGC06787, UGC08699, UGC11914
constant 4: NGC0891, NGC5005, NGC6195, NGC6946, UGC02885, UGC03205, UGC03580, UGC05253
linear 1: NGC4138, UGC02953, UGC06614, UGC06973
linear 2: NGC5033, NGC7814, UGC02487, UGC09133
linear 3: NGC2841, NGC5005, NGC5985, UGC02885
linear 4: NGC4013, NGC6946, UGC03546, UGC08699
linear 5: NGC0891, NGC6195, UGC03205, UGC05253
linear 6: NGC2683, NGC4217, NGC7331, UGC06787
linear 7: NGC2955, NGC6674, UGC02916, UGC06786
linear 8: IC4202, NGC4157, UGC03580, UGC11914
qp 1: NGC4138, UGC02953, UGC06614, UGC06973
qp 2: NGC5033, NGC7814, UGC02487, UGC09133
qp 3: NGC2841, NGC5005, NGC5985, UGC02885
qp 4: NGC4013, NGC6946, UGC03546, UGC08699
qp 5: NGC0891, NGC6195, UGC03205, UGC05253
qp 6: NGC2683, NGC4217, NGC7331, UGC06787
qp 7: NGC2955, NGC6674, UGC02916, UGC06786
qp 8: IC4202, NGC4157, UGC03580, UGC11914
gp 1: UGC02953, UGC06614
gp 2: NGC4138, UGC09133
gp 3: NGC2841, UGC06973
gp 4: NGC6946, UGC02487
gp 5: NGC0891, UGC05253
gp 6: NGC6674, UGC06787
gp 7: NGC4013, UGC02916
gp 8: NGC2683, UGC11914
gp 9: UGC02885, UGC03580
gp 10: NGC4157, UGC06786
gp 11: NGC7331, NGC7814
gp 12: NGC6195, UGC03205
gp 13: NGC4217, UGC08699
gp 14: NGC2955, NGC5985
gp 15: NGC5033, UGC03546
gp 16: IC4202, NGC5005

Preflight: python3 -B /mnt/users/koe/RAR_mono_grid/run_bulge_reruns_20261015.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/bulge_normalized_20261015 --dry-run
After each grid family completes: use merge_run_shards.py on its isolated family directory; do not merge old/new normalization shards.

## Merge completed corrected grid shards
```bash
python3 -B merge_run_shards.py --profile constant --output-root /mnt/users/koe/RAR_mono_grid/outputs/bulge_normalized_20261015/constant --shard-ids 1 2 3 4
python3 -B merge_run_shards.py --profile linear --output-root /mnt/users/koe/RAR_mono_grid/outputs/bulge_normalized_20261015/linear --shard-ids 1 2 3 4 5 6 7 8
```
Never mix old-normalization shards into these directories. GP JSONs need no CSV merger; replacement remains manual and QC-dependent.
