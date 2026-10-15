# GP full-sample submission review

PREPARED ONLY. 25 jobs x 7 galaxies = 175.
1 CPU and 8 GB per job; up to 25 CPUs / 200 GB concurrently.
Scheduling comment 5 hrs is not a wall-time limit. Runtime is galaxy-dependent.
5000 directions per each of 4 ell/Rd, 16 amplitudes, unchanged physics/probability controls.
Reduced retention; shared galaxies/ outputs, per-shard reports and per-galaxy logs.
First finite-search survey; frontier convergence remains unestablished.

```bash
addqueue -q cmb -n 1 -m 8 -c "5 hrs" /mnt/users/koe/RAR_mono_grid/run_gp_rar_retained_batch.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/gp_rar_full_retained_20261014_d5000_a16 --shard 1
addqueue -q cmb -n 1 -m 8 -c "5 hrs" /mnt/users/koe/RAR_mono_grid/run_gp_rar_retained_batch.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/gp_rar_full_retained_20261014_d5000_a16 --shard 2
addqueue -q cmb -n 1 -m 8 -c "5 hrs" /mnt/users/koe/RAR_mono_grid/run_gp_rar_retained_batch.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/gp_rar_full_retained_20261014_d5000_a16 --shard 3
addqueue -q cmb -n 1 -m 8 -c "5 hrs" /mnt/users/koe/RAR_mono_grid/run_gp_rar_retained_batch.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/gp_rar_full_retained_20261014_d5000_a16 --shard 4
addqueue -q cmb -n 1 -m 8 -c "5 hrs" /mnt/users/koe/RAR_mono_grid/run_gp_rar_retained_batch.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/gp_rar_full_retained_20261014_d5000_a16 --shard 5
addqueue -q cmb -n 1 -m 8 -c "5 hrs" /mnt/users/koe/RAR_mono_grid/run_gp_rar_retained_batch.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/gp_rar_full_retained_20261014_d5000_a16 --shard 6
addqueue -q cmb -n 1 -m 8 -c "5 hrs" /mnt/users/koe/RAR_mono_grid/run_gp_rar_retained_batch.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/gp_rar_full_retained_20261014_d5000_a16 --shard 7
addqueue -q cmb -n 1 -m 8 -c "5 hrs" /mnt/users/koe/RAR_mono_grid/run_gp_rar_retained_batch.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/gp_rar_full_retained_20261014_d5000_a16 --shard 8
addqueue -q cmb -n 1 -m 8 -c "5 hrs" /mnt/users/koe/RAR_mono_grid/run_gp_rar_retained_batch.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/gp_rar_full_retained_20261014_d5000_a16 --shard 9
addqueue -q cmb -n 1 -m 8 -c "5 hrs" /mnt/users/koe/RAR_mono_grid/run_gp_rar_retained_batch.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/gp_rar_full_retained_20261014_d5000_a16 --shard 10
addqueue -q cmb -n 1 -m 8 -c "5 hrs" /mnt/users/koe/RAR_mono_grid/run_gp_rar_retained_batch.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/gp_rar_full_retained_20261014_d5000_a16 --shard 11
addqueue -q cmb -n 1 -m 8 -c "5 hrs" /mnt/users/koe/RAR_mono_grid/run_gp_rar_retained_batch.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/gp_rar_full_retained_20261014_d5000_a16 --shard 12
addqueue -q cmb -n 1 -m 8 -c "5 hrs" /mnt/users/koe/RAR_mono_grid/run_gp_rar_retained_batch.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/gp_rar_full_retained_20261014_d5000_a16 --shard 13
addqueue -q cmb -n 1 -m 8 -c "5 hrs" /mnt/users/koe/RAR_mono_grid/run_gp_rar_retained_batch.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/gp_rar_full_retained_20261014_d5000_a16 --shard 14
addqueue -q cmb -n 1 -m 8 -c "5 hrs" /mnt/users/koe/RAR_mono_grid/run_gp_rar_retained_batch.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/gp_rar_full_retained_20261014_d5000_a16 --shard 15
addqueue -q cmb -n 1 -m 8 -c "5 hrs" /mnt/users/koe/RAR_mono_grid/run_gp_rar_retained_batch.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/gp_rar_full_retained_20261014_d5000_a16 --shard 16
addqueue -q cmb -n 1 -m 8 -c "5 hrs" /mnt/users/koe/RAR_mono_grid/run_gp_rar_retained_batch.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/gp_rar_full_retained_20261014_d5000_a16 --shard 17
addqueue -q cmb -n 1 -m 8 -c "5 hrs" /mnt/users/koe/RAR_mono_grid/run_gp_rar_retained_batch.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/gp_rar_full_retained_20261014_d5000_a16 --shard 18
addqueue -q cmb -n 1 -m 8 -c "5 hrs" /mnt/users/koe/RAR_mono_grid/run_gp_rar_retained_batch.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/gp_rar_full_retained_20261014_d5000_a16 --shard 19
addqueue -q cmb -n 1 -m 8 -c "5 hrs" /mnt/users/koe/RAR_mono_grid/run_gp_rar_retained_batch.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/gp_rar_full_retained_20261014_d5000_a16 --shard 20
addqueue -q cmb -n 1 -m 8 -c "5 hrs" /mnt/users/koe/RAR_mono_grid/run_gp_rar_retained_batch.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/gp_rar_full_retained_20261014_d5000_a16 --shard 21
addqueue -q cmb -n 1 -m 8 -c "5 hrs" /mnt/users/koe/RAR_mono_grid/run_gp_rar_retained_batch.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/gp_rar_full_retained_20261014_d5000_a16 --shard 22
addqueue -q cmb -n 1 -m 8 -c "5 hrs" /mnt/users/koe/RAR_mono_grid/run_gp_rar_retained_batch.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/gp_rar_full_retained_20261014_d5000_a16 --shard 23
addqueue -q cmb -n 1 -m 8 -c "5 hrs" /mnt/users/koe/RAR_mono_grid/run_gp_rar_retained_batch.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/gp_rar_full_retained_20261014_d5000_a16 --shard 24
addqueue -q cmb -n 1 -m 8 -c "5 hrs" /mnt/users/koe/RAR_mono_grid/run_gp_rar_retained_batch.py --output-root /mnt/users/koe/RAR_mono_grid/outputs/gp_rar_full_retained_20261014_d5000_a16 --shard 25
```

## Galaxy allocation

Shard 01: D512-2, F563-V2, F565-V2, F567-2, UGC02953, UGC05750, UGC05918
Shard 02: F571-V1, F583-4, KK98-251, NGC2403, UGC05999, UGC06667, UGCA442
Shard 03: CamB, F583-1, UGC02455, UGC02487, UGC04499, UGC07232, UGC09133
Shard 04: D564-8, F574-1, NGC0055, NGC2841, UGC06446, UGC07577, UGCA444
Shard 05: F568-1, F571-8, NGC5585, NGC6946, UGC05716, UGC07089, UGC07690
Shard 06: D631-7, DDO064, NGC0891, NGC3917, UGC00634, UGC02259, UGC05253
Shard 07: DDO168, F574-2, NGC3769, NGC4217, UGC05414, UGC06787, UGC06923
Shard 08: ESO444-G084, UGC02916, UGC03546, UGC06399, UGC06917, UGC07603, UGC10310
Shard 09: F561-1, NGC4559, UGC07151, UGC07399, UGC08837, UGC11914, UGC12732
Shard 10: F563-V1, NGC0801, NGC2683, NGC3198, NGC3992, NGC4100, UGC06628
Shard 11: NGC1003, NGC4010, PGC51017, UGC01281, UGC06614, UGC09037, UGC11820
Shard 12: DDO154, IC4202, NGC0300, NGC4389, UGC05829, UGC06973, UGC08550
Shard 13: F563-1, NGC3741, NGC4068, NGC4157, NGC5907, UGC07323, UGC07559
Shard 14: F568-V1, NGC0100, NGC5371, UGC03580, UGC05721, UGC05764, UGC12506
Shard 15: ESO563-G021, NGC3109, NGC3953, NGC4013, NGC6195, NGC7793, UGC02023
Shard 16: F568-3, NGC3893, NGC4138, NGC5005, NGC6015, UGC04305, UGC11455
Shard 17: NGC1090, NGC2903, NGC2976, UGC00191, UGC00891, UGC07125, UGC12632
Shard 18: ESO079-G014, NGC2915, NGC3949, NGC4051, NGC5033, UGC00128, UGC00731
Shard 19: NGC0024, NGC0247, NGC4085, NGC4183, NGC6503, NGC6674, UGC05986
Shard 20: NGC1705, NGC2366, NGC3521, NGC5985, UGC04325, UGC06983, UGC09992
Shard 21: ESO116-G012, NGC3726, UGC04278, UGC06786, UGC07261, UGC07524, UGCA281
Shard 22: F579-V1, NGC2955, NGC3877, NGC6789, NGC7331, UGC06930, UGC08699
Shard 23: NGC2998, NGC5055, UGC01230, UGC02885, UGC03205, UGC05005, UGC06818
Shard 24: DDO161, NGC0289, NGC3972, NGC4088, UGC07608, UGC07866, UGC08286
Shard 25: DDO170, IC2574, NGC4214, NGC7814, UGC04483, UGC08490, UGC11557

Read-only preflight: python3 -B run_gp_rar_retained_batch.py --dry-run

After completion: python3 -B run_gp_rar_retained_batch.py --audit
No CSV merger is needed for these selected-profile JSON products.
