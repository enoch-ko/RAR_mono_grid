#!/bin/bash -l
echo =========================================================
echo Job submitted date = $(date)
date_start=`date +%s`
echo "$SLURM_JOB_NUM_NODES nodes ($SLURM_CPUS_ON_NODE processes per node)"
echo $SLURM_JOB_NUM_NODES hosts used: $SLURM_JOB_NODELIST
export OMPI_MCA_pml=^cm
export OMP_NUM_THREADS=1
cd /mnt/users/koe/RAR_mono_grid
echo Job output begins
echo -----------------
echo
ulimit -l unlimited
/usr/local/shared/slurm/bin/srun -u -n 1 --mpi=pmix --mem-per-cpu=8192 nice -n 10 python3 -u /mnt/users/koe/RAR_mono_grid/mono_SPARC_RAR.py
echo ---------------
echo Job output ends
date_end=`date +%s`
seconds=$((date_end-date_start))
minutes=$((seconds/60))
seconds=$((seconds-60*minutes))
hours=$((minutes/60))
minutes=$((minutes-60*hours))
echo =========================================================
echo Slurm job finished date = `date`
echo Total run time : $hours Hours $minutes Minutes $seconds Seconds
echo =========================================================
