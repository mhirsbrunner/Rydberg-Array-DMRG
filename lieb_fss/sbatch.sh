#!/bin/bash
#SBATCH -N 1
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -J 8_threads
#SBATCH --mail-user=mark@mhirsbrunner.com
#SBATCH --mail-type=ALL
#SBATCH -t 1:00:00
#SBATCH --account=m4370
#SBATCH --reservation=fss_23x12_test

module load julia
module load python

julia --sysimage ../sys_rydberg.so fss.jl $N_Y $N_WORKERS $N_THREADS $CONFIG_DIR $DATA_DIR --write_dir $WRITE_DIR --outputlevel 1

mv slurm-${SLURM_JOB_ID}.out $TOP_DIR