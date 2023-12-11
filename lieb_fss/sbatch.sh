#!/bin/bash
#SBATCH -N 1
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -J 23x12_x_1500
#SBATCH --mail-user=mark@mhirsbrunner.com
#SBATCH --mail-type=ALL
#SBATCH -t 24:00:00
#SBATCH --account=m4370
##SBATCH --reservation=fss_23x12

module load julia
module load python

srun -n 1 -u julia --sysimage ../sys_rydberg.so fss.jl $N_Y $N_WORKERS $N_THREADS $CONFIG_DIR $DATA_DIR --write_dir $WRITE_DIR --outputlevel 1

mv slurm-${SLURM_JOB_ID}.out $TOP_DIR