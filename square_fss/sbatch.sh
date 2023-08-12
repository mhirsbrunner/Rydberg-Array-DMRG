#!/bin/bash
#SBATCH -N 1
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -J FSS_DMRG
#SBATCH --mail-user=markhirsbrunner@lbl.gov
#SBATCH --mail-type=ALL
#SBATCH -t 12:00:00
#SBATCH --account=m888

module load julia
module load python

julia --sysimage ../sys_rydberg.so fss.jl $N_Y $N_WORKERS $N_THREADS $CONFIG_DIR $DATA_DIR --write_dir $WRITE_DIR

mv slurm-${SLURM_JOB_ID}.out $TOP_DIR