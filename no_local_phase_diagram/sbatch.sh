#!/bin/bash
#SBATCH -N 1
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -J 8x4_pbc_pd
#SBATCH --mail-user=hrsbrnn2@illinois.edu
#SBATCH --mail-type=ALL
#SBATCH -t 24:00:00
#SBATCH --account=m888

module load julia
module load python

ml use /global/common/software/nersc/n9/julia/modules

srun -n 1 -u julia -J../sys_rydberg.so phase_diagram.jl $N_WORKERS $N_THREADS $CONFIG_DIR $DATA_DIR --outputlevel 1

mv slurm-${SLURM_JOB_ID}.out $TOP_DIR