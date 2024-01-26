#!/bin/bash
#SBATCH -N 1
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -J lieb_pd_obc
#SBATCH --mail-user=hrsbrnn2@illinois.edu
#SBATCH --mail-type=ALL
#SBATCH -t 24:00:00
#SBATCH --account=m888

module load julia
module load python

# srun -n 1 -u julia -J../sys_rydberg.so phase_diagram.jl $N_WORKERS $N_THREADS $CONFIG_DIR $DATA_DIR --write_dir $WRITE_DIR
srun -n 1 -u julia -J../sys_rydberg.so phase_diagram.jl $N_WORKERS $N_THREADS $CONFIG_DIR $DATA_DIR

# python ../data_analysis/compile_results.py $TOP_DIR

# rm -rf $TOP_DIR/phase_diagram_data

mv slurm-${SLURM_JOB_ID}.out $TOP_DIR