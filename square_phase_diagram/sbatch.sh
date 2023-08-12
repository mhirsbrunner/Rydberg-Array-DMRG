#!/bin/bash
#SBATCH -N 1
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -J square_pd_DMRG
#SBATCH --mail-user=markhirsbrunner@lbl.gov
#SBATCH --mail-type=ALL
#SBATCH -t 12:00:00
#SBATCH --account=m888

module load julia
module load python

julia -Jsys_rydberg.so square_phase_diagram.jl $N_WORKERS $N_THREADS $CONFIG_DIR $DATA_DIR --write_dir $WRITE_DIR

python ../data_analysis/compile_results.py $TOP_DIR

rm -rf $TOP_DIR/phase_diagram_data

mv slurm-${SLURM_JOB_ID}.out $TOP_DIR