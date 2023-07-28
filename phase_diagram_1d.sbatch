#!/bin/bash
#SBATCH -N 1
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -J DMRG
#SBATCH --mail-user=markhirsbrunner@lbl.gov
#SBATCH --mail-type=ALL
#SBATCH -t 6:00:00
#SBATCH --account=m888

module load julia
module load python

julia -Jsys_rydberg.so -p$N_WORKERS generate_phase_diagram.jl $CONFIG_DIR $DATA_DIR --write_dir $WRITE_DIR

python compile_results.py $TOP_DIR

rm -rf $TOP_DIR/phase_diagram_data

mv slurm-${SLURM_JOB_ID}.out $TOP_DIR