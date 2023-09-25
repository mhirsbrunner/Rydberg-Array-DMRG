#!/bin/bash
#SBATCH -N 1
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -J Lieb_Local_PD
#SBATCH --mail-user=mark@mhirsbrunner.com
#SBATCH --mail-type=ALL
#SBATCH -t 24:00:00
#SBATCH --account=m888

module load julia
module load python

julia -J../sys_rydberg.so local_detuning_phase_diagram.jl $N_WORKERS $N_THREADS $CONFIG_DIR $DATA_DIR --write_dir $WRITE_DIR --outputlevel 1

# python ../data_analysis/compile_results.py $TOP_DIR

# rm -rf $TOP_DIR/phase_diagram_data

mv slurm-${SLURM_JOB_ID}.out $TOP_DIR