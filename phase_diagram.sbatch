#!/bin/bash
#SBATCH -N 1
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -J 1d_z2
#SBATCH --mail-user=markhirsbrunner@lbl.gov
#SBATCH --mail-type=ALL
#SBATCH -t 06:00:00
#SBATCH --account=m888

podman-hpc run --rm -it \
--mount type=bind,src=/global/homes/m/mhirsbru/Rydberg-Array-DMRG/$DATA_DIR/config,dst=/config,readonly \
--mount type=bind,src=/global/homes/m/mhirsbru/Rydberg-Array-DMRG/$DATA_DIR/phase_diagram_data,dst=/data_dir \
--mount type=bind,src=/pscratch/sd/m/mhirsbru,dst=/SCRATCH \
mhirsbru/rydberg_dmrg julia -JRydbergArrayDMRGSysimage.so -p20 generate_phase_diagram_1d.jl