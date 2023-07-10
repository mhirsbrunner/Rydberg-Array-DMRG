mkdir temp_data_dir

podman-hpc run --rm -it \
--mount type=bind,src=/global/homes/m/mhirsbru/Rydberg-Array-DMRG/config,dst=/config,readonly \
--mount type=bind,src=/global/homes/m/mhirsbru/Rydberg-Array-DMRG/temp_data_dir,dst=/data_dir \
--mount type=bind,src=/pscratch/sd/m/mhirsbru,dst=/SCRATCH \
mhirsbru/rydberg_dmrg julia -JRydbergArrayDMRGSysimage.so