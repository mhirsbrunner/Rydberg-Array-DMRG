#!/bin/bash

# docker run --rm \
# --mount type=bind,src=/home/mark/Dropbox/LBL/Rydberg_Array_DMRG/config,dst=/config,readonly \
# --mount type=bind,src=/home/mark/Dropbox/LBL/Rydberg_Array_DMRG/test_data,dst=/data \
# mhirsbrunner/rydberg_array_dmrg 4 1.5 4.0 data test

# /global/homes/m/mhirsbru/Rydberg-Array-DMRG/config
# /global/homes/m/mhirsbru/Rydberg-Array-DMRG/test_data

CURRENTTIME=`date +"%T"`
echo "Calling container_payload at ${CURRENTTIME} for params ${6}!"

podman-hpc run --rm \
--mount type=bind,src=$1,dst=/config,readonly \
--mount type=bind,src=$2,dst=/data_dir \
mhirsbru/rydberg_image $3 $4 $5 data_dir $6

CURRENTTIME=`date +"%T"`
echo "Finished container_payload call at ${CURRENTTIME} for params ${6}!"