docker run --rm \
--mount type=bind,src=/home/mark/Dropbox/LBL/Rydberg_Array_DMRG/config,dst=/config,readonly \
--mount type=bind,src=/home/mark/Dropbox/LBL/Rydberg_Array_DMRG/test_data,dst=/data \
mhirsbrunner/rydberg_array_dmrg 4 1.5 4.0 data test