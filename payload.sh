#!/bin/bash
julia --project=. --sysimage ~/.julia/sysimages/sys_itensors.so -t $1 1d_array_dmrg.jl $2 $3 --data_dir $4 --fname $5