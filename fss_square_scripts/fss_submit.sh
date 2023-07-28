#!/bin/bash
N_WORKERS=20
N_THREADS=12

declare -a n_y_vals=(4 6, 8)

TIMESTAMP=$(date +%y-%m-%d_%H-%M-%S)
TOP_DIR=/global/homes/m/mhirsbru/Rydberg-Array-DMRG/data/square/${TIMESTAMP}

if [ -d $TOP_DIR ] 
then
    echo "Data directory already exists" 
    exit 1
else
    mkdir $TOP_DIR
fi

CONFIG_DIR=$TOP_DIR/config
cp -r /global/homes/m/mhirsbru/Rydberg-Array-DMRG/configs/square $CONFIG_DIR

DATA_DIR=$TOP_DIR/phase_diagram_data
mkdir $DATA_DIR

WRITE_DIR=/pscratch/sd/m/mhirsbru/${TIMESTAMP}
mkdir $WRITE_DIR


for n_y in "${n_y_vals[@]}"
do
   sbatch --export=ALL,TOP_DIR=$TOP_DIR,CONFIG_DIR=$CONFIG_DIR,DATA_DIR=$DATA_DIR,WRITE_DIR=$WRITE_DIR,N_WORKERS=$N_WORKERS,N_THREADS=$N_THREADS,N_Y=$n_y fss.sbatch
done