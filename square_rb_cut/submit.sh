#!/bin/bash
N_WORKERS=31
N_THREADS=8

TIMESTAMP=$(date +%y-%m-%d_%H-%M-%S)
# TOP_DIR=/global/homes/m/mhirsbru/Rydberg-Array-DMRG/data/square/${TIMESTAMP}
TOP_DIR=/global/homes/m/mhirsbru/Rydberg-Array-DMRG/data/square/fss/16x8
CONFIG_DIR=$TOP_DIR/config
DATA_DIR=$TOP_DIR/phase_diagram_data

if [ -d $TOP_DIR ] 
then
    echo "Data directory already exists"
else
    mkdir $TOP_DIR
    cp -r /global/homes/m/mhirsbru/Rydberg-Array-DMRG/configs/square $CONFIG_DIR
    mkdir $DATA_DIR
fi

WRITE_DIR=/pscratch/sd/m/mhirsbru/${TIMESTAMP}
mkdir $WRITE_DIR

sbatch --export=ALL,TOP_DIR=$TOP_DIR,CONFIG_DIR=$CONFIG_DIR,DATA_DIR=$DATA_DIR,WRITE_DIR=$WRITE_DIR,N_WORKERS=$N_WORKERS,N_THREADS=$N_THREADS sbatch.sh