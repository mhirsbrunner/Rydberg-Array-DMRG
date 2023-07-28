#!/bin/bash
N_WORKERS=128

TIMESTAMP=$(date +%y-%m-%d_%H-%M-%S)
TOP_DIR=data/1d/${TIMESTAMP}

if [ -d $TOP_DIR ] 
then
    echo "Data directory already exists" 
    exit 1
else
    mkdir $TOP_DIR
fi

CONFIG_DIR=$TOP_DIR/config
cp -r configs/1d $CONFIG_DIR

DATA_DIR=$TOP_DIR/phase_diagram_data
mkdir $DATA_DIR

WRITE_DIR=/pscratch/sd/m/mhirsbru/${TIMESTAMP}
mkdir $WRITE_DIR

sbatch --export=ALL,TOP_DIR=$TOP_DIR,CONFIG_DIR=$CONFIG_DIR,DATA_DIR=$DATA_DIR,WRITE_DIR=$WRITE_DIR,N_WORKERS=$N_WORKERS phase_diagram_1d.sbatch