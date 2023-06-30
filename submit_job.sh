#!/bin/bash
TIMESTAMP=$(date +%y-%m-%d_%H-%M-%S)
DATA_DIR="data/${TIMESTAMP}"

if [ -d $DATA_DIR ] 
then
    echo "Data directory already exists" 
    exit 1
else
    mkdir $DATA_DIR
fi

cp -r config $DATA_DIR/config
mkdir $DATA_DIR/phase_diagram_data

sbatch --export=ALL,DATA_DIR=$DATA_DIR phase_diagram.sbatch