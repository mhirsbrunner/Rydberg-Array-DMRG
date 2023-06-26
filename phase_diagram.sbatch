#!/bin/bash
#SBATCH -N 1
#SBATCH -C cpu
#SBATCH -q debug
#SBATCH -J 1d_z2
#SBATCH --mail-user=markhirsbrunner@lbl.gov
#SBATCH --mail-type=ALL
#SBATCH -t 00:20:00
#SBATCH --account=m888

module load parallel

N_CORES_TOT=120 #leave some cores idle, change this to threads
N_THREADS_PER_PROC=8
let N_PARALLEL=2*$N_CORES_TOT/$N_THREADS_PER_PROC

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

RB_OVER_A_MIN=1.0
RB_OVER_A_MAX=4.0
RB_OVER_A_STEP=0.25
RB_OVER_A_AX=($(seq $RB_OVER_A_MIN $RB_OVER_A_STEP $RB_OVER_A_MAX))

DELTA_OVER_OMEGA_MIN=0.0
DELTA_OVER_OMEGA_MAX=4.0
DELTA_OVER_OMEGA_STEP=0.5
DELTA_OVER_OMEGA_AX=($(seq $DELTA_OVER_OMEGA_MIN $DELTA_OVER_OMEGA_STEP $DELTA_OVER_OMEGA_MAX))


INPUT_FNAME="${DATA_DIR}/inputs.txt"

touch $INPUT_FNAME

for RB_OVER_A in "${RB_OVER_A_AX[@]}"
do
for DELTA_OVER_OMEGA in "${DELTA_OVER_OMEGA_AX[@]}"
do
    echo "${PWD}/${DATA_DIR}/config ${PWD}/${DATA_DIR}/phase_diagram_data $N_THREADS_PER_PROC $RB_OVER_A $DELTA_OVER_OMEGA ${RB_OVER_A}_${DELTA_OVER_OMEGA}" >> $INPUT_FNAME
done
done

echo "Got to srun"

srun -n 1 -c $N_THREADS_PER_PROC --cpu_bind=cores \
parallel  --colsep ' ' --jobs $N_PARALLEL \
./container_payload.sh {} :::: $INPUT_FNAME



# Need to ask Johannes why this takes so long to get started, if it starts at all


# try fireworks rather than GNU parallel, or parcel, or snake_make**, Distributed and ClusterManagers within Julia
# podman hpc on nersc documentation
# use pkg to make your code into packages