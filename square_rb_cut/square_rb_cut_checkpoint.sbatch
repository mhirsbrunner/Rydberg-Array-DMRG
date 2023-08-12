#!/bin/bash
#SBATCH -N 1
#SBATCH -C cpu
#SBATCH -q debug
#SBATCH -J sq_rb_cut_DMRG
#SBATCH --account=m4370
#SBATCH --mail-user=markhirsbrunner@lbl.gov
#SBATCH --mail-type=ALL
#SBATCH -t 00:3:00
#SBATCH --time-min=00:02:00
#SBATCH --comment=01:00:00
#SBATCH --signal=B:USR1@60
#SBATCH --requeue
#SBATCH --open-mode=append


module load julia
module load python
module load dmtcp nersc_cr

export LD_DEBUG=all

# source $HOME/dmtcp/env_setup_2.sh

start_coordinator

#c/r jobs
if [[ $(restart_count) == 0 ]]; then

    #user setting
#     $HOME/dmtcp/bin/dmtcp_launch -j julia --sysimage ../sys_rydberg.so square_rb_cut.jl $N_WORKERS $N_THREADS $CONFIG_DIR $DATA_DIR --write_dir $WRITE_DIR &
    dmtcp_launch -j julia --sysimage ../sys_rydberg.so square_rb_cut.jl $N_WORKERS $N_THREADS $CONFIG_DIR $DATA_DIR --write_dir $WRITE_DIR &
elif [[ $(restart_count) > 0 ]] && [[ -e dmtcp_restart_script.sh ]]; then

    ./dmtcp_restart_script.sh &
else

    echo "Failed to restart the job, exit"; exit
fi

# requeueing the job if remaining time >0
ckpt_command=ckpt_dmtcp    #additional checkpointing right before the job hits the wall limit 
requeue_job func_trap USR1

wait

# mv slurm-${SLURM_JOB_ID}.out $TOP_DIR