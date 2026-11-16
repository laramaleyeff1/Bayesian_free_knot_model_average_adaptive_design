#!/bin/bash
#SBATCH --account=def-laramal
#SBATCH --nodes=1                # number of node MUST be 1
#SBATCH --cpus-per-task=4        # number of processes
#SBATCH --mem-per-cpu=2048M      # memory; default unit is megabytes
#SBATCH --time=0-02:30           # time (DD-HH:MM)
#SBATCH --array=1-100

module load StdEnv/2020 gcc/9.3.0 r/4.3.1
Rscript simulation_study_i.R $1 $2 $3 $4 $5 $6 $7 $8 $9 ${SLURM_ARRAY_TASK_ID}
