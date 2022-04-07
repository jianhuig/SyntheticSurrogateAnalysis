#!/bin/bash
#SBATCH --nodes=20
#SBATCH --ntasks-per-node=40
#SBATCH --time=0-23:00           # time (DD-HH:MM)
#module load NiaEnv/2018a
module purge
module load CCEnv
module load StdEnv/2020
module load gcc/9.3.0 openmpi/4.0.3 r/4.0.2

R_PROFILE=${HOME}/R/x86_64-pc-linux-gnu-library/4.0/snow/RMPISNOWprofile; export R_PROFILE
mpirun -np 800 -bind-to core:overload-allowed R CMD BATCH --no-save ${name}.R

# sbatch --export=name="mean_comparision" submit.sh