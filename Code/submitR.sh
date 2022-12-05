#!/bin/bash
#SBATCH --nodes=20
#SBATCH --ntasks-per-node=40
#SBATCH --time=0-02:30           # time (DD-HH:MM)
#module load NiaEnv/2018a
module purge
module load CCEnv
module load StdEnv/2020
module load gcc/9.3.0 openmpi/4.0.3 r/4.0.2

R_PROFILE=${HOME}/R/x86_64-pc-linux-gnu-library/4.0/snow/RMPISNOWprofile; export R_PROFILE
mpirun -np 800 -bind-to core:overload-allowed R CMD BATCH --no-save "--args $pheno_id" ${name}.R
${name}_pheno${pheno_id}.Rout

# for pheno_id in f.23248.2.0 f.23260.2.0 f.23265.2.0 f.23283.2.0 f.23277.2.0 f.23287.2.0
# do
# 	sbatch --export=pheno_id=$pheno_id,name="data_analysis" submitR.sh
# done
