#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=0-01:00           # time (DD-HH:MM)
module purge
module load CCEnv
module load StdEnv/2020
module load plink/2.00-10252019-avx2

plink2 --bfile allchromosome --pheno height_permuted.txt --covar covariates.txt --covar-variance-standardize --linear