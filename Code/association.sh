#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=0-00:15           # time (DD-HH:MM)
module purge
module load CCEnv
module load StdEnv/2020
module load plink/2.00-10252019-avx2

plink2 --bfile allchromosome --pheno Data/height_imputed.txt --pheno-name int_oracle --covar Data/height_covariate.txt --glm hide-covar