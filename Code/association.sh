#!/bin/bash
#SBATCH --nodes=20
#SBATCH --ntasks-per-node=1
#SBATCH --time=0-05:00           # time (DD-HH:MM)
module purge
module load CCEnv
module load StdEnv/2020
module load plink2

./plink2 --bfile allchromosome  \
	  --keep id.txt  \
	  --pheno phenotype.txt \
	  --covar covariate.txt \
	  --covar-variance-standardize \
	  --linear