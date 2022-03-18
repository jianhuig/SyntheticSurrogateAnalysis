#!/bin/bash
#SBATCH --nodes=20
#SBATCH --ntasks-per-node=1
#SBATCH --time=0-05:00           # time (DD-HH:MM)
module purge
module load CCEnv
module load StdEnv/2020

cd $SCRATCH
# Obtain 10 genetic PC by randomized algorithm 
# https://pubmed.ncbi.nlm.nih.gov/26924531/
./plink2 --bfile allchromosome --pca approx 10 