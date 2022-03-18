#!/bin/bash
#SBATCH --nodes=20
#SBATCH --ntasks-per-node=1
#SBATCH --time=0-05:00

module load plink
module load gnu-parallel/20191122
scontrol show hostname ${SLURM_JOB_NODELIST} > ./node_list_${SLURM_JOB_ID}

seq 20 | parallel --sshloginfile ./node_list_${SLURM_JOB_ID}
-j 20 /scratch/l/leisun/jianhuig/plink2 # 20 parallel process
-bfile /scratch/l/leisun/jianhuig/plink
--make-king triangle bin4 # kingship matrix
-parallel {} 20 # i out of 20 small batch file
--out /scratch/l/leisun/jianhuig/split