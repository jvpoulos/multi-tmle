#!/bin/bash
#SBATCH --job-name="misTreat"
#SBATCH -e misTreat.err
#SBATCH -p mpi
#SBATCH -n 160
#SBATCH --mincpus=18
#SBATCH --mem-per-cpu=4G
#SBATCH -t 5-00:00:00
#SBATCH --array=18

ulimit -l unlimited

module load gcc/6.2.0 R/4.0.1
module load openmpi/4.1.1

#Run the job
mpirun --mca btl tcp,self Rscript tmle_MultinomialTrts.R ${SLURM_ARRAY_TASK_ID} 'binomial' 'TRUE' 'TRUE' 'FALSE' 'FALSE' 'TRUE' 'FALSE' 'FALSE'