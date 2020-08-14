#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=16:00:00
#SBATCH --mem=20GB
#SBATCH --job-name=myTest
#SBATCH --mail-type=END
#SBATCH --mail-user=jll616@nyu.edu
#SBATCH --output=slurm_%j.out

module purge
module load matlab/2020a

cd /scratch/jll616/matlab-scripts

cat<<EOF | matlab -nodisplay 
maxNumCompThreads($SLURM_CPUS_PER_TASK);
A_modelfitting($SLURM_ARRAY_TASK_ID)
exit
EOF





