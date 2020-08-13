#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --time=16:00:00
#SBATCH --mem=20GB
#SBATCH --job-name=myTest
#SBATCH --mail-type=END
#SBATCH --mail-user=jll616@nyu.edu
#SBATCH --output=slurm_%j.out

module purge
module load matlab/2017a

cd
cd scratch/jll616/matlab-scripts
cat A_modelfitting.m | srun matlab -nodisplay