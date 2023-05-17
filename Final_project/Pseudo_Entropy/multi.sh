#!/bin/bash

# Runtime and memory

#SBATCH --mem-per-cpu=8G 
#SBATCH --time=144:00:00

# For parallel jobs

#SBATCH --cpus-per-task=3
#SBATCH --ntasks=1

#SBATCH -p mem

module load python
source activate compAstro

which python

srun python QMulti2D.py
