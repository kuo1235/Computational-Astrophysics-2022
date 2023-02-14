#!/bin/bash
#SBATCH --ntasks=1 
#SBATCH --time=144:00:00
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH -p mem

module load python
source activate compAstro

which python

srun python PE_3.py
