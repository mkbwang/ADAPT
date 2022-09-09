#!/bin/sh

#SBATCH --job-name=zoib_model
#SBATCH --time=2:00:00
#SBATCH --mail-user=wangmk@umich.edu
#SBATCH --mail-type=END,FAIL,BEGIN
#SBATCH --mem=1g
#SBATCH --cpus-per-task=4


Rscript /home/wangmk/MDAWG/DiffRatio/ZOIB_trial.R 

