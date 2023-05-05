#!/bin/sh

#SBATCH --job-name=sim_ANCOMBC
#SBATCH --time=1:00:00
#SBATCH --mail-user=wangmk@umich.edu
#SBATCH --mail-type=END,FAIL,BEGIN
#SBATCH --array=1-100
#SBATCH --mem=3g
#SBATCH --cpus-per-task=12
#SBATCH --output=/home/wangmk/MDAWG/POLDA/simulation/ANCOMBC/slurm-sim-ANCOMBC.out

module load R/4.1.2-gcc8.3.0
Rscript --vanilla /home/wangmk/MDAWG/POLDA/simulation/ANCOMBC/ancombc_simulation.R

