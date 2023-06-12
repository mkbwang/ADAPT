#!/bin/sh

#SBATCH --job-name=sim_corncob
#SBATCH --time=3:00:00
#SBATCH --mail-user=wangmk@umich.edu
#SBATCH --mail-type=END,FAIL,BEGIN
#SBATCH --array=1-100
#SBATCH --mem=5g
#SBATCH --cpus-per-task=6
#SBATCH --output=/home/wangmk/MDAWG/POLDA/simulation/corncob/slurm-sim-corncob.out

module load R/4.2.2
Rscript --vanilla /home/wangmk/MDAWG/POLDA/simulation/corncob/corncob_simulation.R

