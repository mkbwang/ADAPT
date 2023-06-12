#!/bin/sh

#SBATCH --job-name=sim_ALDEx2
#SBATCH --time=2:00:00
#SBATCH --mail-user=wangmk@umich.edu
#SBATCH --mail-type=END,FAIL,BEGIN
#SBATCH --array=1-100
#SBATCH --mem=4g
#SBATCH --cpus-per-task=6
#SBATCH --output=/home/wangmk/MDAWG/POLDA/simulation/Aldex2/slurm-sim-Aldex2.out

module load R/4.2.2
Rscript --vanilla /home/wangmk/MDAWG/POLDA/simulation/Aldex2/Aldex_simulation.R

