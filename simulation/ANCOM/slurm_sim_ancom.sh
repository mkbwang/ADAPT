#!/bin/sh

#SBATCH --job-name=sim_ANCOM
#SBATCH --time=8:00:00
#SBATCH --mail-user=wangmk@umich.edu
#SBATCH --mail-type=END,FAIL,BEGIN
#SBATCH --array=1-100
#SBATCH --mem=6g
#SBATCH --cpus-per-task=12
#SBATCH --output=/home/wangmk/MDAWG/POLDA/simulation/ANCOM/slurm-sim-ANCOM.out

module load R/4.1.2-gcc8.3.0
Rscript --vanilla /home/wangmk/MDAWG/POLDA/simulation/ANCOM/ancom_simulation.R

