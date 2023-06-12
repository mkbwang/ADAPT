#!/bin/sh

#SBATCH --job-name=sim_DACOMP
#SBATCH --time=3:00:00
#SBATCH --mail-user=wangmk@umich.edu
#SBATCH --mail-type=END,FAIL,BEGIN
#SBATCH --array=1-100
#SBATCH --mem=7g
#SBATCH --cpus-per-task=12
#SBATCH --output=/home/wangmk/MDAWG/POLDA/simulation/DACOMP/slurm-sim-DACOMP.out

module load R/4.2.2
Rscript --vanilla /home/wangmk/MDAWG/POLDA/simulation/DACOMP/DACOMP_Simulation.R

