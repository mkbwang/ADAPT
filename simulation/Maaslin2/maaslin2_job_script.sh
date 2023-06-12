#!/bin/sh

#SBATCH --job-name=sim_maaslin2
#SBATCH --time=3:00:00
#SBATCH --mail-user=wangmk@umich.edu
#SBATCH --mail-type=END,FAIL,BEGIN
#SBATCH --array=8,16,21,43,61,64
#SBATCH --mem=7g
#SBATCH --cpus-per-task=12
#SBATCH --output=/home/wangmk/MDAWG/POLDA/simulation/Maaslin2/slurm-sim-Maaslin2.out

module load R/4.2.2
Rscript --vanilla /home/wangmk/MDAWG/POLDA/simulation/Maaslin2/maaslin2_simulation.R

