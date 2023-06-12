#!/bin/sh

#SBATCH --job-name=DirGamma_generation
#SBATCH --time=00:40:00
#SBATCH --mail-user=wangmk@umich.edu
#SBATCH --mail-type=END,FAIL,BEGIN
#SBATCH --array=1-100
#SBATCH --mem=2g
#SBATCH --cpus-per-task=2
#SBATCH --output=/home/wangmk/MDAWG/POLDA/simulation/simulate_DirMultinom/slurm-sim.out

module load R/4.2.2
Rscript --vanilla /home/wangmk/MDAWG/POLDA/simulation/simulate_DirMultinom/simulate_DirMultinom.R
