#!/bin/sh

#SBATCH --job-name=sim_ebay
#SBATCH --time=3:00:00
#SBATCH --mail-user=wangmk@umich.edu
#SBATCH --mail-type=END,FAIL,BEGIN
#SBATCH --array=1-100
#SBATCH --mem=6g
#SBATCH --cpus-per-task=6
#SBATCH --output=/home/wangmk/MDAWG/POLDA/simulation/eBay/slurm-sim-eBay.out

module load R/4.2.2
Rscript --vanilla /home/wangmk/MDAWG/POLDA/simulation/eBay/ebay_simulation.R

