#!/bin/bash
# Download MACA monthly projections.
# Submit from the repository root:  sbatch hpc/download_maca.sh
#SBATCH --job-name=maca_download
#SBATCH --account=swbsc
#SBATCH --partition=cpu
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/maca_download_%j.out

module load cray-R-spatial/4.6

Rscript Analysis/Cover/DataPrep/00_download_maca_rasters.R