#!/bin/bash
# gridMET 1991-2020 climate normals: aggregate daily gridMET to monthly,
# then summarise to _CLIM normals.
# Submit from the repository root:   sbatch hpc/summarise_gridmet.sh
# Both R scripts skip finished years, so a failed job can be resubmitted.

#SBATCH --job-name=gridmet_summarise
#SBATCH --account=swbsc
#SBATCH --partition=cpu
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=logs/gridmet_summarise_%j.out

set -e # don't run the normals step if the monthly step fails

module load cray-R-spatial/4.6

Rscript Analysis/Cover/DataPrep/00_gridmet_daily_to_monthly.R
Rscript Analysis/Cover/DataPrep/01_summarise_gridmet_climate_data.R
