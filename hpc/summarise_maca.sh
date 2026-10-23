#!/bin/bash
# Summarise MACA projections to _CLIM normals: one array task per
# GCM x scenario-period (20 x 5 = 100 tasks).
# Submit from the repository root:   sbatch hpc/summarise_maca.sh
# Test a single task (BNU-ESM historical 1971-2000) first:
#                                     sbatch --array=0 hpc/summarise_maca.sh
# this script works b/ slurm 'loops' over the array (0-99) 
# with SLURM_ARRAY_TASK_ID being like 'i' in the loop

#SBATCH --job-name=maca_summarise
#SBATCH --account=swbsc
#SBATCH --partition=cpu
#SBATCH --array=0-99
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=logs/maca_summarise_%A_%a.out

gcms=(
  BNU-ESM CCSM4 CNRM-CM5 CSIRO-Mk3-6-0 CanESM2
  GFDL-ESM2G GFDL-ESM2M HadGEM2-CC365 HadGEM2-ES365
  IPSL-CM5A-LR IPSL-CM5A-MR IPSL-CM5B-LR MIROC-ESM
  MIROC-ESM-CHEM MIROC5 MRI-CGCM3 NorESM1-M
  bcc-csm1-1 bcc-csm1-1-m inmcm4
)
scenarios=(historical rcp45     rcp45     rcp85     rcp85)
periods=(  1971-2000  2040-2069 2070-2099 2040-2069 2070-2099)

n=${#scenarios[@]}
# Note: ${array[index]} returns the value at that index
gcm=${gcms[$((SLURM_ARRAY_TASK_ID / n))]} # integer devision e.g. $((2/5)) = 0 
scenario=${scenarios[$((SLURM_ARRAY_TASK_ID % n))]} # subsetting by the remainder, e.g. 7%5 = 2
period=${periods[$((SLURM_ARRAY_TASK_ID % n))]}

module load cray-R-spatial/4.6

Rscript Analysis/Cover/DataPrep/01_summarise_projected_climate_data.R --gcm="$gcm" --scenario="$scenario" --period="$period"
