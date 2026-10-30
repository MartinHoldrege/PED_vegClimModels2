#!/bin/bash
# Fit cover models and/or predict them onto the CONUS rasters: one array task
# per row of the runs table in Analysis/Cover/main_cover.R (cover_runs and
# vc_model). The evaluation reports are run locally.
#
# 1. Set cover_runs / vc_model in main_cover.R, push, and git pull here.
# 2. Copy the inputs (locally):      bash hpc/sync_cover_to_hpc.sh
# 3. Submit from the repository root, one task per row of runs (main_cover.R
#    prints them with their run_index; an index past the last row is an
#    error). Fit only (the default):
#                                    sbatch --array=1-6 hpc/run_cover.sh
#    Fit and predict:                sbatch --array=1-6 hpc/run_cover.sh --run_predict=TRUE
#    Predict only (fits already here):
#                                    sbatch --array=1-6 hpc/run_cover.sh --run_fit=FALSE --run_predict=TRUE
#    Test one row first:             sbatch --array=1 hpc/run_cover.sh
# 4. Copy the outputs back (locally): bash hpc/sync_cover_from_hpc.sh
# logs/ must exist (mkdir -p logs).
#
# 10 cores: grplasso fits its 10 CV folds in parallel (.grpl_folds() uses
# future::availableCores(), which follows SLURM_CPUS_PER_TASK), and ranger
# uses them as threads (R_RANGER_NUM_THREADS below). grpreg, glmnet and the
# prediction run on one core, so for predict-only jobs add
# --cpus-per-task=1 to sbatch.
#
# Combined rows: the combining scripts need the component predictions, which
# the other tasks of the same array may not have written yet, so submit them
# in a later array (or run them locally).

#SBATCH --job-name=run_cover
#SBATCH --account=swbsc
#SBATCH --partition=cpu
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=64G
#SBATCH --output=logs/run_cover_%A_%a.out

module load cray-R-spatial/4.6

# ranger uses 2 threads unless told otherwise
export R_RANGER_NUM_THREADS=$SLURM_CPUS_PER_TASK

# arguments to this script ("$@") come after the defaults, so they replace
# them (the last copy of an option wins); reports are never run here
Rscript Analysis/Cover/main_cover.R --run_index="$SLURM_ARRAY_TASK_ID" \
  --run_fit=TRUE --run_predict=FALSE "$@" --run_evaluate=FALSE
