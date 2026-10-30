#!/bin/bash
# Copy the inputs of the cover model fits and predictions (hpc/run_cover.sh)
# from the local large-file drive to Hovenweep with rsync. Only new or changed
# files are transferred.
#
# Fitting (02_fit_classification.R, 02_fit_cover.R) reads only the training
# data (read_cover_training()) and the scaling parameters
# (read_scale_params()); CV folds are made at fit time (make_cover_folds()).
# Prediction (03_predict_raster.R) reads the fitted model, the snap mask
# (read_mask()), and the climate and soils rasters (read_climate_raster()).
#
# Run locally from Git Bash:  bash hpc/sync_cover_to_hpc.sh
# Preview, copying nothing:   bash hpc/sync_cover_to_hpc.sh --dry-run
# Turn sections on/off with the true/false switches below.
#
# Requires rsync, the `hwdtn` alias in ~/.ssh/config, and calderaProject
# exported in ~/.bashrc.

set -e

if [ -z "$calderaProject" ]; then
  echo "calderaProject is not set; export it in ~/.bashrc" >&2
  exit 1
fi

remote="hwdtn:${calderaProject%/}/PED_vegClimModels2"
local_root="/d/USGS/large_files/PED_vegClimModels2"  # paths$large

# which parts to run
sync_fit_data=true      # training data and scaling parameters
sync_predict_data=true  # mask, climate and soils rasters (~2.3 GB)
sync_fits=false         # models fit locally, to predict them on the HPC

# cover data versions (training data and anomaly scaling parameters)
vcs=(c04 c05 c06)
# fitted models to copy, if sync_fits (rsync pattern)
fits_include='*.rds'

dry_run=""
if [ "$1" = "--dry-run" ]; then
  dry_run="--dry-run"
fi

cd "$local_root"

files=()

# fitting inputs ----------------------------------------------------------
if $sync_fit_data; then
  files+=(Data_processed/scale_params/scale_params_climate_soils.csv)
  for vc in "${vcs[@]}"; do
    files+=("Data_processed/cover_combined/cover_clim_soils_${vc}.csv"
            "Data_processed/scale_params/scale_params_anomalies_${vc}.csv")
  done
fi

# prediction inputs -------------------------------------------------------
# the MACA _CLIM rasters made on the HPC (summarise_maca.sh) are already there
if $sync_predict_data; then
  files+=(Data_processed/masks/daymet_conus_snap_1000m.tif
          Data_processed/soils/soil_covariates_solus100_1000m.tif
          Data_processed/WallToWallClimateData/DaymetClimateData_1991-2020_CLIM.tif
          Data_processed/WallToWallClimateData/ForecastedClimateData_BNU-ESM_rcp85_CLIM.tif
          Data_processed/WallToWallClimateData/ForecastedClimateData_IPSL-CM5A-MR_rcp85_CLIM.tif)
fi

if [ ${#files[@]} -gt 0 ]; then
  for f in "${files[@]}"; do
    if [ ! -f "$f" ]; then
      echo "missing: $local_root/$f" >&2
      exit 1
    fi
  done
  # -R (--relative) recreates the Data_processed/... paths under the repo
  rsync -avPR $dry_run "${files[@]}" "$remote/"
fi

# fitted models -----------------------------------------------------------
# --update: files that are newer on the HPC are kept
if $sync_fits; then
  fits=(Data_processed/CoverData/Fit/$fits_include)
  if [ ! -f "${fits[0]}" ]; then
    echo "no fitted models match: $fits_include" >&2
    exit 1
  fi
  rsync -avPRu $dry_run "${fits[@]}" "$remote/"
fi
