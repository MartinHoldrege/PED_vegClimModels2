#!/bin/bash
# Copy cover model outputs (hpc/run_cover.sh) from Hovenweep to the local
# large-file drive with rsync. Only new or changed files are transferred, and
# files that are newer locally are kept (--update).
#
# Run locally from Git Bash:  bash hpc/sync_cover_from_hpc.sh
# Preview, copying nothing:   bash hpc/sync_cover_from_hpc.sh --dry-run
# Turn sections on/off with the true/false switches below.
#
# Requires rsync, the `hwdtn` alias in ~/.ssh/config, and calderaProject
# exported in ~/.bashrc

set -e

if [ -z "$calderaProject" ]; then
  echo "calderaProject is not set; export it in ~/.bashrc" >&2
  exit 1
fi

remote="hwdtn:${calderaProject%/}/PED_vegClimModels2"
local_root="/d/USGS/large_files/PED_vegClimModels2"  # paths$large

# which parts to run
sync_fits=true
sync_predictions=true

# files to copy (rsync patterns), e.g. '*_c05-m01.51*' for one version
fits_include='*.rds'
predictions_include='*.tif'

dry_run=""
if [ "$1" = "--dry-run" ]; then
  dry_run="--dry-run"
fi

# fitted models (02_fit_*.R) ----------------------------------------------
# <cover_type>_<cover_model>_<vc>-<vmc>.rds (cover_model_path())
if $sync_fits; then
  sub=Data_processed/CoverData/Fit
  mkdir -p "$local_root/$sub"
  rsync -avPu $dry_run --include="$fits_include" --exclude='*' \
    "$remote/$sub/" "$local_root/$sub/"
fi

# prediction rasters (03_predict_raster.R) --------------------------------
# <cover_type>_<cover_model>_<vc>-<vmc>_<scenario>.tif (cover_pred_path())
if $sync_predictions; then
  sub=Data_processed/CoverData/Predictions
  mkdir -p "$local_root/$sub"
  rsync -avPu $dry_run --include="$predictions_include" --exclude='*' \
    "$remote/$sub/" "$local_root/$sub/"
fi
