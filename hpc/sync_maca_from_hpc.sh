#!/bin/bash
# Copy the final MACA climate-normal rasters (not the per-year intermediate
# files) from Hovenweep to the local paths$large directory.
#
# Run locally from Git Bash (any directory):  bash hpc/sync_maca_from_hpc.sh
# Re-running copies only new or changed files.
#
# Requires:
#   - rsync, and ssh access to the HPC through the `hwdtn` alias in ~/.ssh/config
#   - calderaProject exported in ~/.bashrc (the team project directory on the
#     HPC)

set -e

if [ -z "$calderaProject" ]; then
  echo "calderaProject is not set; export it in ~/.bashrc" >&2
  exit 1
fi

# local paths$large (same logic as config/paths.R)
if [ -d /e/USGS ]; then
  local_base=/e/USGS
elif [ -d /d/USGS ]; then
  local_base=/d/USGS
else
  echo "local large-file directory not found" >&2
  exit 1
fi

sub=Data_processed/WallToWallClimateData
local_dir=$local_base/large_files/PED_vegClimModels2/$sub
mkdir -p "$local_dir"

rsync -avP \
  --include='*maca*_CLIM.tif' --exclude='*' \
  "hwdtn:${calderaProject%/}/PED_vegClimModels2/$sub/" \
  "$local_dir/"