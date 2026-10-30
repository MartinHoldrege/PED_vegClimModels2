#!/bin/bash
# Copy outputs from Hovenweep to the local large-file drive with rsync.
# Only new or changed files are transferred.
#
# Run locally from Git Bash:  bash hpc/sync_from_hpc.sh
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
sync_maca_climate=true

# MACA climate normals (01_summarise_projected_climate_data.R) ------------
# final rasters only; the per-year intermediate folders are skipped
if $sync_maca_climate; then
  sub=Data_processed/WallToWallClimateData/maca
  mkdir -p "$local_root/$sub"
  rsync -avP --include='*maca*_CLIM.tif' --exclude='*' \
    "$remote/$sub/" "$local_root/$sub/"
fi
