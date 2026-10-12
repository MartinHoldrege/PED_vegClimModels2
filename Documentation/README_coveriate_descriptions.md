# Layer definitions: soils and Daymet climate rasters

Both rasters are on the CONUS Daymet 1 km grid, masked to the snap raster 
daymet_conus_snap_1000m.tif which other downstream products will be aligned with 
(created in `Analysis/Cover/00_create_snap_raster.R`)

## `soil_covariates_solus100_1000m.tif`

Derived from SOLUS100 (Nauman et al. 2024) by `01_soils_aggregate.R` and
`02_soils_calculate_variables.R`. "Profile mean" is a thickness-weighted mean
over the soil profile, truncated at the depth to bedrock.

| Layer | Definition | Units |
|---|---|---|
| `soilDepth` | Depth to a root-restricting layer, capped at 201 cm | cm |
| `clay_surface` | Clay content, 0–3 cm | % of fine earth |
| `clay` | Clay content, profile mean | % of fine earth |
| `sand` | Sand content, profile mean | % of fine earth |
| `coarse` | Coarse fragments, profile mean | % by volume |
| `carbon` | Soil organic carbon, 0–3 cm | % |
| `AWHC` | Available water holding capacity, summed over the profile | cm |

## `DaymetClimateData_1991-2020_CLIM.tif`

Annual metrics from monthly Daymet v4 data, calculated per year and then
reduced across years by `01_summarise_daymet_climate_data.R`. The suffix gives
the across-year reduction: `_meanAnnAvg_CLIM` and `_meanAnnTotal_CLIM` are
means, `_95percentile_CLIM` and `_5percentile_CLIM` are percentiles. Layers are
listed here.

| Layer | Definition | Units |
|---|---|---|
| `tmin_meanAnnAvg_CLIM` | Annual mean minimum temperature | °C |
| `tmax_meanAnnAvg_CLIM` | Annual mean maximum temperature | °C |
| `tmean_meanAnnAvg_CLIM` | Annual mean temperature | °C |
| `prcp_meanAnnTotal_CLIM` | Total annual precipitation | mm |
| `T_warmestMonth_meanAnnAvg_CLIM` | Mean daily maximum temperature of the warmest month | °C |
| `T_coldestMonth_meanAnnAvg_CLIM` | Mean daily minimum temperature of the coldest month | °C |
| `precip_wettestMonth_meanAnnAvg_CLIM` | Precipitation of the wettest month | mm |
| `precip_driestMonth_meanAnnAvg_CLIM` | Precipitation of the driest month | mm |
| `precip_Seasonality_meanAnnAvg_CLIM` | Precipitation seasonality, as a coefficient of variation among months | unitless |
| `PrecipTempCorr_meanAnnAvg_CLIM` | Correlation between monthly precipitation and temperature | correlation |
| `aboveFreezing_month_meanAnnAvg_CLIM` | First month with minimum temperature above freezing | month (1–12) |
| `isothermality_meanAnnAvg_CLIM` | Isothermality (monthly temperature range relative to annual range) | % |
| `annWaterDeficit_meanAnnAvg_CLIM` | Annual water deficit | mm |
| `annWetDegDays_meanAnnAvg_CLIM` | Annual wet degree days | degree days |
| `annVPD_mean_meanAnnAvg_CLIM` | Annual mean vapor pressure deficit | Pa |
| `annVPD_max_meanAnnAvg_CLIM` | Vapor pressure deficit of the highest month | Pa |
| `annVPD_min_meanAnnAvg_CLIM` | Vapor pressure deficit of the lowest month | Pa |
| `annVPD_max_95percentile_CLIM` | Vapor pressure deficit of the highest month, 95th percentile across years | Pa |
| `annWaterDeficit_95percentile_CLIM` | Annual water deficit, 95th percentile across years | mm |
| `annWetDegDays_5percentile_CLIM` | Annual wet degree days, 5th percentile across years | degree days |
| `durationFrostFreeDays_5percentile_CLIM` | Frost-free period, 5th percentile across years | days |
| `durationFrostFreeDays_meanAnnAvg_CLIM` | Frost-free period | days |


