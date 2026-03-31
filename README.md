# PED_vegClimModels

This version of the repository is for the part of the PED project focusing on modeling biomass and phenology

# Versioning naming used

## Data versions (`vd`)

### Simulated data (for testing)

- **s01** -- early version, can't be recreated without older commit
- **s02** -- 10x range in intercepts (20, 148, 148, 12, 7, 4); multiplicative per-PFT noise (0.2), observation noise (0.4), region bias (0.4)
- **s03** -- 100x range in intercepts (2, 148, 148, 1.2, 0.7, 0.4); same noise structure as s02
- **s04** -- 100x range; tree cover < 1% trimmed to zero
- **s05** -- 100x range; tree cover < 10% trimmed to zero
- **s06** -- 100x range; tree cover < 10% trimmed to zero; correctly specified (no per-PFT or region noise, sigma_obs = 0.5)
- **s07** -- 100x range; tree cover < 10% and shrub cover < 10% trimmed to zero; correctly specified (sigma_obs = 0.5)

### Real data

- **d01** -- no trimming of tree cover
- **d02** -- tree cover < 1% trimmed to zero
- **d03** -- tree cover < 10% trimmed to zero
- **d04** -- tree cover < 10% and shrub cover < 10% trimmed to zero

## Purer pixel selection versions (`vp`)

- **p01** -- q = 0.9, min_raw_cover = 0.05 (moderate purity)
- **p02** -- q = 0.98, min_raw_cover = 0.05 (high purity)
- **p03** -- q = 0 (no purity filtering)

## Model versions (`vm`)

- **m01** -- 6 PFTs | standard LL | 4 pred | no pre-fit alpha. Predictors: tmean, precip, PrecipTempCorr, sand.
- **m02** -- 6 PFTs | standard LL | 3 pred | no pre-fit alpha. Predictors: tmean, precip, PrecipTempCorr.
- **m03** -- 6 PFTs | log(y+1) LL | 4 pred | no pre-fit alpha. Predictors: tmean, precip, PrecipTempCorr, sand.
- **m04** -- 6 PFTs | log(y+1) LL | 3 pred | no pre-fit alpha. Predictors: tmean, precip, PrecipTempCorr.
- **m05** -- 3 PFTs | log(y+1) LL | 3 pred | no pre-fit alpha. PFTs: shrub, totalTree, totalHerbaceous.
- **m06** -- 3 PFTs | alpha² LL | 3 pred | no pre-fit alpha. Alpha constrained >= 0 via alpha² parameterization.
- **m07** -- 3 PFTs | log(y+1) LL | 3 pred | pre-fit alpha (totalHerbaceous). Tree/shrub cover < 10% filter for alpha subset.
- **m08** -- 6 PFTs | log(y+1) LL | 4 pred | pre-fit alpha (C3Gram, C4Gram, Forb). Tree/shrub cover < 10% filter for alpha subset.

