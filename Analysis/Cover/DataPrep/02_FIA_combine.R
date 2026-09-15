# 02_FIA_combine.R
#
# Combine the FIA cover and basal-area tables into one plot-year table
# matching the LDC and LFRDB outputs.
#
# Cover in P2VEG is percent of the CONDITION, not the plot, so conditions are
# averaged weighted by CONDPROP_UNADJ. 01_ already dropped water, disturbed,
# treated and agricultural/developed conditions, so the plot value describes
# the retained portion only - prop_retained records how much of the plot that
# is, and plots with little retained area should be treated cautiously.
#
# Needle/broad is a basal-area fraction, not a cover fraction: FIA records no
# species-level tree cover. Total tree cover comes from P2VEG (TT + NT).
#
# FIA has no usable bare ground (see 01_) for our (above canopy) perspective,
# so cov_bare is absent.
#
# Coordinates are FIA's fuzzed LAT/LON, displaced up to ~1.6 km and
# in some cases swapped between plots within a county.
#
# Inputs:  vegetationComposition<suffix>.csv, TREEtable<suffix>.csv
# Output:  fia_cover_by_plot_year.csv
#
# September, 2026

source('Functions/init.R')

out_dir <- file.path(paths$large, "Data_processed/FIA/")
suffix  <- '_fia1'

veg  <- read_csv(file.path(out_dir, paste0("vegetationComposition", suffix, ".csv")),
                 guess_max = Inf)
tree <- read_csv(file.path(out_dir, paste0("TREEtable_", suffix, ".csv")),
                 guess_max = Inf)

cond_keys <- c("PLT_CN", "INVYR", "STATECD", "UNITCD", "COUNTYCD", "PLOT", "CONDID")
plot_keys <- setdiff(cond_keys, "CONDID")

# ---- condition-level cover, in the LDC/LFRDB groups ----------------------
# A growth habit absent from a condition means zero cover, not unassessed:
# P2VEG records all five habits where it was collected.
cond_cover <- veg |>
  mutate(across(ends_with("_AerialCover"), \(z) replace_na(z, 0)),
         cov_tree       = TallyTree_AerialCover + NonTallyTree_AerialCover,
         cov_tree       = pmin(cov_tree, 100),
         cov_shrub      = Shrub_AerialCover,
         cov_herbaceous = Forbs_AerialCover + Graminoid_AerialCover,
         cov_herbaceous = pmin(cov_herbaceous, 100),
         cov_forb       = Forbs_AerialCover,
         cov_graminoid  = Graminoid_AerialCover) |>
  select(all_of(cond_keys), CONDPROP_UNADJ, MEASYEAR, LAT, LON, STATENAME,
         starts_with("cov_"))

# ---- condition-level needle/broad basal-area fractions -------------------
cond_frac <- tree |>
  transmute(across(all_of(cond_keys)),
            ba_total  = basalArea_allGroups_in2,
            ba_needle = basalArea_needle_in2,
            ba_broad  = basalArea_broad_in2)

# ---- average conditions to plot, weighted by condition area --------------
#' Area-weighted mean over retained conditions
#' @param x Values to average.
#' @param w Condition proportions.
wmean <- function(x, w) {
  ok <- !is.na(x) & !is.na(w)
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

plot_year <- cond_cover |>
  left_join(cond_frac, by = cond_keys) |>
  summarise(
    across(starts_with("cov_"), \(z) wmean(z, CONDPROP_UNADJ)),
    # basal area is an amount, not a percentage, so sum weighted by area
    ba_needle     = sum(ba_needle * CONDPROP_UNADJ, na.rm = TRUE),
    ba_broad      = sum(ba_broad  * CONDPROP_UNADJ, na.rm = TRUE),
    ba_total      = sum(ba_total  * CONDPROP_UNADJ, na.rm = TRUE),
    prop_retained = sum(CONDPROP_UNADJ, na.rm = TRUE),
    n_conditions  = n(),
    year          = first(MEASYEAR),
    Latitude_NAD83  = first(LAT),
    Longitude_NAD83 = first(LON),
    STATENAME     = first(STATENAME),
    .by = all_of(plot_keys)
  ) |>
  mutate(
    ba_known = ba_needle + ba_broad, # some of ba_total is unknown whether broad or needle
    frac_needle = if_else(ba_known > 0, ba_needle / ba_known, NA_real_),
    frac_broad  = if_else(ba_known > 0, ba_broad  / ba_known, NA_real_)
  ) |> 
  select(-ba_known)

min_prop_retained <- 0.5
mean(plot_year$prop_retained < min_prop_retained)*100 # % of plots with less than 50% of the plot retained
plot_year <- plot_year |> filter(prop_retained >= min_prop_retained)

# ---- checks --------------------------------------------------------------
message("plot-years: ", nrow(plot_year),
        " | conditions in: ", nrow(cond_cover),
        " | with >1 condition: ", sum(plot_year$n_conditions > 1))

print("\nprop_retained:")
summary(plot_year$prop_retained) |> print()

print("\nmedian cover:")
plot_year |>
  summarise(across(c(starts_with("cov_"), starts_with("frac_")),
                   \(z) round(median(z, na.rm = TRUE), 2))) |>
  glimpse()

# cover is a percentage of condition area, so plot values must stay in [0, 100]
stopifnot(all(plot_year$cov_tree       <= 100 + 1e-6, na.rm = TRUE),
          all(plot_year$cov_shrub      <= 100 + 1e-6, na.rm = TRUE),
          all(plot_year$cov_herbaceous <= 100 + 1e-6, na.rm = TRUE))

# fractions must sum to 1 where any basal area exists
frac_sum <- plot_year$frac_needle + plot_year$frac_broad
stopifnot(all(is.na(frac_sum) | abs(frac_sum - 1) < 0.02))

dup <- plot_year |>
  select(Latitude_NAD83, Longitude_NAD83, year) |>
  duplicated() |> sum()
if (dup > 0) message("WARNING: ", dup, " duplicate location-year rows")

write_csv(plot_year, file.path(out_dir, 
                               paste0("fia_cover_by_plot_year", suffix,
                                      ".csv")))
