#///////////////////////////////////////////////////////////////////////////
# FIA data wrangling, based heavily on Alice Stears code
# Produces plot-level functional-group cover, bare ground, and tree basal area
# for the vegetation-cover modeling pipeline.
#
# Outputs :
#   vegetationComposition_fiaX.csv - plot-avg aerial cover by functional group
#   groundCover_fiaXcsv           - bare ground (averaged across FIA's 2 sources)
#   TREEtable_fiaX.csv             - plot-avg basal area by taxonomic group + props
#
# Analysis decisions (preserved from the original):
#   - filter water / disturbed / treated / nonforest conditions
#   - drop AK, HI, and island territories
#   - aerial cover = P2VEG LAYER 5, summed within subplot by group, averaged
#     across subplots within a plot
#   - angio/conifer split from BASAL AREA proportions (FIA has no cover split)
#   - NA basal-area group = true zero
#

#///////////////////////////////////////////////////////////////////////////

source('Functions/init.R')
library(FIESTA)

# for taxonlookup (used below)
#install.packages("devtools")
# pak::pak("ropenscilabs/datastorr")
# pak::pak("wcornwell/taxonlookup")

# --- params ----------------------------------------------------------
raw_dir <- file.path(paths$large, "Data_raw/FIA/CSV_FIADB_ENTIRE")   # input FIADB CSV tables
out_dir <- file.path(paths$large, "Data_processed/FIA/")              # output directory
dir.create(out_dir, showWarnings = FALSE)
suffix <- '_fia1' # the output version
drop_states <- c("VI","PR","PW","MP","MH","GU","FM","AS","HI","AK")

# see https://research.fs.usda.gov/sites/default/files/2024-05/wo-v9-2_apr2024_ug_fiadb_database_description_nfi.pdf
# COND filter codes (meanings from the FIA manual):
#   COND_STATUS_CD 3,4 = census / noncensus water
#   DSTRBCD 30/31/32 fire, 46 grazing [decided not to include], 21, 22 disease, 80 human, 91-95 mass movement/geology
#   TRTCD   10 cut, 20 site prep, 30/40 regen, 50 other silviculture
#   PRESNFCD 10-12 agricultural, 30/31/32 developed / cultural / right-of-way
disturb_codes   <- c(30, 31, 32, 21, 22, 80, 91:95)
treat_codes     <- c(10, 20, 30, 40, 50)
nonforest_codes <- c(10, 11, 12, 13, 16, 17,31, 30, 32)

missing_flag <- 999   # sentinel for missing bare-ground values

# --- reference: state codes -------------------------------------------------

state_codes <- FIESTA::ref_statecd |>
  rename(STATECD = VALUE, STATENAME = MEANING, STATEABB = ABBR) |>
  select(-c(RS, RSCD, REGION, REGION_SPGRPCD))

# --- PLOT: location + inventory year ----------------------------------------
plot_tbl <- read_csv(file.path(raw_dir, "ENTIRE_PLOT.csv")) |>  # plot table
  filter(!is.na(LAT), !is.na(LON))

# --- COND: filter to usable conditions, attach location ---------------------
cond <- read_csv(file.path(raw_dir, "ENTIRE_COND.csv"), guess_max = Inf)  |> 
  filter(!(COND_STATUS_CD %in% c(3, 4))) |> # remove Noncensus water and Census water
  filter(!(DSTRBCD1 %in% disturb_codes),
         !(DSTRBCD2 %in% disturb_codes),
         !(DSTRBCD3 %in% disturb_codes)) |>
  filter(!(TRTCD1 %in% treat_codes),
         !(TRTCD2 %in% treat_codes),
         !(TRTCD3 %in% treat_codes)) |>
  filter(!(PRESNFCD %in% nonforest_codes)) |>
  inner_join(plot_tbl[, c("STATECD","UNITCD","COUNTYCD","PLOT","INVYR","LAT","LON")],
            by = c("STATECD","UNITCD","COUNTYCD","PLOT","INVYR")) |>
  left_join(state_codes, by = "STATECD") |>
  filter(!(STATEABB %in% drop_states))

# condition-level location + attributes carried onto each plot-level table.
# PCTBARE_RMRS (RMRS interior West bare-ground estimate) rides along here.
cond_loc <- cond |>
  select(PLT_CN, INVYR, STATECD, UNITCD, COUNTYCD, PLOT, CONDID,
         PCTBARE_RMRS, SLOPE, ASPECT, STATENAME, LAT, LON)

# join keys shared by the plot-level tables and cond_loc
plot_keys <- c("PLT_CN","INVYR","STATECD","UNITCD","COUNTYCD","PLOT","CONDID")

# --- Vegetation composition: aerial cover by functional group ---------------
# LAYER 5 = aerial cover (summed across vegetation layers). Sum cover within
# each subplot by growth-habit group, then average across subplots in a plot.
veg_group_keys <- plot_keys

veg_composition <- read.csv(file.path(raw_dir, "ENTIRE_P2VEG_SUBP_STRUCTURE.csv")) |>
  # filtering to get rid of some old protocol codes that may not
  # have split aout shrubs etc. 
  filter(LAYER == 5, GROWTH_HABIT_CD %in% c("FB","GR","SH","TT","NT")) |>
  inner_join(cond_loc, by = plot_keys) |>  # only keeping locations with good condition
  group_by(across(all_of(c(veg_group_keys, "SUBP", "GROWTH_HABIT_CD")))) |> 
  summarise(cover = sum(COVER_PCT, na.rm = TRUE), .groups = "drop") |>   # within subplot
  group_by(across(all_of(c(veg_group_keys, "GROWTH_HABIT_CD")))) |>
  summarise(cover = mean(cover, na.rm = TRUE), .groups = "drop") |>      # across subplots
  pivot_wider(names_from = GROWTH_HABIT_CD, values_from = cover) |>
  # FIA growth-habit codes: FB forb, GR graminoid, SH shrub,
  # TT tally tree, NT non-tally tree
  rename(Forbs_AerialCover        = FB,
         Graminoid_AerialCover    = GR,
         Shrub_AerialCover        = SH,
         TallyTree_AerialCover    = TT,
         NonTallyTree_AerialCover = NT)

write_csv(veg_composition,
          file.path(out_dir, paste0("vegetationComposition", suffix, ".csv")))

# --- Ground cover: bare ground, averaged across FIA's two sources -----------
# FIA measures bare ground two ways, in largely non-overlapping regions:
#   PCTBARE_RMRS  - COND attribute, RMRS interior West (already on cond_loc)
#   GRND_CVR BARE - ground-cover transect, Pacific states (OR/CA/WA)
# missing_flag -> NA in both, then average per plot (na.rm = TRUE coalesces
# where only one source exists, means where both do).

# STOP--this is conceptually wrong (fix)

# in this dataset 'Ground cover items must be in contact with the ground'
bare_transect <- read_csv(file.path(raw_dir, "ENTIRE_GRND_CVR.csv")) |>
  pivot_wider(names_from = GRND_CVR_TYP, values_from = CVR_PCT) |>
  group_by(PLT_CN, INVYR, STATECD, UNITCD, COUNTYCD, PLOT) |>
  summarise(BareGround_transect = mean(BARE, na.rm = TRUE), .groups = "drop") |>
  mutate(BareGround_transect = na_if(BareGround_transect, NaN),
         BareGround_transect = na_if(BareGround_transect, missing_flag))

ground_cover <- cond_loc |>
  mutate(PCTBARE_RMRS = na_if(PCTBARE_RMRS, missing_flag)) |>
  left_join(bare_transect,
            by = c("PLT_CN","INVYR","STATECD","UNITCD","COUNTYCD","PLOT")) |>
  mutate(BareGroundCover = rowMeans(cbind(PCTBARE_RMRS, BareGround_transect),
                                    na.rm = TRUE),
         BareGroundCover = na_if(BareGroundCover, NaN)) |>
  select(PLT_CN, INVYR, STATECD, UNITCD, COUNTYCD, PLOT, CONDID,
         SLOPE, ASPECT, STATENAME, LAT, LON, BareGroundCover) 

stopifnot(all(!is.na(ground_cover$LAT)))

write_csv(ground_cover,
          file.path(out_dir, paste0("groundCover", suffix, ".csv")))

# --- Species -> taxonomic group lookup --------------------------------------
species_groups <- FIESTA::ref_species |>
  left_join(
    {
      tax <- taxonlookup::lookup_table(
        unique(FIESTA::ref_species$SCIENTIFIC_NAME), by_species = TRUE)
      tax$SCIENTIFIC_NAME <- rownames(tax)
      tax
    },
    by = "SCIENTIFIC_NAME"
  ) |>
  # manual group assignments for species taxonlookup missed (kept from original)
  mutate(group = case_when(
    GENUS %in% c("Reynoldsia","Feijoa","Exorrhiza","Gulubia","Trukia",
                 "Neolaugeria","Nesoluma","Hyeronima","Guamia","Carmona",
                 "Munroidendron")                            ~ "Angiosperms",
    GENUS %in% c("Acoelorraphe","Family Arecaceae","Howeia") ~ "Pteridophytes",
    GENUS == "Cupressocyparis"                               ~ "Gymnosperms",
    SCIENTIFIC_NAME == "Tree broadleaf"                      ~ "Angiosperms",
    SCIENTIFIC_NAME == "Tree evergreen"                      ~ "Gymnosperms",
    TRUE                                                     ~ group
  ))

# --- TREE: basal area by taxonomic group ------------------------------------
# Pull the TREE table in non-overlapping 5-state chunks (original hardcoded
# ranges that overlapped -> some states double-counted).
state_ids  <- unique(cond$STATECD)
state_chunks <- split(state_ids, ceiling(seq_along(state_ids) / 5))

tree_basal_area <- map_dfr(state_chunks,
                           ~ FIESTA::DBgetCSV("TREE", states = .x)) |>
  left_join(species_groups[, c("SPCD","SCIENTIFIC_NAME","family","group")],
            by = "SPCD") |>
  # generic broadleaf / evergreen / unknown rows get a group directly
  mutate(group = case_when(
    SCIENTIFIC_NAME == "Tree broadleaf" ~ "Angiosperms",
    SCIENTIFIC_NAME == "Tree evergreen" ~ "Gymnosperms",
    SCIENTIFIC_NAME == "Tree unknown"   ~ "Unknown",
    TRUE ~ group
  )) |>
  filter(STATUSCD == 1) |>                    # live trees only, before aggregation
  mutate(basalArea_in2 = pi * (DIA / 2)^2)     # basal area per tree (sq in)

tree_keys <- c("PLT_CN","INVYR","STATECD","UNITCD","COUNTYCD","PLOT","CONDID","group")

tree_use <- tree_basal_area |>
  group_by(across(all_of(c(tree_keys, "SUBP")))) |>
  summarise(basalAreaSum_in2 = sum(basalArea_in2, na.rm = TRUE), .groups = "drop") |>  # within subplot
  group_by(across(all_of(tree_keys))) |>
  summarise(basalArea_in2 = mean(basalAreaSum_in2, na.rm = TRUE), .groups = "drop") |> # across subplots
  mutate(PLT_CN = as.double(PLT_CN)) |>
  pivot_wider(names_from = group, values_from = basalArea_in2,
              names_glue = "basalArea_{group}_in2") |>
  mutate(across(starts_with("basalArea_"), ~ replace_na(.x, 0))) |>  # absent group = true zero
  left_join(cond_loc, by = plot_keys) |>
  filter(!is.na(LAT))

# total basal area + per-group proportions (%)
ba_cols <- names(tree_use)[str_starts(names(tree_use), "basalArea_") &
                             str_ends(names(tree_use), "_in2")]
tree_use <- tree_use |>
  mutate(basalArea_allGroups_in2 = rowSums(across(all_of(ba_cols))),
         across(all_of(ba_cols),
                ~ .x / basalArea_allGroups_in2 * 100,
                .names = "{str_replace(.col, '_in2', '_perc')}"))

write_csv(tree_use, file.path(out_dir, 
                              paste0("TREEtable_", suffix, ".csv")))
