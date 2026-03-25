# various functions for preparing data for modeling and other tasks


# prep dataframe ----------------------------------------------------------

# for preparing version 1 of the biomass dataset, this 
# isn't the final version of the dataset to use
# but more for development/proof of concept
# trim_tree_cov--provide a proportional cover of trees below it's considered
# 0 (i.e. none zero tree predictions are leaking into areas
# where there probably aren't trees)
prepare_d01 <- function(data, cover_suffix, pfts, trim_tree_cov = NULL) {
  dat2 <- data %>% 
    # for consistent downstream naming
    rename_with(.fn = \(x) str_replace(x, cover_suffix, 'Cov'),
                .cols = ends_with(cover_suffix)) %>% 
    filter(!is.na(biomass_MgPerHect)) %>% 
    rename(totalBio = biomass_MgPerHect)
  
  # set.seed(12)
  # dat2 <- sample_n(dat2, size = n_sample)
  if(!is.null(trim_tree_cov)) {

    tree_cols <- c("totalTreeCov", "needleLeavedTreeCov", "broadLeavedTreeCov")
    stopifnot(trim_tree_cov >=0,
              trim_tree_cov <= 1,
              tree_cols %in% names(dat2))
    for(col in tree_cols) {
      # make small "noise' values 0
      dat2[[col]][dat2[[col]] < trim_tree_cov] <- 0
    }
    
  }
  
  dat3 <- dat2
  # make cover sum to 1 -----------------------------------------------------
  
  cov_cols <- c(paste0(pfts, "Cov"), "bareGroundCov")
  

  
  dat4 <- dat3 %>%
    mutate(.cov_total = rowSums(across(all_of(cov_cols)), na.rm = TRUE)) %>%
    mutate(across(all_of(cov_cols), ~ .x / .cov_total)) %>%
    select(-.cov_total)
  
  
  # normalize potential predictor variables ---------------------------------
  
  # normalize predictors (all potential predictors, not just those in the model)
  # doing this before sim_bio so all data are on standardized scale
  
  vars_to_scale <- names(dat4)[!names(dat4) %in% c(
    # skip non-numeric / identifier columns
    "biomassSource", "biomassType", "year", 
    "burnedMoreThan20YearsAgo", "region",
    # skip cover columns (these are proportions, not climate)
    paste0(pfts, "Cov"), "bareGroundCov",
    "totalTreeCov", "totalHerbaceousCov",
    # skip response
    "totalBio",
    # skip coordinates (keep raw for spatial analyses)
    "x", "y"
  )]
  
  scale_df <- tibble::tibble(
    variable = vars_to_scale,
    mean = sapply(dat4[vars_to_scale], mean, na.rm = TRUE),
    sd = sapply(dat4[vars_to_scale], sd, na.rm = TRUE)
  )
  
  dat4 <- dat4 |>
    mutate(across(all_of(vars_to_scale), ~ as.numeric(scale(.x))))
  
  # add artificial regions --------------------------------------------------
  # based on climate clustering
  
  region_vars <- c("tmean_CLIM", 
                   "precip_CLIM", 
                   "x",
                   "y")
  
  reg <- make_region_kmeans(
    dat = dat4,
    vars = region_vars,
    nstart = 1,
    k = 20,
    seed = 42,
    max_iter = 10000
  )
  
  dat4$region <- reg$data
  
  # 'cleaning' 

  stopifnot(mean(dat4$totalBio == 0) < 0.01) # if there are many zeros that's suspect
  
  # making it a very small number so log is possible
  dat4$totalBio[dat4$totalBio == 0] <- min(dat4$totalBio[dat4$totalBio > 0])/2 
  
  list(data = dat4,
       scale = scale_df)
}
