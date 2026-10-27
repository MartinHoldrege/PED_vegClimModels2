# Purpose:
# Dot plots of projected climate change by GCM (MACAv2-METDATA). One PDF per
# climate variable, one page per scenario-period (RCP 4.5/8.5 x 2040-2069/
# 2070-2099). Each page has one panel for CONUS and one per National Climate
# Assessment region. Per GCM: median (point), interquartile range (line) and mean
# (x) of the area-weighted change from 1971-2000. GCMs are ordered by their
# CONUS median for that page, in every panel. Axis limits are shared within
# a page but differ between pages.
#
# Input: maca_delta_summaries.csv (02_summarise_maca_deltas.R)
# Needs `gcm_colors` in Functions/constants.R (printed by that script).
#
# October 2026

# dependencies ------------------------------------------------------------

source("Functions/init.R")
source_functions()
source('Functions/viz/fig_params.R')

# params ------------------------------------------------------------------

test_run <- FALSE # if TRUE, only the first variable

# every variable is plotted as absolute change; these are also plotted as
# percent change, on extra pages after the absolute ones
pct_vars <- pct_vars <- c("MAP", "P_wettestMonth", "P_driestMonth",
                          "WDD_mean", "WDD_p05",
                          "VPD_mean", "VPD_max", "VPD_min", "VPD_max_p95")

p_in <- file.path(paths$large, "Data_processed/WallToWallClimateData/maca",
                  "summaries", "maca_delta_summaries.csv")
out_dir <- "Figures/Climate/maca_gcm_dotplots"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fig_width <- 11
fig_height <- 8.5

# read --------------------------------------------------------------------

stopifnot(file.exists(p_in))
dat <- readr::read_csv(p_in, show_col_types = FALSE)

if (!exists("gcm_colors")) {
  stop("gcm_colors not found; paste the vector printed by ",
       "02_summarise_maca_deltas.R into Functions/viz/fig_params.R")
}

stopifnot(all(unique(dat$gcm) %in% names(gcm_colors)))

# absolute change for all variables, plus percent change for pct_vars
dat <- dat |>
  filter(delta_type == "abs" | variable %in% pct_vars) |>
  mutate(plot_type = delta_type) |>
  mutate(scen_label = paste0(toupper(sub("rcp", "RCP ", scenario)), ", ",
                             period),
         scen_label = sub("RCP (\\d)(\\d)", "RCP \\1.\\2", scen_label))

# CONUS first, then NCA regions by code (roughly west to east)
region_levels <- dat |>
  distinct(region, region_code) |>
  arrange(!is.na(region_code), region_code) |>
  pull(region)

vars <- sort(unique(dat$variable))
if (test_run) vars <- vars[1]

# plotting function -------------------------------------------------------

#' Dot plot of GCM changes for one variable and scenario-period
#'
#' @param df Rows of `dat` for one variable and one scenario-period.
#' @return ggplot object.
plot_gcm_deltas <- function(df) {
  gcm_order <- df |>
    filter(region == "CONUS") |>
    arrange(median) |>
    pull(gcm)
  df <- mutate(df, gcm = factor(gcm, levels = gcm_order),
               region = factor(region, levels = region_levels))
  
  type <- unique(df$plot_type)
  ylab <- paste0("Change in ", unique(df$variable),
                 if (type == "pct") " (%)" else "")
  
  ggplot(df, aes(x = gcm, color = gcm)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
    geom_linerange(aes(ymin = q25, ymax = q75)) +
    geom_point(aes(y = median), size = 1.5) +
    geom_point(aes(y = mean), shape = 4, size = 1.5) +
    scale_color_manual(values = gcm_colors, guide = "none") +
    facet_wrap(~region, labeller = label_wrap_gen(width = 25)) +
    labs(x = NULL, y = ylab,
         title = paste0(unique(df$variable), ": ", unique(df$scen_label)),
         subtitle = paste0("Change from 1971-2000. Point = median, line = ",
                           "interquartile range, x = mean (area-weighted). ",
                           "GCMs ordered by CONUS median.")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                     size = 6))
}

# figures -----------------------------------------------------------------

for (v in vars) {
  dv <- filter(dat, variable == v)
  # absolute-change pages first, then percent-change pages (if any)
  dv$page <- paste(dv$plot_type, dv$scen_label)
  page_order <- dv |>
    distinct(page, plot_type, scen_label) |>
    arrange(plot_type, scen_label) |>
    pull(page)
  pages <- split(dv, dv$page)[page_order]
  
  out_file <- file.path(out_dir, paste0("maca_gcm_delta_", v, ".pdf"))
  pdf(out_file, width = fig_width, height = fig_height)
  for (pg in pages) print(plot_gcm_deltas(pg))
  dev.off()
  message("Wrote ", out_file)
}