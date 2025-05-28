# ============================================================================
#  HARD-CLAM SCLEROCHRONOLOGY ‒ MASTER ANALYSIS SCRIPT           (2025-05-28)
#  --------------------------------------------------------------------------
#  Reproduces every data product and statistical result reported in:
#      “Regional spring–summer SST, not seagrass recovery, governs hard-clam
#       growth” — Limnology & Oceanography Letters (in review)

# Repo   : https://github.com/mrcornish/VCR_hard_clam_chron

#     Authors: Michael R. Cornish1* 
#         and Max C.N. Castorani1

#Affiliations:
# Department of Environmental Sciences, University of Virginia, 
#  Charlottesville, Virginia, USA. 
# *Corresponding author: mcornish@virginia.edu 	


#Author Contribution Statement: This research was supported by Virginia 
#Sea Grant through a research fellowship to MRC and the U.S. National 
#Science Foundation through sustained support of the Virginia Coast Reserve 
#Long-Term Ecological Research project (award no. 1832221). 
#MRC conducted the field and laboratory activities, analyzed the data, 
#and wrote the manuscript, with guidance from MCNC. Both authors contributed 
#to designing the study, interpreting the results, and revising the 
#manuscript. Both authors gave final approval for publication and declare 
#they have no competing interests.  







#  ── Inputs ────────────────────────────────────────────────────────────────
#  • data/
#      ├─ VIMS_DAWT.csv           : hind-cast daily SST 1963-1981
#      ├─ vcr_envirn_data.csv     : in-situ logger SST 1982-2020
#      ├─ site_chrons3.csv        : shell-level SGI (cross-dated, detrended)
#      ├─ master_chron_neg_exp.csv (Population-level chrinology, used only for window search)
#      └─ seagrass/ <*.shp>       : Annual VIMS seagrass polygons 2001-2020 
#
#  ── Outputs ───────────────────────────────────────────────────────────────
#  • site_df.rds         : row-level DF with SGI, SST windows, distance
#  • site_means.csv      : site-level means for SAR analysis
#  • models_tbl.csv/html : model-comparison table S2
#  • best_mixed.rds      : selected LME4 fit
#  • best_sar.rds        : selected SAR-error fit
#
#  HOW TO REPRODUCE
#  --------------------------------------------------------------------------
#  1. Clone repository  ▸ renv::restore()   # recreates package versions
#  2. source("clam_sclero_master_script.R") # <-- this file
#  3. Inspect the console for ✔ success messages.
#

## ───────────────────────── 00) SET-UP ───────────────────────────────────────
need <- c("tidyverse",  "lubridate",   "janitor",  "dplR",      "zoo",
          "sf",         "spdep",       "spatialreg",
          "lme4",       "lmerTest",    "broom.mixed",
          "gt",         "xtable",      "here")

new <- setdiff(need, rownames(installed.packages()))
if (length(new)) install.packages(new, quiet = TRUE)

suppressPackageStartupMessages(lapply(need, library, character.only = TRUE))

# Ensure dplyr verbs outrank generic ones from sf/raster
conflicts_prefer(dplyr::select, dplyr::filter, dplyr::rename,
                 dplyr::mutate, dplyr::summarise, dplyr::arrange)

DIR_DATA     <- here("data")
FILE_SITES   <- file.path(DIR_DATA, "site_chrons3.csv")
FILE_MASTER  <- file.path(DIR_DATA, "master_chron_neg_exp.csv")
DIR_SEAGRASS <- file.path(DIR_DATA, "seagrass")

CRS_UTM  <- 32618
BBOX_POLY <- st_as_sfc(st_bbox(c(xmin = 421575, xmax = 429970.3,
                                 ymin = 4116561, ymax = 4127426),
                               crs = CRS_UTM))

# ----- helper functions -----------------------------------------------------
se <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

safe_arima <- function(x){
  if (length(na.omit(x)) < 3) return(rep(NA_real_, length(x)))
  fit <- try(arima(x, order = c(1, 0, 0)), silent = TRUE)
  if (inherits(fit, "try-error")) return(rep(NA_real_, length(x)))
  residuals(fit)
}

run_pw <- function(df){
  # pre-whiten SGI & SST (AR1) before correlating
  if (!all(c("anom", "SGI") %in% names(df)) || nrow(df) < 5) return(NULL)
  df <- arrange(df, Year)
  list(x = safe_arima(df$anom), y = safe_arima(df$SGI))
}

message("✔ libraries loaded & helpers ready")



## ───────────────────────── 01) MONTHLY SST 1963-2020 ───────────────────────
DAWT_raw <- read_csv(file.path(DIR_DATA, "VIMS_DAWT.csv"),
                     show_col_types = FALSE) %>%
  rename(Year = Y, Month = M, Day = D) %>%
  mutate(across(DAWT, as.numeric))

pre1982 <- DAWT_raw %>%
  filter(between(Year, 1963, 1981)) %>%
  group_by(Year, Month) %>%
  summarise(bay_temp_mean = mean(DAWT, na.rm = TRUE), .groups = "drop")

logger <- read_csv(file.path(DIR_DATA, "vcr_envirn_data.csv"),
                   show_col_types = FALSE) %>%
  clean_names() %>%
  rename(bay_temp = cbw_meas) %>%
  mutate(date  = parse_date_time(date, orders = c("ymd", "dmy", "mdy")),
         Year  = year(date),
         Month = month(date)) %>%
  filter(between(Year, 1982, 2020)) %>%
  group_by(Year, Month) %>%
  summarise(bay_temp_mean = mean(bay_temp, na.rm = TRUE), .groups = "drop")

monthly_sst <- bind_rows(pre1982, logger) %>% arrange(Year, Month)

message("✔ monthly SST assembled (1963-2020)")



## ───────────────────────── 02) CLIMATE-WINDOW SEARCH ───────────────────────
#  (Find the SST month set most strongly correlated with the master SGI)

window_sweep <- function(mons){
  clim <- monthly_sst %>% filter(Month %in% mons) %>%
    group_by(Month) %>%
    summarise(clim = mean(bay_temp_mean), .groups = "drop")
  
  monthly_sst %>% filter(Month %in% mons) %>%
    left_join(clim, by = "Month") %>%
    mutate(anom = bay_temp_mean - clim) %>%
    group_by(Year) %>%
    summarise(anom = mean(anom), .groups = "drop")
}

if (file.exists(FILE_MASTER)) {
  sgimaster <- read_csv(FILE_MASTER, show_col_types = FALSE)
  
  win_grid <- tidyr::expand_grid(start = 1:12, len = 3:6) %>%
    mutate(months = purrr::map2(start, len,
                                ~ ((.x - 1 + 0:(.y - 1)) %% 12) + 1))
  
  best <- win_grid %>%
    mutate(r = purrr::map_dbl(months, ~{
      ser <- window_sweep(.x) %>% inner_join(sgimaster, by = "Year")
      pw  <- run_pw(ser)
      if (is.null(pw)) return(NA_real_)
      cor(pw$x, pw$y, use = "complete.obs")
    })) %>%
    filter(!is.na(r)) %>%
    slice_max(abs(r), n = 1)
  
  WIN_BEST_MONTHS <- best$months[[1]]
  BEST_LABEL      <- paste(month.abb[WIN_BEST_MONTHS], collapse = " ")
  message(glue::glue("✔ best window = {BEST_LABEL} ({length(WIN_BEST_MONTHS)} mo)"))
} else {
  WIN_BEST_MONTHS <- 3:5
  BEST_LABEL      <- "Mar Apr May"
  warning("master chronology missing – default Mar-May window used")
}

get_window <- function(mons, lbl){
  window_sweep(mons) %>%
    mutate(window = lbl,
           anom_c = scale(anom, center = TRUE, scale = FALSE)[, 1])
}

clim_tbl <- bind_rows(
  get_window(2:6,              "Feb–Jun (5 mo)"),
  get_window(WIN_BEST_MONTHS,  glue::glue("{BEST_LABEL} ({length(WIN_BEST_MONTHS)} mo)"))
)



## ───────────────────────── 03) SGI + DISTANCE-TO-EDGE ──────────────────────
#  (Add Euclidean distance of each clam site to meadow edge for each year)

site_sgi <- read_csv(FILE_SITES, show_col_types = FALSE) %>% 
  select(-starts_with("..."))  # drop unnamed cols from Excel

yrs_layers <- c("2020","2019","2018","2017","2015","2013",
                "2011","2010","2009","2008","2007","2003","2001")

read_bed <- function(layer){
  st_read(DIR_SEAGRASS, layer = layer, quiet = TRUE) %>%
    st_transform(CRS_UTM) %>%
    st_make_valid() %>%
    st_crop(BBOX_POLY) %>%
    filter(st_dimension(geometry) == 2) %>%
    mutate(Year = as.integer(layer)) %>%
    select(Year, geometry)
}

sg_all <- purrr::map_dfr(yrs_layers, read_bed)

site_pts <- site_sgi %>%
  st_as_sf(coords = c("Long", "Lat"), crs = 4326, remove = FALSE) %>%
  st_transform(CRS_UTM)

edge_dist <- function(yr, pts){
  beds <- sg_all %>% filter(Year == yr)
  if (!nrow(beds))
    return(tibble(Site = pts$Site, Year = yr, dist_edge = NA_real_))
  
  border <- st_boundary(st_union(beds))
  d      <- st_distance(pts, border)
  inside <- st_within(pts, st_union(beds), sparse = FALSE)[, 1]
  
  tibble(Site = pts$Site, Year = yr,
         dist_edge = ifelse(inside, -as.numeric(d), as.numeric(d)))
}

dist_df <- purrr::map_dfr(unique(site_pts$Year),
                          ~ edge_dist(.x, filter(site_pts, Year == .x)))

site_sgi <- site_sgi %>%
  left_join(dist_df, by = c("Site", "Year")) %>%
  mutate(dist_edge = ifelse(Year < 2001, 0, dist_edge)) %>%   # pre-restoration yrs
  arrange(Site, Year) %>%
  group_by(Site) %>%
  mutate(dist_edge = zoo::na.approx(dist_edge, x = Year,
                                    na.rm = FALSE, rule = 2)) %>%
  ungroup()



## ───────────────────────── 04) BUILD & SAVE ANALYSIS FRAME ─────────────────
site_df <- tidyr::expand_grid(site_sgi,
                              window = unique(clim_tbl$window)) %>%
  left_join(clim_tbl, by = c("Year", "window")) %>%
  tidyr::drop_na(SGI, anom_c, dist_edge) %>%
  mutate(dist_edge      = as.numeric(dist_edge),
         SST_def_anom   = if_else(window == "Feb–Jun (5 mo)", anom_c, NA_real_),
         SST_best_anom  = if_else(window != "Feb–Jun (5 mo)", anom_c, NA_real_))

saveRDS(site_df, "site_df.rds")

site_means <- site_df %>%
  group_by(Site, Long, Lat) %>%
  summarise(
    SGI           = mean(SGI),
    dist_edge     = mean(dist_edge),
    SST_def_anom  = mean(SST_def_anom,  na.rm = TRUE),
    SST_best_anom = mean(SST_best_anom, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(site_means, "site_means.csv")
message("✔ site_df.rds and site_means.csv written")



## ───────────────────────── 05) MIXED-EFFECTS MODELS ────────────────────────
form_list <- list(
  M1 = SGI ~ SST_def_anom  + dist_edge + (1              | Site),
  M2 = SGI ~ SST_def_anom  + dist_edge + (1 + SST_def_anom  | Site),
  M3 = SGI ~ SST_best_anom + dist_edge + (1              | Site),
  M4 = SGI ~ SST_best_anom + dist_edge + (1 + SST_best_anom | Site)
)

mix_mods <- purrr::map(form_list, \(f)
                       lmer(f, data = site_df, REML = FALSE,
                            control = lmerControl(optimizer = "bobyqa",
                                                  optCtrl   = list(maxfun = 2e5)))
)

message("✔ mixed-effects candidate models fitted")



## ───────────────────────── 06) SAR ERROR-MODELS ────────────────────────────
#   Build k-nearest-neighbour spatial weights (min(4, n-1))
k_nn   <- max(1, min(4, nrow(site_means) - 1))
coords <- as.matrix(site_means[, c("Long", "Lat")])
nb     <- knn2nb(knearneigh(coords, k = k_nn))
lw     <- nb2listw(nb, style = "W")

sar_mods <- list(
  SAR_def  = errorsarlm(SGI ~ SST_def_anom  + dist_edge,
                        data = site_means, listw = lw),
  SAR_best = errorsarlm(SGI ~ SST_best_anom + dist_edge,
                        data = site_means, listw = lw)
)

message("✔ SAR error-lag models fitted")



## ───────────────────────── 07) MODEL-COMPARISON TABLE ──────────────────────
get_mixed_row <- function(mod, label){
  broom.mixed::tidy(mod, effects = "fixed") %>%
    mutate(Model = label,
           rho   = NA_real_,
           AIC   = AIC(mod)) %>%
    rename(Beta = estimate, SE = std.error, `t/z` = statistic, p = p.value) %>%
    select(Model, term, Beta, SE, `t/z`, p, rho, AIC)
}

get_sar_row <- function(mod, label){
  co <- summary(mod)$Coef
  tibble(
    Model = label,
    term  = rownames(co),
    Beta  = co[, "Estimate"],
    SE    = co[, "Std. Error"],
    `t/z` = co[, "z value"],
    p     = co[, "Pr(>|z|)"],
    rho   = summary(mod)$rho,
    AIC   = AIC(mod)
  )
}

tbl <- bind_rows(
  purrr::imap_dfr(mix_mods, get_mixed_row),
  purrr::imap_dfr(sar_mods,  get_sar_row)
) %>% arrange(AIC)

write_csv(tbl, "models_tbl.csv")

gt(tbl) %>%
  tab_header(title = md("**Table S2. Mixed-effects vs SAR model comparison**")) %>%
  fmt_number(columns = c(Beta, SE, `t/z`, p, rho, AIC), decimals = 3) %>%
  cols_label(p = md("*p*")) %>%
  gtsave("models_tbl.html")

message("✔ models_tbl.csv and .html exported")



## ───────────────────────── 08) SAVE “BEST” MODELS ──────────────────────────
best_mix_name <- tbl %>% filter(Model %in% names(mix_mods),
                                term == "(Intercept)") %>%
  slice_min(AIC) %>% pull(Model)

best_sar_name <- tbl %>% filter(Model %in% names(sar_mods),
                                term == "(Intercept)") %>%
  slice_min(AIC) %>% pull(Model)

saveRDS(mix_mods[[best_mix_name]], "best_mixed.rds")
saveRDS(sar_mods [[best_sar_name]], "best_sar.rds")

message(glue::glue("✔ best models saved (LME = {best_mix_name}; SAR = {best_sar_name})"))



## ───────────────────────── 09) REPRODUCIBILITY FOOTER ──────────────────────
cat("\n=========== SESSION INFO ===========\n")
print(sessionInfo(), locale = FALSE)
cat("\n====================================\n")

message("\n🎉  ALL DONE — analysis complete\n")
