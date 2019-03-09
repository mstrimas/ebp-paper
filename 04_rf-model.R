# random forest model for encounter rate

library(auk)
library(sf)
library(raster)
library(ranger)
library(scam) 
library(dggridR)
library(PresenceAbsence)
library(verification)
library(edarf)
library(viridis)
library(fields)
library(dplyr)
library(stringr)
library(purrr)
library(tidyr)
library(tibble)
library(readr)
library(lubridate)
walk(list.files("R", full.names = TRUE), source)
# resolve namespace conflicts
select <- dplyr::select
projection <- raster::projection

set.seed(1)
# set species for analysis
species <- "Wood Thrush"
sp_code <- ebird_species(species, "code")
# setup spatial sampling regime
sample_regime <- "both"
sample_spacing <- 5


# load data ----

# ebird data
ebird <- read_csv("data/ebd_june_bcr27_zf.csv", na = "") %>% 
  filter(species_code == sp_code)

# modis covariates
habitat <- read_csv("data/modis_pland_checklists.csv", 
                    col_types = cols(
                      .default = col_double(),
                      checklist_id = col_character()))
pred_surface <- read_csv("data/modis_pland_prediction-surface.csv", 
                         col_types = cols(
                           .default = col_double(),
                           id = col_integer(),
                           year = col_integer()))

# combine modis and ebird data
ebird_habitat <- inner_join(ebird, habitat, by = "checklist_id")

# optional checklist calibration index
cci_file <- "data/cci_june_bcr27.csv"
if (file.exists(cci_file)) {
  cci <- read_csv(cci_file)
  ebird_habitat <- inner_join(ebird_habitat, cci, by = "checklist_id") %>% 
    filter(!is.na(checklist_calibration_index))
}


# map data ----

map_proj <- st_crs(102003)
# borders
f_gpkg <- "data/gis-data.gpkg"
ne_land <- read_sf(f_gpkg, "ne_land") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
ne_country_lines <- read_sf(f_gpkg, "ne_country_lines") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
ne_state_lines <- read_sf(f_gpkg, "ne_state_lines") %>% 
  filter(country_code == "US") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
bcr <- read_sf(f_gpkg, "bcr") %>% 
  filter(bcr_code == 27) %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()


# spatial subsampling ----

ebird_ss <- hex_sample(ebird_habitat, spacing = sample_spacing,
                       regime = sample_regime)


# test/train split ----

# prepare data
ebird_split <- ebird_habitat %>% 
  # random forests requires an integer repsonse
  mutate(species_observed = as.integer(species_observed)) %>% 
  # select only the columns to be used in the model
  select(species_observed,
         observation_date, 
         time_observations_started, duration_minutes,
         effort_distance_km, number_observers, 
         contains("checklist_calibration_index"),
         starts_with("pland_")) %>% 
  # test/train split, 80/20
  split(if_else(runif(nrow(.)) <= 0.8, "train", "test"))
# ensure there are no rows with NAs
stopifnot(all(complete.cases(ebird_split$test)))
stopifnot(all(complete.cases(ebird_split$train)))

# fit rf model ----

rf <- ranger(formula =  species_observed ~ ., 
             num.trees = 1000, mtry = sqrt(ncol(ebird_split$train) - 1),
             importance = "impurity",
             data = ebird_split$train)


# model calibration ----

# predict on training partition
rf_pred_train <- predict(rf, data = ebird_split$train, type = "response")
rf_pred_train <- data.frame(id = nrow(ebird_split$train),
                            obs = ebird_split$train$species_observed,
                            pred = rf_pred_train$predictions) %>% 
  drop_na()
# calibration model
calibration_model <- scam(obs ~ s(pred, k = 5, bs = "mpi"), 
                          data = rf_pred_train, 
                          family = binomial)


# validation ----

# predict on test partition
rf_pred_test <- predict(rf, data = ebird_split$test, type = "response")
rf_pred_test <- data.frame(id = nrow(ebird_split$test),
                           obs = ebird_split$test$species_observed,
                           pred = rf_pred_test$predictions) %>% 
  drop_na()

# mean squared error (mse)
mse <- mean((rf_pred_test$obs - rf_pred_test$pred)^2, na.rm = TRUE)
# pick threshold to maximize kappa
opt_thresh <- optimal.thresholds(rf_pred_test, opt.methods = "MaxKappa")[1, 2]
# calculate accuracy metrics: auc, kappa, sensitivity, specificity, brier
pa_metrics <- presence.absence.accuracy(rf_pred_test, threshold = opt_thresh, 
                                        na.rm = TRUE, st.dev = FALSE)
# choose the threshold that maximizes kappa
rf_assessment <- c(
  mse = mse,
  auc = pa_metrics$AUC,
  brier = brier(rf_pred_test$obs, rf_pred_test$pred)$bs,
  # threshold required
  thresh = opt_thresh,
  kappa = pa_metrics$Kappa,
  sensitivity = pa_metrics$sensitivity,
  specificity = pa_metrics$specificity
) %>% 
  enframe("metric")
str_glue("output/04_rf-model_assessment_{sp_code}.csv") %>% 
  write_csv(rf_assessment, .)


# predict ----

# find optimal time of day to maximimise detectability
# constrain to times with at least 1% of checklists
pred_t <- partial_dependence(rf, vars = "time_observations_started",
                             n = c(24 * 6, nrow(ebird_split$train)),
                             data = ebird_split$train)

# hours with at least 1% of checklists
search_hours <- ebird_split$train %>%
  mutate(hour = floor(time_observations_started)) %>%
  count(hour) %>%
  mutate(pct = n / sum(n)) %>%
  filter(pct >= 0.01)

# constrained peak time
t_max <- pred_t %>%
  filter(floor(time_observations_started) %in% search_hours$hour) %>%
  top_n(1, wt = species_observed) %>%
  pull(time_observations_started)

# add effort covariates to prediction 
pred_surface_eff <- pred_surface %>% 
  mutate(observation_date = ymd("2016-06-15"),
         time_observations_started = t_max,
         duration_minutes = 60,
         effort_distance_km = 1,
         number_observers = 1,
         checklist_calibration_index = 2)
# predict
rf_pred <- predict(rf, data = pred_surface_eff, type = "response")
# calibrate
rf_pred_calibrate <- predict(calibration_model, 
                             data.frame(pred = rf_pred$predictions), 
                             type = "response")
rf_pred <- bind_cols(pred_surface_eff, prob = rf_pred_calibrate)

# rasterize
r <- raster("data/modis_5xagg.tif")
r_pred <- rf_pred %>% 
  select(prob, latitude, longitude) %>% 
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(crs = projection(r)) %>% 
  rasterize(r)
r_pred <- r_pred[[-1]]
str_glue("output/04_rf-model_predictions_{sp_code}.tif") %>% 
  writeRaster(r_pred, ., overwrite = TRUE)


# map predictions ----

r_pred_proj <- projectRaster(r_pred, crs = map_proj$proj4string, method = "ngb")

str_glue("figures/04_rf-model_predictions_{sp_code}.png") %>% 
  png(width = 2400, height = 1800, res = 300)
par(mar = c(4, 0.5, 0.5, 0.5))

# set up plot area
plot(bcr, col = NA, border = NA)
plot(ne_land, col = "#dddddd", border = "#888888", lwd = 0.5, add = TRUE)


# modified plasma palette
plasma_rev <- rev(plasma(25, end = 0.9))
gray_int <- colorRampPalette(c("#dddddd", plasma_rev[1]))
pal <- c(gray_int(4)[2], plasma_rev)

# probability of detection
mx <- ceiling(100 * max(r_pred_proj[], na.rm = TRUE)) / 100
brks <- seq(0, mx, length.out = length(pal) + 1)
plot(r_pred_proj, col = pal, breaks = brks, maxpixels = ncell(r_pred_proj),
     legend = FALSE, add = TRUE)

# borders
plot(bcr, col = NA, border = "#000000", lwd = 1, add = TRUE)
plot(ne_state_lines, col = "#ffffff", lwd = 0.75, add = TRUE)
plot(ne_country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
box()

# legend
par(new = TRUE, mar = c(0, 0, 0, 0))
title <- str_glue("{species} Encounter Rate")
lbl_brks <- seq(0, mx, by = 0.1)
image.plot(zlim = range(brks), legend.only = TRUE, col = pal,
           smallplot = c(0.25, 0.75, 0.05, 0.08),
           horizontal = TRUE,
           axis.args = list(at = lbl_brks, labels = lbl_brks,
                            fg = "black", col.axis = "black",
                            cex.axis = 0.5, lwd.ticks = 0.5,
                            padj = -2.5),
           legend.args = list(text = title,
                              side = 3, col = "black",
                              cex = 1, line = 0))
dev.off()