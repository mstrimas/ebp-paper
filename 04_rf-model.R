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

# resolve namespace conflicts
select <- dplyr::select
projection <- raster::projection

# read in functions
walk(list.files("R", full.names = TRUE), source)

set.seed(1)
# set species for analysis
species <- "Wood Thrush"
sp_code <- ebird_species(species, "code")

# setup spatial sampling regime
sample_regime <- "both"
sample_spacing <- 5


# load data ----

data_folder <- "data/"
data_tag <- "mayjune_201718_bcr27"

# ebird data
# ebird <- read_csv(paste0(data_folder, "ebd_", data_tag, "_zf.csv"), na = "") %>% 
#   filter(species_code == sp_code)

ebird <- read_csv(paste0(data_folder, "data_4_models_", data_tag, ".csv"), na = "") 
species_count <- ebird[,which(colnames(ebird)==sp_code)] %>%
                  as.matrix() %>% as.vector() %>% as.numeric()
species_binary <- ifelse(is.na(species_count), 1, ifelse(species_count==0, 0, 1))
ebird$species_observed <- species_binary
ebird <- ebird %>%
        select(checklist_id, sampling_event_identifier, species_observed, latitude, longitude,
                protocol_type, all_species_reported, observation_date, time_observations_started,
                duration_minutes, effort_distance_km, number_observers, type)

# modis covariates
habitat <- read_csv(paste0(data_folder, "modis_pland_checklists_", data_tag, ".csv"), 
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
cci_file <- paste0("data/cci_", data_tag, ".csv")
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


# spatial subsampling for train and test_2017 ----

set.seed(1)
ebird_habitat <- ebird_habitat %>%
        mutate(type_week = paste(type, lubridate::week(observation_date), sep="_"))
ebird_ss <- hex_sample(ebird_habitat, spacing = sample_spacing,
                       regime = sample_regime, byvar = "type_week") %>%
        select(checklist_id, sampling_event_identifier) %>%
        mutate(selected = 1) %>% 
        right_join(ebird_habitat) %>%
        select(-type_week) %>%
        mutate(selected = ifelse(is.na(selected), 0, 1)) %>%
        mutate(selected = ifelse(type == "test_bbs", 1, ifelse(type == "other", 0, selected))) %>%
        filter(selected == 1)


# prepare data for modelling
ebird_split <- ebird_ss %>% 
  # random forests requires an integer repsonse
  mutate(species_observed = as.integer(species_observed)) %>% 
  mutate(species_observed = as.factor(species_observed)) %>%
  # select only the columns to be used in the model
  select(species_observed,
         observation_date, 
         time_observations_started, duration_minutes,
         effort_distance_km, number_observers, 
         contains("checklist_calibration_index"),
         starts_with("pland_")) %>% 
  split(ebird_ss$type)
# ensure there are no rows with NAs
stopifnot(all(complete.cases(ebird_split$test_2017)))
stopifnot(all(complete.cases(ebird_split$test_bbs)))
stopifnot(all(complete.cases(ebird_split$train)))


# fit rf model ----
pos_prop <- mean(as.numeric(as.character(ebird_split$train$species_observed)))

set.seed(1)
rf <- ranger(formula =  species_observed ~ ., 
             data = ebird_split$train,
             num.trees = 1000, 
             importance = "impurity",
             probability = TRUE,
             mtry = sqrt(ncol(ebird_split$train) - 1),
             replace = TRUE,
             sample.fraction = c(pos_prop, pos_prop),
             keep.inbag = TRUE)


# variable importance

data.frame(var =names(rf$variable.importance), imp = as.numeric(rf$variable.importance)) %>%
            mutate(imp = imp/max(imp)) %>%
            mutate(imp = round(imp, 3)) %>% 
            arrange(-imp)

#                            var   imp
# 1             duration_minutes 1.000
# 2    time_observations_started 0.913
# 3           effort_distance_km 0.788
# 4             observation_date 0.725
# 5                     pland_08 0.580
# 6  checklist_calibration_index 0.510
# 7                     pland_09 0.398
# 8                     pland_04 0.377
# 9                     pland_13 0.359
# 10                    pland_05 0.299

# model calibration ----

# predict on training partition
logit <- function(p) {log(p/(1-p))}
pred_train <- predict(rf, data = ebird_split$train, type = "response")
pred_train <- data.frame(obs = as.numeric(as.character(ebird_split$train$species_observed)),
                            pred = as.numeric(as.vector(pred_train$predictions[,c("1")]))) %>% 
                  mutate(logit_pred = logit(pred)) %>%
  drop_na()
# calibration model\
calibration_model <- gam(obs ~ s(pred, k = 3), 
                          data = pred_train, 
                          family = binomial, gamma = 2)
pred_train$pred_calibrated <- predict(calibration_model, newdata=pred_train, type="response")

# plot the model calibration
plot(pred_train$pred, jitter(pred_train$obs) + 0.2 - 0.4*pred_train$obs )
sq <- seq(0.01, 0.99, by=0.01)
lines(sq,  
  predict(calibration_model, data.frame(pred = sq), type = "response"), 
  col="red", lwd=2)
abline(0, 1, col="grey80")


# # prevalence adjustment ----

# # step 1. estimate the optimal time of day when the species is observed



# # step 2. estimate the target prevalence
# # predict to all training data, but set standardised effort variables
# train_high_effort <- ebird_split$train %>%
#           mutate(effort_distance_km = 2, duration_minutes = 120, 
#               checklist_calibration_index = 1.85, time_observations_started = optimal_time,
#               number_observers = 1)


# validation ----

# predict on test partition 2017
pred_test_2017_raw <- predict(rf, data = ebird_split$test_2017, type = "response")
pred_test_2017_cal <- predict(calibration_model, data.frame(pred = pred_test_2017_raw$predictions[,"1"]), type="response")
# pred_test_2017_adj <- 
pred_test_2017 <- data.frame(id = 1:nrow(ebird_split$test_2017),
                           obs = as.numeric(as.character(ebird_split$test_2017$species_observed)),
                           pred = pred_test_2017_cal) %>% 
  drop_na()


# mean squared error (mse)
mse <- mean((pred_test_2017$obs - pred_test_2017$pred)^2, na.rm = TRUE)
# pick threshold to maximize kappa
opt_thresh <- optimal.thresholds(pred_test_2017, opt.methods = "MaxKappa")[1, 2]
# calculate accuracy metrics: auc, kappa, sensitivity, specificity, brier
pa_metrics <- presence.absence.accuracy(pred_test_2017, threshold = opt_thresh, 
                                        na.rm = TRUE, st.dev = FALSE)

# choose the threshold that maximizes kappa
rf_assessment <- c(
  mse = mse,
  auc = pa_metrics$AUC,
  brier = brier(pred_test_2017$obs, pred_test_2017$pred)$bs,
  # threshold required
  thresh = opt_thresh,
  kappa = pa_metrics$Kappa,
  sensitivity = pa_metrics$sensitivity,
  specificity = pa_metrics$specificity
) %>% 
  enframe("metric")
str_glue("output/04_rf-model_assessment_test_2017_{sp_code}.csv") %>% 
  write_csv(rf_assessment, .)



# predict on test partition bbs
pred_test_bbs_raw <- predict(rf, data = ebird_split$test_bbs, type = "response")
pred_test_bbs_cal <- predict(calibration_model, data.frame(pred = pred_test_bbs_raw$predictions[,"1"]), type="response")
# pred_test_bbs_adj <- 
pred_test_bbs <- data.frame(id = 1:nrow(ebird_split$test_bbs),
                           obs = as.numeric(as.character(ebird_split$test_bbs$species_observed)),
                           pred = pred_test_bbs_cal) %>% 
  drop_na()


# mean squared error (mse)
mse <- mean((pred_test_bbs$obs - pred_test_bbs$pred)^2, na.rm = TRUE)
# pick threshold to maximize kappa
opt_thresh <- optimal.thresholds(pred_test_bbs, opt.methods = "MaxKappa")[1, 2]
# calculate accuracy metrics: auc, kappa, sensitivity, specificity, brier
pa_metrics <- presence.absence.accuracy(pred_test_bbs, threshold = opt_thresh, 
                                        na.rm = TRUE, st.dev = FALSE)

# choose the threshold that maximizes kappa
rf_assessment <- c(
  mse = mse,
  auc = pa_metrics$AUC,
  brier = brier(pred_test_bbs$obs, pred_test_bbs$pred)$bs,
  # threshold required
  thresh = opt_thresh,
  kappa = pa_metrics$Kappa,
  sensitivity = pa_metrics$sensitivity,
  specificity = pa_metrics$specificity
) %>% 
  enframe("metric")
str_glue("output/04_rf-model_assessment_test_bbs_{sp_code}.csv") %>% 
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