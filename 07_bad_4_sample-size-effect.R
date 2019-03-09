# bad practice random forest model for encounter rate

library(auk)
library(sf)
library(raster)
library(dggridR)
library(ranger)
library(scam) 
library(PresenceAbsence)
library(verification)
library(edarf)
library(fields)
library(dplyr)
library(purrr)
library(tidyr)
library(readr)
library(stringr)
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
ebird <- read_csv("data/ebd_june_bcr27_bad_zf.csv", na = "") %>% 
  filter(species_code == sp_code) %>% 
  mutate(species_observed = as.integer(species_observed))

# modis covariates
habitat <- read_csv("data/modis_pland_checklists.csv", 
                    col_types = cols(
                      .default = col_double(),
                      checklist_id = col_character()))

# combine modis and ebird data
ebird_habitat <- inner_join(ebird, habitat, by = "checklist_id")

# optional checklist calibration index
cci_file <- "data/cci_june_bcr27.csv"
if (file.exists(cci_file)) {
  cci <- read_csv(cci_file)
  ebird_habitat <- left_join(ebird_habitat, cci, by = "checklist_id")
}


# test dataset ----

# standardized test dataset for all models
ebird_test <- ebird_habitat %>%
  filter(all_species_reported,
         protocol_type %in% c("Stationary", "Traveling"),
         effort_distance_km <= 5,
         duration_minutes <= 5 * 60,
         number_observers <= 10)

# spatial subsampling
test_data_ss <- hex_sample(ebird_test, spacing = sample_spacing,
                           regime = sample_regime)

# take a 20% testing dataset
test_data <- sample_frac(test_data_ss, 0.2)

# remove these checklists from training data
model_data <- filter(ebird_habitat, !checklist_id %in% test_data$checklist_id)


# setup bad practice combinations ----

n_runs <- 20
bp_runs <- cross_df(list(good = c(TRUE, FALSE), 
                         prop_data = seq(0.1, 1, by = 0.1),
                         run = seq.int(n_runs)))


# fit bad practice model ----

fit_bp_model <- function(good, prop_data,
                         data, test_data,
                         spacing, regime, ...) {
  if (good) {
    # complete checklists only
    data <- filter(data, as.logical(all_species_reported))
    
    # filter on effort covariates
    data <- data %>%
      filter(protocol_type %in% c("Stationary", "Traveling"),
             effort_distance_km <= 5,
             duration_minutes <= 5 * 60,
             number_observers <= 10)
    
    # spatial subsample
    data <- hex_sample(data, spacing = spacing, regime = regime)
    
    # select covariates
    data <- data %>% 
      mutate(protocol_type = factor(protocol_type, 
                                    levels = c("Stationary" , 
                                               "Traveling"))) %>% 
      select(species_observed,
             observation_date, 
             time_observations_started, duration_minutes,
             effort_distance_km, number_observers, protocol_type,
             contains("checklist_calibration_index"),
             starts_with("pland_"))
  } else {
    data <- select(data, species_observed, starts_with("pland_"))
  }
  stopifnot(all(complete.cases(data)))
  
  # subsample some data
  data <- sample_frac(data, size = prop_data)
  
  # fit model
  mod <- ranger(formula =  species_observed ~ ., 
                num.trees = 1000, mtry = sqrt(ncol(data) - 1),
                importance = "impurity",
                data = data)
  
  # calibrate
  mod_pred_train <- predict(mod, data = data, type = "response")
  mod_pred_train <- data.frame(id = nrow(data),
                               obs = data$species_observed,
                               pred = mod_pred_train$predictions) %>% 
    drop_na()
  mod_cal <- scam(obs ~ s(pred, k = 5, bs = "mpi"), 
                  data = mod_pred_train, 
                  family = binomial)
  
  # validation
  # predict on test data set
  pred_rf <- predict(mod, data = test_data, type = "response")
  pred_cal <- predict(mod_cal, 
                      newdata = data.frame(pred = pred_rf$predictions), 
                      type = "response")
  pred <- data.frame(id = nrow(test_data),
                     obs = test_data$species_observed,
                     pred = pred_cal) %>% 
    drop_na()
  
  # validation metrics
  # mean squared error (mse)
  mse <- mean((pred$obs - pred$pred)^2, na.rm = TRUE)
  # pick threshold to maximize kappa
  thresh <- optimal.thresholds(pred, opt.methods = "MaxKappa")[1, 2]
  # calculate accuracy metrics: auc, kappa, sensitivity, specificity, brier
  pa_metrics <- presence.absence.accuracy(pred, threshold = thresh, 
                                          na.rm = TRUE, st.dev = FALSE)
  
  # summarise the performance of this model
  data.frame(n_checklists = nrow(data),
             n_sightings = sum(data$species_observed),
             mse = round(mse, 5),
             auc = round(pa_metrics$AUC, 5),
             # threshold required
             threshold = round(thresh, 5),
             kappa = round(pa_metrics$Kappa, 5), 
             sensitivity = round(pa_metrics$sensitivity, 5),
             specificity = round(pa_metrics$specificity, 5)
  )
}
bp_runs$ppms <- pmap(bp_runs, fit_bp_model, 
                     data = model_data, test_data = test_data,
                     spacing = sample_spacing, regime = sample_regime)


# validation comparison ----

ppm <- unnest(bp_runs) %>% 
  mutate(n_checklists_log = log10(n_checklists),
         n_sightings_log = log10(n_sightings))
str_glue("output/07_bad_4_sample-size-effect_assessment_{sp_code}.csv") %>% 
  write_csv(ppm, .)

# plot ppms comparison
metrics <- tibble(metric = c("mse", "auc",
                             "kappa", "sensitivity", "specificity"),
                  label = c("Mean Squared Error (MSE)", "AUC",
                            "Kappa", "Sensitivity", "Specificity"))
ppm_best <- filter(ppm, good)
ppm_bad <- filter(ppm, !good)

str_glue("figures/07_bad_4_sample-size-effect_assessment_{sp_code}.png") %>% 
  png(width = 2200, height = 1500, res = 300)

par(mfrow = c(2, 3), mar = c(1, 5, 1, 1), oma = c(4, 1, 1, 1))
cols <- c(bad = "grey65", best = "black")

# performance metrics
for (i in seq.int(nrow(metrics))) {
  m <- metrics$metric[i]
  brks <- pretty(c(0, ppm[[m]]), n = 4)
  boxplot(as.formula(paste(m,  "~ prop_data")), 
          data = ppm_bad,
          range = 0, boxwex = 0.2, lty = 1, staplewex = 0,
          boxcol = "transparent", 
          border = cols["bad"], 
          xlim = c(0.5, 11), ylim = range(brks), 
          xaxt = "n", yaxt = "n", 
          xlab = "", ylab = metrics$label[i])
  par(new = TRUE)
  boxplot(as.formula(paste(m,  "~ prop_data")), 
          data = ppm_best,
          range = 0, boxwex = 0.2, lty = 1, staplewex = 0,
          boxcol = "transparent", 
          border = cols["best"], 
          xlim = c(0.15, 10.65), ylim = range(brks), 
          xaxt = "n", yaxt = "n",
          xlab = "", ylab = "")
  par(new = TRUE)
  axis(side = 2, at = brks)
  if (i > 3) {
    axis(side = 1, at = c(1, 5, 10), labels = c(0, 0.5, 1))
  }
  text(x = 0.5, y = max(brks) - 0.02 * diff(range(brks)), 
       labels = LETTERS[i], 
       font = 2)
}

# number checklists
brks <- pretty(c(ppm$n_checklists_log, ppm$n_sightings_log), n = 4)

boxplot(n_sightings_log ~ prop_data, 
        data = ppm_bad,
        range = 0, boxwex = 0.2, lty = 1, staplewex = 0,
        boxcol = "transparent", 
        border = cols["bad"], 
        xlim = c(0.5, 11), ylim = range(brks), 
        xaxt = "n", yaxt = "n", 
        xlab = "", ylab = "No. checklists")
l <- distinct(ppm_bad, prop_data, n_checklists_log)
lines(l$prop_data * 10, l$n_checklists_log, col = cols["bad"])
par(new = TRUE)
boxplot(n_sightings_log ~ prop_data, 
        data = ppm_best,
        range = 0, boxwex = 0.2, lty = 1, staplewex = 0,
        boxcol = "transparent", 
        border = cols["best"], 
        xlim = c(0.15, 10.65), ylim = range(brks), 
        xaxt = "n", yaxt = "n",
        xlab = "", ylab = "")
l <- distinct(ppm_best, prop_data, n_checklists_log)
lines(l$prop_data * 10, l$n_checklists_log, col = cols["best"])

# axes
axis(side = 1, at = c(0, 5, 10), labels = c(0, 0.5, 1))
lbl_brks <- (1:10) %>% 
  keep(~ . <= max(brks)) %>% 
  keep(~ . >= min(brks))
lbl <- sci_notation(10^lbl_brks)
axis(side = 2, at = lbl_brks, labels = lbl)
text(x = 0.5, y = max(brks) - 0.02 * diff(range(brks)), 
     labels = LETTERS[6], 
     font = 2)

# add legend
legend(x = 4, y = lbl_brks[2], 
       lwd = 2, col = cols, 
       legend = c("Model 2", "Best practice"), cex = 0.8)

par(mfrow = c(1, 1), mar = c(5, 0, 0, 0), oma = c(0, 0, 0, 0))
title(xlab = "Proportion of total data", cex.lab = 0.75)

dev.off()