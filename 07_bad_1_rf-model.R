# bad practice random forest model for encounter rate

library(auk)
library(sf)
library(raster)
library(dggridR)
library(ranger)
library(maxnet)
library(scam) 
library(PresenceAbsence)
library(verification)
library(edarf)
library(ggplot2)
library(ggthemes)
library(hexbin)
library(viridis)
library(fields)
library(dplyr)
library(purrr)
library(tidyr)
library(readr)
library(stringr)
library(lubridate)
library(forcats)
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
  ebird_habitat <- left_join(ebird_habitat, cci, by = "checklist_id")
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

bp_runs <- tibble(run_name = c("maxent", "bad_practice", "complete", 
                               "sss", "effort", "best_practice"),
                  maxnet = c(1, 0, 0, 0, 0, 0),
                  complete = c(0, 0, 1, 1, 1, 1),
                  spatial_subsample = c(0, 0, 0, 1, 1, 1),
                  effort_filter = c(0, 0, 0, 0, 1, 1),
                  effort_covs = c(0, 0, 0, 0, 0, 1)) %>% 
  mutate_if(is.numeric, as.logical) %>% 
  mutate(run_id = row_number()) %>% 
  select(run_id, everything())


# fit bad practice model ----

fit_bp_model <- function(maxnet, 
                         complete, spatial_subsample, 
                         effort_filter, effort_covs,
                         data,
                         spacing, regime, ...) {
  if (maxnet) {
    present <- data %>%
      filter(as.logical(species_observed)) %>%
      select(starts_with("pland"))
    # randomly select background points
    bg_n <- 10000
    background <- pred_surface %>%
      sample_n(size = bg_n, replace = FALSE) %>%
      select(starts_with("pland"))
    p <- c(rep(1, nrow(present)), rep(0, nrow(background)))
    combined <- bind_rows(present, background)
    mod <- maxnet(p, combined, maxnet.formula(p, combined, classes = "lq"))
    return(list(model = mod, t_max_det = 7,
                n_checklists = nrow(present), 
                n_sightings = nrow(present)))
  } else {
    # complete checklists only
    if (complete) {
      data <- filter(data, as.logical(all_species_reported))
    }
    # filter on effort covariates
    if (effort_filter){
      data <- data %>%
        filter(protocol_type %in% c("Stationary", "Traveling"),
               effort_distance_km <= 5,
               duration_minutes <= 5 * 60,
               number_observers <= 10)
    }
    # spatial subsample
    if (spatial_subsample) {
      data <- hex_sample(data, spacing = spacing, regime = regime)
    }
    
    # select covariates
    if (effort_covs) {
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
      
      # ensure there are no rows with NAs
      stopifnot(all(complete.cases(data)))
    } else {
      data <- select(data, species_observed, starts_with("pland_"))
      
      # ensure there are no rows with NAs
      stopifnot(all(complete.cases(data)))
    }
    
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
    
    # maximum time of day for detection
    if (effort_covs) {
      pred_t <- partial_dependence(mod, vars = "time_observations_started", 
                                   n = c(24 * 6, nrow(data)), data = data) 
      
      # hours with at least 1% of checklists
      search_hours <- data %>%
        mutate(hour = floor(time_observations_started)) %>%
        count(hour) %>%
        mutate(pct = n / sum(n)) %>%
        filter(pct >= 0.01)
      
      # constrained peak time
      t_max_det <- pred_t %>%
        filter(floor(time_observations_started) %in% search_hours$hour) %>%
        top_n(1, wt = species_observed) %>%
        pull(time_observations_started)
    } else {
      t_max_det <- 7
    }
    
    return(list(model = mod, calibration = mod_cal, t_max_det = t_max_det,
                n_checklists = nrow(data), 
                n_sightings = sum(data$species_observed)))
  }
}
bp_runs$models <- pmap(bp_runs, fit_bp_model, data = model_data,
                       spacing = sample_spacing, regime = sample_regime)
# amount of data in each run
run_counts <- bp_runs %>% 
  mutate(n_checklists = map_int(models, "n_checklists"),
         n_sightings = map_int(models, "n_sightings")) %>% 
  select(-models)
str_glue("output/07_bad_1_rf-model_counts_{sp_code}.csv") %>% 
  write_csv(run_counts, .)

# validation ----

validate <- function(model, data) {
  # predict on test data set
  if (inherits(model$model, "maxnet")) {
    # predict on standardised fixed test dataset
    pred <- predict(model$model, newdata = data, 
                    type = "logistic", clamp = FALSE)
    pred <- data.frame(id = nrow(data), 
                       obs = data$species_observed,
                       pred = pred)
  } else {
    pred_rf <- predict(model$model, data = data, type = "response")
    pred_cal <- predict(model$calibration, 
                        newdata = data.frame(pred = pred_rf$predictions), 
                        type = "response")
    pred <- data.frame(id = nrow(data),
                       obs = data$species_observed,
                       pred = pred_cal)
  }
  pred <- drop_na(pred)
  
  # validation metrics
  # mean squared error (mse)
  mse <- mean((pred$obs - pred$pred)^2, na.rm = TRUE)
  # pick threshold to maximize kappa
  thresh <- optimal.thresholds(pred, opt.methods = "MaxKappa")[1, 2]
  # calculate accuracy metrics: auc, kappa, sensitivity, specificity, brier
  pa_metrics <- presence.absence.accuracy(pred, threshold = thresh, 
                                          na.rm = TRUE, st.dev = FALSE)
  
  # summarise the performance of this model
  tibble(mse = round(mse, 5),
         auc = round(pa_metrics$AUC, 5),
         # threshold required
         threshold = round(thresh, 5),
         kappa = round(pa_metrics$Kappa, 5), 
         sensitivity = round(pa_metrics$sensitivity, 5),
         specificity = round(pa_metrics$specificity, 5)
  )
}
bp_runs <- mutate(bp_runs, ppms = map(models, validate, data = test_data))
ppm <- bp_runs %>% 
  select(-models) %>% 
  unnest()
str_glue("output/07_bad_1_rf-model_assessment_{sp_code}.csv") %>% 
  write_csv(ppm, .)

# plot comparing ppms
ppm_plot <- ppm %>% 
  select_if(~ !is.logical(.)) %>% 
  gather("metric", "value", -run_id, -run_name) %>%
  filter(metric != "threshold") %>% 
  arrange(run_id) %>% 
  mutate(metric = factor(metric, 
                         levels = c("mse", "auc",
                                    "kappa", "sensitivity", "specificity"),
                         labels = c("Mean Squared Error (MSE)", "AUC",
                                    "Kappa", "Sensitivity", "Specificity")),
         run = if_else(run_id == 6, "Best practice", paste("Model", run_id)),
         run = as_factor(run),
         start = if_else(metric == "AUC", 0.5, 0))
g_ppm <- ggplot(ppm_plot) +
  aes(x = run, y = value) +
  geom_point() +
  geom_point(aes(y = start), color = "transparent") +
  facet_wrap(~ metric, nrow = 2, scales = "free_y") +
  scale_y_continuous(expand = expand_scale(mult = c(0, 0.2))) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text = element_text(size = 12, hjust = 0),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = "transparent"),
        axis.ticks.y = element_line(),
        panel.grid = element_blank())
str_glue("figures/07_bad_1_rf-model_assessment_{sp_code}.png")  %>% 
  ggsave(g_ppm, width = 20, height = 20, units = "cm", dpi = 300)


# prediction ----

predict_raster <- function(model, data, template) {
  # add effort covariates to prediction surface
  data <- data %>% 
    mutate(observation_date = ymd("2016-06-15"),
           time_observations_started = model$t_max_det,
           duration_minutes = 60,
           effort_distance_km = 1,
           number_observers = 1, 
           checklist_calibration_index = 2,
           protocol_type = "Traveling")
  
  # predict
  if (inherits(model$model, "maxnet")) {
    pred <- predict(model$model, newdata = data, 
                    type = "logistic", clamp = FALSE) %>% 
      as.vector()
  } else {
    pred_rf <- predict(model$model, data = data, type = "response")
    pred <- predict(model$calibration, 
                    newdata = data.frame(pred = pred_rf$predictions), 
                    type = "response") %>% 
      as.vector()
  }
  pred_df <- bind_cols(data, prob = pred)
  
  # rasterize
  r_pred <- pred_df %>% 
    select(prob, latitude, longitude) %>% 
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
    st_transform(crs = projection(template)) %>% 
    rasterize(template)
  r_pred[[-1]]
}
r <- raster("data/modis_5xagg.tif")
r_pred <- map(bp_runs$models, predict_raster, 
              data = pred_surface, template = r) %>% 
  stack()
r_pred <- str_glue("output/07_bad_1_rf-model_predictions_{sp_code}.tif") %>% 
  writeRaster(r_pred, ., overwrite = TRUE) %>% 
  setNames(bp_runs$run_name)


# density plot ----

# prepare data
pred_compare <- rasterToPoints(r_pred) %>% 
  as_tibble() %>% 
  drop_na() %>% 
  select(-x, -y) %>% 
  gather("run_name", "bad_practice", -best_practice) %>% 
  inner_join(bp_runs %>% select(run_id, run_name), by = "run_name") %>% 
  select(run_id, bad_practice, best_practice)
pred_compare <- pred_compare %>% 
  filter(run_id == 1) %>% 
  mutate(bad_practice = best_practice,
         run_id = 6) %>% 
  bind_rows(pred_compare, .) %>% 
  arrange(run_id) %>% 
  mutate(run = if_else(run_id == 6, "Best practice", paste("Model", run_id)),
         run = as_factor(run))

# plot
g_density <- ggplot(pred_compare) + 
  aes(x = best_practice, y = bad_practice) +
  geom_hex(aes(fill = stat(count))) + 
  scale_x_continuous(breaks = c(0, 0.5, 1)) + 
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  coord_equal() +
  scale_fill_viridis_c(trans = "log10", labels = scales::comma) + 
  facet_wrap(~ run, nrow = 2) + 
  labs(x = "best practice predictions", y = "bad practice predictions") + 
  guides(fill = guide_colorbar(title = "# predictions", 
                               title.position = "left",
                               barwidth = 0.5, barheight = 12)) +
  theme_few() +
  theme(legend.position = "right",
        legend.title = element_text(angle = 90, hjust = 0.5))
str_glue("figures/07_bad_1_rf-model_density_{sp_code}.png")  %>% 
  ggsave(g_density, width = 20, height = 12, units = "cm", dpi = 300)


# maps ----

r_pred_proj <- projectRaster(r_pred, crs = map_proj$proj4string, method = "ngb")

# breaks and palette
plasma_rev <- rev(plasma(25, end = 0.9))
gray_int <- colorRampPalette(c("#dddddd", plasma_rev[1]))
pal <- c(gray_int(4)[2], plasma_rev)
mx <- ceiling(100 * max(cellStats(r_pred_proj, max))) / 100
brks <- seq(0, mx, length.out = length(pal) + 1)

str_glue("figures/07_bad_1_rf-model_predictions_{sp_code}.png") %>% 
  png(width = 2400, height = 3000, res = 300)

par(mfrow = c(3, 2), mar = c(0.5, 0.5, 0.5, 0.5), omi = c(0.6, 0, 0, 0))
for (i in seq.int(nlayers(r_pred_proj))) {
  r_plot <- r_pred_proj[[i]]
  name <- if_else(i == 6, "Best practice", paste("Model", i))
  
  # set up plot area
  plot(bcr, col = NA, border = NA)
  plot(ne_land, col = "#dddddd", border = "#888888", lwd = 0.5, add = TRUE)
  
  # probability of detection
  plot(r_plot, col = pal, breaks = brks, maxpixels = ncell(r_plot),
       legend = FALSE, add = TRUE)
  
  # borders
  plot(bcr, col = NA, border = "#000000", lwd = 1, add = TRUE)
  plot(ne_state_lines, col = "#ffffff", lwd = 0.75, add = TRUE)
  plot(ne_country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
  
  title(main = name, line = -2, cex.main = 1.5, font.main = 1)
  box()
}

# legend
par(new = TRUE, mfrow = c(1, 1), mar = c(0, 0, 0, 0), omi = c(0, 0, 0, 0))
lbl_brks <- seq(0, mx, by = 0.1)
image.plot(zlim = range(brks), legend.only = TRUE, col = pal,
           smallplot = c(0.25, 0.75, 0.035, 0.055),
           horizontal = TRUE,
           axis.args = list(at = lbl_brks, labels = lbl_brks,
                            fg = "black", col.axis = "black",
                            cex.axis = 0.75, lwd.ticks = 0.5,
                            padj = -1.8),
           legend.args = list(text = NULL,
                              side = 3, col = "black",
                              cex = 1, line = 0))
dev.off()