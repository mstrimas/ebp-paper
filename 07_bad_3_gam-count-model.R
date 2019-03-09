# bad practice random forest model for encounter rate

library(auk)
library(sf)
library(raster)
library(dggridR)
library(mgcv)
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
# gam parameters
# degrees of freedom for smoothing
k <- 5
# degrees of freedom for cyclic time of day smooth
k_time <- 7 

# load data ----

# ebird data
ebird <- read_csv("data/ebd_june_bcr27_bad_zf.csv", na = "") %>% 
  filter(species_code == sp_code,
         !is.na(observation_count)) %>% 
  mutate(species_observed = as.integer(species_observed),
         day_of_year = yday(observation_date))

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

bp_runs <- tibble(run_name = c("pres_only", "bad_practice", "complete", 
                               "sss", "effort", "best_practice"),
                  pres_only = c(1, 0, 0, 0, 0, 0),
                  complete = c(0, 0, 1, 1, 1, 1),
                  spatial_subsample = c(0, 0, 0, 1, 1, 1),
                  effort_filter = c(0, 0, 0, 0, 1, 1),
                  effort_covs = c(0, 0, 0, 0, 0, 1)) %>% 
  mutate_if(is.numeric, as.logical) %>% 
  mutate(run_id = row_number()) %>% 
  select(run_id, everything())


# fit bad practice model ----

fit_bp_model <- function(pres_only, 
                         complete, spatial_subsample, 
                         effort_filter, effort_covs,
                         data, test_data,
                         spacing, regime, ...) {
  if (pres_only) {
    present <- data %>%
      filter(as.logical(species_observed)) %>%
      select(observation_count, starts_with("pland"))
    # randomly select background points
    bg_n <- 10000
    background <- pred_surface %>%
      sample_n(size = bg_n, replace = FALSE) %>%
      mutate(observation_count = 0) %>% 
      select(observation_count, starts_with("pland"))
    data <- bind_rows(present, background)
  }
  
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
  
  # pland covariates
  # choose habitat covariates based on prior knowledge of species
  # for wood thrush:
  #   - deciduous (pland_04) and mixed (pland_05) forest; known nesting habitat
  #   - cropland (pland_12) and urban (pland_13); avoided
  # for northern bobwhite:
  #   - savanna (pland_09) and woody savanna (pland_08); preferred
  #   - deciduous forest (pland_04) and urban (pland_13); avoided
  # for white ibis:
  #   - fresh water (pland_00) and wetlands (pland_11); preferred
  #   - savanna (pland_09) and woody savanna (pland_08); avoided
  if (sp_code == "woothr") {
    plands <- paste0("pland_", c("04", "05", "12", "13"))
  } else if (sp_code == "norbob") {
    plands <- paste0("pland_", c("04", "08", "09", "13"))
  } else if (sp_code == "whiibi") {
    plands <- paste0("pland_", c("00", "08", "09", "11"))
  } else {
    stop("species code not valid.")
  }
  
  # select covariates
  if (effort_covs) {
    data <- data %>% 
      mutate(protocol_type = factor(protocol_type, 
                                    levels = c("Stationary" , "Traveling"))) %>%
      # select only the columns to be used in the model
      select(observation_count,
             day_of_year, time_observations_started, duration_minutes,
             effort_distance_km, number_observers, protocol_type,
             contains("checklist_calibration_index"),
             one_of(plands))
    
    # continuous predictors
    # hold out time to treat seperately since it's cyclic
    continuous_covs <- data %>% 
      select(-observation_count, -protocol_type, -time_observations_started) %>% 
      names()
    
    # create model formula for predictors
    gam_formula_rhs <- str_glue("s({var}, k = {k})", 
                                var = continuous_covs, k = k) %>% 
      str_flatten(collapse = " + ") %>% 
      str_glue(" ~ ", .,
               " + protocol_type + ",
               "s(time_observations_started, bs = \"cc\", k = {k})", 
               k = k_time) %>% 
      as.formula()
    
    # ensure there are no rows with NAs
    stopifnot(all(complete.cases(data)))
  } else {
    data <- select(data, observation_count, starts_with("pland_"))
    
    # create model formula for predictors
    gam_formula_rhs <- str_glue("s({var}, k = {k})", 
                                var = plands, k = k) %>% 
      str_flatten(collapse = " + ") %>% 
      paste("~", .) %>% 
      as.formula()
    
    # ensure there are no rows with NAs
    stopifnot(all(complete.cases(data)))
  }
  
  # fit models
  # model formula including response
  gam_formula <- update.formula(observation_count ~ ., gam_formula_rhs)
  
  # explicitly specify where knots should occur for time_observations_started
  # this ensures that the cyclic spline joins the variable at midnight
  # this won't happen by default if there are no data near midnight
  time_knots <- list(time_observations_started = seq(0, 24, 
                                                     length.out = k_time))
  
  # zero-inflated poisson
  m_ziplss <- gam(list(gam_formula, # count model
                       gam_formula_rhs), # presence model
                  data = data, 
                  family = "ziplss", 
                  knots = time_knots)
  
  # negative binomial
  m_nb <- gam(gam_formula,
              data = data, 
              family = "nb",
              knots = time_knots)
  
  # tweedie distribution
  m_tw <- gam(gam_formula,
              data = data, 
              family = "tw",
              knots = time_knots)
  
  # predict on test
  obs_count <- select(test_data, obs = observation_count)
  
  inv_link <- binomial(link = "cloglog")$linkinv
  m_ziplss_pred <- predict(m_ziplss, test_data, type = "link") %>% 
    as.data.frame() %>%
    transmute(family = "Zero-inflated Poisson",
              pred = inv_link(V2) * exp(V1)) %>% 
    bind_cols(obs_count)
  
  m_nb_pred <- predict(m_nb, test_data, type = "response") %>% 
    tibble(family = "Negative Binomial", pred = .) %>% 
    bind_cols(obs_count)
  
  m_tw_pred <- predict(m_tw, test_data, type = "response") %>% 
    tibble(family = "Tweedie", pred = .) %>% 
    bind_cols(obs_count)
  
  test_pred <- bind_rows(m_ziplss_pred, m_nb_pred, m_tw_pred) %>% 
    mutate(family = as_factor(family))
  
  # deviance explained
  de <- c(summary(m_ziplss)$dev.expl,
          summary(m_nb)$dev.expl,
          summary(m_tw)$dev.expl)
  
  # mse and rank correlation
  gam_assessment <- test_pred %>% 
    group_by(family) %>% 
    summarise(mse = mean((obs - pred)^2, na.rm = TRUE),
              rank_cor = cor.test(obs, pred, 
                                  method = "spearman", 
                                  exact = FALSE)$estimate) %>% 
    ungroup() %>% 
    mutate(dev_exp = de)
  
  # select model
  best_mod <- gam_assessment %>% 
    arrange(desc(rank_cor), mse, dev_exp) %>% 
    slice(1) %>% 
    pull(family) %>% 
    as.character()
  if (best_mod == "Zero-inflated Poisson") {
    pred_model <- m_ziplss
  } else if (best_mod == "Negative Binomial") {
    pred_model <- m_nb
  } else {
    pred_model <- m_tw
  }
  
  # maximum time of day for detection
  if (effort_covs) {
    seq_tod <- seq(0, 24, length.out = 300)
    tod_df <- data %>% 
      # find average pland habitat covariates
      select(starts_with("pland")) %>% 
      summarize_all(mean, na.rm = TRUE) %>% 
      ungroup() %>% 
      # use standard checklist
      mutate(day_of_year = yday(ymd("2016-06-15")),
             duration_minutes = 60,
             effort_distance_km = 1,
             number_observers = 1,
             checklist_calibration_index = 2,
             protocol_type = "Traveling") %>% 
      cbind(time_observations_started = seq_tod)
    
    # predict at different start times
    pred_tod <- predict(pred_model, newdata = tod_df, 
                        type = "link", 
                        se.fit = TRUE)
    if (best_mod == "Zero-inflated Poisson") {
      pred_tod <- pred_tod %>% 
        as.data.frame() %>% 
        setNames(c("fit_count", "fit_zero", "se_count", "se_zero")) %>% 
        mutate(time_observations_started = seq_tod,
               pred_lcl = inv_link(fit_zero - 1.96 * se_zero * 
                                     exp(fit_count - 1.96 * se_count)))
    } else {
      pred_tod <- pred_tod %>% 
        as_tibble() %>% 
        # calculate backtransformed confidence limits
        transmute(time_observations_started = seq_tod,
                  pred_lcl = pred_model$family$linkinv(fit - 1.96 * se.fit))
    }
    
    # select detectability peak
    t_max <- pred_tod$time_observations_started[which.max(pred_tod$pred_lcl)]
  } else {
    t_max <- 7
  }
  
  return(list(best_model = best_mod,
              model = pred_model, 
              ppms = gam_assessment,
              t_max_det = t_max))
}
bp_runs$models <- pmap(bp_runs, fit_bp_model, 
                       data = model_data, test_data = test_data,
                       spacing = sample_spacing, regime = sample_regime)


# validation ----

ppm <- bp_runs %>% 
  mutate(ppms = map(models, "ppms"),
         best_model = map_chr(models, "best_model")) %>% 
  select(-models) %>% 
  unnest()
str_glue("output/07_bad_3_gam-count-model_assessment_{sp_code}.csv") %>% 
  write_csv(ppm, .)

# plot comparing ppms
ppm_plot <- ppm %>% 
  filter(best_model == family) %>% 
  select(run_id, run_name, family, mse, rank_cor) %>% 
  gather("metric", "value", -run_id, -run_name, -family) %>%
  arrange(run_id) %>% 
  mutate(metric = factor(metric, 
                         levels = c("mse", "rank_cor"),
                         labels = c("Mean Squared Error (MSE)", 
                                    "Spearman's Rank Correlation")),
         run = if_else(run_id == 6, "Best practice", paste("Model", run_id)),
         run = as_factor(run))
g_ppm <- ggplot(ppm_plot) +
  aes(x = run, y = value) +
  geom_point() +
  facet_wrap(~ metric, nrow = 1, scales = "free_y") +
  scale_y_continuous(expand = expand_scale(mult = c(0, 0.2)), 
                     limits = c(0, NA)) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text = element_text(size = 12, hjust = 0),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = "transparent"),
        axis.ticks.y = element_line(),
        panel.grid = element_blank())
str_glue("figures/07_bad_3_gam-count-model_assessment_{sp_code}.png")  %>% 
  ggsave(g_ppm, width = 15, height = 7, units = "cm", dpi = 300)


# prediction ----

predict_raster <- function(model, data, template) {
  # add effort covariates to prediction surface
  data <- data %>% 
    mutate(observation_date = ymd("2016-06-15"),
           day_of_year = yday(observation_date),
           time_observations_started = model$t_max_det,
           duration_minutes = 60,
           effort_distance_km = 1,
           number_observers = 1, 
           checklist_calibration_index = 2,
           protocol_type = "Traveling")
  
  # predict
  if (model$best_model == "Zero-inflated Poisson") {
    inv_link <- binomial(link = "cloglog")$linkinv
    pred <- predict(model$model, data, type = "link") %>% 
      as.data.frame() %>%
      mutate(pred = inv_link(V2) * exp(V1)) %>% 
      pull(pred)
  } else {
    pred <- predict(model$model, data, type = "response") %>% 
      as.vector()
  }
  pred_df <- bind_cols(data, abd = pred)
  
  # rasterize
  r_pred <- pred_df %>% 
    select(abd, latitude, longitude) %>% 
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
    st_transform(crs = projection(template)) %>% 
    rasterize(template)
  r_pred[[-1]]
}
r <- raster("data/modis_5xagg.tif")
r_pred <- map(bp_runs$models, predict_raster, 
              data = pred_surface, template = r) %>% 
  stack()
r_pred <- str_glue("output/07_bad_3_gam-count-model_predictions_{sp_code}.tif") %>% 
  writeRaster(r_pred, ., overwrite = TRUE) %>% 
  setNames(bp_runs$run_name)


# maps ----

r_pred_proj <- projectRaster(r_pred, crs = map_proj$proj4string, method = "ngb")

# breaks and palette
plasma_rev <- rev(plasma(25, end = 0.9))
gray_int <- colorRampPalette(c("#dddddd", plasma_rev[1]))
pal <- c(gray_int(4)[2], plasma_rev)
mx <- ceiling(100 * max(cellStats(r_pred_proj, max))) / 100
mx <- ceiling(100 * max(cellStats(r_pred_proj, max))) / 100
brks <- seq(0, mx, length.out = length(pal) + 1)

str_glue("figures/07_bad_3_gam-count-model_predictions_{sp_code}.png") %>% 
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
lbl_brks <- seq(0, mx, by = 0.5)
lbl_brks[length(lbl_brks)] <- mx
image.plot(zlim = range(brks), legend.only = TRUE, col = pal,
           smallplot = c(0.25, 0.75, 0.03, 0.05),
           horizontal = TRUE,
           axis.args = list(at = lbl_brks, labels = lbl_brks,
                            fg = "black", col.axis = "black",
                            cex.axis = 0.75, lwd.ticks = 0.5,
                            padj = -1.8),
           legend.args = list(text = NULL,
                              side = 3, col = "black",
                              cex = 1, line = 0))
dev.off()