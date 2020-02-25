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
sample_regime <- "together"
sample_spacing <- 5
calibrate <- TRUE

# --------------------------------------------------------------------
# set paths
figure_folder <- "figures/"
output_folder <- "output/"


date <- Sys.Date()
run_name <- paste0("ss_together_rf_val_both_", date)

# load data ----

data_folder <- "data/"
data_tag <- "mayjune_201718_bcr27"

# ebird data
ebird <- read_csv(paste0(data_folder, "data_bad_4_models_", data_tag, ".csv"), na = "") 
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
  ebird_habitat <- left_join(ebird_habitat, cci, by = "checklist_id")
  figure_folder <- paste0(figure_folder, "with_cci/encounter/", run_name, "/")
  output_folder <- paste0(output_folder, "with_cci/encounter/", run_name, "/")
} else {
  figure_folder <- paste0(figure_folder, "without_cci/encounter/", run_name, "/")
  output_folder <- paste0(output_folder, "without_cci/encounter/", run_name, "/")  
}

dir.create(figure_folder, recursive = TRUE)
dir.create(output_folder, recursive = TRUE)


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

# standardized test dataset for all models
test_bbs_id <- ebird_habitat %>%
  filter(type == "test_bbs") %>%
  mutate(week = lubridate::week(observation_date)) %>%
  hex_sample(spacing = sample_spacing,
                       regime = "both", byvar = "week") %>%
  select(checklist_id, sampling_event_identifier, type) %>%
  mutate(selected = 1)

test_2017_id <- ebird_habitat %>%
  filter(type == "test_2017") %>%
  filter(all_species_reported) %>%
  drop_na() %>%
  mutate(week = lubridate::week(observation_date)) %>%
  hex_sample(spacing = sample_spacing,
                       regime = "both", byvar = "week") %>%
  select(checklist_id, sampling_event_identifier, type) %>%
  mutate(selected = 1)

train_2018_id <- ebird_habitat %>%
  filter(type == "train") %>% 
  mutate(week = lubridate::week(observation_date)) %>%
  hex_sample(spacing = sample_spacing, regime = sample_regime, byvar = "week") %>%
  select(checklist_id, sampling_event_identifier, type) %>%
  mutate(selected = 1)

ebird_ss <- rbind(test_bbs_id, test_2017_id, train_2018_id) %>%
        right_join(ebird_habitat) %>%
        mutate(selected = ifelse(is.na(selected), 0, selected)) %>%
        mutate(protocol_traveling = ifelse(protocol_type == "Traveling", 1, 0)) %>%
        mutate(time_observations_started = as.numeric(as.character(time_observations_started))) %>%
        mutate(number_observers = as.numeric(as.character(number_observers))) %>%
        mutate(day_of_year = yday(observation_date))


# test dataset ----

# standardized test dataset for all models
ebird_test_bbs <- ebird_ss %>%
  filter(type == "test_bbs")

ebird_test_2017 <- ebird_ss %>%
  filter(type == "test_2017", selected == 1)

# define the training data
model_data <- filter(ebird_ss, type == "train")


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
  select(run_id, everything()) %>%
  mutate(calibrate_plot_name = paste0(figure_folder, "calibration_model_", run_id, "_", date, ".png"))


# fit bad practice model ----
bp_runs$models <- pmap(bp_runs, fit_bp_model, data = model_data,
                       spacing = sample_spacing, regime = sample_regime, 
                       calibrate = calibrate, calibrate_plot = TRUE)

# amount of data in each run
run_counts <- bp_runs %>% 
  mutate(n_checklists = unlist(map(models, "n_checklists")),
         n_sightings = unlist(map(models, "n_sightings"))) %>% 
  select(-models)
str_glue("{output_folder}/07_bad_1_rf-model_counts_{sp_code}.csv") %>% 
  write_csv(run_counts, .)


# ####################################################################
# validation ----
bp_runs_bbs <- mutate(bp_runs, ppms = map(models, validate, data = ebird_test_bbs))
ppm_bbs <- bp_runs_bbs %>% 
  select(-models) %>% 
  unnest(cols = c(ppms))
str_glue("{output_folder}/07_bad_1_rf-model_assessment_test_bbs_{sp_code}_{date}.csv") %>% 
  write_csv(ppm_bbs, .)

bp_runs_2017 <- mutate(bp_runs, ppms = map(models, validate, data = ebird_test_2017))
ppm_2017 <- bp_runs_2017 %>% 
  select(-models) %>% 
  unnest(cols = c(ppms))
str_glue("{output_folder}/07_bad_1_rf-model_assessment_test_2017_{sp_code}_{date}.csv") %>% 
  write_csv(ppm_2017, .)

for(i in 1:2){

  ppm <- ppm_bbs
  if(i==2) ppm <- ppm_2017
  val_type <- c("bbs", "2017")[i]

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
  str_glue("{figure_folder}/07_bad_1_rf-model_assessment_{val_type}_{sp_code}.png")  %>% 
    ggsave(g_ppm, width = 20, height = 20, units = "cm", dpi = 300)

}

# prediction ----

predict_raster <- function(model, data, template) {
  # add effort covariates to prediction surface
  data <- data %>% 
    mutate(observation_date = ymd("2018-06-15"),
           day_of_year = yday(observation_date),
           time_observations_started = model$t_max_det,
           duration_minutes = 60,
           effort_distance_km = 1,
           number_observers = 1, 
           checklist_calibration_index = 2,
           protocol_type = "Traveling", 
           protocol_traveling = 1)
  
  # predict
  pred <- predict_bp_model(model, data)
  pred_df <- bind_cols(data, prob = pred)
  
  # rasterize
  print("rasterize")
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
r_pred <- str_glue("{output_folder}/07_bad_1_rf-model_predictions_{sp_code}_{date}.tif") %>% 
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
str_glue("{figure_folder}/07_bad_1_rf-model_density_{sp_code}.png")  %>% 
  ggsave(g_density, width = 20, height = 12, units = "cm", dpi = 300)


# maps ----

r_pred_proj <- projectRaster(r_pred, crs = map_proj$proj4string, method = "ngb")

# breaks and palette
plasma_rev <- rev(plasma(25, end = 0.9))
gray_int <- colorRampPalette(c("#dddddd", plasma_rev[1]))
pal <- c(gray_int(4)[2], plasma_rev)
# mx <- ceiling(100 * max(cellStats(r_pred_proj, max))) / 100
mx <- 1
brks <- seq(0, mx, length.out = length(pal) + 1)

str_glue("{figure_folder}/07_bad_1_rf-model_predictions_{sp_code}_{date}.png") %>% 
  png(width = 2400, height = 3000, res = 300)

  par(mfrow = c(3, 2), mar = c(0.5, 0.5, 0.5, 0.5), omi = c(0.6, 0, 0, 0))
  for (i in seq.int(nlayers(r_pred_proj))) {
    r_plot <- r_pred_proj[[i]]
    name <- paste("Model", i)
    
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



