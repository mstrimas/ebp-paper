# occupancy model

library(auk)
library(sf)
library(raster)
library(dggridR)
library(unmarked)
library(MuMIn)
library(AICcmodavg)
library(fields)
library(viridis)
library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(lubridate)
# resolve namespace conflicts
select <- dplyr::select
projection <- raster::projection

set.seed(1)
# set species for analysis
species <- "Wood Thrush"
sp_code <- ebird_species(species, "code")
# setup spatial sampling regime
sample_spacing <- 5


# load data ----


data_folder <- "data/"
data_tag <- "mayjune_201718_bcr27"

# ebird data
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
  ebird_habitat <- left_join(ebird_habitat, cci, by = "checklist_id")
}


ebird_habitat <- ebird_habitat %>%
        mutate(type_week = paste(type, lubridate::week(observation_date), sep="_"))

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


# setup bad practice combinations ----

bp_runs <- tibble(run_name = c("bad_practice", "complete", 
                               "sss", "effort", "best_practice"),
                  complete = c(0, 1, 1, 1, 1),
                  spatial_subsample = c(0, 0, 1, 1, 1),
                  effort_filter = c(0, 0, 0, 1, 1),
                  effort_covs = c(0, 0, 0, 0, 1)) %>% 
  mutate_if(is.numeric, as.logical) %>% 
  mutate(run_id = row_number() + 1) %>% 
  select(run_id, everything())


# fit bad practice model ----
set.seed(1)

  print("")
  print("=======================================")
  print(paste("running for sim", i))

  # reduce to 75% of all data (in all subsets)
  ebird_habitat_sample <- ebird_habitat %>%
      sample_frac(0.75)


fit_bp_model_occu <- function(complete, spatial_subsample, 
                         effort_filter, effort_covs,
                         data, spacing, ...) {


bp_runs$models <- pmap(bp_runs, fit_bp_model, data = ebird_habitat,
                       spacing = sample_spacing)


# predict occupancy ----

# prediction ----

predict_raster <- function(model, data, template) {
  pred <- predict(model,newdata = as.data.frame(data), type = "state") %>% 
    select(occ_prob = Predicted, occ_se = SE) %>% 
    bind_cols(pred_surface, .)
  
  # rasterize
  r_pred <- pred %>% 
    select(occ_prob, occ_se, latitude, longitude) %>% 
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
    st_transform(crs = projection(template)) %>% 
    rasterize(template)
  r_pred[[-1]]
}
r <- raster("data/modis_5xagg.tif")
bp_runs <- bp_runs %>% 
  mutate(r_pred = map(models, predict_raster, data = pred_surface, 
                      template = r))

# occupancy probability
f_prob <- str_glue("output/07_bad_2_occupancy-model_pred-occ_{sp_code}.tif")
r_occ_prob <- map(bp_runs$r_pred, ~ .[["occ_prob"]]) %>% 
  stack() %>% 
  writeRaster(f_prob, overwrite = TRUE) %>% 
  setNames(bp_runs$run_name)

# occupancy se
f_se <- str_glue("output/07_bad_2_occupancy-model_pred-se_{sp_code}.tif")
r_occ_se <- map(bp_runs$r_pred, ~ .[["occ_se"]]) %>% 
  stack() %>% 
  writeRaster(f_se, overwrite = TRUE) %>% 
  setNames(bp_runs$run_name)


# maps ----

r_pred_proj <- projectRaster(r_occ_prob, crs = map_proj$proj4string,
                             method = "ngb")

# breaks and palette
plasma_rev <- rev(plasma(25, end = 0.9))
gray_int <- colorRampPalette(c("#dddddd", plasma_rev[1]))
pal <- c(gray_int(4)[2], plasma_rev)
mx <- ceiling(100 * max(cellStats(r_pred_proj, max))) / 100
brks <- seq(0, mx, length.out = length(pal) + 1)

str_glue("figures/07_bad_2_occupancy-model_predictions_{sp_code}.png") %>% 
  png(width = 2400, height = 3000, res = 300)

par(mfrow = c(3, 2), mar = c(0.5, 0.5, 0.5, 0.5), omi = c(0.6, 0, 0, 0))
for (i in seq.int(nlayers(r_pred_proj) + 1)) {
  # keep first panel blank
  if (i == 1) {
    plot.new()
    next()
  }
  
  r_plot <- r_pred_proj[[i - 1]]
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