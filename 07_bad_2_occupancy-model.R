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

# ebird data
ebird <- read_csv("data/ebd_june_bcr27_bad_zf.csv", na = "") %>% 
  filter(species_code == sp_code) %>% 
  mutate(species_observed = as.integer(species_observed),
         day_of_month = day(observation_date))

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

fit_bp_model <- function(complete, spatial_subsample, 
                         effort_filter, effort_covs,
                         data, spacing, ...) {
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
             number_observers <= 5) %>% 
      mutate(protocol_type = factor(protocol_type, 
                                    levels = c("Stationary" , 
                                               "Traveling")))
  }
  
  # prepare for unmarked
  # require a minimum of 2 observations per site (location-observers combination)
  # period of closure = june of each year
  occ <- filter_repeat_visits(data, min_obs = 2, max_obs = 10,
                              annual_closure = TRUE,
                              date_var = "observation_date",
                              site_vars = c("locality_id", "observer_id"))
  
  # convert to a wide format for use with unmarked
  # include all detection covariates
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
  o_covs <- c("day_of_month",
              "time_observations_started", 
              "duration_minutes", 
              "effort_distance_km", 
              "number_observers", 
              "protocol_type",
              "checklist_calibration_index") %>% 
    intersect(names(occ))
  occ_wide <- format_unmarked_occu(occ, site_id = "site", 
                                   response = "species_observed",
                                   site_covs = c("n_observations", 
                                                 "latitude", 
                                                 "longitude",
                                                 # habitat covariates
                                                 plands),
                                   obs_covs = o_covs)
  
  # spatial subsampling
  if (spatial_subsample) {
    # generate hexagonal grid
    dggs <- dgconstruct(spacing = spacing)
    occ_wide_cell <- occ_wide %>% 
      mutate(cell = dgGEO_to_SEQNUM(dggs, longitude, latitude)$seqnum)
    # sample one record (set of repeated observations at a site) per grid cell
    occ_wide <- occ_wide_cell %>% 
      group_by(cell) %>% 
      sample_n(size = 1) %>% 
      ungroup() %>% 
      select(-cell)
  }
  
  # make unmarked object from the spatially subsampled data
  occ_um <- formatWide(occ_wide, type = "unmarkedFrameOccu")
  
  # fit an occupancy model
  # four habitat types and seven detection covariates
  # specify detection model
  if (effort_covs) {
    model_formula <- o_covs %>% 
      paste(collapse = " + ") %>% 
      paste("~", ., "~ .") %>% 
      as.formula()
  } else {
    model_formula <- ~ 1 ~ .
  }
  
  # add in species specific occupancy formula
  occ_formula <- keep(names(occ_wide), ~ str_detect(., "^pland")) %>% 
    paste(collapse = " + ") %>% 
    paste("~ ", .) %>% 
    as.formula()
  model_formula <- update(model_formula, occ_formula)
  
  # fit model
  occu(model_formula, data = occ_um)
}
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