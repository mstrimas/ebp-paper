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

# read in functions
walk(list.files("R", full.names = TRUE), source)

set.seed(1)
# set species for analysis
species <- "Wood Thrush"
sp_code <- ebird_species(species, "code")

# setup spatial sampling regime
sample_spacing <- 5

# set paths ----
figure_folder <- "figures/"
output_folder <- "output/"


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
        select(checklist_id, sampling_event_identifier, observer_id, species_observed, latitude, longitude,
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
  figure_folder <- paste0(figure_folder, "with_cci/occupancy/")
  output_folder <- paste0(output_folder, "with_cci/occupancy/")
} else {
  figure_folder <- paste0(figure_folder, "without_cci/occupancy/")
  output_folder <- paste0(output_folder, "without_cci/occupancy/")  
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

# formatting data columns
ebird_habitat <- ebird_habitat %>%
        mutate(protocol_traveling = ifelse(protocol_type == "Traveling", 1, 0)) %>%
        mutate(time_observations_started = as.numeric(as.character(time_observations_started))) %>%
        mutate(number_observers = as.numeric(as.character(number_observers))) %>%
        mutate(duration_minutes = as.numeric(as.character(duration_minutes))) %>%
        mutate(day_of_year = yday(observation_date))


# assign observations to hexagonal grid with cell size ~ 5 km 
dggs <- dgconstruct(spacing = sample_spacing)
ebird_habitat_cell <- ebird_habitat %>% 
  mutate(cell = dgGEO_to_SEQNUM(dggs, longitude, latitude)$seqnum) %>%
  mutate(locality_id = paste(round(longitude, 7), round(latitude, 7), sep="_"))
#  mutate(locality_id = cell)

cell_lat_lon <- ebird_habitat_cell %>%
                filter(type == "train") %>%
                select(cell, latitude, longitude, starts_with("pland_")) %>%
                group_by(cell) %>%
                summarise_all(median) %>%
                ungroup()

# ebird_habitat_cell <- ebird_habitat_cell %>%
#                 select(-latitude, -longitude, -starts_with("pland_")) %>%
#                 left_join(cell_lat_lon)

# prepare for unmarked ----

# require a minimum of 2 observations per site (location-observers combination)
# period of closure = june of each year
occ <- ebird_habitat_cell %>%
        filter(type == "train") %>%
        filter_repeat_visits(min_obs = 2, max_obs = 10,
                            annual_closure = TRUE,
                            date_var = "observation_date",
                            site_vars = c("locality_id")) # , "observer_id"))

occ_test_2017 <- ebird_habitat_cell %>%
        filter(type == "test_2017") %>%
        filter_repeat_visits(min_obs = 2, max_obs = 10,
                            annual_closure = TRUE,
                            date_var = "observation_date",
                            site_vars = c("locality_id")) #, "observer_id"))

occ_test_bbs <- ebird_habitat_cell %>%
        filter(type == "test_bbs") %>%
        filter_repeat_visits(min_obs = 2, max_obs = 10,
                            annual_closure = TRUE,
                            date_var = "observation_date",
                            site_vars = c("locality_id")) #, "observer_id"))


nrow(ebird_habitat_cell)
# 35446

nrow(occ)
# 8378  # lat, lon, observer
# 8659  # lat, lon

length(table(occ$locality_id))
# 2115  # lat, lon, observer
# 2300  # lat, lon


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
o_covs <- c("day_of_year",
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
                                               "cell",
                                               # habitat covariates
                                               plands),
                                 obs_covs = o_covs)


# spatial subsampling ----
# sample one record (set of repeated observations at a site) per grid cell

seed <- 1
set.seed(seed)
occ_sss <- occ_wide %>% 
  sample_frac(0.75) %>%
  group_by(cell) %>% 
  sample_n(size = 1) %>% 
  ungroup() %>% 
  select(-cell)

# number of 'sites'
nrow(occ_sss)
# 1571  # lat, lon, observer

# average number of visits per 'site'
10 - round(mean(is.na(occ_sss[,paste0("y.", 1:10)])),2)*10
# 3.2  # lat, lon, observer


# fit model ----

# make unmarked object from the spatially subsampled data
occ_um <- formatWide(occ_sss, type = "unmarkedFrameOccu")

# fit a global occupancy model, including all covariates: 
# four habitat types and seven detection covariates
# specify detection model
model_formula <- o_covs %>% 
  paste(collapse = " + ") %>% 
  paste("~", ., "~ .") %>% 
  as.formula()
# add in species specific occupancy formula
occ_formula <- keep(names(occ_wide), ~ str_detect(., "^pland")) %>% 
  paste(collapse = " + ") %>% 
  paste("~ ", .) %>% 
  as.formula()
model_formula <- update(model_formula, occ_formula)

# fit model
occ_model <- occu(model_formula, data = occ_um)


# model assessment and selection ----

# occ_gof <- mb.gof.test(occ_model, nsim = 10, plot.hist = FALSE)
# str_glue("output/05_occupancy-model_gof_{sp_code}.rds") %>% 
#   saveRDS(occ_gof, .)


# model selection ----

# try all possible permutations of the occupancy covariates
det_terms <- getAllTerms(occ_model) %>% 
  keep(str_detect, pattern = "^p\\(")
occ_dredge <- dredge(occ_model, fixed = det_terms)
# alternatively, dredge on all covariates
# later prediction step will take many hours
# occ_dredge <- dredge(occ_model)

# subset to those with the most suport for model averaging
occ_dredge_95 <- get.models(occ_dredge, subset = cumsum(weight) <= 0.95)

# average models based on model weights 
if (length(occ_dredge_95) > 1) {
  occ_avg <- model.avg(occ_dredge_95, fit = TRUE)
} else {
  occ_avg <- occ_model
}


# predict occupancy ----

# predicting on a model vs a model average returns different results
occ_pred <- predict(occ_avg,
                    newdata = as.data.frame(pred_surface), 
                    type = "state")
if (inherits(occ_pred, "data.frame")) {
  occ_pred <- occ_pred %>% 
    select(occ_prob = Predicted, occ_se = SE) %>% 
    bind_cols(pred_surface, .)
} else {
  occ_pred <- as.data.frame(occ_pred) %>% 
    select(occ_prob = fit, occ_se = se.fit) %>% 
    bind_cols(pred_surface, .)
}

# rasterize
r <- raster("data/modis_5xagg.tif")
r_pred <- occ_pred %>% 
  select(occ_prob, occ_se, latitude, longitude) %>% 
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(crs = projection(r)) %>% 
  rasterize(r)
r_pred <- r_pred[[-1]]

str_glue("output/05_occupancy-model_pred-occ_{sp_code}_seed{seed}.tif") %>% 
  writeRaster(r_pred[["occ_prob"]], ., overwrite = TRUE)
str_glue("output/05_occupancy-model_pred-se_{sp_code}_seed{seed}.tif") %>% 
  writeRaster(r_pred[["occ_se"]], ., overwrite = TRUE)


# map predictions ----

r_pred_proj <- projectRaster(r_pred[["occ_prob"]], crs = map_proj$proj4string, 
                             method = "ngb")

str_glue("figures/05_occupancy-model_predictions_{sp_code}_seed{seed}.png") %>% 
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



# predict on test data ----

# predicting on a model vs a model average returns different results
occ_pred_test_bbs <- predict(occ_avg,
                    newdata = as.data.frame(occ_test_bbs),
                    type = "state")
if (inherits(occ_pred_test_bbs, "data.frame")) {
  occ_pred_test_bbs <- occ_pred_test_bbs %>% 
    select(occ_prob = Predicted, occ_se = SE) %>% 
    bind_cols(occ_test_bbs, .)
} else {
  occ_pred_test_bbs <- as.data.frame(occ_pred_test_bbs) %>% 
    select(occ_prob = fit, occ_se = se.fit) %>% 
    bind_cols(occ_test_bbs, .)
}

# predicting detectability
occ_test_bbs$protocol_type <- factor(occ_test_bbs$protocol_type, levels=c("Stationary", "Traveling"))

occ_pred_test_bbs_det <- predict(occ_avg,
                    newdata = as.data.frame(occ_test_bbs),
                    type = "det")
if (inherits(occ_pred_test_bbs_det, "data.frame")) {
  occ_pred_test_bbs_det <- occ_pred_test_bbs_det %>% 
    select(det_prob = Predicted, det_se = SE) %>% 
    bind_cols(occ_pred_test_bbs, .)
} else {
  occ_pred_test_bbs_det <- as.data.frame(occ_pred_test_bbs_det) %>% 
    select(det_prob = fit, det_se = se.fit) %>% 
    bind_cols(occ_pred_test_bbs, .)
}

occ_pred_test_bbs_det$pred_reporting <- occ_pred_test_bbs_det$occ_prob * occ_pred_test_bbs_det$det_prob

str_glue("output/05_occupancy-model_pred-occ_{sp_code}_seed{seed}.tif") %>% 
  writeRaster(r_pred[["occ_prob"]], ., overwrite = TRUE)
str_glue("output/05_occupancy-model_pred-se_{sp_code}_seed{seed}.tif") %>% 
  writeRaster(r_pred[["occ_se"]], ., overwrite = TRUE)


compare_df <- occ_pred_test_bbs_det %>%
  mutate(id = 1:nrow(occ_pred_test_bbs_det), obs = species_observed, pred = pred_reporting) %>%
  select(obs, pred, occ_prob, observer_id, observation_date) %>%
  group_by(observer_id, observation_date) %>%
  summarise(obs = mean(obs), pred = mean(pred), occ_prob = mean(occ_prob), count = n()) %>%
  ungroup() %>%
  filter(count>20)

# visualise model fit
plot(compare_df$pred, compare_df$obs, xlim=c(0, 1), ylim=c(0, 1)); abline(0, 1, col = "grey")
plot(compare_df$occ_prob, compare_df$obs, xlim=c(0, 1), ylim=c(0, 1)); abline(0, 1, col = "grey")




# predict on test data from 2017----

# predicting on a model vs a model average returns different results
occ_pred_test_2017 <- predict(occ_avg,
                    newdata = as.data.frame(occ_test_2017),
                    type = "state")
if (inherits(occ_pred_test_2017, "data.frame")) {
  occ_pred_test_2017 <- occ_pred_test_2017 %>% 
    select(occ_prob = Predicted, occ_se = SE) %>% 
    bind_cols(occ_test_2017, .)
} else {
  occ_pred_test_2017 <- as.data.frame(occ_pred_test_2017) %>% 
    select(occ_prob = fit, occ_se = se.fit) %>% 
    bind_cols(occ_test_2017, .)
}


occ_pred_test_2017_det <- predict(occ_avg,
                    newdata = as.data.frame(occ_test_2017),
                    type = "det")
if (inherits(occ_pred_test_2017_det, "data.frame")) {
  occ_pred_test_2017_det <- occ_pred_test_2017_det %>% 
    select(det_prob = Predicted, det_se = SE) %>% 
    bind_cols(occ_pred_test_2017, .)
} else {
  occ_pred_test_2017_det <- as.data.frame(occ_pred_test_2017_det) %>% 
    select(det_prob = fit, det_se = se.fit) %>% 
    bind_cols(occ_pred_test_2017, .)
}

occ_pred_test_2017_det$pred_reporting <- occ_pred_test_2017_det$occ_prob * occ_pred_test_2017_det$det_prob

str_glue("output/05_occupancy-model_pred-test_2017_{sp_code}_seed{seed}.tif") %>% 
  writeRaster(occ_pred_test_2017[["occ_prob"]], ., overwrite = TRUE)
str_glue("output/05_occupancy-model_pred-test_2017_se_{sp_code}_seed{seed}.tif") %>% 
  writeRaster(occ_pred_test_2017[["occ_se"]], ., overwrite = TRUE)


compare_df <- occ_pred_test_2017_det %>%
  mutate(id = 1:nrow(occ_pred_test_2017_det), obs = species_observed, pred = pred_reporting) %>%
  select(obs, pred, occ_prob, cell) %>%
  group_by(cell) %>%
  summarise(obs = mean(obs), pred = mean(pred), occ_prob = mean(occ_prob), count = n()) %>%
  ungroup() %>%
  filter(count>10)

# visualise model fit
plot(compare_df$pred, compare_df$obs, xlim=c(0, 1), ylim=c(0, 1)); abline(0, 1, col = "grey")
plot(compare_df$occ_prob, compare_df$obs, xlim=c(0, 1), ylim=c(0, 1)); abline(0, 1, col = "grey")

g <- ggplot(data = compare_df, aes(x = pred, y = obs)) + geom_point() + 
  gg_theme()





# occupancy vs encounter rate ----

r_comp <- c("output/04_rf-model_predictions_{sp_code}.tif",
            "output/05_occupancy-model_pred-occ_{sp_code}_seed{seed}.tif") %>% 
  map_chr(str_glue) %>% 
  stack() %>% 
  projectRaster(crs = map_proj$proj4string, method = "ngb") %>% 
  setNames(c("rf", "occ"))

# breaks and palette
plasma_rev <- rev(plasma(25, end = 0.9))
gray_int <- colorRampPalette(c("#dddddd", plasma_rev[1]))
pal <- c(gray_int(4)[2], plasma_rev)
mx <- ceiling(100 * max(cellStats(r_comp, max))) / 100
brks <- seq(0, mx, length.out = length(pal) + 1)

str_glue("figures/05_comparison_rf-occ_{sp_code}.png") %>% 
  png(width = 2400, height = 3000, res = 300)

par(mfrow = c(2, 1), mar = c(0.5, 0.5, 0.5, 0.5), omi = c(0.6, 0, 0, 0))
for (i in names(r_comp)) {
  r_plot <- r_comp[[i]]
  name <- if_else(i == "rf", "Encounter Rate (RF)", "Occupancy Probability")
  
  # set up plot area
  plot(bcr, col = NA, border = NA)
  plot(ne_land, col = "#dddddd", border = "#888888", lwd = 0.5, add = TRUE)
  
  # probability of detection
  plot(r_comp[[i]], col = pal, breaks = brks, maxpixels = ncell(r_comp),
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