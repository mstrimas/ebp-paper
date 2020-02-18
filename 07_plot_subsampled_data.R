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
sample_regime <- "both"
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



set.seed(1)


  # reduce to 75% of all data (in all subsets)
  ebird_habitat_sample <- ebird_habitat %>%
      sample_frac(0.75)

  print("spatial subsampling")
  ebird_ss <- hex_sample(ebird_habitat_sample, spacing = sample_spacing,
                         regime = sample_regime, byvar = "type_week") %>%
          select(checklist_id, sampling_event_identifier) %>%
          mutate(selected = 1) %>% 
          left_join(ebird_habitat) %>%
          select(-type_week) %>%
          mutate(selected = ifelse(is.na(selected), 0, 1)) %>%
          mutate(selected = ifelse(type == "test_bbs", 1, ifelse(type == "other", 0, selected)))

  # test dataset ----

  # standardized test dataset for all models
  ebird_test_bbs <- ebird_ss %>%
    filter(type == "test_bbs")

  ebird_test_2017 <- ebird_ss %>%
    filter(type == "test_2017", selected == 1)

  # define the training data
  model_data <- filter(ebird_ss, type == "train")



# --------------------------------------------------------------------
# project the datasets

train_data_proj <- as.matrix(cbind(model_data$longitude, model_data$latitude)) %>%
        sp::SpatialPoints(proj4string = CRS("+init=epsg:4326")) %>% 
        st_as_sf() %>%
        st_transform(crs = map_proj)


test_data_bbs_proj <- as.matrix(cbind(ebird_test_bbs$longitude, ebird_test_bbs$latitude)) %>%
        sp::SpatialPoints(proj4string = CRS("+init=epsg:4326")) %>% 
        st_as_sf() %>%
        st_transform(crs = map_proj)

test_data_2017_proj <- as.matrix(cbind(ebird_test_2017$longitude, ebird_test_2017$latitude)) %>%
        sp::SpatialPoints(proj4string = CRS("+init=epsg:4326")) %>% 
        st_as_sf() %>%
        st_transform(crs = map_proj)


figure_folder <- "figures/"
plot_name <- paste0(figure_folder, "train_test_subsampled_maps_", Sys.Date(), ".png")
png(plot_name, width = 14, height = 12, units="cm", pointsize=9, res=300)

  par(mfrow=c(2,2), mar = c(0.5, 0.5, 0.5, 0.5))

  # set up plot area
  plot(bcr, col = NA, border = NA)
  plot(ne_land, col = "#dddddd", border = "#888888", lwd = 0.5, add = TRUE)

  # borders
  plot(bcr, col = NA, border = "#000000", lwd = 1, add = TRUE)
  plot(ne_state_lines, col = "#ffffff", lwd = 0.75, add = TRUE)
  plot(ne_country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
  box()

  # add the data!
  plot(train_data_proj, add = TRUE, pch=16, cex=0.2)


  # set up plot area
  plot(bcr, col = NA, border = NA)
  plot(ne_land, col = "#dddddd", border = "#888888", lwd = 0.5, add = TRUE)

  # borders
  plot(bcr, col = NA, border = "#000000", lwd = 1, add = TRUE)
  plot(ne_state_lines, col = "#ffffff", lwd = 0.75, add = TRUE)
  plot(ne_country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
  box()

  # add the data!
  plot(test_data_bbs_proj, add = TRUE, pch=16, cex=0.2)


  # set up plot area
  plot(bcr, col = NA, border = NA)
  plot(ne_land, col = "#dddddd", border = "#888888", lwd = 0.5, add = TRUE)

  # borders
  plot(bcr, col = NA, border = "#000000", lwd = 1, add = TRUE)
  plot(ne_state_lines, col = "#ffffff", lwd = 0.75, add = TRUE)
  plot(ne_country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
  box()

  # add the data!
  plot(test_data_2017_proj, add = TRUE, pch=16, cex=0.2)

dev.off()


