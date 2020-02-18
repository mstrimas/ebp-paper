# occupancy model

library(auk)
library(sf)
library(raster)
library(dggridR)
library(unmarked)
library(MuMIn)
library(AICcmodavg)
library(PresenceAbsence)
library(fields)
library(viridis)
library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(tidyr)
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
ebird <- read_csv(paste0(data_folder, "data_bad_4_models_", data_tag, ".csv"), na = "") 
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

# formatting data columns
ebird_habitat <- ebird_habitat %>%
        mutate(protocol_traveling = ifelse(protocol_type == "Traveling", 1, 0)) %>%
        mutate(time_observations_started = as.numeric(as.character(time_observations_started))) %>%
        mutate(number_observers = as.numeric(as.character(number_observers))) %>%
        mutate(duration_minutes = as.numeric(as.character(duration_minutes))) %>%
        mutate(day_of_year = yday(observation_date)) %>%
        mutate(locality_id = paste0(round(longitude, 7), "_", round(latitude, 7)))

# prepare for unmarked
# require a minimum of 2 observations per site (location-observers combination)
# period of closure = june of each year

number_sites <- number_visits <- vector()

# --------------------------------------------------------------------
# Model 2 - all checklists
paste("Data for model 2")

occ <- filter_repeat_visits(ebird_habitat, min_obs = 2, max_obs = 10,
                            annual_closure = TRUE,
                            date_var = "observation_date",
                            site_vars = c("locality_id")) #, "observer_id"))

str_glue("{data_folder}occ_data_{data_tag}_mod2.csv") %>%
  write_csv(occ, path=.)

number_visits[1] <- nrow(occ)
number_sites[1] <- length(table(occ$locality_id))

# --------------------------------------------------------------------
# Model 3 - complete checklists
paste("Data for model 3")

data <- filter(ebird_habitat, as.logical(all_species_reported))
occ <- filter_repeat_visits(data, min_obs = 2, max_obs = 10,
                            annual_closure = TRUE,
                            date_var = "observation_date",
                            site_vars = c("locality_id")) #, "observer_id"))

str_glue("{data_folder}occ_data_{data_tag}_mod3.csv") %>%
  write_csv(occ, path=.)

number_visits[2] <- nrow(occ)
number_sites[2] <- length(table(occ$locality_id))


# --------------------------------------------------------------------
# Model 4 - complete checklists and spatial subsampling
paste("Data for model 4")

occ_site_only <- occ %>% select(locality_id, longitude, latitude)

# generate hexagonal grid
dggs <- dgconstruct(spacing = 5)
occ_site_subsamp <- occ_site_only %>% 
  mutate(cell = dgGEO_to_SEQNUM(dggs, longitude, latitude)$seqnum) %>%
  group_by(cell) %>% 
  sample_n(size = 1) %>% 
  ungroup() %>% 
  select(-cell) %>%
  select(locality_id) %>%
  left_join(occ)

str_glue("{data_folder}occ_data_{data_tag}_mod4.csv") %>%
  write_csv(occ_site_subsamp, path=.)

number_visits[3] <- nrow(occ_site_subsamp)
number_sites[3] <- length(table(occ_site_subsamp$locality_id))

# --------------------------------------------------------------------
# Model 5 - complete checklists, effort filters, and spatial subsampling
paste("Data for model 5")

data <- data %>%
  mutate(duration_minutes = as.numeric(as.character(duration_minutes))) %>%      
  filter(protocol_type %in% c("Stationary", "Traveling"),
         effort_distance_km <= 5,
         duration_minutes <= 5 * 60,
         number_observers <= 5) %>% 
  mutate(protocol_type = factor(protocol_type, 
                                levels = c("Stationary" , 
                                           "Traveling")))

occ <- filter_repeat_visits(data, min_obs = 2, max_obs = 10,
                            annual_closure = TRUE,
                            date_var = "observation_date",
                            site_vars = c("locality_id")) #, "observer_id"))

occ_site_only <- occ %>% select(locality_id, longitude, latitude)

# generate hexagonal grid
occ_site_subsamp <- occ_site_only %>% 
  mutate(cell = dgGEO_to_SEQNUM(dggs, longitude, latitude)$seqnum) %>%
  group_by(cell) %>% 
  sample_n(size = 1) %>% 
  ungroup() %>% 
  select(-cell) %>%
  select(locality_id) %>%
  left_join(occ)

str_glue("{data_folder}occ_data_{data_tag}_mod5.csv") %>%
  write_csv(occ_site_subsamp, path = .)


number_visits[4] <- nrow(occ_site_subsamp)
number_sites[4] <- length(table(occ_site_subsamp$locality_id))



