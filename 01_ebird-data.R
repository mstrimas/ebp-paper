# extract and zero-fill ebird data
# three species are included, although only wood thrush is explored in paper
# a incorrectly processed dataset is also prepared for the "bad" analysis
# checklist species richness is also calculate to be used for cci

library(auk)
library(sf)
library(dplyr)
library(purrr)
library(lubridate)
library(readr)
# custom functions
walk(list.files("R", full.names = TRUE), source)


# zero-filled ebird data ----

# ebd extraction
f_ebd <- "data/ebd_june_bcr27.txt"
f_sampling <- "data/ebd_june_bcr27_sampling.txt"
if (!file.exists(f_ebd)) {
  ebd_filtered <- auk_ebd("ebd_relAug-2018.txt", 
                          file_sampling = "ebd_sampling_relAug-2018.txt") %>% 
    auk_species(c("Wood Thrush", "Northern Bobwhite", 
                  "White Ibis")) %>% 
    # southeastern coastal plain bcr
    auk_bcr(bcr = 27) %>% 
    # june, any year
    auk_date(date = c("*-06-01", "*-06-30")) %>% 
    # exclude incidental and other non-standard protocols
    auk_protocol(protocol = c("Stationary", "Traveling")) %>% 
    auk_complete() %>% 
    auk_filter(file = f_ebd, file_sampling = f_sampling)
}
# zero fill
ebd_zf <- auk_zerofill(f_ebd, f_sampling, collapse = TRUE)

# clean up variables
ebd_zf <- ebd_zf %>% 
  mutate(
    # use species code
    species_code = ebird_species(scientific_name, "code"),
    # convert X to NA
    observation_count = if_else(observation_count == "X", 
                                NA_character_, observation_count),
    observation_count = as.integer(observation_count),
    # effort_distance_km to 0 for non-travelling counts
    effort_distance_km = if_else(protocol_type != "Traveling", 
                                 0, effort_distance_km),
    # convert time to decimal hours since midnight
    time_observations_started = time_to_decimal(time_observations_started))

# additional filtering
ebd_zf <- ebd_zf %>% 
  filter(
    # effort filters
    duration_minutes <= 5 * 60,
    effort_distance_km <= 5,
    # last 10 years of data
    year(observation_date) >= 2008,
    !is.na(time_observations_started),
    # 10 or fewer observers
    number_observers <= 10)

# output
ebd_zf <- ebd_zf %>% 
  select(checklist_id, observer_id, sampling_event_identifier,
         species_code,
         observation_count, species_observed, 
         state_code, locality_id, latitude, longitude,
         protocol_type, all_species_reported,
         observation_date, time_observations_started, 
         duration_minutes, effort_distance_km,
         number_observers)
write_csv(ebd_zf, "data/ebd_june_bcr27_zf.csv", na = "")


# bad data ----

# data for bad practices
# intentionally keep incomplete checklists and all protocols
f_ebd_bad <- "data/ebd_june_bcr27_bad.txt"
f_sampling_bad <- "data/ebd_june_bcr27_bad_sampling.txt"
if (!file.exists(f_ebd_bad)) {
  ebd_filtered <- auk_ebd("ebd_relAug-2018.txt", 
                          file_sampling = "ebd_sampling_relAug-2018.txt") %>% 
    auk_species(c("Wood Thrush", "Northern Bobwhite", 
                  "White Ibis")) %>% 
    # southeastern coastal plain bcr
    auk_bcr(bcr = 27) %>% 
    # june, any year
    auk_date(date = c("*-06-01", "*-06-30")) %>% 
    auk_filter(file = f_ebd_bad, file_sampling = f_sampling_bad)
}
# zero fill
# intential bad practice: zero-filling incomplete checklists
ebd_zf_bad <- auk_zerofill(f_ebd_bad, f_sampling_bad, 
                           collapse = TRUE, complete = FALSE)

# clean up variables
ebd_zf_bad <- ebd_zf_bad %>% 
  mutate(
    # use species code
    species_code = ebird_species(scientific_name, "code"),
    # convert X to NA
    observation_count = if_else(observation_count == "X", 
                                NA_character_, observation_count),
    observation_count = as.integer(observation_count),
    # effort_distance_km to 0 for non-travelling counts
    effort_distance_km = if_else(protocol_type != "Traveling", 
                                 0, effort_distance_km),
    # convert time to decimal hours since midnight
    time_observations_started = time_to_decimal(time_observations_started))

# additional filtering - last 10 years of data
# additional effort filters intentionally not applied
ebd_zf_bad <- filter(ebd_zf_bad, year(observation_date) >= 2008)

# output
ebd_zf_bad <- ebd_zf_bad %>% 
  select(checklist_id, observer_id, sampling_event_identifier,
         species_code,
         observation_count, species_observed, 
         state_code, locality_id, latitude, longitude,
         protocol_type, all_species_reported,
         observation_date, time_observations_started, 
         duration_minutes, effort_distance_km,
         number_observers)
write_csv(ebd_zf_bad, "data/ebd_june_bcr27_bad_zf.csv", na = "")


# richness ----

# all species data to calculate richness for cci
f_all <- "data/ebd_june_bcr27_all-species.txt"
if (!file.exists(f_all)) {
  cols <- c("sampling_event_identifier", "group_identifier", "observer_id",
            "taxonomic_order", "scientific_name", "subspecies_scientific_name",
            "protocol_type",
            "observation_date", "time_observations_started",
            "duration_minutes", "effort_distance_km", "number_observers",
            "observation count")
  auk_ebd("ebd_relAug-2018.txt") %>% 
    # southeastern coastal plain bcr
    auk_bcr(bcr = 27) %>% 
    # june, any year
    auk_date(date = c("*-06-01", "*-06-30")) %>% 
    # exclude incidental and other non-standard protocols
    auk_protocol(protocol = c("Stationary", "Traveling")) %>% 
    auk_complete() %>% 
    auk_filter(file = f_all, keep = cols)
}
ebd_all_zf <- auk_zerofill(f_all, f_sampling)
richness <- ebd_all_zf$observations %>% 
  group_by(checklist_id) %>% 
  summarize(n_species = sum(species_observed))
write_csv(richness, "data/richness_june_bcr27.csv")