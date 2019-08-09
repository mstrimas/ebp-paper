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

# location of large unzipped eBird data
ebd_data_folder <- "/Volumes/ebird_data/EBD/"
ebd_release <- "May-2019"

# where to save the processed datasets
ebd_save <- "data/"

# name for this specific data extraction
data_tag <- "mayjune_201718_bcr27"

# zero-filled ebird data ----

# ebd extraction
orig_ebd <- paste0(ebd_data_folder, "ebd_rel", ebd_release, "/ebd_rel", ebd_release, ".txt")
orig_sampling <- paste0(ebd_data_folder, "ebd_sampling_rel", ebd_release, "/ebd_sampling_rel", ebd_release, ".txt")

# names of the processed ebird files
f_ebd <- paste0(ebd_save, "ebd_", data_tag, ".txt")
f_sampling <- paste0(ebd_save, "ebd_sampling_", data_tag, ".txt")

if (!file.exists(f_ebd)) {
  ebd_filtered <- auk_ebd(orig_ebd, 
                          file_sampling = orig_sampling) %>% 
    auk_species(c("Wood Thrush", "Northern Bobwhite", 
                  "White Ibis")) %>% 
    # southeastern coastal plain bcr
    auk_bcr(bcr = 27) %>% 
    # june, any year
    auk_date(date = c("*-05-01", "*-06-30")) %>% 
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
    year(observation_date) > 2016,
    year(observation_date) < 2019,
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
write_csv(ebd_zf, paste0(ebd_save, "ebd_", data_tag, "_zf.csv"), na = "")


# bad data ----

# data for bad practices
# intentionally keep incomplete checklists and all protocols

f_ebd_bad <- paste0(ebd_save, "ebd_", data_tag, "_bad.txt")
f_sampling_bad <- paste0(ebd_save, "ebd_sampling_", data_tag, "_bad.txt")

if (!file.exists(f_ebd_bad)) {
  ebd_filtered <- auk_ebd(orig_ebd, 
                          file_sampling = orig_sampling) %>% 
    auk_species(c("Wood Thrush", "Northern Bobwhite", 
                  "White Ibis")) %>% 
    # southeastern coastal plain bcr
    auk_bcr(bcr = 27) %>% 
    # june, any year
    auk_date(date = c("*-05-01", "*-06-30")) %>% 
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
ebd_zf_bad <- filter(ebd_zf_bad, 
    year(observation_date) > 2016,
    year(observation_date) < 2019)

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
write_csv(ebd_zf_bad, paste0(ebd_save, "ebd_", data_tag, "_zf_bad.csv"), na = "")


# richness ----

# all species data to calculate richness for cci
f_all <- paste0(ebd_save, "ebd_", data_tag, "_zf_speciesrichness.csv")
if (!file.exists(f_all)) {
  cols <- c("sampling_event_identifier", "group_identifier", "observer_id",
            "taxonomic_order", "scientific_name", "subspecies_scientific_name",
            "protocol_type",
            "observation_date", "time_observations_started",
            "duration_minutes", "effort_distance_km", "number_observers",
            "observation count")
  auk_ebd(orig_ebd) %>% 
    # southeastern coastal plain bcr
    auk_bcr(bcr = 27) %>% 
    # june, any year
    auk_date(date = c("*-05-01", "*-06-30")) %>% 
    # exclude incidental and other non-standard protocols
    auk_protocol(protocol = c("Stationary", "Traveling")) %>% 
    auk_complete() %>% 
    auk_filter(file = f_all, keep = cols)

  ebd_all_zf <- auk_zerofill(f_all, f_sampling)
  richness <- ebd_all_zf$observations %>% 
    group_by(checklist_id) %>% 
    summarize(n_species = sum(species_observed))
  write_csv(richness, paste0(ebd_save, "ebd_", data_tag, "_zf_speciesrichness.csv"))

}
