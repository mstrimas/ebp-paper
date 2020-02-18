library(tidyverse)
library(auk)
library(lubridate)
library(mgcv)
list.files("R", full.names = TRUE) %>% 
  walk(source)


# data location and name
data_folder <- "data/"
data_tag <- "mayjune_201718_bcr27"


# sampling event information
checklist_ids <- read_sampling(paste0(data_folder, "ebd_sampling_", data_tag, "_bad.txt")) %>%
  filter(yday(observation_date)>134) %>%
  select(checklist_id, group_identifier, protocol_code)

# select only for checklists still present after additional effort filters
sampling_events <- read_csv(paste0(data_folder, "ebd_", data_tag, "_zf_bad.csv"), na = "") %>%
  filter(yday(observation_date)>134) %>%
  select(-observation_count, -species_observed) %>%
  distinct() %>%
  inner_join(checklist_ids) %>%
  select(-species_code) %>%
  distinct()

# modis covariates
habitat <- read_csv(paste0(data_folder, "modis_pland_checklists_", data_tag, ".csv"), 
                    col_types = cols(
                      .default = col_double(),
                      checklist_id = col_character())) %>%
          distinct()

# add richness and modis covariates to sampling event data
richness <- read_csv(paste0(data_folder, "ebd_", data_tag, "_zf_speciesrichness.csv")) %>% 
  right_join(sampling_events, by = "checklist_id") %>% 
  left_join(habitat, by = "checklist_id")

# check that all checklists have richness and habitat info
stopifnot(all(!is.na(richness$n_species)))
stopifnot(all(!is.na(richness$pland_00)))

# run the model to calculate the "checklist_calibration_index"


# --------------------------------------------------------------------
# data manipulation and preparation

# calculate habitat diversity
gs <- function(p) {1 - sum(p^2)}
richness$hab_diversity <- richness %>%
                        select(starts_with("pland")) %>%
                        apply(1, FUN = gs)

# convert observer_id to factor and sort
richness$observer_id <- str_split(richness$observer_id, ",") %>% 
  map_chr(~ paste(sort(.), collapse = ","))
richness$observer_id <- factor(richness$observer_id)

# create checklist_number and log(checklist_number) variable
richness2 <- richness %>%
              arrange(observation_date, time_observations_started) %>%
              group_by(observer_id) %>%
              mutate(checklist_number = row_number()) %>%
              ungroup() %>%
              mutate(log_checklist_number = log10(checklist_number)) %>%

              # sort out some of the effort variables
              mutate(effort_hrs = as.numeric(as.character(duration_minutes)) / 60) %>%
              mutate(sqrt_effort_hrs = sqrt(effort_hrs)) %>%
              mutate(effort_distance_km = ifelse(protocol_type == "Stationary", 0, effort_distance_km)) %>%
              mutate(i_stationary = ifelse(protocol_type == "Stationary", 1, 0)) %>%

              # manipulate day and time to continuous covariates
              mutate(day = yday(ymd(observation_date))) %>%
              mutate(time = time_observations_started)


# --------------------------------------------------------------------
# create prediction file with average landcover

prediction_data <- select(richness2, contains("pland")) %>%

                  # take the mean of each column
                  summarise_all(mean) %>%
                  mutate(hab_diversity = gs(pland_00:pland_15)) %>%
                  mutate(day = 165, time = 7, effort_hrs = 1, sqrt_effort_hrs = 1, 
                         effort_distance_km = 1, number_observers = 1, i_stationary = 0) %>%

                  # merge with observers 
                  cbind(richness2 %>% select(observer_id, log_checklist_number, checklist_id))


# --------------------------------------------------------------------
# break data up into groups of up to 300 observers

group_size <- 300

observer_names <- levels(richness2$observer_id)
no_observers <- length(observer_names)
no_groups <- ceiling(no_observers / group_size)

# split observers into n groups/batches
obs_cat <- rep(1:no_groups, group_size)[1:no_observers]
obs_group_match <- data.frame(observer_id = observer_names, group = obs_cat)


# --------------------------------------------------------------------
# loop through each group of observers and run the model

# set up list to collect the models for each group
# and to keep the factor levels for the observers 
group_model_set <- observer_levels <- list(length = no_groups)

for (i in seq_len(no_groups)){

  # ####################################################################
  # MODEL FITTING STEP
  
  # select only the observers in this group
  obs_this_group <- obs_group_match$observer_id[obs_group_match$group==i]
  sub_richness.data <- filter(richness2, observer_id %in% as.character(obs_this_group))

  # remove unused levels of the observer_id factor
  sub_richness.data$observer_id <- sub_richness.data$observer_id[drop = TRUE]
  
  write("length of obs_this_group", stderr())
  write(length(obs_this_group), stderr())

  if (length(obs_this_group) < 10) {
    gam_group <- "no_model"
  } else  {
    # run the model to estimate observer expertise
    gam_group <- try({bam(n_species ~ effort_hrs + sqrt_effort_hrs
                        + pland_00 + pland_01 + pland_02 + pland_03
                        + pland_04 + pland_05 + pland_06
                        + pland_08 + pland_09 + pland_10 + pland_11
                        + pland_12 + pland_13 + pland_14 + pland_15
                        + hab_diversity
                        + i_stationary + effort_distance_km + number_observers
                        + day + s(time, bs = "cc", k = 5) 
                        + log_checklist_number
                        + s(observer_id, bs = "re") 
                        + s(observer_id, effort_hrs, bs = "re") 
                        + s(observer_id, log_checklist_number, bs = "re"),
                        data = sub_richness.data, family = "poisson", discrete = TRUE)}, silent = TRUE)

    if(inherits(gam_group, "try-error")) gam_group <- "no_model"
  }

  # ####################################################################
  # PREDICTION STEP

  # --------------------------------------------------------------------
  # predict for each group using the bam model

  prediction_data$observer_id <- as.character(prediction_data$observer_id)    
  sub_prediction.data <- filter(prediction_data, observer_id %in% obs_this_group)

  # convert the OBSERVER_ID to a factor with levels derived from the model
  sub_prediction.data$observer_id <- factor(as.character(sub_prediction.data$observer_id), levels = levels(sub_richness.data$observer_id))

  if (inherits(gam_group, "bam")){
    # predict from the model to estimate log expected number of species
    pred_group <- predict(gam_group, newdata = sub_prediction.data, type = "link")
  } else {
    # fill prediction column with zeros, if the model didn't run
    pred_group <- rep(NA, nrow(sub_prediction.data))
  }

  # combine with the newdata
  pred_data <- select(sub_prediction.data, 
                      checklist_id, observer_id) %>%
              cbind(exp_species = pred_group) 

  # convert OBSERVER_ID back to character
  pred_data$observer_id <- as.character(pred_data$observer_id)    

  # join together predictions from all the groups
  if (i == 1) {
    preds_all_groups <- pred_data
  } else {
    preds_all_groups <- rbind(preds_all_groups, pred_data)
  }

} # close i

# convert to the standardised index
results <- preds_all_groups %>%
          mutate(cci1 = (exp_species - mean(exp_species)) / sd(exp_species)) %>%
          mutate(cci2 = ifelse(cci1 > 4.5, 4.5, ifelse(cci1 < -4.5, -4.5, cci1))) %>%
          mutate(checklist_calibration_index = round(cci2, 0)) %>%
          select(checklist_id, checklist_calibration_index)

write_csv(results, paste0("data/cci_", data_tag, ".csv"))

