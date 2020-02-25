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
library(tibble)
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
anchor_model <- 6


# set paths ----
figure_folder <- "figures/"
output_folder <- "output/"

date <- Sys.Date()
date <- "2020-02-19"
run_name <- paste0("ss_", sample_regime, "_rf_val_both_", date)

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

ebird_habitat <- ebird_habitat %>%
        mutate(type_week = paste(type, lubridate::week(observation_date), sep="_")) %>%
        mutate(protocol_traveling = ifelse(protocol_type == "Traveling", 1, 0)) %>%
        mutate(time_observations_started = as.numeric(as.character(time_observations_started))) %>%
        mutate(number_observers = as.numeric(as.character(number_observers))) %>%
        mutate(day_of_year = yday(observation_date))


# fit bad practice model ----


# setup bad practice combinations ----

bp_runs_master <- tibble(run_name = c("maxent", "bad_practice", "complete", 
                               "sss", "effort", "best_practice"),
                  maxnet = c(1, 0, 0, 0, 0, 0),
                  complete = c(0, 0, 1, 1, 1, 1),
                  spatial_subsample = c(0, 0, 0, 1, 1, 1),
                  effort_filter = c(0, 0, 0, 0, 1, 1),
                  effort_covs = c(0, 0, 0, 0, 0, 1)) %>% 
  mutate_if(is.numeric, as.logical) %>% 
  mutate(run_id = row_number()) %>% 
  select(run_id, everything())


nsim <- 25

for(i in 1:nsim){

  print("")
  print("=======================================")
  print(paste("running for sim", i))
  set.seed(i)

  # reduce to 75% of all data (in all subsets)
  ebird_habitat_sample <- ebird_habitat %>%
      sample_frac(0.75)


  # standardized test dataset for all models
  test_bbs_id <- ebird_habitat_sample %>%
    filter(type == "test_bbs") %>%
    mutate(week = lubridate::week(observation_date)) %>%
    hex_sample(spacing = sample_spacing,
                         regime = "both", byvar = "week") %>%
    select(checklist_id, sampling_event_identifier, type) %>%
    mutate(selected = 1)

  test_2017_id <- ebird_habitat_sample %>%
    filter(type == "test_2017") %>%
    filter(all_species_reported) %>%
    drop_na() %>%
    mutate(week = lubridate::week(observation_date)) %>%
    hex_sample(spacing = sample_spacing,
                         regime = "both", byvar = "week") %>%
    select(checklist_id, sampling_event_identifier, type) %>%
    mutate(selected = 1)

  train_2018_id <- ebird_habitat_sample %>%
    filter(type == "train") %>% 
    mutate(week = lubridate::week(observation_date)) %>%
    hex_sample(spacing = sample_spacing, regime = sample_regime, byvar = "week") %>%
    select(checklist_id, sampling_event_identifier, type) %>%
    mutate(selected = 1)

  ebird_ss <- rbind(test_bbs_id, test_2017_id, train_2018_id) %>%
          right_join(ebird_habitat_sample) %>%
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

  test_2017_pos <- ebird_test_2017 %>%
    filter(species_observed==1) 
  test_2017_neg <- ebird_test_2017 %>%
    sample_n(size = nrow(test_2017_pos))

  ebird_test_2017_bal <- rbind(test_2017_pos, test_2017_neg)

  # define the training data
  model_data <- filter(ebird_ss, type == "train")

  # if(validation=="8020"){
  #   s_train <- rbinom(nrow(model_data), 1, 0.8) %>% as.logical()
  #   ebird_test_2017 <- model_data[!s_train,]  %>% filter(selected == 1)
  #   model_data <- model_data[s_train,]
  # }

  print("running models")
  bp_runs <- bp_runs_master

  bp_runs$models <- pmap(bp_runs, fit_bp_model, data = model_data,
                       spacing = sample_spacing, regime = sample_regime,
                       calibrate = calibrate, calibrate_plot = FALSE, 
                       subsample_seed = i)

  # validation ----
  print("validation - bbs data")

  bp_runs_bbs <- mutate(bp_runs, ppms = map(models, validate, data = ebird_test_bbs))
  ppm_bbs <- bp_runs_bbs %>% 
    select(-models) %>% 
    unnest(cols = c(ppms)) %>%
    mutate(tss = sensitivity + specificity - 1) %>%
    # select(run_id, mse:tss) %>%
    # rename_all(function(x){paste0(x, "_bbs")}) %>%
    # rename("run_id" = "run_id_bbs") %>%
    mutate(val_type = "bbs")

  print("validation - 2017 data")

  bp_runs_2017 <- mutate(bp_runs, ppms = map(models, validate, data = ebird_test_2017))
  ppm_bbs_2017 <- bp_runs_2017 %>% 
    select(-models) %>% 
    unnest(cols = c(ppms)) %>%
    mutate(tss = sensitivity + specificity - 1) %>%
    mutate(val_type = "2017") %>%
    rbind(ppm_bbs)

  print("validation - 2017 balanced data")

  bp_runs_2017_bal <- mutate(bp_runs, ppms = map(models, validate, data = ebird_test_2017_bal))
  ppms <- bp_runs_2017_bal %>% 
    select(-models) %>% 
    unnest(cols = c(ppms)) %>%
    mutate(tss = sensitivity + specificity - 1) %>%
    mutate(val_type = "2017_bal") %>%
    rbind(ppm_bbs_2017) %>%
    mutate(sim_id = i)

  if(i==1) all_ppms <- ppms
  if(i>1) all_ppms <- rbind(all_ppms, ppms)

}

str_glue("{output_folder}/07_bad_1_rf-model_assessment_multi_{sp_code}_{date}.png") %>% 
  write_csv(all_ppms, .)

# all_ppms <- str_glue("{output_folder}/07_bad_1_rf-model_assessment_multi_{sp_code}_{date}.png") %>% read_csv()


# 2017 validation plot

  # plot comparing ppms
  ppm_plot <- all_ppms %>% 
    select_if(~ !is.logical(.)) %>% 
    gather("metric", "value", -run_id, -run_name, -sim_id, -val_type) %>%
    filter(metric != "threshold", metric != "n_checklists", metric != "n_pos") %>% 
    arrange(val_type, sim_id, run_id) %>% 
    mutate(metric_label = ifelse(metric %in% c("auc", "mse", "tss"), str_to_upper(metric), str_to_title(metric))) %>%
    mutate(metric_label = factor(metric_label, levels = c("MSE", "AUC", "Kappa", "Sensitivity", "Specificity", "TSS"))) %>%
    mutate(run = paste("Model", run_id),
           run = as_factor(run),
           start = if_else(metric == "AUC", 0.5, 0)) %>%
    filter(!is.na(metric))

for(i in 1:3){

  val_type_plot <- c("bbs", "2017", "2017_bal")[i]
  plot_data <- filter(ppm_plot, val_type==val_type_plot)

  g_ppm <- ggplot(plot_data) +
    aes(x = run, y = value) +
    geom_boxplot(coef=5) +
    facet_wrap(~ metric_label, nrow = 2, scales = "free_y") +
    scale_y_continuous(breaks = c(0, 0.5, 1), limits=c(0, 1)) +
    labs(x = NULL, y = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          strip.text = element_text(size = 12, hjust = 0),
          panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white"),
          panel.border = element_rect(color = "black", fill = "transparent"),
          axis.ticks.y = element_line(),
          panel.grid = element_blank())
  plotname <- str_glue("{figure_folder}/07_bad_1_rf-model_assessment_multi_{val_type_plot}_{sp_code}.png") %>%
    ggsave(g_ppm, width = 20, height = 20, units = "cm", dpi = 300)

}

for(i in 1:3){

  val_type_plot <- c("bbs", "2017", "2017_bal")[i]
  plot_data <- filter(ppm_plot, val_type==val_type_plot)

  g_ppm <- ggplot(plot_data) +
    aes(x = run, y = value) +
    geom_boxplot(coef=5) +
    facet_wrap(~ metric_label, nrow = 2, scales = "free_y") +
    labs(x = NULL, y = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          strip.text = element_text(size = 12, hjust = 0),
          panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white"),
          panel.border = element_rect(color = "black", fill = "transparent"),
          axis.ticks.y = element_line(),
          panel.grid = element_blank())
  plotname <- str_glue("{figure_folder}/07_bad_1_rf-model_assessment_multi_{val_type_plot}_{sp_code}_unscaled.png") %>%
    ggsave(g_ppm, width = 20, height = 20, units = "cm", dpi = 300)

}

# --------------------------------------------------------------------
# plot differences

ppms_gather <- all_ppms %>%
  select_if(~ !is.logical(.)) %>% 
  select(-threshold) %>%
  gather("metric", "value", -run_id, -run_name, -sim_id, -val_type)

diff_plot <- ppms_gather %>%
  filter(run_id == anchor_model) %>%
  select(metric, sim_id, value, val_type) %>%
  rename(value_anchor = value) %>%
  right_join(ppms_gather) %>%
  mutate(diff = value - value_anchor) %>%
  arrange(val_type, sim_id, run_id) %>% 
  mutate(run = paste("Model", run_id),
         run = as_factor(run),
         start = if_else(metric == "AUC", 0.5, 0)) %>%
  filter(!is.na(metric))

## ggplot is the worst

add_grey <- TRUE

for(i in 1:3){

  val_type_plot <- c("bbs", "2017", "2017_bal")[i]
  plot_data <- diff_plot %>%
      filter(val_type==val_type_plot) %>%
      filter(! metric %in% c("n_checklists", "n_pos")) %>%
      mutate(metric_short = metric) %>%
      mutate(metric = ifelse(metric %in% c("auc", "mse", "tss"), str_to_upper(metric), str_to_title(metric))) %>%
      mutate(metric = factor(metric, levels = c("MSE", "AUC", "Kappa", "Sensitivity", "Specificity", "TSS")))

  maxy <- plot_data %>% 
            select(metric, diff) %>% 
            group_by(metric) %>%
            summarise(max_abs = max(abs(diff))) %>%
            ungroup()

  str_glue("{figure_folder}/07_bad_1_rf-model_assessment_multi_DIFF_{val_type_plot}_anchor{anchor_model}_{sp_code}_grey{add_grey}.png") %>%
    png(width = 21, height = 17, units="cm", res = 600)

        par(mfrow = c(2, 3), mar = c(1, 5, 1, 1), oma = c(4, 1, 1, 1))

        # performance metrics
        for (j in 1:6) {
          m <- levels(plot_data$metric)[j]
          maxyy <- maxy$max_abs[j]
          xnames <- rep("", 6)
          if(j>3) xnames <- paste("Model", 1:6)           
          boxplot(as.formula("diff ~ run"), 
                  data = plot_data[plot_data$metric==m,],
                  range = 0, boxwex = 0.8, lty = 1, staplewex = 0,
                  boxcol = "white", 
                  col = "white",
                  xlab = "", ylab = "",
                  ylim=c(-1*maxyy, maxyy), las = 2, names = xnames)
          if(!add_grey) abline(h=0, lwd=2, col="grey70")
          if(add_grey) {
            ymin <- -1
            ymax <- 0
            if(j==1) {ymin <- 0; ymax <- 1 }
            polygon(x = c(-1, 10, 10, -1, -1), y = c(ymin, ymin, ymax, ymax, ymin), col="grey78", border = alpha("white", 0))
          }
          par(new=TRUE)
          boxplot(as.formula("diff ~ run"), 
                  data = plot_data[plot_data$metric==m,],
                  range = 0, boxwex = 0.8, lty = 1, staplewex = 0,
                  xlab = "", ylim = c(-1*maxyy, maxyy), names = rep("", 6), 
                  col = alpha("white", 0.4), 
                  xaxt="n", yaxt="n", 
                  ylab = bquote(Delta~.(levels(plot_data$metric)[j])))
          text(x = 0.5, y = maxyy*0.95, 
               labels = LETTERS[j], 
               font = 2, pos=4)
        }

    dev.off()

}
