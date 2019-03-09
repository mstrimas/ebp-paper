# gam abundance models

library(auk)
library(sf)
library(raster)
library(dggridR)
library(mgcv)
library(viridis)
library(fields)
library(dplyr)
library(purrr)
library(readr)
library(ggplot2)
library(stringr)
library(lubridate)
library(forcats)
walk(list.files("R", full.names = TRUE), source)
# resolve namespace conflicts
select <- dplyr::select
projection <- raster::projection

set.seed(1)
# set species for analysis
species <- "Wood Thrush"
sp_code <- ebird_species(species, "code")
# setup spatial sampling regime
sample_regime <- "both"
sample_spacing <- 5


# load data ----

# ebird data
ebird <- read_csv("data/ebd_june_bcr27_zf.csv", na = "") %>% 
  filter(species_code == sp_code,
         !is.na(observation_count))

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
  ebird_habitat <- inner_join(ebird_habitat, cci, by = "checklist_id") %>% 
    filter(!is.na(checklist_calibration_index))
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


# spatial subsampling ----

ebird_ss <- hex_sample(ebird_habitat, spacing = sample_spacing,
                       regime = sample_regime)


# test/train split ----

# prepare data
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
ebird_split <- ebird_ss %>% 
  mutate(day_of_year = yday(observation_date),
         protocol_type = factor(protocol_type, 
                                levels = c("Stationary" , "Traveling"))) %>%
  # select only the columns to be used in the model
  select(observation_count,
         day_of_year, time_observations_started, duration_minutes,
         effort_distance_km, number_observers, protocol_type,
         contains("checklist_calibration_index"),
         one_of(plands)) %>% 
  # test/train split, 80/20
  split(if_else(runif(nrow(.)) <= 0.8, "train", "test"))


# fit gam models ----

# gam parameters
# degrees of freedom for smoothing
k <- 5
# degrees of freedom for cyclic time of day smooth
k_time <- 7 

# continuous predictors
# hold out time to treat seperately since it's cyclic
continuous_covs <- ebird_split$train %>% 
  select(-observation_count, -protocol_type, -time_observations_started) %>% 
  names()

# create model formula for predictors
gam_formula_rhs <- str_glue("s({var}, k = {k})", 
                            var = continuous_covs, k = k) %>% 
  str_flatten(collapse = " + ") %>% 
  str_glue(" ~ ", .,
           " + protocol_type + ",
           "s(time_observations_started, bs = \"cc\", k = {k})", 
           k = k_time) %>% 
  as.formula()

# model formula including response
gam_formula <- update.formula(observation_count ~ ., gam_formula_rhs)

# explicitly specify where the knots should occur for time_observations_started
# this ensures that the cyclic spline joins the variable at midnight
# this won't happen by default if there are no data near midnight
time_knots <- list(time_observations_started = seq(0, 24, length.out = k_time))

# zero-inflated poisson
m_ziplss <- gam(list(gam_formula, # count model
                     gam_formula_rhs), # presence model
                data = ebird_split$train, 
                family = "ziplss", 
                knots = time_knots)

# negative binomial
m_nb <- gam(gam_formula,
            data = ebird_split$train, 
            family = "nb",
            knots = time_knots)

# tweedie distribution
m_tw <- gam(gam_formula,
            data = ebird_split$train, 
            family = "tw",
            knots = time_knots)

# validation ----

obs_count <- select(ebird_split$test, obs = observation_count)

# presence probability is on the complimentary log-log scale
# when can get the inverse link function with
inv_link <- binomial(link = "cloglog")$linkinv
m_ziplss_pred <- predict(m_ziplss, ebird_split$test, type = "link") %>% 
  as.data.frame() %>%
  transmute(family = "Zero-inflated Poisson",
            pred = inv_link(V2) * exp(V1)) %>% 
  bind_cols(obs_count)

m_nb_pred <- predict(m_nb, ebird_split$test, type = "response") %>% 
  tibble(family = "Negative Binomial", pred = .) %>% 
  bind_cols(obs_count)

m_tw_pred <- predict(m_tw, ebird_split$test, type = "response") %>% 
  tibble(family = "Tweedie", pred = .) %>% 
  bind_cols(obs_count)

# combine predictions
test_pred <- bind_rows(m_ziplss_pred, m_nb_pred, m_tw_pred) %>% 
  mutate(family = as_factor(family))

# deviance explained
de <- c(summary(m_ziplss)$dev.expl,
        summary(m_nb)$dev.expl,
        summary(m_tw)$dev.expl)

# mse and rank correlation
gam_assessment <- test_pred %>% 
  group_by(family) %>% 
  summarise(pct_prob_pred = mean(obs / pred > 10),
            mse = mean((obs - pred)^2, na.rm = TRUE),
            rank_cor = cor.test(obs, pred, 
                                method = "spearman", 
                                exact = FALSE)$estimate) %>% 
  ungroup() %>% 
  mutate(dev_exp = de)
str_glue("output/06_gam-count-model_assessment_{sp_code}.csv") %>% 
  write_csv(gam_assessment, .)


# model selection ----

# plot predicted vs. observed
ticks <- c(0, 1, 10, 100, 1000)
mx <- round(max(test_pred$obs))
g <- ggplot(test_pred) +
  aes(x = log10(obs + 1), 
      y = log10(pred + 1)) +
  geom_jitter(alpha = 0.2, height = 0) +
  # y = x line
  geom_abline(slope = 1, intercept = 0, alpha = 0.5) +
  # area where counts off by a factor of 10
  geom_area(data = tibble(x = log10(seq(0, mx - 1) + 1), 
                          y = log10(seq(0, mx - 1) / 10 + 1)),
            mapping = aes(x = x, y = y),
            fill = "red", alpha = 0.2) +
  # loess fit
  geom_smooth(method = "loess", 
              method.args = list(span = 2 / 3, degree = 1)) +
  scale_x_continuous(breaks = log10(ticks + 1), labels = ticks) +
  scale_y_continuous(breaks = log10(ticks + 1), labels = ticks) +
  labs(x = "Observed count",
       y = "Predicted count") +
  facet_wrap(~ family, nrow = 1)
str_glue("figures/06_gam-count-model_model-comparison_{sp_code}.png") %>% 
  ggsave(g)

# select model
if (sp_code == "woothr") {
  # negative binomial has highest deviance explained and equal lowest mse 
  # visually some predictions are a bit higher for abundance values
  # tweedie and negative binomial have low number of problematic errors 
  # (true count / predicted count ) > 10
  pred_model <- m_nb
} else if (sp_code == "norbob") {
  # tweedie has lowest mse, although negative binomial is very close
  # negative binomial has the highest proportion of deviance explained
  # tweedie and negative binomial have low number of problematic errors 
  # (true count / predicted count ) > 10
  # many real zero counts are modelled to have abundance. 
  pred_model <- m_nb
} else if (sp_code == "whiibi") {
  # tweedie model has lowest mse, equal highest deviance explained and 
  # low number of problematic errors: (true count/predicted count) > 10
  # however, high counts are consistently underestimated 
  # many real zero counts are predicted with relatively high abundance
  # none of the models have high percent deviance explained, suggesting the 
  # high counts are hard to predict
  # visually negative binomial and tweedie look like a similar fit.
  pred_model <- m_tw
} else {
  stop("species code not valid.")
}


# predict ----

# find optimal time of day to maximimise detectability
# use lower confidence limit of prediction to avoid spurious relationship

# create a dataframe of covariates with a range of start times
seq_tod <- seq(0, 24, length.out = 300)
tod_df <- ebird_split$train %>% 
  # find average pland habitat covariates
  select(starts_with("pland")) %>% 
  summarize_all(mean, na.rm = TRUE) %>% 
  ungroup() %>% 
  # use standard checklist
  mutate(day_of_year = yday(ymd("2016-06-15")),
         duration_minutes = 60,
         effort_distance_km = 1,
         number_observers = 1,
         checklist_calibration_index = 2,
         protocol_type = "Traveling") %>% 
  cbind(time_observations_started = seq_tod)

# predict at different start times
pred_tod <- predict(pred_model, newdata = tod_df, 
                    type = "link", 
                    se.fit = TRUE) %>% 
  as_tibble() %>% 
  # calculate backtransformed confidence limits
  transmute(time_observations_started = seq_tod,
            pred = pred_model$family$linkinv(fit),
            pred_lcl = pred_model$family$linkinv(fit - 1.96 * se.fit),
            pred_ucl = pred_model$family$linkinv(fit + 1.96 * se.fit))

# select detectability peak
t_max <- pred_tod$time_observations_started[which.max(pred_tod$pred_lcl)]

# add effort covariates to prediction surface
pred_surface_eff <- pred_surface %>% 
  mutate(day_of_year = yday(ymd("2016-06-15")),
         time_observations_started = t_max,
         duration_minutes = 60,
         effort_distance_km = 1,
         number_observers = 1,
         checklist_calibration_index = 2,
         protocol_type = "Traveling") %>% 
  as.data.frame()

# predict
abd_pred <- predict(pred_model, 
                    newdata = pred_surface_eff, 
                    type = "link", 
                    se.fit = TRUE) %>% 
  as_tibble() %>% 
  # calculate confidence limits and back transform
  transmute(abd = pred_model$family$linkinv(fit),
            abd_se = pred_model$family$linkinv(se.fit),
            abd_lcl = pred_model$family$linkinv(fit - 1.96 * se.fit),
            abd_ucl = pred_model$family$linkinv(fit + 1.96 * se.fit)) %>%
  # add to prediction surface
  bind_cols(pred_surface, .) %>% 
  select(latitude, longitude, abd, abd_se, abd_lcl, abd_ucl)

# rasterize
r <- raster("data/modis_5xagg.tif")
r_pred <- abd_pred %>% 
  select(abd, abd_se, latitude, longitude) %>% 
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(crs = projection(r)) %>% 
  rasterize(r)
r_pred <- r_pred[[-1]]

str_glue("output/06_gam-count-model_pred-abd_{sp_code}.tif") %>% 
  writeRaster(r_pred[["abd"]], ., overwrite = TRUE)
str_glue("output/06_gam-count-model_pred-se_{sp_code}.tif") %>% 
  writeRaster(r_pred[["abd_se"]], ., overwrite = TRUE)


# map predictions ----

zero_threshold <- 0.05

# project predictions
r_pred_proj <- projectRaster(r_pred[["abd"]], 
                             crs = map_proj$proj4string,
                             method = "ngb")

str_glue("figures/06_gam-count-model_predictions_{sp_code}.png") %>% 
  png(width = 2400, height = 1800, res = 300)
par(mar = c(4, 0.25, 0.25, 0.25))

# set up plot area
plot(bcr, col = NA, border = NA)
plot(ne_land, col = "#dddddd", border = "#888888", lwd = 0.5, add = TRUE)

# modified plasma palette
plasma_rev <- rev(plasma(25, end = 0.9))
gray_int <- colorRampPalette(c("#dddddd", plasma_rev[1]))
pal <- c(gray_int(4)[2], plasma_rev)

# set very low values to zero
r_pred_proj[r_pred_proj <= zero_threshold] <- NA
# log transform
r_pred_proj <- log10(r_pred_proj)
# breaks and legend
mx <- ceiling(100 * cellStats(r_pred_proj, max)) / 100
mn <- floor(100 * cellStats(r_pred_proj, min)) / 100
brks <- seq(mn, mx, length.out = length(pal) + 1)
lbl_brks <- sort(c(-2:2, mn, mx))
lbls <- round(10^lbl_brks, 2)

# abundance
plot(r_pred_proj, 
     col = pal, breaks = brks, 
     maxpixels = ncell(r_pred_proj),
     legend = FALSE, add = TRUE)

# borders
plot(bcr, border = "#000000", col = NA, lwd = 1, add = TRUE)
plot(ne_state_lines, col = "#ffffff", lwd = 0.75, add = TRUE)
plot(ne_country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
box()

# legend
par(new = TRUE, mar = c(0, 0, 0, 0))
image.plot(zlim = range(brks), legend.only = TRUE, col = pal,
           smallplot = c(0.25, 0.75, 0.06, 0.09),
           horizontal = TRUE,
           axis.args = list(at = lbl_brks, 
                            labels = lbls,
                            fg = "black", col.axis = "black",
                            cex.axis = 0.75, lwd.ticks = 0.5,
                            padj = -1.5),
           legend.args = list(text = str_glue("{species} Relative Abundance"),
                              side = 3, col = "black",
                              cex = 1, line = 0))
dev.off()