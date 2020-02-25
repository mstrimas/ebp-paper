# take the processed ebird data and separate into train, test1, test2
# train = 2018 non-bbs data
# test1 = 2018 bbs data
# test2 = 2017 non-bbs data


library(tidyverse)
library(lubridate)
library(rgdal)
library(sf)

# custom functions
walk(list.files("R", full.names = TRUE), source)

# folder to save training and validation datasets
ebd_save <- "data/"

# name of data extraction
data_tag <- "mayjune_201718_bcr27"

# ####################################################################
# READ IN PROCESSED EBIRD DATA

eb_zf_loc <- paste0(ebd_save, "ebd_", data_tag, "_zf.csv")
eb_zf <- read.csv(eb_zf_loc) %>%
          select(-species_observed) %>%
          spread(species_code, observation_count) %>%
          filter(yday(observation_date)>134)

# read in badly processed data
eb_bad_zf_loc <- paste0(ebd_save, "ebd_", data_tag, "_zf_bad.csv")
eb_bad_zf <- read.csv(eb_bad_zf_loc) %>%
          select(-species_observed) %>%
          spread(species_code, observation_count) %>%
          filter(yday(observation_date)>134)



# ####################################################################
# FIND AND FLAG SUSPECTED BBS COUNTS

# Some BBS observers add their observations into the eBird database

# identified by many 3min point counts by the same person on the same day
# also identified by several 3min point counts in the proximity

bbs_obs_date <- eb_zf %>%
        filter(protocol_type == "Stationary") %>%
        filter(duration_minutes == 3) %>% 
        group_by(observation_date, observer_id) %>%
        summarise(no_cl = n()) %>%
        ungroup() %>%
        filter(no_cl >= 40)

bbs_obs_date_loc <- eb_zf %>%
        filter(paste0(observation_date, observer_id) %in% paste0(bbs_obs_date$observation_date, bbs_obs_date$observer_id)) %>%
        filter(duration_minutes==3, all_species_reported==TRUE)

# split by observer and date
sp <- split(bbs_obs_date_loc, list(bbs_obs_date_loc$observer_id, bbs_obs_date_loc$observation_date), drop = TRUE)


# --------------------------------------------------------------------
# find nearest neighbour for each point (with the same observer and date)
# also find the number of neighbours within 2 miles, 4 miles, 6 miles of 
# the focal point

find_nearest_neighbour <- function(df){

  locs <- df[,c("longitude", "latitude")] %>%
            as.matrix()

  df_sf <- st_as_sf(sp::SpatialPoints(coords=locs, proj4string = CRS("+init=epsg:4326")))

  dist_mat <- st_distance(x=df_sf, y=df_sf)
  min_dist <- apply(dist_mat, 1, function(x){min(x[x>0])})
  number_near1 <- apply(dist_mat, 1, function(x){sum((x[x>0])<3300)})
  number_near2 <- apply(dist_mat, 1, function(x){sum((x[x>0])<6700)})
  number_near3 <- apply(dist_mat, 1, function(x){sum((x[x>0])<10000)})

  df$min_dist <- round(min_dist)
  df$number_near1 <- number_near1
  df$number_near2 <- number_near2
  df$number_near3 <- number_near3

  return(df)
}

nn <- lapply(sp, FUN=find_nearest_neighbour)
obs_date_nn <- bind_rows(nn)


# --------------------------------------------------------------------
# define which are bbs ones by min_dist and number of points within threshold

threshold1 <- 2
threshold2 <- 5
threshold3 <- 9

obs_date_nn$bbs <- ifelse(obs_date_nn$min_dist<2000 & obs_date_nn$number_near1 > (threshold1 - 1)  & obs_date_nn$number_near3 > (threshold3 - 1)  & obs_date_nn$number_near3 > (threshold3 - 1), 1, 0)

bbs <- obs_date_nn %>% filter(bbs==1)
not_bbs <- obs_date_nn %>% filter(bbs==0)

# --------------------------------------------------------------------
# plot out to visualise the points selected as bbs and not bbs

par(mfrow=c(2,1))

# estimated 3min point counts from bbs routes shown in red
plot(obs_date_nn$longitude, obs_date_nn$latitude, pch=16, cex=0.5,
#  xlim=c(-77, -76), ylim=c(36, 38))
  xlim=c(-86, -85), ylim=c(30, 32))
points(bbs$longitude, bbs$latitude, pch=16, col="red", cex=0.5)
points(not_bbs$longitude, not_bbs$latitude, pch=16, col="black", cex=0.5)

# estimated 3min point counts NOT from bbs routes shown in red
plot(obs_date_nn$longitude, obs_date_nn$latitude, pch=".",
  xlim=c(-86, -85), ylim=c(30, 32))
points(bbs$longitude, bbs$latitude, pch=16, col="black", cex=0.5)
points(not_bbs$longitude, not_bbs$latitude, pch=16, col="red", cex=0.5)


# filter to only estimated bbs stops 
bbs_for_merge <- obs_date_nn %>% 
          filter(bbs==1) %>%
          mutate(checklist_id = as.character(checklist_id)) %>%
          select(checklist_id, bbs) #%>%



# --------------------------------------------------------------------
# merge back in with main dataset

eb_zf_all <- eb_zf %>%
            mutate(checklist_id = as.character(checklist_id)) %>%
            left_join(bbs_for_merge) %>%
            mutate(bbs = ifelse(is.na(bbs), 0, 1)) %>%
            mutate(type = case_when(
                  bbs==0 & year(observation_date)==2018 ~ "train",
                  bbs==0 & year(observation_date)==2017 ~ "test_2017",
                  bbs==1 & year(observation_date)==2018 ~ "test_bbs",
                  TRUE                                   ~ "other"))



# ####################################################################
# FIND AND FLAG SUSPECTED BBS COUNTS WITHIN EBIRD DATASET

# Some BBS observers add their observations into the eBird database

# identified by many 3min point counts by the same person on the same day
# also identified by several 3min point counts in the proximity

eb_bad_zf_all <- eb_bad_zf %>%
        mutate(checklist_id = as.character(checklist_id)) %>%
        left_join(select(eb_zf_all, checklist_id, type)) %>%
        mutate(type = ifelse(is.na(type), 
                            ifelse(year(observation_date)==2018, "train", "other"), 
                            type))


# ####################################################################
# WRITE TO FILES

data_loc <- paste0(ebd_save, "data_4_models_", data_tag, ".csv")
write_csv(eb_zf_all, data_loc)

data_bad_loc <- paste0(ebd_save, "data_bad_4_models_", data_tag, ".csv")
write_csv(eb_bad_zf_all, data_bad_loc)



# ####################################################################
# SUMMARISE NUMBERS

# all data
nrow(eb_bad_zf_all)

# split by year
table(year(eb_bad_zf_all$observation_date))

# split by year and bbs
eb_bad_zf %>%
        mutate(checklist_id = as.character(checklist_id)) %>%
        left_join(bbs_for_merge) %>%
        mutate(bbs = ifelse(is.na(bbs), 0, 1)) %>%
        mutate(yr = year(observation_date)) %>%
        select(bbs, yr) %>%
        table()

test <- eb_bad_zf %>%
        mutate(checklist_id = as.character(checklist_id)) %>%
        left_join(bbs_for_merge) %>%
        mutate(bbs = ifelse(is.na(bbs), 0, 1)) %>%
        mutate(yr = year(observation_date)) %>%
        filter(yr == 2017, bbs == 0)



# ####################################################################
# PLOT DIFFERENT DATASETS

map_proj <- st_crs(102003)
# borders
f_gpkg <- paste0(ebd_save, "gis-data.gpkg")
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


# --------------------------------------------------------------------
# project the datasets

train_data <- eb_zf_all %>%
        filter(type=="train")

train_data_proj <- as.matrix(cbind(train_data$longitude, train_data$latitude)) %>%
        sp::SpatialPoints(proj4string = CRS("+init=epsg:4326")) %>% 
        st_as_sf() %>%
        st_transform(crs = map_proj)


test_data_bbs <- eb_zf_all %>%
        filter(type=="test_bbs")

test_data_bbs_proj <- as.matrix(cbind(test_data_bbs$longitude, test_data_bbs$latitude)) %>%
        sp::SpatialPoints(proj4string = CRS("+init=epsg:4326")) %>% 
        st_as_sf() %>%
        st_transform(crs = map_proj)


test_data_2017 <- eb_zf_all %>%
        filter(type=="test_2017")

test_data_2017_proj <- as.matrix(cbind(test_data_2017$longitude, test_data_2017$latitude)) %>%
        sp::SpatialPoints(proj4string = CRS("+init=epsg:4326")) %>% 
        st_as_sf() %>%
        st_transform(crs = map_proj)


figure_folder <- "figures/"
plot_name <- paste0(figure_folder, "train_test_maps_", Sys.Date(), ".png")
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

# ####################################################################
# FIND THE DATES OF BBS SURVEYS

hist(yday(test_data_bbs$observation_date))
min(yday(test_data_bbs$observation_date))
# day 138 = may 18

min(yday(train_data$observation_date))








