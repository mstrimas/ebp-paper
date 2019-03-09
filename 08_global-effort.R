# create map of global ebird effort
# for figure 1

library(auk)
library(data.table)
library(sf)
library(raster)
library(janitor)
library(viridis)
library(fields)
library(dplyr)
library(purrr)
library(stringr)
# resolve namespace conflicts
select <- dplyr::select
# custom functions
walk(list.files("R", full.names = TRUE), source)
# script variables
proj <- "+proj=eqearth +lon_0=-85 +wktext"
cell_res <- 25000


# ebird data ----

# global checklist location from ebd
f <- "data/ebd_global-checklists.txt"
if (!file.exists(f)) {
  ebd_locs <- auk_sampling("ebd_sampling_relDec-2018.txt") %>% 
    auk_select(file = f, overwrite = TRUE,
               select = c("sampling_event_identifier",
                          "locality_id", "latitude", "longitude",
                          "all_species_reported", "duration_minutes",
                          "protocol_type"))
}
ebd_locs <- fread(f) %>% 
  clean_names()


# calculate global effort ----

# summarize to effort at each unique location
effort <- ebd_locs[, 
                   .(n_all = .N, 
                     n_complete = sum(all_species_reported, na.rm = TRUE),
                     effort_hours = sum(all_species_reported * duration_minutes, 
                                        na.rm = TRUE) / 60),
                   by = .(locality_id, latitude, longitude)]

# convert to spatial
effort_sf <- effort %>% 
  select(-locality_id) %>% 
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(crs = proj)
saveRDS(effort_sf, "output/08_global-effort_sf.rds")

# create a raster template
r <- raster(xmn = -180, xmx = 180, ymn = -90, ymx = 90) %>% 
  projectRaster(crs = proj, res  = cell_res)
r_effort <- rasterize(effort_sf, r, fun = sum)
r_effort <- r_effort[[-1]]
for (i in 1:nlayers(r_effort)) {
  f <- names(r_effort)[i] %>% 
    str_replace("_", "-") %>% 
    sprintf("output/08_global-effort_%s.tif", .)
  writeRaster(r_effort[[i]], filename = f, overwrite = TRUE)
}


# map ----

# load map data
f_gpkg <- "data/gis-data.gpkg"
ne_bbox <- read_sf(f_gpkg, "ne_bbox") %>% 
  recenter_sf(crs = proj) %>% 
  st_geometry() %>% 
  st_set_precision(1e2) %>% 
  st_union()
ne_graticules <- read_sf(f_gpkg, "ne_graticules") %>% 
  recenter_sf(crs = proj) %>% 
  st_geometry()
ne_land <- read_sf(f_gpkg, "ne_land") %>% 
  recenter_sf(crs = proj) %>% 
  st_geometry()
ne_country_lines <- read_sf(f_gpkg, "ne_country_lines") %>% 
  recenter_sf(crs = proj) %>% 
  st_geometry()

# plot
layers <- c(n_all = "# of eBird checklists per 25 km grid cell",
            n_complete = "# of complete eBird checklists per 25 km grid cell",
            effort_hours = paste("eBird effort hours per 25 km grid cell",
                                 "(complete checklists)"))
for (i in seq_along(layers)) {
  x <- names(layers)[i] %>% 
    str_replace("_", "-") %>% 
    sprintf("output/08_global-effort_%s.tif", .) %>% 
    raster()
  names(layers)[i] %>% 
    str_replace("_", "-") %>% 
    sprintf("figures/08_global-effort_%s.png", .) %>% 
    png(width = 2200, height = 1200)
  par(mfrow = c(1, 1), mar = c(6, 0, 0, 0), bg = "#ffffff")
  
  # land background
  plot(ne_bbox, col = "#ffffff", border = "black", lwd = 3, 
       xaxt = "n", yaxt = "n", bty = "n")
  plot(ne_graticules, col = "#888888", lwd = 0.5, add = TRUE)
  plot(ne_land, col = "#dddddd", border = "#888888", lwd = 1, add = TRUE)
  
  # data
  pal <- viridis(25)
  mx <- cellStats(x, max) %>% ceiling()
  brks <- 10^seq(0, log10(mx), length.out = length(pal) + 1)
  lbl_brks <- seq(0, log10(mx), length.out = 5)
  if (names(layers)[i] == "effort_hours") {
    lbl_lbl <- label_hours(10^lbl_brks)
  } else {
    lbl_lbl <- signif(10^lbl_brks, 3) %>% 
      round() %>% 
      format(big.mark = ",") %>% 
      trimws()
  }
  plot(x, col = pal, breaks = brks, maxpixels = ncell(x),
       legend = FALSE, axes = FALSE, add = TRUE)
  
  # lines
  plot(ne_country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
  # bounding box
  plot(ne_bbox, col = NA, border = "black", lwd = 3, add = TRUE)
  
  # legend
  image.plot(zlim = range(lbl_brks), legend.only = TRUE, col = pal,
             smallplot = c(0.35, 0.65, 0.03, 0.05),
             horizontal = TRUE,
             axis.args = list(at = lbl_brks, labels = lbl_lbl,
                              fg = "black", col.axis = "black",
                              cex.axis = 2, lwd.ticks = 2),
             legend.args = list(text = layers[i],
                                side = 3, col = "black",
                                cex = 3, line = 1))
  dev.off()
}
