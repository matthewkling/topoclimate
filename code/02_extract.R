
# This script adds climate and topographic variables for each FIA subplot.
# inputs: fia_subplots.rds; and fia_plot_coordinates.csv with true plot coordinates
# output: fia_topoclimate.rds

library(tidyverse)
library(raster)
library(geosphere)
library(foreign)
library(furrr)
library(parallel)
library(conflicted)

conflict_prefer("select", "dplyr", quiet = T)
conflict_prefer("filter", "dplyr", quiet = T)
conflict_prefer("extract", "raster", quiet = T)
conflict_prefer("lag", "stats", quiet = T)

source("code/02_extract_functions.R")

coordinates_file <- "data/derived/fia_plot_coordinates.csv"


# LOAD FIA DATA #######################

sp <- readRDS("data/derived/fia_subplots.rds") # generated in 01_assemble.R



# PLOT COORDINATES #####################################

# stash a copy of fuzzed plot coordinates
fuzz <- sp %>% select(subplot_id, lon, lat)

# join to in true plot coordinates
sp <- sp %>% select(-lon, -lat)
sp <- read_csv(coordinates_file) %>%
  inner_join(sp) %>%
  filter(is.finite(lon))

# remove merging variables and select unique locations
sp <- sp %>% 
  select(plot_id, subplot_id, lon, lat, slope, aspect) %>%
  distinct()


# derive subplot coordinates 
sp <- sp %>% 
  select(plot_id, subplot_id, lon, lat) %>%
  distinct() %>%
  group_by(plot_id) %>% 
  filter(n() == 4) %>%
  arrange(subplot_id) %>%
  mutate(bearing = c(0, 0, 120, -120), # in degrees from N
         distance = c(0, 120, 120, 120) * .3048) %>% # feet to meters
  rowwise() %>%
  mutate(ll = destPoint(c(lon, lat), bearing, distance),
         lon_sp = ll[,1],
         lat_sp = ll[,2]) %>% 
  select(-ll, -lon, -lat, -bearing, -distance) %>%
  ungroup() %>%
  right_join(sp) %>%
  select(-lon, -lat) %>%
  rename(lon = lon_sp, lat = lat_sp) %>%
  filter(is.finite(lon))

# coordinate reference system for FIA
fia_crs <- "+init=epsg:4269 +proj=longlat +ellps=GRS80 +datum=NAD83"


# TOPOGRAPHY ###########################################

topo_ak <- sp %>%
  topography(veg_file = "data/raw/landfire/LF2016_EVT_200_AK/Tif/LA16_EVT_200.tif", 
             veg_dbf_file = "data/raw/landfire/LF2016_EVT_200_AK/Tif/LA16_EVT_200.tif.vat.dbf", 
             slope_file = "data/raw/landfire/LF2016_Slp_200_AK/Tif/LA16_Slp_200.tif", 
             aspect_file = "data/raw/landfire/LF2016_Asp_200_AK/Tif/LA16_Asp_200.tif", 
             elev_file = "data/raw/landfire/LF2016_Elev_200_AK/Tif/LA16_Elev_200.tif")

topo_conus <- sp %>%
  topography(veg_file = "data/raw/landfire/LF2016_EVT_200_CONUS/Tif/LC16_EVT_200.tif",
             veg_dbf_file = "data/raw/landfire/LF2016_EVT_200_CONUS/Tif/LC16_EVT_200.tif.vat.dbf",
             slope_file = "data/raw/landfire/LF2016_Slp_200_CONUS/Tif/LC16_Slp_200.tif",
             aspect_file = "data/raw/landfire/LF2016_Asp_200_CONUS/Tif/LC16_Asp_200.tif",
             elev_file = "data/raw/landfire/LF2016_Elev_200_CONUS/Tif/LC16_Elev_200.tif")

topo <- bind_rows(topo_conus, topo_ak)
# saveRDS(topo, "data/derived/topography.rds")

# merge with master subplot table
sp <- topo %>% distinct() %>% select(-lon, -lat) %>%
  left_join(sp, .) %>%
  filter(is.finite(northness), # drop plots without topography data (mainly wrong veg class)
         is.finite(tpi))



# MACROCLIMATE ##############

sp <- add_macroclimate(sp)



# WIND ###################################

sp <- add_wind(sp)



# FINAL STEPS ###############

# average over the multi-inventory variance in subplot topographic measurements
sp <- sp %>%
  select(-aspect_deg, -aspect_rad, -wdir) %>% # can't be averaged except by circular mean, and not currently used downstream
  group_by(plot_id, subplot_id) %>%
  summarize_all(mean, na.rm = T) %>%
  ungroup()

# swap out exact coordinates with fuzzed coordinates, for confidentiality
sp <- sp %>%
  select(-lon, -lat) %>%
  left_join(fuzz)

# export data
saveRDS(op, "data/derived/fia_topoclimate.rds")
