
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

coordinates_file <- "G:/KatieB/Collaborations/Berkeley Tree Distributions & Climate Change/M_Kling/Actual_Coordinates_Request/Coords_from_Andy/AllPlots_Min2.csv"


# LOAD FIA DATA #######################

sp <- readRDS("data/derived/fia_subplots.rds") # generated in 01_assemble.R



# PLOT COORDINATES #####################################

# stash a copy of fuzzed plot coordinates
fuzz <- sp %>% select(subplot_id, lon, lat)

# I don't have (and I don't think we actually need) the UNITCD, so I am 
# going to overwrite the plot_id variable without UNITCD so I can match 
# my coordinates set to sp by the updated plot_id variable. I'll replace 
# with the original plot id values after completing the merge.
plot_id_components<-strsplit(sp$plot_id, split=" ")
pid_1<-sapply(plot_id_components, "[[",1)
pid_3<-sapply(plot_id_components, "[[",3)
pid_4<-sapply(plot_id_components, "[[",4)
pid.new<-paste(pid_1, pid_3, pid_4)

# sub in the new plot id values for the original ones 
# (to allow join operation below) but retain originals
# for reinsertion later:
sp$plot_id.orig<-sp$plot_id
sp$plot_id<-pid.new

# join to in true plot coordinates
sp <- sp %>% select(-lon, -lat)
sp <- read_csv(coordinates_file) %>%
        inner_join(sp) %>%
        filter(is.finite(lon))

# remove merging variables and select unique locations
# (some plots w multiple slope/aspect measurements from resurveys are included twice)
sp <- sp %>% 
        select(plot_id.orig, subplot_id, lon, lat, slope, aspect) %>%
        distinct()
colnames(sp)[1]<-'plot_id'

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
# saveRDS(sp, "data/derived/fia_subplots_coords.rds")

# coordinate reference system for FIA
fia_crs <- "+init=epsg:4269 +proj=longlat +ellps=GRS80 +datum=NAD83"


# TOPOGRAPHY ###########################################

############################################
##
# Calculate/Extract AK Topography Data:
### veg type ###

# (used to select non-sampled FIA sites that are not forest or developed)
veg <- raster("G:/KatieB/Collaborations/Berkeley Tree Distributions & Climate Change/M_Kling/Spatial_Data/TPI_Calc/LA16_EVT_200.tif")

# prepare subplot locations for extraction
low_lat <- grepl("CONUS", "G:/KatieB/Collaborations/Berkeley Tree Distributions & Climate Change/M_Kling/Spatial_Data/TPI_Calc/LA16_EVT_200.tif")
e <- filter(sp, (lat > 50), is.finite(lat))
coordinates(e) <- c("lon", "lat")
crs(e) <- fia_crs
e <- spTransform(e, crs(veg))

# extract vegetation type at subplots (SLOW)
message("... extracting vegetation type ...")
e$veg <- extract(veg, e)


### slope and aspect for non-forest sites ###

# load LF slope/aspect data (slope and aspect are both in degrees)
slope <- raster("G:/KatieB/Collaborations/Berkeley Tree Distributions & Climate Change/M_Kling/Spatial_Data/TPI_Calc/LA16_Slp_200.tif")
aspect <- raster("G:/KatieB/Collaborations/Berkeley Tree Distributions & Climate Change/M_Kling/Spatial_Data/TPI_Calc/LA16_Asp_200.tif")
NAvalue(slope) <- -9999
NAvalue(aspect) <- -9999

# identify sites with wild, non-forest vegetation
nonforest_veg <- c("Herb", "Shrub", "Sparse", "Barren", "Snow-Ice")
nonforest <- read.dbf("G:/KatieB/Collaborations/Berkeley Tree Distributions & Climate Change/M_Kling/Spatial_Data/TPI_Calc/LA16_Slp_200.tif") %>%
  as_tibble() %>%
  filter(EVT_LF %in% nonforest_veg)
nf <- e[e$veg %in% nonforest$Value,]

# extract slope and aspect
message("... extracting slope ...")
nf$slope <- extract(slope, nf)
message("... extracting aspect ...")
nf$aspect <- extract(aspect, nf)

# for LF variables, deal with LF and FIA encoding quirks for shallow slopes
nft <- nf %>% as.data.frame() %>% as_tibble() %>%
  select(subplot_id, lon, lat, 
         slope_deg_lf = slope, 
         aspect_deg_lf = aspect) %>%
  mutate(# in LF, "non-defined aspect (slope is less than or =2) are assigned a value of -1"
    slope_deg_lf = ifelse(slope_deg_lf <= 2, 0, slope_deg_lf),
    aspect_deg_lf = ifelse(slope_deg_lf == 0, 0, aspect_deg_lf),
    
    # replicate FIA encoding quirks for LF variables
    slope_pct_lf = tan(slope_deg_lf * pi / 180) * 100, 
    slope_deg_lf = ifelse(slope_pct_lf < 5, 0, slope_deg_lf)) %>%
  select(-slope_pct_lf) %>%
  mutate(nonforest = TRUE)

# merge FIA topography with landfire topography
topo <- sp %>%
  filter((lat < 50) == low_lat, is.finite(lat)) %>%
  # select(-lon, -lat) %>% # this is introducing an error downstream
  full_join(nft) %>%
  rename(slope_pct = slope,
         aspect_deg = aspect) %>%
  mutate(# slopes below 5% have aspect set to 360; set these to zero to avoid biasing exposure
    slope_pct = ifelse(slope_pct < 5, 0, slope_pct),
    slope_rad = atan(slope_pct/100),
    slope_deg = slope_rad * 180 / pi,
    
    # merge with nonforest topo, at stage when units match
    slope_deg = ifelse(is.na(slope_deg), slope_deg_lf, slope_deg),
    aspect_deg = ifelse(is.na(aspect_deg), aspect_deg_lf, aspect_deg),
    aspect_rad = aspect_deg / 180 * pi,
    slope_rad = slope_deg / 180 * pi) # not redundant -- have to recalculate post-merge

# calculate northness and eastness
topo <- topo %>% mutate(northness = cos(aspect_rad) * sin(slope_rad),
                        eastness = sin(aspect_rad) * sin(slope_rad)) %>%
  filter(is.finite(northness),
         is.finite(lon)) %>%
  select(-slope_deg_lf, -aspect_deg_lf)


### TPI ###
message("... computing TPI ...")

# mTPI neighborhood radii, in meters
radii <- c(100, 225, 500)

# load elevation data
dem <- raster("G:/KatieB/Collaborations/Berkeley Tree Distributions & Climate Change/M_Kling/Spatial_Data/TPI_Calc/LA16_Elev_200.tif")# %>% extend(50)
NAvalue(dem) <- -9999

# precalculate focal matrices for each radius
fws <- list()
for(r in radii) fws[[paste0("r", r)]] <- focalWeight(dem, r, "circle")

# coordinates for extraction
e <- sp %>% select(subplot_id, lon, lat)
coordinates(e) <- c("lon", "lat")
crs(e) <- fia_crs
e <- spTransform(e, crs(dem))
e <- crop(e, dem) %>% as.data.frame()

# compute mTPI in parallel
options(future.rng.onMisuse="ignore")
plan(multisession, workers = detectCores() - 1)
tpi <- e %>%
  select(x = lon, y = lat) %>%
  future_pmap(., possibly(mtpi, NA_real_), r = radii, dem = dem, fws = fws)
e$tpi <- map_dbl(tpi, "mtpi")
e$tpis <- map_dbl(tpi, "mtpis")
topo_ak <- e %>% select(subplot_id, tpi, tpis) %>% right_join(topo)

# saveRDS(topo_ak, "data/derived/topography_ak.rds")


############################################
##
# Calculate/Extract CONUS Topography Data:

### veg type ###

# (used to select non-sampled FIA sites that are not forest or developed)
veg <- raster("G:/KatieB/Collaborations/Berkeley Tree Distributions & Climate Change/M_Kling/Spatial_Data/TPI_Calc/LC16_EVT_200.tif")

# prepare subplot locations for extraction
low_lat <- grepl("CONUS", "G:/KatieB/Collaborations/Berkeley Tree Distributions & Climate Change/M_Kling/Spatial_Data/TPI_Calc/LC16_EVT_200.tif")
e <- filter(sp, lat <= 50, is.finite(lat))
coordinates(e) <- c("lon", "lat")
crs(e) <- fia_crs
e <- spTransform(e, crs(veg))

# extract vegetation type at subplots (SLOW)
message("... extracting vegetation type ...")
e$veg <- extract(veg, e)


### slope and aspect for non-forest sites ###

# load LF slope/aspect data (slope and aspect are both in degrees)
slope <- raster("G:/KatieB/Collaborations/Berkeley Tree Distributions & Climate Change/M_Kling/Spatial_Data/TPI_Calc/LC16_Slp_200.tif")
aspect <- raster("G:/KatieB/Collaborations/Berkeley Tree Distributions & Climate Change/M_Kling/Spatial_Data/TPI_Calc/LC16_Asp_200.tif")
NAvalue(slope) <- -9999
NAvalue(aspect) <- -9999

# identify sites with wild, non-forest vegetation
nonforest_veg <- c("Herb", "Shrub", "Sparse", "Barren", "Snow-Ice")
nonforest <- read.dbf("G:/KatieB/Collaborations/Berkeley Tree Distributions & Climate Change/M_Kling/Spatial_Data/TPI_Calc/LC16_EVT_200.tif.vat.dbf") %>%
  as_tibble() %>%
  filter(EVT_LF %in% nonforest_veg)
nf <- e[e$veg %in% nonforest$Value,]

# extract slope and aspect
message("... extracting slope ...")
nf$slope <- extract(slope, nf)
message("... extracting aspect ...")
nf$aspect <- extract(aspect, nf)

# for LF variables, deal with LF and FIA encoding quirks for shallow slopes
nft <- nf %>% as.data.frame() %>% as_tibble() %>%
  select(subplot_id, lon, lat, 
         slope_deg_lf = slope, 
         aspect_deg_lf = aspect) %>%
  mutate(# in LF, "non-defined aspect (slope is less than or =2) are assigned a value of -1"
    slope_deg_lf = ifelse(slope_deg_lf <= 2, 0, slope_deg_lf),
    aspect_deg_lf = ifelse(slope_deg_lf == 0, 0, aspect_deg_lf),
    
    # replicate FIA encoding quirks for LF variables
    slope_pct_lf = tan(slope_deg_lf * pi / 180) * 100, 
    slope_deg_lf = ifelse(slope_pct_lf < 5, 0, slope_deg_lf)) %>%
  select(-slope_pct_lf) %>%
  mutate(nonforest = TRUE)

# merge FIA topography with landfire topography
topo <- sp %>%
  filter(lat < 50, is.finite(lat))%>%
  # select(-lon, -lat) %>% # this is introducing an error downstream
  full_join(nft) %>%
  rename(slope_pct = slope,
         aspect_deg = aspect) %>%
  mutate(# slopes below 5% have aspect set to 360; set these to zero to avoid biasing exposure
    slope_pct = ifelse(slope_pct < 5, 0, slope_pct),
    slope_rad = atan(slope_pct/100),
    slope_deg = slope_rad * 180 / pi,
    
    # merge with nonforest topo, at stage when units match
    slope_deg = ifelse(is.na(slope_deg), slope_deg_lf, slope_deg),
    aspect_deg = ifelse(is.na(aspect_deg), aspect_deg_lf, aspect_deg),
    aspect_rad = aspect_deg / 180 * pi,
    slope_rad = slope_deg / 180 * pi) # not redundant -- have to recalculate post-merge

# calculate northness and eastness
topo <- topo %>% mutate(northness = cos(aspect_rad) * sin(slope_rad),
                        eastness = sin(aspect_rad) * sin(slope_rad)) %>%
  filter(is.finite(northness),
         is.finite(lon)) %>%
  select(-slope_deg_lf, -aspect_deg_lf)

# saveRDS(topo, "data/derived/topography_conus_pre_tpi_calc.rds")

### TPI ###
message("... computing TPI ...")

# mTPI neighborhood radii, in meters
radii <- c(100, 225, 500)

# load elevation data
dem <- raster("G:/KatieB/Collaborations/Berkeley Tree Distributions & Climate Change/M_Kling/Spatial_Data/TPI_Calc/LC16_Elev_200.tif")# %>% extend(50)
NAvalue(dem) <- -9999

# precalculate focal matrices for each radius
fws <- list()
for(r in radii) fws[[paste0("r", r)]] <- focalWeight(dem, r, "circle")

# coordinates for extraction
e <- sp %>% select(subplot_id, lon, lat)
coordinates(e) <- c("lon", "lat")
crs(e) <- fia_crs
e <- spTransform(e, crs(dem))
e <- crop(e, dem) %>% as.data.frame()

# compute mTPI in parallel
options(future.rng.onMisuse="ignore")
plan(multisession, workers = detectCores() - 1)
tpi <- e %>%
  select(x = lon, y = lat) %>%
  future_pmap(., possibly(mtpi, NA_real_), r = radii, dem = dem, fws = fws)
e$tpi <- map_dbl(tpi, "mtpi")
e$tpis <- map_dbl(tpi, "mtpis")
topo_conus <- e %>% select(subplot_id, tpi, tpis) %>% right_join(topo)
saveRDS(topo_conus, "data/derived/topography_conus.rds")

# Bind the AK & CONUS datasets together:
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
saveRDS(sp, "data/derived/fia_topoclimate.rds")
