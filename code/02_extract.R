
# This script adds climate and topographic variables for each FIA subplot.
# inputs: fia_subplots.rds, and a TBD dataset with true plot coordinates
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



# LOAD FIA DATA #######################

sp <- readRDS("data/derived/fia_subplots.rds") # generated in fia_assemble.R



# PLOT COORDINATES #####################################

# stash a copy of fuzzed plot coordinates
fuzz <- sp %>% select(subplot_id, lon, lat)

# ** TO DO: swap in true plot coordinates here **


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

topography <- function(sp, veg_file, veg_dbf_file, slope_file, aspect_file, elev_file){
  
  ### veg type ###
  
  # (used to select non-sampled FIA sites that are not forest or developed)
  veg <- raster(veg_file)
  
  # prepare subplot locations for extraction
  low_lat <- grepl("CONUS", veg_file)
  e <- filter(sp, (lat < 50) == low_lat, is.finite(lat))
  coordinates(e) <- c("lon", "lat")
  crs(e) <- fia_crs
  e <- spTransform(e, crs(veg))
  
  # extract vegetation type at subplots (SLOW)
  message("... extracting vegetation type ...")
  e$veg <- extract(veg, e)
  
  
  ### slope and aspect for non-forest sites ###
  
  # load lLF slope/aspect data (slope and aspect are both in degrees)
  slope <- raster(slope_file)
  aspect <- raster(aspect_file)
  NAvalue(slope) <- -9999
  NAvalue(aspect) <- -9999
  
  # identify sites with wild, non-forest vegetation
  nonforest_veg <- c("Herb", "Shrub", "Sparse", "Barren", "Snow-Ice")
  nonforest <- read.dbf(veg_dbf_file) %>%
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
    select(-slope_pct_lf)
  
  # merge FIA topography with landfire topography
  topo <- sp %>%
    filter((lat < 50) == low_lat, is.finite(lat)) %>%
    select(-lon, -lat) %>%
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
      aspect_rad = aspect_deg/360*2*pi,
      slope_rad = slope_deg/360*2*pi) # not redundant -- have to recalculate post-merge
  
  # calculate northness and eastness
  topo <- topo %>% mutate(northness = cos(aspect_rad) * sin(slope_rad),
                          eastness = sin(aspect_rad) * sin(slope_rad)) %>%
    filter(is.finite(northness),
           is.finite(lon)) %>%
    select(-slope_deg_lf, -aspect_deg_lf)
  
  
  ### TPI ###
  # this variable is not available elsewhere at 30m that I can find,
  # and is too computationally burdensome to compute wall-to-wall,
  # so rather than extracting it from a raster we are computing it here
  # for just the neighborhood around each plot, based on an elevation raster.
  message("... computing TPI ...")
  
  # mTPI neighborhood radii, in meters
  radii <- c(100, 225, 500)
  
  # load elevation data
  dem <- raster(elev_file)# %>% extend(50)
  NAvalue(dem) <- -9999
  
  # precalculate focal matrices for each radius
  fws <- list()
  for(r in radii) fws[[paste0("r", r)]] <- focalWeight(dem, r, "circle")
  
  # compute topographic position index (TPI) for a given point location
  # (standardized multi-scale TPI -- mTPIs sensu Theobald et al. (2015))
  mtpi <- function(x, y, radii, dem){
    
    # neighborhood elevation data for the largest radius
    rr <- rev(sort(radii))
    ext <- extent(x - rr[1], x + rr[1], y - rr[1], y + rr[1]) %>% extend(30)
    dc <- crop(dem, ext)
    e <- as.matrix(dc)
    
    ids <- dc
    ids[] <- 1:ncell(ids)
    id <- extract(ids, matrix(c(x, y), 1))
    cl <- id %% ncol(dc)
    rw <- (id - cl) / ncol(dc) + 1
    
    tpi <- rep(NA, length(rr))
    for(i in 1:length(rr)){
      r <- rr[i]
      fw <- fws[[paste0("r", r)]]
      rad <- (dim(fw)[1] - 1) / 2
      
      er <- e[(rw - rad):(rw + rad), 
              (cl - rad):(cl + rad)]
      # if(extract(dc, matrix(c(x, y), 1)) != er[rad + 1, rad + 1]) stop("alignment issues 1")
      er <- as.vector(er)
      
      # if(length(er) != length(fw)) stop("alignment issues 2")
      er[fw == 0] <- NA
      
      # TPIs
      eo <- er[(length(er) + 1) / 2] # elev of central cell
      er <- na.omit(er[er != 0])
      le <- length(er)
      en <- sum(er) / le # mean elevation of neighborhood
      es <- sqrt(sum((er - en)^2) / le) # standard deviation of elevation
      tpi[i] <- (eo - en) / es # z-score
    }
    mean(tpi, na.rm = T) # average TPI across radii
  }
  
  # coordinates for extraction
  e <- topo %>% select(subplot_id, lon, lat)
  coordinates(e) <- c("lon", "lat")
  crs(e) <- crs(dem)
  e <- crop(e, dem) %>% as.data.frame()
  
  # compute mTPI in parallel
  plan(multisession, workers = detectCores() - 1)
  e$tpi <- e %>%
    select(x = lon, y = lat) %>%
    future_pmap_dbl(., possibly(mtpi, NA_real_), r = radii, dem = dem)
  topo <- e %>% select(subplot_id, tpi) %>% right_join(topo)
  
  return(topo)
}

topo_conus <- sp %>%
  topography(veg_file = "data/landfire/LF2016_EVT_200_CONUS/Tif/LC16_EVT_200.tif", 
             veg_dbf_file = "data/landfire/LF2016_EVT_200_CONUS/Tif/LC16_EVT_200.tif.vat.dbf", 
             slope_file = "data/landfire/LF2016_Slp_200_CONUS/Tif/LC16_Slp_200.tif", 
             aspect_file = "data/landfire/LF2016_Asp_200_CONUS/Tif/LC16_Asp_200.tif", 
             elev_file = "data/landfire/LF2016_Elev_200_CONUS/Tif/LC16_Elev_200.tif")

topo_ak <- sp %>%
  topography(veg_file = "data/landfire/LF2016_EVT_200_AK/Tif/LA16_EVT_200.tif", 
             veg_dbf_file = "data/landfire/LF2016_EVT_200_AK/Tif/LA16_EVT_200.tif.vat.dbf", 
             slope_file = "data/landfire/LF2016_Slp_200_AK/Tif/LA16_Slp_200.tif", 
             aspect_file = "data/landfire/LF2016_Asp_200_AK/Tif/LA16_Asp_200.tif", 
             elev_file = "data/landfire/LF2016_Elev_200_AK/Tif/LA16_Elev_200.tif")

topo <- bind_rows(topo_conus, topo_ak)
# saveRDS(topo, "data/derived/topography.rds")

# merge with master subplot table
sp <- topo %>% distinct() %>% select(-lon, -lat) %>%
  left_join(sp, .) %>%
  filter(is.finite(northness), # drop plots without topography data (mainly wrong veg class)
         is.finite(tpi))



# MACROCLIMATE ##############

# load bioclimatic variables from CHELSA, extract, and add to dataset
clim <- list.files("data/raw/chelsa", full.names = T) %>%
  stack() %>%
  setNames(paste0("bio", c(1, 12, 5, 6)))

s <- sp %>% select(subplot_id, lon, lat) %>% as.data.frame()
coordinates(s) <- c("lon", "lat")
crs(s) <- fia_crs
s <- spTransform(s, crs(clim))
sp <- cbind(sp, extract(clim, s))

# log precip for normality
log10p1 <- function(x) log10(x+1)
sp <- sp %>% 
  filter(is.finite(bio1)) %>%
  mutate_at(vars(bio12), log10p1)



# WIND ###################################

# load wind data, extract, and add to dataset
# (data is from Kling & Ackerly 2020 Nature Climate Change)
w <- stack("data/raw/wind/regimes.tif")
names(w) <- c("wspeed", "wdir", "waniso", "wu", "wv")
w <- rotate(w)
s <- sp %>% select(subplot_id, lon, lat)
coordinates(s) <- c("lon", "lat")
crs(s) <- CRS("+proj=longlat +datum=NAD27")
s <- spTransform(s, crs(w))
w <- raster::extract(w, s)
sp <- cbind(sp, w)

# calculate "windwardness": how strongly plot faces into prevailing wind
# (depends on plot slope, plot aspect vs wind direction, wind anisotropy)
anglediff <- function(x, y){
  z <- abs(x - y)
  ifelse(z > pi, abs(2*pi - z), z)
}
sp <- sp %>%
  ungroup() %>%
  mutate(windiff = pi - anglediff(atan2(wv, wu), # downwind direction, 
                                  atan2(northness, eastness)), # downslope direction
         windward = cos(windiff) * sin(slope_rad) * waniso / mean(waniso))



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
