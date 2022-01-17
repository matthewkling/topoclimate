

# This script adds climate and topographic variables for each CNFI plot.

library(tidyverse)
library(raster)
# library(geosphere)
# library(foreign)
library(furrr)
library(parallel)
library(conflicted)

conflict_prefer("select", "dplyr", quiet = T)
conflict_prefer("filter", "dplyr", quiet = T)
conflict_prefer("extract", "raster", quiet = T)
# conflict_prefer("lag", "stats", quiet = T)

source("code/02_extract_functions.R")



# LOAD CNFI DATA #######################

sp <- readRDS("data/derived/cnfi_plots.rds") # generated in 01_assemble_canada.R
fia_crs <- "+init=epsg:4269 +proj=longlat +ellps=GRS80 +datum=NAD83"



# TOPOGRAPHY ###########################################

add_topo <- function(sp){
        
        # CNFI native topo fomatting:
        # slope is percent, aspect is degrees.
        # -1 means NA for both variables
        # 999 aspect reported for slopes <= 2%
        topo <- sp %>% 
                rename(slope_pct = slope,
                       aspect_deg = aspect) %>%
                mutate(slope_pct = ifelse(slope_pct == -1, NA, slope_pct),
                       aspect_deg = ifelse(aspect_deg == -1, NA, aspect_deg)) %>%
                mutate(slope_pct = ifelse(slope_pct <= 2, 0, slope_pct),
                       aspect_deg = ifelse(slope_pct == 0, 0, aspect_deg)) %>%
                mutate(slope_rad = atan(slope_pct / 100),
                       slope_deg = slope_rad * 180 / pi,
                       aspect_rad = aspect_deg / 360*2*pi) %>%
                select(-slope_pct)
        
        # calculate northness and eastness
        topo <- topo %>% mutate(northness = cos(aspect_rad) * sin(slope_rad),
                                eastness = sin(aspect_rad) * sin(slope_rad)) %>%
                filter(is.finite(northness),
                       is.finite(lon)) %>%
                select(-slope_deg, -aspect_deg)
        
        
        ## TPI ##
        
        # download elevation tile and compute mTPI for focal site
        get_tpi <- function(x = topo$lon[1],
                            y = topo$lat[1],
                            radii = c(100, 225, 500)){
                
                dem <- microclima::get_dem(lat = y, long = x, resolution = 30)
                fws <- list()
                for(r in radii) fws[[paste0("r", r)]] <- focalWeight(dem, r, "circle")
                
                p <- data.frame(x = x, y = y)
                coordinates(p) <- c("x", "y")
                crs(p) <- "+proj=longlat +ellps=GRS80 +datum=NAD83"
                p <- spTransform(p, crs(dem))
                p <- coordinates(p)
                
                mtpi(p[,1], p[,2], r, dem, fws)
        }
        
        options(future.rng.onMisuse="ignore")
        plan(multisession, workers = detectCores() - 1)
        tpi <- topo %>%
                select(x = lon, y = lat) %>%
                future_pmap(., possibly(get_tpi, 
                                        c(mtpi = NA_real_, mtpis = NA_real_)))
        topo$tpi <- map_dbl(tpi, "mtpi")
        topo$tpis <- map_dbl(tpi, "mtpis")
        
        return(topo)
}

sp <- sp %>% add_topo()



# MACROCLIMATE ##############

sp <- sp %>% 
        mutate(subplot_id = nfi_plot) %>%
        add_macroclimate()



# WIND ######################

sp <- sp %>% add_wind()



# EXPORT #####################

# select variables to export
sp <- sp %>% select(nfi_plot, lon:tpis, bio1:windward)

# obfuscate coordinates for privacy
sp <- sp %>%
        mutate(lon = lon + runif(nrow(.), -.005, .005),
               lat = lat + runif(nrow(.), -.005, .005))

# export
saveRDS(sp, "data/derived/cnfi_topoclimate.rds")
