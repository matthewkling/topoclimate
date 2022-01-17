
# these functions are shared by the US and Canada extraction scripts



# compute topographic position index (TPI) for a given point location
# (standardized multi-scale TPI -- mTPIs sensu Theobald et al. (2015))
mtpi <- function(x, y, radii, dem, fws){
        
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
        tpis <- rep(NA, length(rr))
        for(i in 1:length(rr)){
                r <- rr[i]
                fw <- fws[[paste0("r", r)]]
                rad <- (dim(fw)[1] - 1) / 2
                
                er <- e[(rw - rad):(rw + rad), 
                        (cl - rad):(cl + rad)]
                er <- as.vector(er)
                
                er[fw == 0] <- NA
                
                # mTPI
                eo <- er[(length(er) + 1) / 2] # elev of central cell
                er <- na.omit(er[er != 0])
                le <- length(er)
                en <- sum(er) / le # mean elevation of neighborhood
                tpi[i] <- eo - en
                
                # mTPIs
                es <- sqrt(sum((er - en)^2) / le) # standard deviation of elevation
                tpis[i] <- tpi[i] / es # z-score
                
        }
        
        # average TPI across radii
        c(mtpi = mean(tpi, na.rm = T), 
          mtpis = mean(tpis, na.rm = T))
}



add_macroclimate <- function(sp){
        
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
        
        return(sp)
}


add_wind <- function(sp){
        
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
        
        return(sp)
}