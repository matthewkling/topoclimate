
library(tidyverse)
library(doParallel)
library(raster)
library(furrr)

select <- dplyr::select


# Katie: change this file path to point to your coordinates data.
# This should be a csv with 3 variables: PLT_CN, LON, LAT
coordinates_file <- "tpi_check/coordinates.csv"



## prep FIA data ##########

# load preassembled data with plot coordinates and physiography
d <- read_csv("tpi_check/plot_physiography.csv")

# add plot coordinates
d <- read_csv(coordinates_file) %>%
        inner_join(d) %>%
        filter(is.finite(LON))

# select a random subset of plots
d <- d %>% sample_n(2000)



## load elevation data ############

elev_file <- "data/raw/landfire/LF2016_Elev_200_CONUS/Tif/LC16_Elev_200.tif"
dem <- raster(elev_file)
NAvalue(dem) <- -9999



## compute TPI variables ###############

# mTPI neighborhood radii, in meters
radii <- c(100, 225, 500)

# precalculate focal matrices for each radius
fws <- list()
for(r in radii) fws[[paste0("r", r)]] <- focalWeight(dem, r, "circle")

# compute topographic position index (TPI) for a given point location
# (standardized multi-scale TPI -- mTPIs sensu Theobald et al. (2015))
tpi <- function(x, y, radii, dem, std = T){
        
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
                er <- as.vector(er)
                
                er[fw == 0] <- NA
                
                # mTPI
                eo <- er[(length(er) + 1) / 2] # elev of central cell
                er <- na.omit(er[er != 0])
                le <- length(er)
                en <- sum(er) / le # mean elevation of neighborhood
                tpi[i] <- eo - en
                
                # mTPIs
                if(std){
                        es <- sqrt(sum((er - en)^2) / le) # standard deviation of elevation
                        tpi[i] <- tpi[i] / es # z-score
                }
        }
        mean(tpi, na.rm = T) # average TPI across radii
}

# coordinates for extraction
e <- d
coordinates(e) <- c("LON", "LAT")
fia_crs <- "+init=epsg:4269 +proj=longlat +ellps=GRS80 +datum=NAD83"
crs(e) <- fia_crs
e <- spTransform(e, crs(dem))
e <- crop(e, dem) %>% as.data.frame()

# test for a single site
tpi(x = e$LON[1], y = e$LAT[1], radii, dem)

# compute TPI for all sites (in parallel for speed)
plan(multisession, workers = detectCores() - 1)
e$mtpis <- e %>%
        select(x = LON, y = LAT) %>%
        future_pmap_dbl(., possibly(tpi, NA_real_), r = radii, dem = dem, std = T)
e$mtpi <- e %>%
        select(x = LON, y = LAT) %>%
        future_pmap_dbl(., possibly(tpi, NA_real_), r = radii, dem = dem, std = F)




## plot results #########

e <- e %>% mutate(phys = recode(PHYSCLCD,
                                "11" = "xeric: dry tops",
                                "12" = "xeric: dry slopes",
                                "13" = "xeric: deep sands",
                                "19" = "xeric: other",
                                "21" = "mesic: flatwoods",
                                "22" = "mesic: rolling uplands",
                                "23" = "mesic: moist slopes and coves",
                                "24" = "mesic: narrow floodplains/bottomlands",
                                "25" = "mesic: broad floodplains/bottomlands",
                                "29" = "mesic: other",
                                "31" = "hydric: swamps/bogs",
                                "32" = "hydric: small drains",
                                "33" = "hydric: bays and wet pocosins",
                                "34" = "hydric: beaver ponds",
                                "35" = "hydric: cypress ponds",
                                "39" = "hydric: other"))


pd <- e %>%
        pivot_longer(c(mtpi, mtpis), 
                     names_to = "variable", values_to = "raw") %>%
        group_by(variable) %>%
        mutate(quantile = ecdf(raw)(raw)) %>%
        pivot_longer(c(raw, quantile), 
                     names_to = "stat", values_to = "value")

# version with all plots
p <- pd %>%
        ggplot(aes(value, phys)) +
        facet_grid(. ~ variable + stat, scales = "free") +
        geom_vline(xintercept = 0) +
        geom_violin(color = "red", fill = "red", alpha = .5) +
        geom_boxplot(fill = NA) +
        theme_bw() +
        theme(axis.title.y = element_blank())
ggsave("tpi_check/PHYSCLCD_TPI_boxplots.png",
       p, width = 12, height = 8, units = "in")

p <- pd %>%
        ggplot(aes(phys)) +
        facet_grid(. ~ variable + stat, scales = "free") +
        geom_histogram(stat = "count") +
        theme_bw() +
        theme(axis.title.y = element_blank()) +
        coord_flip() +
        scale_y_log10()
ggsave("tpi_check/PHYSCLCD_TPI_histograms.png",
       p, width = 12, height = 8, units = "in")


# version separated by UNITCD (crude control for regional climate)
p <- pd %>%
        ggplot(aes(value, phys)) +
        facet_grid(UNITCD ~ variable + stat, scales = "free") +
        geom_vline(xintercept = 0) +
        geom_violin(color = "red", fill = "red", alpha = .5) +
        geom_boxplot(fill = NA) +
        theme_bw() +
        theme(axis.title.y = element_blank())
ggsave("tpi_check/PHYSCLCD_TPI_UNITCD_boxplots.png",
       p, width = 12, height = 16, units = "in")

p <- pd %>%
        ggplot(aes(phys)) +
        facet_grid(UNITCD ~ variable + stat, scales = "free") +
        geom_histogram(stat = "count") +
        theme_bw() +
        theme(axis.title.y = element_blank()) +
        coord_flip() +
        scale_y_log10()
ggsave("tpi_check/PHYSCLCD_TPI_UNITCD_histograms.png",
       p, width = 12, height = 16, units = "in")
