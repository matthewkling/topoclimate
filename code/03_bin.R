
# This script summarizes subplot-level topoclimate data into bins.
# inputs: fia_occurrences.rds, fia_topoclimate.rds
# outputs: various versions of binned data

library(tidyverse)
library(raster)
library(furrr)
library(parallel)
library(conflicted)

conflict_prefer("select", "dplyr", quiet = T)
conflict_prefer("filter", "dplyr", quiet = T)
conflict_prefer("extract", "raster", quiet = T)
conflict_prefer("lag", "stats", quiet = T)


# input file paths
topoclimate_file <- "data/derived/fia_topoclimate.rds"
occurrence_file <- "data/derived/fia_occurrences.rds"


# function to generate binned data for all species
bin_data <- function(nspecies = NULL, randspecies = F, topspecies = T, species = NULL, # used to select subsets of species
                     bbox = 2, # lat/lon buffer around species range
                     nbins = 5, # number of summary bins in each dimension
                     vars = c("bio5", "bio6", "bio12"), # climate variables modified by topography
                     topo_vars = c("northness", "eastness"), # topographic variables
                     avars = "bio1", # additional variables that alter topoclimate effects
                     ncores = 5, # number of cores for parallel processing
                     outfile){ # path to save output file
        
        # load and filter input data
        d <- left_join(readRDS(topoclimate_file),
                       readRDS(occurrence_file)) %>%
                filter(!is.na(northness), !is.na(tpi),
                       !is.na(windward),
                       is.finite(bio1), is.finite(bio12)) %>%
                rename(species = gs)
        focal_species <- unique(d$species[!is.na(d$focal_spp) & d$focal_spp == T])
        
        # rescale predictors
        svars <- c(vars, avars, topo_vars, "lon", "lat")
        mv <- d %>% # FIA-wide mean and sd for every variable, for later use
                group_by(subplot_id) %>% slice(1) %>% ungroup() %>%
                select(all_of(c(topo_vars, vars, avars)), lon, lat) %>%
                summarize_all(list(mean = mean, sd = sd))
        dc <- d %>%
                select(all_of(c(topo_vars, vars, avars)), lon, lat,
                       plot_id, subplot_id, species)
        for(v in svars) dc[[v]] <- (dc[[v]] - mv[[paste0(v, "_mean")]]) / mv[[paste0(v, "_sd")]]
        
        # build a focal species list, and potentially subset it
        spp <- count(dc, species) %>% 
                filter(species %in% focal_species, !is.na(species)) %>% 
                arrange(desc(n))
        if(!topspecies) spp <- spp %>% arrange(n)
        if(randspecies) spp <- spp %>% sample_n(nrow(.))
        if(!is.null(nspecies)) spp <- slice(spp, 1:nspecies)
        if(!is.null(species)) spp <- spp[spp$species %in% species,]
        
        # species-wise binning function
        sp_bins <- function(s, data = dc, bbox, topo_vars, clim_vars, avars = NULL,
                            buffer = .1, trim = .001, # topoclimate bbox params
                            pb = NULL, nbins = 5){
                
                require(tidyverse)
                select <- dplyr::select
                # if(!is.null(pb)) pb$tick()$print()
                
                # limit to geographic bounding box;
                # classify subplot-level presences and absences
                ll <- 110/85 # lat/lon ratio at ~40N
                sd <- data %>%
                        mutate(pres = species == s & !is.na(species)) %>%
                        select(pres, plot_id, subplot_id, lon, lat, 
                               all_of(c(clim_vars, topo_vars, avars))) %>%
                        filter(between(lon, min(lon[pres]) - bbox*ll, max(lon[pres]) + bbox*ll),
                               between(lat, min(lat[pres]) - bbox, max(lat[pres]) + bbox)) %>%
                        group_by(subplot_id) %>%
                        mutate(species = s,
                               pres = any(pres)) %>%
                        slice(1) %>%
                        ungroup()
                
                # limit to climatic bounding box
                buffrange <- function(x, buffer = .1, trim = .001){
                        x <- quantile(x, c(trim, 1-trim)) # trim extreme outliers
                        x + c(-1, 1) * diff(x) * buffer
                }
                for(v in clim_vars){
                        vrange <- buffrange(sd[[v]][sd$pres], buffer = buffer, trim = trim)
                        sd <- sd %>% filter(between(sd[[v]], vrange[1], vrange[2]))
                }
                
                # binning function
                bin <- function(x, n = 5){
                        cuts <- seq(min(x) - .00001, max(x) + .00001, length.out = n+1)
                        y <- cut(x, cuts, labels = F)
                        midpoints <- ((cuts + lag(cuts))/2)[2:(n+1)]
                        y <- midpoints[y]
                        y
                }
                
                # assign subplots to bins
                binvars <- c(clim_vars, topo_vars)
                x <- sd %>%
                        mutate_at(all_of(binvars), funs(b = bin), n = nbins) %>% 
                        group_by(across(all_of(paste0(binvars, "_b"))))
                
                # count subplots and species presences per bin
                x1 <- x %>% summarize(n = n(),
                                      npres = sum(pres),
                                      species = s) %>% ungroup()
                
                # calculate bin topoclimate (mean of subplots in bin)
                x2 <- x %>% summarize_at(c(binvars, avars), mean) %>% ungroup()
                
                # join and return
                left_join(x1, x2) %>% select(-paste0(binvars, "_b"))
        }
        
        # bin data for every species
        plan(multisession, workers = ncores)
        md <- future_map_dfr(spp$species, sp_bins, nbins = nbins,
                             data = dc, bbox = bbox,
                             topo_vars = topo_vars, clim_vars = vars, avars = avars)
        
        # metadata on species ID numbers
        spid <- tibble(id = as.integer(factor(md$species)),
                       species = md$species) %>%
                distinct()
        
        # package and export data
        y <- list(md = md,
                  mv = mv,
                  spid = spid,
                  dc = dc) #### fixme: don't export this--look at downstream scripts to see what is needed that CAN be exported.
        saveRDS(y, outfile)
}



# generate several variants of binned data, with different variable sets

bin_data(nspecies = NULL, bbox = 10, nbins = 4,
         vars = c("bio5", "bio6", "bio12"),
         avars = c("bio1", "wspeed"),
         topo_vars = c("northness", "eastness", "windward"),
         outfile = "data/derived/binned/4bn_wind.rds")

bin_data(nspecies = NULL, bbox = 10, nbins = 4,
         vars = c("bio5", "bio6", "bio12"),
         avars = c("bio1"),
         topo_vars = c("northness", "eastness"),
         outfile = "data/derived/binned/4bn_nowind.rds")

bin_data(nspecies = NULL, bbox = 10, nbins = 4,
         vars = c("bio5", "bio6", "bio12"),
         avars = c("bio1"),
         topo_vars = c("northness", "eastness", "tpi"),
         outfile = "data/derived/binned/4bn_tpi_nowind.rds")

bin_data(nspecies = NULL, bbox = 10, nbins = 4,
         vars = c("bio5", "bio6", "bio12"),
         avars = c("bio1", "wspeed"),
         topo_vars = c("northness", "eastness", "tpi", "windward"),
         outfile = "data/derived/binned/4bn_tpi_wind.rds")
