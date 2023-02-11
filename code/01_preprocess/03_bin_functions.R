
# function to generate binned data for all species
bin_data <- function(spp,
                     subplots = 1:4, # indices of subplots to use
                     bbox = 2, # lat/lon buffer around species range
                     nbins = 5, # number of summary bins in each dimension
                     vars = c("bio5", "bio6", "bio12"), # climate variables modified by topography
                     topo_vars = c("northness", "eastness"), # topographic variables
                     avars = "bio1", # additional variables that alter topoclimate effects
                     ncores = 5, # number of cores for parallel processing
                     file_path = NULL){ # path to save output file
        
        message(paste("generating", file_path))
        
        # rescale predictors
        svars <- c(vars, avars, topo_vars, "lon", "lat")
        mv <- d %>% # dataset-wide mean and sd for every variable, for later use
                group_by(subplot_id) %>% slice(1) %>% ungroup() %>%
                filter(subplot %in% subplots) %>%
                select(all_of(c(topo_vars, vars, avars)), lon, lat) %>%
                summarize_all(list(mean = mean, sd = sd))
        dc <- d %>%
                filter(subplot %in% subplots) %>%
                select(all_of(c(topo_vars, vars, avars)), lon, lat,
                       plot_id, subplot_id, species)
        for(v in svars) dc[[v]] <- (dc[[v]] - mv[[paste0(v, "_mean")]]) / mv[[paste0(v, "_sd")]]
        
        # species-wise binning function
        sp_bins <- function(s, data = dc, bbox, topo_vars, clim_vars, avars = NULL,
                            buffer = .1, trim = .001, # topoclimate bbox params
                            pb = NULL, nbins = 5){
                
                require(tidyverse)
                select <- dplyr::select
                
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
        md <- future_map_dfr(spp, sp_bins, 
                             nbins = nbins, data = dc, bbox = bbox,
                             topo_vars = topo_vars, clim_vars = vars, avars = avars)
        
        # metadata on species ID numbers
        spid <- tibble(id = as.integer(factor(md$species)),
                       species = md$species) %>%
                distinct()
        
        # package and export data
        y <- list(md = md,
                  mv = mv,
                  spid = spid,
                  # dc = dc, # not exporting this, for locational privacy
                  subplots = subplots, bbox = bbox, nbins = nbins,
                  vars = vars, avars = avars, topo_vars = topo_vars)
        
        if(!is.null(file_path)) saveRDS(y, file_path)
        if(is.null(file_path)) return(y)
}
