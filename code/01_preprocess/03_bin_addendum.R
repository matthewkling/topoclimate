
# This script is an addendum to 03_bin.R, 
# generating some additional binning formulations



library(tidyverse)
library(raster)
library(furrr)
library(parallel)
library(conflicted)

conflict_prefer("select", "dplyr", quiet = T)
conflict_prefer("filter", "dplyr", quiet = T)
conflict_prefer("extract", "raster", quiet = T)
conflict_prefer("lag", "stats", quiet = T)

source("code/01_preprocess/03_bin_functions.r")


## load and merge data ##################

species <- readRDS("data/derived/species.rds")

fia <- left_join(readRDS("data/derived/fia_topoclimate.rds"),
                 readRDS("data/derived/fia_occurrences.rds")) %>%
        mutate(wua = wu * waniso,
               wva = wv * waniso) %>%
        select(plot_id, subplot_id, lon, lat, gs,
               bio1, bio12, bio5, bio6, wspeed,
               northness, eastness, tpi, tpis, windward,
               wua, wva)

nfi <- left_join(readRDS("data/derived/cnfi_topoclimate.rds"),
                 readRDS("data/derived/cnfi_occurrences.rds")) %>%
        left_join(species %>% select(gs, genus, species) %>% na.omit()) %>%
        select(-genus, -species) %>%
        mutate(wua = wu * waniso,
               wva = wv * waniso) %>%
        select(plot_id = nfi_plot, lon, lat, gs,
               bio1, bio12, bio5, bio6, wspeed,
               northness, eastness, tpi, tpis, windward,
               wua, wva) %>%
        mutate(plot_id = as.character(plot_id))

d <- bind_rows(fia, nfi) %>%
        filter(!is.na(northness), !is.na(tpi),
               !is.na(windward),
               is.finite(bio1), is.finite(bio12)) %>%
        rename(species = gs)

# add subplot ID for CNFI plots, and subplot index for all plots
d <- d %>% mutate(subplot_id = ifelse(is.na(subplot_id), 
                                      paste(plot_id, "1"), subplot_id),
                  subplot = str_sub(subplot_id, -1, -1) %>% as.integer())



# identify focal species -- those occurring on 100+ plots
exclude <- species$gs[!is.na(species$exclude) & species$exclude]
species <- d %>%
        filter(! is.na(species),
               ! species %in% exclude) %>%
        group_by(species) %>%
        summarize(n_plots = length(unique(plot_id)),
                  n_subplots = length(unique(subplot_id))) %>%
        mutate(focal = n_plots >= 100)
focal_species <- species %>% filter(focal) %>% pull(species)





## binning ###########


### define several variants of binned data, with different variable sets ###

# parameter combinations
b <- expand_grid(
        spp = list(focal_species), 
        bbox = 10, 
        nbins = 4,
        vars = list(c("bio5", "bio6", "bio12", "wua", "wva", "wspeed"),
                    c("bio5", "bio6", "bio12", "wua", "wva", "wspeed", "lat"),
                    c("bio5", "bio6", "bio12", "lat")),
        avars = list(c("bio1")),
        topo_vars = list(c("northness", "eastness", "tpi")),
        subplots = list(1,
                        1:4))

# file paths for binned data
b <- b %>% mutate(file_path = paste0("data/derived/binned/param_A", 1:nrow(.), ".rds"))

# save parameter metadata
saveRDS(b, "data/derived/binned/metadata_A.rds")


### perform binning ###

b %>% pmap(bin_data, ncores = 1)
