
# This script merges various raw FIA inputs to generate two output tables:
# - subplot attributes (one row per subplot)
# - species occurrences per subplot (one row per unique species-subplot)


library(data.table)
library(tidyverse)
library(raster)
library(geosphere)
library(furrr)
library(parallel)
library(foreign)
library(readxl)
library(janitor)
library(conflicted)

conflict_prefer("select", "dplyr", quiet = T)
conflict_prefer("filter", "dplyr", quiet = T)
conflict_prefer("between", "dplyr", quiet = T)
conflict_prefer("extract", "raster", quiet = T)

## interior Alaska -- a separate SQLite database ##########

library(dbplyr)
library(RSQLite)
ak <- dbConnect(SQLite(), "data/raw/fia/interior_AK/SQLite-FIADB-INTAK-TANANA/SQLite_FIADB_INTAK_TANANA.db")
src_dbi(ak)

plot_AK <- tbl(ak, "PLOT") %>% 
        select(PLT_CN = CN, MANUAL, DESIGNCD, PLOT_STATUS_CD,
               STATECD, INVYR, UNITCD, COUNTYCD, PLOT, 
               LAT, LON, ELEV) %>%
        as_tibble()

subplot_AK <- tbl(ak, "SUBPLOT") %>% 
        select(PLT_CN, SUBP, SUBP_CN = CN, PREV_SBP_CN, POINT_NONSAMPLE_REASN_CD,
               SLOPE, ASPECT) %>%
        as_tibble()

cond_AK <- tbl(ak, "COND") %>% 
        select(PLT_CN, CONDID, PROP_BASIS, STDORGCD, PHYSCLCD, COND_STATUS_CD, CONDPROP_UNADJ,
               SLOPE_COND = SLOPE, ASPECT_COND = ASPECT) %>%
        as_tibble()

tree_AK <- tbl(ak, "TREE") %>%
        select(PLT_CN, PLOT, SUBP, TREE, SPCD, DIA, HT) %>%
        as_tibble()

species_AK <- tbl(ak, "REF_SPECIES") %>%
        select(SPCD, COMMON_NAME, GENUS, SPECIES) %>%
        as_tibble()


## primary FIA database ############

# PLOT table
# primary key: CN
# unique key: STATECD, INVYR, UNITCD, COUNTYCD, PLOT
plot <- fread("data/raw/fia/ENTIRE/PLOT.csv", stringsAsFactors = F, 
              colClasses = c(CN = "character")) %>% 
        select(PLT_CN = CN, MANUAL, DESIGNCD, PLOT_STATUS_CD,
               STATECD, INVYR, UNITCD, COUNTYCD, PLOT, 
               LAT, LON, ELEV) %>% as_tibble() %>%
        bind_rows(plot_AK)

# SUBPLOT table
# primary key: CN
# unique key: PLT_CN, SUBP
# When PLOT.MANUAL <1.0, the field crew measured slope at a condition level
# When PLOT.MANUAL >=1.0, slope is collected at a subplot level, and then the aspect 
# from the subplot representing the greatest proportion of the condition is assigned as a surrogate
subplot <- fread("data/raw/fia/ENTIRE/SUBPLOT.csv", stringsAsFactors = F, 
                 colClasses = c(CN = "character", PLT_CN = "character")) %>% 
        select(PLT_CN, SUBP, SUBP_CN = CN, PREV_SBP_CN, POINT_NONSAMPLE_REASN_CD,
               SLOPE, ASPECT) %>% as_tibble() %>%
        bind_rows(subplot_AK %>% mutate(PREV_SBP_CN = as.integer(PREV_SBP_CN)))

# CONDITION table (used to identify and exclude planted plots)
# primary key: CN
# unique key: PLT_CN, CONDID
cond <- fread("data/raw/fia/ENTIRE/COND.csv", stringsAsFactors = F, 
              colClasses = c(CN = "character", PLT_CN = "character")) %>% 
        select(PLT_CN, CONDID, PROP_BASIS, STDORGCD, PHYSCLCD, COND_STATUS_CD, CONDPROP_UNADJ,
               SLOPE_COND = SLOPE, ASPECT_COND = ASPECT) %>% as_tibble() %>%
        bind_rows(cond_AK)
planted <- cond %>%
        filter(STDORGCD == 1) %>%
        select(PLT_CN) %>%
        distinct() %>%
        mutate(planted = TRUE)

# TREE table
# primary key: CN
# unique key: PLT_CN, SUBP, TREE
tree <- fread("data/raw/fia/ENTIRE/TREE.csv", 
              stringsAsFactors = F, 
              select = c("PLT_CN", "PLOT", "SUBP", "TREE", "SPCD"),
              colClasses = c(PLT_CN = "character")) %>% as_tibble() %>%
        bind_rows(tree_AK)

# SPECIES table
species <- fread("data/raw/fia/REF_SPECIES.csv", stringsAsFactors = F) %>%
        select(SPCD, COMMON_NAME, GENUS, SPECIES) %>% as_tibble() %>%
        bind_rows(species_AK)



## join tables ########

d <- plot %>%
        left_join(planted) %>%
        left_join(subplot) %>%
        left_join(tree) %>%
        left_join(species) %>% 
        mutate(PLOT_ID = paste(STATECD, UNITCD, COUNTYCD, PLOT),
               SUBPLOT_ID = paste(PLOT_ID, SUBP)) %>%
        filter(is.na(planted),
               DESIGNCD %in% c(1, 111:118, 230:242, 311:323, 328, 501:506), # 4-subplot design
               SUBP <= 4, # remove extraneous subplots
               MANUAL >= 1) %>% # for versions < 1, slope and aspect can't be attributed to subplots
        select(PLT_CN, STATECD, UNITCD, COUNTYCD, PLOT, SUBP, INVYR, DESIGNCD,
               PLOT_ID, SUBPLOT_ID, LON, LAT, ELEV, PLOT_STATUS_CD, POINT_NONSAMPLE_REASN_CD,
               SLOPE, ASPECT, GENUS, SPECIES) %>%
        clean_names() %>%
        mutate(lon = ifelse(lon>0, lon-360, lon)) %>%
        filter(lat > 24)

# subplots to exclude -- those marked as unsampled for ALL inventories
exclude <- d %>% 
        group_by(subplot_id) %>%
        summarize(exclude = !any(is.na(point_nonsample_reasn_cd))) %>%
        filter(exclude)
d <- filter(d, ! subplot_id %in% exclude$subplot_id)
d <- select(d, -point_nonsample_reasn_cd)

# subplots, including those without trees
sp <- d %>%
        select(plt_cn, statecd, unitcd, countycd, plot, invyr,
               plot_id, subplot_id, lon, lat, slope, aspect) %>%
        distinct()

# tree occurrences
o <- d %>%
        filter(!is.na(elev), !is.na(slope), !is.na(lon), !is.na(lat), 
               !is.na(genus), genus != "Tree", species != "spp.") %>%
        mutate(gs = paste(genus, species)) %>%
        select(plot_id, subplot_id, gs) %>%
        distinct()


## focal species ######

# criteria from Katie:
# "...exclude any that had exclusions (columns A-D) previously, just for the sake of continuity of data collection among regions. 
# We may also want to exclude any species that are purely 'urban' (have only a "U" in columns L-R)"
fia_spp <- read_xlsx("data/raw/fia/2019 Master Species FGver90_10_01_2019_rev_2_10_2020.xlsx") %>%
        clean_names() %>%
        filter(non_urban_exclusion == ".",
               nfs_region_exclusion == ".",
               mainland_exclusion == ".") %>%
        mutate(urbanity = paste(past_nrs_tally, past_pnw_tally, past_rmrs_tally, past_srs_tally,
                                caribbean, pacific, urban),
               urban = str_detect(urbanity, "U") & !str_detect(urbanity, "P|R|C|K")) %>%
        filter(!urban) %>%
        mutate(gs = paste(genus, species))

spp <- o %>%
        group_by(gs) %>%
        summarize(n_plot = length(unique(plot_id)),
                  n_subp = length(unique(subplot_id)),
                  exclude = ! gs %in% fia_spp$gs) %>%
        distinct()


## export ############

o %>% saveRDS("data/derived/fia_occurrences.rds")
sp %>% saveRDS("data/derived/fia_subplots.rds")
spp %>% saveRDS("data/derived/fia_species.rds")

sp %>% select(statecd, unitcd, countycd, plot, invyr, lon, lat) %>%
        distinct() %>%
        write_csv("data/derived/fia_plot_coordinates.csv")


# cleanup
ls() %>% rm()
gc()
