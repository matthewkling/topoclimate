
library(tidyverse)
library(raster)
select <- dplyr::select


# load and combine data #######################


# plot attributes, exact coords
ct <- cols_only(nfi_plot = "d", land_cover = "c", veg_type = "c", 
                plot_completion = "c", incomp_reason = "c",
                utm_zone = "d", utm_e = "d", utm_n = "d", province = "c",
                slope = "d", aspect = "d", elevation = "d")
d <- list.files("data/raw/cnfi/2021_05_03_Kling_GP_exactCoords_QC_NS_NB_NT_ON", 
                pattern = "gp_site_info", recursive = T, full.names = T)
d <- d[str_detect(d, "climate")]
d <- d %>%
        map(read_csv, col_types = ct) %>%
        bind_rows()

# plot attributes, mixed coords
d <- read_csv("data/raw/cnfi/2021_05_03_Kling_GP_mixed/all_gp_climate_mixed/all_gp_site_info.csv",
              col_types = ct) %>%
        filter(!province %in% d$province) %>%
        bind_rows(d)

# remove planted plots 
regen <- list.files("data/raw/cnfi/2021_05_03_Kling_GP_exactCoords_QC_NS_NB_NT_ON", 
                    pattern = "origin", recursive = T, full.names = T) %>%
        map(read_csv) %>%
        bind_rows()
regen <- read_csv("data/raw/cnfi/2021_05_03_Kling_GP_mixed/all_gp_origin_mixed/all_gp_origin.csv") %>%
        filter(! nfi_plot %in% regen$nfi_plot) %>%
        bind_rows(regen) %>%
        filter(regen_type %in% c("PLA", "SUP")) %>%
        select(nfi_plot) %>%
        distinct() %>%
        mutate(planted = "planted")
d <- d %>%
        left_join(regen) %>%
        filter(is.na(planted))

# standardize coord projections
fia_crs <- CRS("+proj=longlat +datum=WGS84")
utm2ll <- function(x){
        x <- split(x, x$utm_zone) %>%
                map(function(x){
                        x$x <- x$utm_e
                        x$y <- x$utm_n
                        coordinates(x) <- c("utm_e", "utm_n")
                        crs(x) <- paste0("+proj=utm +zone=", x$utm_zone[1], 
                                         " +datum=WGS84 +units=m +ellps=WGS84")
                        spTransform(x, fia_crs)
                })
        x <- do.call("rbind", x)
        bind_cols(x@data, as.data.frame(coordinates(x))) %>%
                rename(lon = utm_e, lat = utm_n)
}
d <- utm2ll(d) 

# remove plots with inexact coordinates
d <- d %>%
        mutate(exact = paste0(x %>% as.character() %>% substr(3, 6),
                              y %>% as.character() %>% substr(4, 7)) != "00000000") %>%
        filter(exact) %>%
        select(-exact)

# remove duplicates
d <- d %>% group_by(nfi_plot) %>%
        mutate(keep = length(nfi_plot) == 1 | plot_completion == "F") %>%
        filter(keep) %>%
        select(-keep)

# remove plots with suspicious province code (a couple near MT/SK border)
provinces <- getData("GADM", country = "CAN", level = 1) %>%
        spTransform(fia_crs)
prov <- extract(provinces, d %>% ungroup() %>% select(lon, lat)) %>%
        separate(HASC_1, c("country", "province"))
d <- d %>%
        ungroup() %>%
        mutate(prov = ifelse(prov$province != "NF", prov$province, "NL")) %>%
        filter(prov == province) %>%
        select(-prov)


# tree occurrences, exact coords
o <- list.files("data/raw/cnfi/2021_05_03_Kling_GP_exactCoords_QC_NS_NB_NT_ON", 
                pattern = "ltp_tree_species_comp", recursive = T, full.names = T)%>%
        map(read_csv) %>%
        bind_rows() %>%
        select(nfi_plot, genus, species, variety) %>%
        distinct()

# tree occurrences, mixed coords
o <- read_csv("data/raw/cnfi/2021_05_03_Kling_GP_mixed/all_gp_trees_mixed/all_gp_ltp_tree_species_comp.csv") %>%
        select(nfi_plot, genus, species, variety) %>%
        filter(! nfi_plot %in% unique(o$nfi_plot)) %>%
        distinct() %>%
        bind_rows(o) %>%
        filter(nfi_plot %in% unique(d$nfi_plot))

# species list
spp <- o %>% count(genus, species, variety, name = "n_plots")


## export ############

o %>% saveRDS("data/derived/cnfi_occurrences.rds")
d %>% saveRDS("data/derived/cnfi_plots.rds")
spp %>% saveRDS("data/derived/cnfi_species.rds")
