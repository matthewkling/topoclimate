
# reconcile taxonomy between FIA and NFI

library(tidyverse)

# load species lists
fia <- readRDS("data/derived/fia_species.rds")
nfi <- readRDS("data/derived/cnfi_species.rds")

# format
d <- fia %>%
        separate(gs, c("g", "s"), sep = " ", remove = F) %>%
        mutate(genus = str_sub(g, 1, 4) %>% str_to_upper(),
               species = str_sub(s, 1, 3) %>% str_to_upper()) %>%
        rename(fia_plots = n_plot,
               fia_subplots = n_subp) %>%
        full_join(nfi %>% rename(nfi_plots = n_plots)) %>%
        arrange(genus, species) %>%
        filter(species != "SPP")

# map NFI anomalies to FIA taxonomy
d <- d %>%
        mutate(gs = case_when(genus == "ACER" & species == "SAH" ~ "Acer saccharum",
                              genus == "ALNU" & species == "INC" ~ "Alnus rugosa",
                              genus == "ALNU" & species == "RUG" ~ "Alnus rugosa",
                              genus == "POPU" & species == "TRI" ~ "Populus balsamifera",
                              genus == "RHAM" & species == "CAT" ~ "Rhamnus cathartica",
                              genus == "SALI" & species == "DIS" ~ "Salix discolor",
                              genus == "SALI" & species == "HUM" ~ NA_character_, # not in CNFI docs
                              genus == "SALI" & species == "PLA" ~ NA_character_, # not in CNFI docs
                              genus == "SALI" & species == "SCO" ~ "Salix scouleriana",
                              TRUE ~ gs))

# export
d %>% saveRDS("data/derived/species.rds")
