

library(tidyverse)
library(data.table)

select <- dplyr::select


## assemble dataset of plot coordinates and physiography ###########

# PLOT table
# primary key: CN
# unique key: STATECD, INVYR, UNITCD, COUNTYCD, PLOT
plot <- fread("data/raw/fia/ENTIRE/PLOT.csv", stringsAsFactors = F, 
              colClasses = c(CN = "character")) %>% 
        select(PLT_CN = CN, MANUAL, DESIGNCD, PLOT_STATUS_CD,
               STATECD, INVYR, UNITCD, COUNTYCD, PLOT, 
               LAT, LON, ELEV) %>% as_tibble() %>%
        filter(DESIGNCD %in% c(1, 506),
               MANUAL >= 1,
               between(LAT, 25, 50),
               between(LON, -125, -65))

# CONDITION table
# primary key: CN
# unique key: PLT_CN, CONDID
cond <- fread("data/raw/fia/ENTIRE/COND.csv", stringsAsFactors = F, 
              colClasses = c(CN = "character", PLT_CN = "character")) %>% 
        select(PLT_CN, CONDID, PROP_BASIS, STDORGCD, PHYSCLCD, COND_STATUS_CD, CONDPROP_UNADJ,
               SLOPE_COND = SLOPE, ASPECT_COND = ASPECT) %>% as_tibble()

# physiographic class most representative of each plot
# (have to do this plot-level, because conditions aren't defined at the subplot level)
physiog <- cond %>%
        group_by(PLT_CN) %>%
        filter(CONDPROP_UNADJ == max(CONDPROP_UNADJ, na.rm = T)) %>% # the most abundant condition
        arrange(CONDID) %>% slice(1) %>% # if there's an abundance tie, condition at center of plot
        ungroup() %>%
        select(PLT_CN, PHYSCLCD, CONDPROP_UNADJ, COND_STATUS_CD)

d <- plot %>%
        left_join(physiog) %>%
        filter(CONDPROP_UNADJ == 1,
               !is.na(PHYSCLCD)) %>%
        select(UNITCD, PLT_CN, LON, LAT, PHYSCLCD)

d %>%
        select(UNITCD, PLT_CN, PHYSCLCD) %>%
        write_csv("tpi_check/plot_physiography.csv")

d %>%
        select(PLT_CN, LON, LAT) %>%
        write_csv("tpi_check/coordinates.csv")
