
# figure 4 panels a and b 

library(tidyverse)
library(raster)
library(patchwork)
select <- dplyr::select


nfi <- readRDS("data/derived/cnfi_plots.rds")

fia <- readRDS("data/derived/fia_subplots_coords.rds") %>%
        group_by(plot_id) %>%
        sample_n(1) %>%
        ungroup()

plots <- bind_rows(select(nfi, lon, lat) %>% mutate(forest = T),
                   select(fia, lon, lat, forest = fia_occ))

coordinates(plots) <- c("lon", "lat")

climate <- list.files("data/raw/chelsa", pattern = "bio1", full.names = T) %>%
        stack() %>%
        extract(plots)

d <- bind_cols(plots@data, 
               as.data.frame(coordinates(plots)) %>% setNames(c("lon", "lat")), 
               as.data.frame(climate) %>% setNames(c("bio1", "bio12")))

f <- filter(d, forest)
nf <- filter(d, !forest)

f$color <- colormap::colors2d(select(f, bio1, bio12) %>% as.data.frame() %>%
                                      mutate(bio12 = log10(bio12)),
                              colors = c("forestgreen", "gold", "red", "cornflowerblue"),
                              xtrans = "rank", ytrans = "rank")
nf$color <- "gray40"

d <- bind_rows(f, nf) %>%
        arrange(forest)

md <- map_data("state")
wd <- map_data("world") %>%
        filter(region %in% c("Canada", "Mexico", "USA"))

style <- theme_bw() + 
        theme(strip.text = element_text(color = "white"),
              strip.background = element_rect(fill = "black", color = "black"),
              legend.position = "bottom")

map <- ggplot(d, aes(lon, lat, color = color)) +
        geom_polygon(data = wd, aes(long, lat, group = group), 
                     size = .1, fill = "gray80", color = "gray40") +
        geom_point(size = .001) +
        geom_path(data = md, aes(long, lat, group = group), 
                  size = .25, alpha = .5, color = "white") +
        scale_color_identity() +
        theme_void() +
        style + 
        theme(axis.text = element_blank(), 
              axis.title = element_blank(),
              axis.ticks = element_blank(), 
              panel.grid = element_blank()) +
        coord_map("azequalarea", orientation = c(65, -100, 0),
                  xlim = c(-131, -69),
                  ylim = c(27, 70))

gd <- f %>%
        select(bio1, bio12, color) %>%
        mutate(bio12 = log10(bio12+1),
               bio1r = plyr::round_any(bio1, 5),
               bio12r = plyr::round_any(bio12, .33),
               diff = sqrt((bio1r-bio1)^2 + (bio12r-bio12)^2),
               bio12r = (10^bio12r)-1) %>%
        group_by(bio1r, bio12r) %>%
        filter(diff == min(diff),
               between(bio1r, min(f$bio1), max(f$bio1)),
               between(bio12r, min(f$bio12), max(f$bio12))) %>%
        ungroup()


sctr <- ggplot(d, aes(bio1, bio12, color = color)) +
        geom_point(size = .25) +
        
        geom_path(data = gd, aes(bio1r, bio12r, group = bio1r), 
                  color = "black", size = .25) +
        geom_path(data = gd, aes(bio1r, bio12r, group = bio12r), 
                  color = "black", size = .25) +
        geom_point(data = gd, aes(bio1r, bio12r, fill = color), 
                   color = "black", shape = 21, size = 3) +
        
        scale_y_log10() +
        scale_color_identity() +
        scale_fill_identity() +
        coord_cartesian(xlim = range(f$bio1),
                        ylim = range(f$bio12),
                        expand = c(0, 0)) +
        theme_classic() +
        theme(panel.background = element_blank(),
              panel.grid = element_blank(),
              plot.background = element_blank(),
              legend.position = "none") +
        labs(x = "macro temperature (°C)",
             y = "macro precipitation (mm)") +
        style

p <- sctr + map + plot_layout(widths = c(1, 1.5))
ggsave("figures/plot_climate_map.png",
       p, width = 8, height = 3.78, units = "in", dpi = 1500)

p <- sctr + map + plot_layout(widths = c(1, 1))
ggsave("figures/plot_climate_map_v2.png",
       p, width = 8, height = 3.3, units = "in", dpi = 1500)

saveRDS(gd, "data/derived/macroclimate_grid.rds")
