
library(tidyverse)
library(raster)
library(patchwork)

select <- dplyr::select




location <- "idaho"

# elevation data ##############

# from https://apps.nationalmap.gov/downloader/#/
files <- list(idaho = "data/DEM/USGS_13_n43w113.tif")
elev <- raster(files[[location]])
names(elev) <- "elevation"

# plot(elev)
# ext <- drawExtent()
extents <- list(idaho = extent(-112.7458, -112.5855, 42.53036, 42.5841))
ext <- extents[[location]]

ep <- as(ext* 1.05, "SpatialPolygons") %>% 
        SpatialPolygonsDataFrame(data = data.frame(x = 4))
rgdal::writeOGR(ep, paste0("data/DEM/", location, ".kml"), "boundary", driver="KML",
                overwrite_layer = T)

e <- crop(elev, ext)
size <- c(geosphere::distm(c(ext@xmin, ext@ymin), c(ext@xmax, ext@ymin)),
          geosphere::distm(c(ext@xmin, ext@ymin), c(ext@xmin, ext@ymax))) / 1000

et <- crop(elev, ext * 1.5)
tpi <- et * 3 - 
        focal(et, focalWeight(et, .5 / 111)) -
        focal(et, focalWeight(et, .15 / 111)) -
        focal(et, focalWeight(et, .05 / 111))
tpi <- crop(tpi, e)

topo <- terrain(e, c("slope", "aspect"))
northness <- sin(topo$slope) * cos(topo$aspect)
eastness <- sin(topo$slope) * sin(topo$aspect)
southwestness <- sin(topo$slope) * sin(topo$aspect + pi * 5/4)

hillshade <- hillShade(topo$slope, topo$aspect, angle = 30, direction = -45)
# plot(hillshade, col = colorRampPalette(c("black", "white"))(100))

topography <- stack(topo, northness, eastness, southwestness, hillshade, tpi, e)
names(topography)[3:7] <- c("northness", "eastness", "southwestness", 
                            "hillshade", "tpi")


# climate data #############

# biovars
clim <- list.files("f:/chelsa/v2/climatology/", full.names=T,
                   pattern = "_bio1_|_bio5_|_bio6_|_bio12_") %>% 
        stack() %>% 
        crop(ext * 3)
names(clim) <- paste0("bio", c(1, 12, 5, 6))

clim[c(1, 3, 4)] <- clim[c(1, 3, 4)] * 10

climate <- clim %>% 
        resample(e) %>%
        crop(e)

climate_nn <- clim %>% 
        resample(e, "ngb") %>%
        crop(e)
names(climate_nn) <- paste0(names(climate_nn), "nn")

# topoclimate parameters ############


fitfile <- "stan/v21_windward/optim_ALLsp_nowind.rds"
datafile <- "stan/data/windex_ALLspp_4bn_nowind.rds"

fit <- readRDS(fitfile)
d <- readRDS(datafile)
md <- d$md
mv <- d$mv
dc <- d$dc


get_deltas <- function(fitfile){
        
        fit <- readRDS(fitfile)
        f <- as.data.frame(fit$par) %>%
                rownames_to_column("param") %>% as_tibble()
        names(f)[2] <- "value"
        
        lookup <- function(var){
                list(bio5 = "heat", 
                     bio6 = "cold", 
                     bio12 = "moisture",
                     bio1 = "tmean")[[var]]
        }
        
        vars = c("bio5", "bio6", "bio12")
        topo_vars = c("northness", "eastness")
        topo_mods = c("bio1", "bio12")
        
        topo_mods1 <- c("int", topo_mods)
        params <- tibble(mod = rep(topo_mods1, length(topo_vars)),
                         topo = rep(topo_vars, each = length(topo_mods1)))
        
        
        deltas <- f %>% 
                filter(str_detect(param, "delta")) %>%
                mutate(param = str_remove(param, "delta\\["),
                       param = str_remove(param, "\\]"),
                       param = str_replace(param, ",", "_")) %>%
                separate(param, c("param", "var")) %>%
                mutate(var = recode(var, 
                                    "1" = lookup(vars[1]), 
                                    "2" = lookup(vars[2]), 
                                    "3" = lookup(vars[3])),
                       topo = params$topo[as.integer(param)],
                       mod = params$mod[as.integer(param)],
                       mod = factor(mod, levels = params$mod[1:(length(topo_vars)+1)]))
        deltas
}

deltas <- get_deltas(fitfile)

b <- deltas %>%
        select(-param) %>%
        mutate(var = recode(var, "heat"  = "h", "cold" = "c", "moisture" = "m"),
               topo = recode(topo, "eastness"  = "e", "northness" = "n")) %>%
        unite(param, var, topo, mod) %>%
        spread(param, value)



# downscale ############

d <- stack(topography, climate, climate_nn) %>% 
        rasterToPoints() %>% as.data.frame() %>% as_tibble() %>%
        filter(is.finite(slope), is.finite(bio1)) %>%
        
        # convert to main units
        mutate(bio1 = bio1 / 10, 
               bio5 = bio5 / 10,
               bio6 = bio6 / 10,
               bio12 = log10(bio12 + 1),
               bio1nn = bio1nn / 10, 
               bio5nn = bio5nn / 10,
               bio6nn = bio6nn / 10,
               bio12nn = log10(bio12nn + 1)) %>%
        
        # standardize
        mutate(bio1 = (bio1 - mv$bio1_mean) / mv$bio1_sd,
               bio5 = (bio5 - mv$bio5_mean) / mv$bio5_sd,
               bio6 = (bio6 - mv$bio6_mean) / mv$bio6_sd,
               bio12 = (bio12 - mv$bio12_mean) / mv$bio12_sd,
               bio1nn = (bio1nn - mv$bio1_mean) / mv$bio1_sd,
               bio5nn = (bio5nn - mv$bio5_mean) / mv$bio5_sd,
               bio6nn = (bio6nn - mv$bio6_mean) / mv$bio6_sd,
               bio12nn = (bio12nn - mv$bio12_mean) / mv$bio12_sd,
               northness = (northness - mv$northness_mean) / mv$northness_sd,
               eastness = (eastness - mv$eastness_mean) / mv$eastness_sd) %>%
        
        # calculate topo effects
        mutate(heat_northness =     b$h_n_int + b$h_n_bio1 * bio1 + b$h_n_bio12 * bio12,
               heat_eastness =      b$h_e_int + b$h_e_bio1 * bio1 + b$h_e_bio12 * bio12,
               cold_northness =     b$c_n_int + b$c_n_bio1 * bio1 + b$c_n_bio12 * bio12,
               cold_eastness =      b$c_e_int + b$c_e_bio1 * bio1 + b$c_e_bio12 * bio12,
               moisture_northness = b$m_n_int + b$m_n_bio1 * bio1 + b$m_n_bio12 * bio12,
               moisture_eastness =  b$m_e_int + b$m_e_bio1 * bio1 + b$m_e_bio12 * bio12) %>%
        
        # calculate microclimate
        mutate(heat = bio5 + heat_northness * northness + heat_eastness * eastness,
               cold = bio6 + cold_northness * northness + cold_eastness * eastness,
               moisture = bio12 + moisture_northness * northness + moisture_eastness * eastness) %>%
        
        # destandardize
        mutate(h = heat * mv$bio5_sd + mv$bio5_mean,
               c = cold * mv$bio6_sd + mv$bio6_mean,
               m = 10 ^ (moisture * mv$bio12_sd + mv$bio12_mean) - 1)


## separate plots ######

pt <- ggplot(d, aes(x, y, fill = hillshade, z = elevation)) +
        geom_raster() +
        stat_contour(color = "black", alpha = .25, size = .25) +
        scale_fill_gradientn(colors = c("black", "white")) +
        labs(fill = "hillshade  ")

p_topo <- d %>%
        filter(is.finite(slope)) %>%
        select(x, y, northness, eastness, southwestness, tpi, elevation) %>%
        gather(var, value, -x, -y, -elevation) %>%
        group_by(var) %>%
        mutate(value = scales::rescale(value)) %>%
        ggplot(aes(x, y, fill = value, z = elevation)) +
        facet_wrap(~var, ncol = 1) +
        geom_raster() +
        stat_contour(color = "black", alpha = .25, size = .25) +
        # scale_fill_gradientn(colors = c("black", "blue", "red", "yellow")) +
        scale_fill_viridis_c() +
        labs(fill = "hillshade  ") +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        theme_classic() + 
        theme(legend.position = "none",
              strip.text = element_blank(),
              axis.text = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank())
ggsave("figures/downscale/topo_vars.png", p_topo,
       width = 12.5, height = 20, units = "in")
ggsave("figures/downscale/topo_vars_bw.png", 
       p_topo + scale_fill_gradientn(colors = c("black", "white")),
       width = 12.5, height = 20, units = "in")

px <- d %>%
        filter(between(x, quantile(x, .2), quantile(x, .4)),
               between(y, quantile(y, .25), quantile(y, .75))) %>%
        select(x, y, northness, eastness, southwestness, tpi, elevation) %>%
        gather(var, value, -x, -y, -elevation) %>%
        group_by(var) %>%
        mutate(value = scales::rescale(value),
               # value = ecdf(value)(value),
               var = factor(var, levels = c("northness", "eastness",
                                            "southwestness", "tpi"))) %>%
        ggplot(aes(x, y, fill = var, alpha = value, z = elevation)) +
        facet_wrap(~var, ncol = 2) +
        geom_raster(fill = "black", alpha = 1) +
        geom_raster() +
        stat_contour(color = "black", alpha = .7, size = .25) +
        # scale_fill_manual(values = c("magenta", "cyan", "green", "yellow")[c(2,3,1,4)]) +
        scale_fill_manual(values = c("magenta", "cyan", "green", "#ffa600")[c(2,3,1,4)]) +
        scale_alpha_continuous(range = 0:1) +
        labs(fill = "hillshade  ") +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        theme_classic() + 
        theme(legend.position = "none",
              strip.text = element_blank(),
              axis.text = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank())
ggsave("figures/downscale/topo_vars_v2.png", px,
       width = 12, height = 12, units = "in")
ggsave("figures/downscale/topo_vars_v2b.png", 
       px + facet_wrap(~var, ncol = 4),
       width = 24, height = 6, units = "in")



px <- d %>%
        filter(between(x, quantile(x, .2), quantile(x, .4)),
               between(y, quantile(y, .25), quantile(y, .75))) %>%
        mutate(eastness = -eastness) %>%
        mutate(rgb = colormap::colors3d(select(., eastness, northness, tpi),
                                        trans = "ecdf")) %>%
        ggplot(aes(x, y, z = elevation)) +
        geom_raster(aes(fill = rgb)) +
        stat_contour(color = "black", alpha = .7, size = .25) +
        scale_fill_identity() +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        theme_classic() + 
        theme(legend.position = "none",
              strip.text = element_blank(),
              axis.text = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank())
ggsave("figures/downscale/topo_vars_rgb.png", px,
       width = 6, height = 6, units = "in")

k <- 1.3
px <- d %>%
        filter(between(x, quantile(x, .2), quantile(x, .4)),
               between(y, quantile(y, .25), quantile(y, .75))) %>%
        mutate(eastness = -eastness) %>%
        mutate(eastness = scales::rescale(eastness)^k,
               northness = scales::rescale(northness)^k,
               tpi = scales::rescale(tpi)^k) %>%
        mutate(rgb = colormap::colors3d(select(., eastness, northness, tpi),
                                        trans = "fit")) %>%
        ggplot(aes(x, y, z = elevation)) +
        geom_raster(aes(fill = rgb)) +
        stat_contour(color = "black", alpha = .7, size = .25) +
        scale_fill_identity() +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        theme_classic() + 
        theme(legend.position = "none",
              strip.text = element_blank(),
              axis.text = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank())
ggsave("figures/downscale/topo_vars_rgb2.png", px,
       width = 6, height = 6, units = "in")



p_topography <- d %>%
        select(x, y, eastness, northness) %>%
        gather(topo, value, -x, -y) %>%
        ggplot(aes(x, y, fill = value)) +
        facet_grid(topo ~ .) +
        geom_raster() +
        scale_fill_gradientn(colors = c("black", "white")) +
        labs(fill = "topography")

p_macroclimate <- d %>%
        select(x, y, 
               bio5_int = bio5, bio6_int = bio6, bio12_int = bio12, 
               bio5_nn = bio5nn, bio6_nn = bio6nn, bio12_nn = bio12nn) %>%
        gather(clim, value, -x, -y) %>%
        separate(clim, c("clim", "stat")) %>%
        mutate(clim = factor(clim, levels = c("bio5", "bio6", "bio12")),
               stat = factor(stat, levels = c("nn", "int"))) %>%
        group_by(clim) %>%
        mutate(value = scale(value)) %>%
        ggplot(aes(x, y, fill = value)) +
        facet_grid(stat ~ clim) +
        geom_raster() +
        scale_fill_gradientn(colors = c("black", "blue", "red", "yellow")) +
        labs(fill = "macroclimate")

p_deltas <- d %>% mutate(h_n = northness * heat_northness,
                         h_e = eastness * heat_eastness,
                         c_n = northness * cold_northness,
                         c_e = eastness * cold_eastness,
                         m_n = northness * moisture_northness,
                         m_e = eastness * moisture_eastness) %>%
        select(x, y, h_n:m_e) %>%
        gather(var, value, -x, -y) %>%
        separate(var, c("clim", "topo")) %>%
        mutate(clim = factor(clim, levels = c("h", "c", "m"))) %>%
        ggplot(aes(x, y, fill = value)) +
        facet_grid(topo ~ clim) +
        geom_raster() +
        scale_fill_gradient2(low = "blue", mid = "gray90", high = "red") +
        labs(fill = "topoclimate\ndelta")

p_topoclimate <- d %>%
        select(x, y, h, c, m) %>%
        gather(clim, value, -x, -y) %>%
        mutate(clim = factor(clim, levels = c("h", "c", "m"))) %>%
        group_by(clim) %>%
        mutate(value = scale(value)) %>%
        ggplot(aes(x, y, fill = value)) +
        facet_grid(. ~ clim) +
        geom_raster() +
        scale_fill_viridis_c() +
        labs(fill = "topoclimate")

p <- guide_area() + plot_spacer() + p_macroclimate +
        plot_spacer() + plot_spacer() + plot_spacer() + 
        p_topography + plot_spacer() + p_deltas +
        plot_spacer() + plot_spacer() + plot_spacer() + 
        plot_spacer() + plot_spacer() + p_topoclimate +
        plot_layout(nrow = 5, 
                    widths = c(1, .02, 3), 
                    heights = c(2, .03, 2, .03, 1),
                    guides = "collect") &
        scale_x_continuous(expand = c(0, 0)) &
        scale_y_continuous(expand = c(0, 0)) &
        theme_void() & 
        theme(legend.direction = "horizontal", 
              legend.box = "vertical",
              legend.box.just = "right",
              strip.text = element_blank()) &
        guides(fill = guide_colorbar(barwidth = 8))
ggsave("figures/downscale/methods_grid.png", p,
       width = 18, height = 10, units = "in")


# ph <- ggplot(d, aes(x, y, fill = h)) +
#         geom_raster() +
#         scale_fill_gradientn(colors = c("darkred", "orange")) +
#         labs(fill = "heat (°C)  ")
# 
# pc <- ggplot(d, aes(x, y, fill = c)) +
#         geom_raster() +
#         scale_fill_gradientn(colors = c("cyan", "darkblue")) +
#         labs(fill = "cold (°C)  ")
# 
# pm <- ggplot(d, aes(x, y, fill = m)) +
#         geom_raster() +
#         scale_fill_gradientn(colors = c("darkkhaki", "forestgreen")) +
#         labs(fill = "moisture (mm)  ")
# 
# p <- pt + ph + pc + pm + 
#         plot_layout(nrow = 2) &
#         scale_x_continuous(expand = c(0, 0)) &
#         scale_y_continuous(expand = c(0, 0)) &
#         guides(fill = guide_colorbar(barwidth = 10)) &
#         theme_void() + theme(legend.position = "top")
# 
# ggsave("figures/downscale/combo.png", p,
#        width = 10, height = 10, units = "in")



## satelite photo ####

# satellite <- png::readPNG("figures/downscale/idaho_satellite_cropped.png")
satellite <- magick::image_read("figures/downscale/idaho_satellite_cropped.png")
aerial <- ggplot() + ggpubr::background_image(satellite)

scat <- magick::image_read("figures/downscale/scatterplot_3d.JPG")
scat <- ggplot() + ggpubr::background_image(scat)



## RGB plots ####

rgb_plot <- function(order, inversion, trans = "ecdf",
                     orientation = "wide",
                     outdir = "figures/downscale/"){
        
        
        pc <- d %>% select(heat, cold, moisture) %>% prcomp(scale. = T)
        y <- predict(pc, d %>% select(heat = bio5nn, cold = bio6nn, moisture = bio12nn))
        
        clr <- rbind(y, pc$x) %>% scale() %>%
                colormap::colors3d(order = order, inversion = inversion, 
                                   trans = "fit")
        
        d$color_macro <- clr[1:(length(clr)/2)]
        d$color_micro <- clr[(length(clr)/2+1):length(clr)]
        
        
        
        # sd <- d %>% sample_n(10000) %>% 
        #         mutate(h = rank(h), c = rank(c), m = rank(m)) %>%
        #         arrange(desc(m)) %>%
        #         mutate(m = plyr::round_any(m, 2000, ceiling))
        # sd %>% ggplot(aes(h, c)) +
        #         facet_grid(.~m, scales = "free", space = "free",
        #                    labeller = label_both) +
        #         geom_point(color = sd$color_micro) +
        #         theme_bw()
        
        
        
        prgb_M <- ggplot(d, aes(x, y)) +
                geom_raster(fill = d$color_macro) +
                theme_void()
        
        prgb_m <- ggplot(d, aes(x, y)) +
                geom_raster(fill = d$color_micro) +
                theme_void()
        
        
        # p <- (pt | aerial) / (prgb_M | scat) / (prgb_m) +
        #         plot_layout(heights = c(1, 1, 2), nrow = 3)
        p <- (pt | prgb_M) / (prgb_m) / (scat | aerial) +
                plot_layout(heights = c(1, 2, 1), nrow = 3)
        
        p <- p &
                scale_x_continuous(expand = c(0, 0)) &
                scale_y_continuous(expand = c(0, 0)) &
                theme_classic() + 
                theme(legend.position = "none",
                      axis.text = element_blank(),
                      axis.title = element_blank(),
                      axis.ticks = element_blank(),
                      axis.line = element_blank())
        
        ggsave(paste0(outdir, location, ".png"), 
               p, width = 25, height = 20, units = "in")
        
        # # zoom in
        # dd <- filter(d, 
        #              between(x, quantile(x, .4), quantile(x, .6)),
        #              between(y, quantile(y, .7), quantile(y, .9)))
        # ggplot(dd, aes(x, y)) +
        #         geom_raster(fill = dd$color_micro) +
        #         theme_void()
}


rgb_plot(4, 5)

# for(order in 1:6) for(inversion in 1:8) rgb_plot(order, inversion, orientation = "wide", 
#                                                  outdir = "figures/downscale/rgb/")




library("plotly")
sd <- d %>% sample_n(10000) #%>% mutate(h = rank(h), c = rank(c), m = rank(m))
fig <- plot_ly(sd, x = ~h, y = ~c, z = ~m,
               marker = list(color = sd$color_micro, size = 3)) %>% 
        add_markers() %>% 
        layout(scene = list(xaxis = list(title = 'heat (°C)'),
                            yaxis = list(title = 'cold (°C)'),
                            zaxis = list(title = 'moisture (mm)')),
               paper_bgcolor="#c2c2c2",
               plot_bgcolor="#c2c2c2")
fig
