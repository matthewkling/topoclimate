
# figure 5

library(tidyverse)
library(raster)
library(patchwork)
select <- dplyr::select


# data ##############

# from https://apps.nationalmap.gov/downloader/#/
# metadata at https://www.sciencebase.gov/catalog/item/5f7784fc82ce1d74e7d6cbec
elev <- raster("data/raw/DEM/USGS_13_n43w113.tif")
names(elev) <- "elevation"
ext <- extent(-112.7458, -112.5855, 42.53036, 42.5841)
e <- crop(elev, ext)

# bounding box as polygon (for import into GE; not used)
ep <- as(ext* 1.05, "SpatialPolygons") %>% 
        SpatialPolygonsDataFrame(data = data.frame(x = 4))
rgdal::writeOGR(ep, paste0("data/raw/DEM/", "bbox", ".kml"), "boundary", driver="KML",
                overwrite_layer = T)

# bounding box size in km (for reference; not used)
size <- c(geosphere::distm(c(ext@xmin, ext@ymin), c(ext@xmax, ext@ymin)),
          geosphere::distm(c(ext@xmin, ext@ymin), c(ext@xmin, ext@ymax))) / 1000

# topoclimate
et <- crop(elev, ext * 1.5) # pad DEM to avoid edge effects
tc <- topoclimate.pred::topoclimate(et, include_inputs = TRUE)
tc <- crop(tc, e)

# additional carto vars
terr <- terrain(et, c("slope", "aspect"))
hill <- hillShade(terr$slope, terr$aspect)
tc$hillshade <- crop(hill, tc)
tc$elevation <- crop(elev, tc)
macro <- topoclimate.pred::macroclimate(e, "ngb")
tc <- stack(tc, setNames(macro, paste0(names(macro), "nn")))

# figure
rgb_plot <- function(order, inversion, trans = "none", data = d,
                     outdir = "figures/downscale/rgb/",
                     test_plot = F, final = F){
        
        require(tidyverse)
        require(patchwork)
        
        data <- as.data.frame(rasterToPoints(tc))
        
        pc <- data %>% select(heat = high_temp, cold = low_temp, moisture) %>% prcomp(scale. = T)
        y <- predict(pc, data %>% select(heat = bio5nn, cold = bio6nn, moisture = bio12nn))
        
        clr <- rbind(y, pc$x) %>% scale() %>%
                colors3d::colors3d(order = order, inversion = inversion, 
                                   trans = trans)
        data$color_macro <- clr[1:(length(clr)/2)]
        data$color_micro <- clr[(length(clr)/2+1):length(clr)]
        
        prgb_M <- ggplot(data, aes(x, y)) +
                geom_raster(fill = data$color_macro) +
                theme_void()
        
        prgb_m <- ggplot(data, aes(x, y)) +
                geom_raster(fill = data$color_micro) +
                theme_void()
        
        if(test_plot){
                p <- prgb_M / prgb_m
                ggsave(paste0(outdir, "RGB_", inversion, "_", order, ".png"), 
                       p, width = 25, height = 20, units = "in")
                return("done")
        }
        
        pt <- ggplot(data, aes(x, y, fill = hillshade, z = elevation)) +
                geom_raster() +
                stat_contour(color = "black", alpha = .25, size = .25) +
                scale_fill_gradientn(colors = c("black", "white"))
        
        satellite <- magick::image_read("figures/downscale/idaho_satellite_cropped.png")
        aerial <- ggplot() + ggpubr::background_image(satellite)
        
        # plot micro v macro
        sd <- data %>% sample_n(100000)
        s1 <- ggplot(sd, 
                     aes(high_temp, moisture, color = color_micro)) +
                geom_point() +
                geom_point(aes(bio5nn, bio12nn), color = "black") +
                scale_color_identity() +
                labs(x = "High temp. (°C)",
                     y = "Moisture (mm)")
        s2 <- ggplot(sd, 
                     aes(high_temp, low_temp, color = color_micro)) +
                geom_point() +
                geom_point(aes(bio5nn, bio6nn), color = "black") +
                scale_color_identity() +
                labs(x = "High temp. (°C)",
                     y = "Low temp. (°C)")
        s3 <- ggplot(sd, 
                     aes(low_temp, moisture, color = color_micro)) +
                geom_point() +
                geom_point(aes(bio6nn, bio12nn), color = "black") +
                scale_color_identity() +
                labs(x = "Low temp. (°C)",
                     y = "Moisture (mm)")
        s <- s1 + s2 + s3 +
                plot_layout(nrow = 1, widths = rep(1, 3)) & 
                theme_minimal(base_size = 25) & 
                theme(plot.margin = unit(c(0, 0, 0, 0), "in"))
        tmp <- "figures/fig5d.png"
        ggsave(tmp, s, width = 12.5, height = 5, units = "in")
        scat <- ggplot() + ggpubr::background_image(magick::image_read(tmp))
        file.remove(tmp)
        
        p <- aerial + pt + prgb_M + prgb_m + scat + 
                plot_layout(ncol = 1)
        p <- p &
                scale_x_continuous(expand = c(0, 0)) &
                scale_y_continuous(expand = c(0, 0)) &
                theme_classic() + 
                theme(legend.position = "none",
                      axis.text = element_blank(),
                      axis.title = element_blank(),
                      axis.ticks = element_blank(),
                      axis.line = element_blank())
        if(!final) ggsave(paste0(outdir, "cmb_", inversion, "_", order, ".png"), 
                          p, width = 10, height = 20, units = "in")
        
        if(final) ggsave("figures/manuscript/downscale.png", 
                         p, width = 10, height = 20, units = "in")
        
        # check how many micro values are outside macro hull
        hull <- data %>% select(macro_bio5, macro_bio6, macro_bio12) %>% geometry::convhulln()
        data$in_hull <- data %>% select(high_temp, low_temp, moisture) %>% 
                as.matrix() %>% geometry::inhulln(hull, .)
        out <- 1 - mean(data$in_hull)
}

rgb_plot(1, 8, final = T)
