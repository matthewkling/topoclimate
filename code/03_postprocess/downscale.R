
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
rgb_plot <- function(order, inversion, trans = "fit", data = d,
                     outdir = "figures/downscale/rgb/",
                     test_plot = F, scatter = F, final = F){
        
        require(tidyverse)
        require(patchwork)
        
        data <- as.data.frame(rasterToPoints(tc))
        
        pc <- data %>% select(heat = high_temp, cold = low_temp, moisture) %>% prcomp(scale. = T)
        y <- predict(pc, data %>% select(heat = bio5nn, cold = bio6nn, moisture = bio12nn))
        
        clr <- rbind(y, pc$x) %>% scale() %>%
                colormap::colors3d(order = order, inversion = inversion, 
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
        
        scat <- magick::image_read("figures/downscale/scatterplot_3d_8_1.png")
        scat <- ggplot() + ggpubr::background_image(scat)
        
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
        ggsave(paste0(outdir, "cmb_", inversion, "_", order, ".png"), 
               p, width = 25, height = 20, units = "in")
        
        if(final) ggsave("figures/manuscript/downscale.png", 
                         p, width = 25, height = 20, units = "in")
        
        if(scatter){
                library("plotly")
                sd <- data %>% sample_n(10000)
                fig <- plotly::plot_ly(sd, x = ~high_temp, y = ~moisture, z = ~low_temp,
                                       marker = list(color = sd$color_micro, size = 3)) %>%
                        plotly::add_markers() %>%
                        plotly::layout(scene = list(xaxis = list(title = 'High temp. (°C)',
                                                                 dtick = 2, 
                                                                 tick0 = 24, 
                                                                 tickmode = "linear"),
                                                    zaxis = list(title = 'Low temp. (°C)   ',
                                                                 dtick = 2, 
                                                                 tick0 = -8, 
                                                                 tickmode = "linear"),
                                                    yaxis = list(title = 'Moisture (mm)',
                                                                 dtick = 100, 
                                                                 tick0 = 450, 
                                                                 tickmode = "linear")),
                                       paper_bgcolor = "#c2c2c2",
                                       plot_bgcolor = "#c2c2c2")
                fig 
        }
}
rgb_plot(1, 8, final = T)
