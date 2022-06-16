
library(tidyverse)
library(patchwork)
library(grid)

source("code/utils.R")


# metadata on pre-processed model inputs
metadata <- readRDS("data/derived/binned/metadata.rds") %>%
        rename(data_file = file_path) %>%
        mutate(mod_vars = map(avars, function(x) c("bio12", x)),
               out_dir = paste0("data/derived/stan/", 
                                str_remove(basename(data_file), "\\.rds"))) %>%
        filter(dir.exists(out_dir)) %>%
        # select(data_file, out_dir, vars, topo_vars, mod_vars) %>%
        mutate(out_dir = str_replace(out_dir, "param", "pbs_d2_k8d3_le5"),
               s_knots = 8, s_degree = 3) %>%
        filter(file.exists(paste0(out_dir, "/pml/fit.rds")))



# load stan samples from MCMC or MLE 
load_samples <- function(fitdir){
        list.files(fitdir, pattern = "csv", full.names = T) %>%
                map_df(read_csv, comment = "#") %>%
                mutate(i = 1:nrow(.)) %>%
                gather(param, value, -i)
}

define <- function(var, vars, define = T){
        def <- function(var){
                list(bio5 = "summer max temp", 
                     bio6 = "winter min temp", 
                     bio12 = "precipitation",
                     bio1 = "tmean")[[var]]
        }
        recode(var,
               "1" = def(vars[1]),
               "2" = def(vars[2]),
               "3" = def(vars[3]))
        y <- recode(var,
                    "1" = vars[1],
                    "2" = vars[2],
                    "3" = vars[3])
        if(define) y <- recode(var,
                               "1" = def(vars[1]),
                               "2" = def(vars[2]),
                               "3" = def(vars[3]))
        return(y)
}

label_vars <- function(x){
        data.frame(z = recode(as.character(x[,1]),
                              "windward" = "windward exposure",
                              "tpi" = "prominence",
                              "tpis" = "prominence",
                              "bio12" = "precipitation",
                              "bio5" = "summer max temp",
                              "bio6" = "winter min temp"))
}

set_order <- function(x) x %>%
        mutate(clim_var = factor(clim_var, 
                                 levels = c("bio12", "bio5", "bio6", "MEAN")),
               topo_var = factor(topo_var, 
                                 levels = c("northness", "eastness", "windward", "tpi", "MEAN")))



GeomArrow2 <- ggplot2::ggproto("GeomArrow2", Geom,
                               required_aes = c("x", "y", "mag", "angle"),
                               default_aes = ggplot2::aes(color = "black", fill = "white", scale = 1),
                               draw_key = draw_key_polygon,
                               draw_panel = function(data, panel_scales, coord, scale = 1, 
                                                     theta = 15, length = .5, alpha = 1) {
                                       coords <- coord$transform(data, panel_scales)
                                       
                                       Mmag <- max(coords$mag)
                                       coords$mag <- with(coords, mag/Mmag*coords$scale)
                                       
                                       coords$dx <- with(coords, cos(angle)*mag)*scale
                                       coords$dy <- with(coords, sin(angle)*mag)*scale
                                       
                                       xx <- unit.c(unit(coords$x, "npc"), 
                                                    unit(coords$x, "npc") + unit(coords$dx, "snpc"))
                                       yy <- unit.c(unit(coords$y, "npc"), 
                                                    unit(coords$y, "npc") + unit(coords$dy, "snpc"))
                                       pol <- grid::polylineGrob(x = xx, y = yy,
                                                                 default.units = "npc",
                                                                 arrow = grid::arrow(angle = theta, type = "closed", 
                                                                                     length = unit(length, "lines")),
                                                                 gp = grid::gpar(col = coords$colour,
                                                                                 fill = coords$fill,
                                                                                 alpha = alpha),
                                                                 id = rep(seq(nrow(coords)), 2))
                                       pol
                                       
                               })

geom_arrow2 <- function(mapping = NULL, data = NULL, stat = "identity",
                        position = "identity", na.rm = FALSE, show.legend = NA,
                        inherit.aes = TRUE, scale = 1, ...) {
        layer(geom = GeomArrow2,
              mapping = mapping,
              data = data,
              stat = stat,
              position = position,
              show.legend = show.legend,
              inherit.aes = inherit.aes,
              params = list(na.rm = na.rm, scale = scale, ...)
        )
}




predict_deltas <- function(md, # model metadata 
                           g, # macroclimate data for which to predict deltas
                           subdir = "hmc"){ # hmc or pml
        
        vars <- md$vars[[1]]
        mod_vars <- md$mod_vars[[1]]
        topo_vars <- md$topo_vars[[1]]
        subplots <- md$subplots[[1]]
        d <- readRDS(md$data_file)
        
        dmd <- d$md
        dmv <- d$mv
        spid <- d$spid
        dmd$sp_id <- as.integer(factor(dmd$species))
        dmd <- left_join(dmd, select(spid, species))
        
        d <- load_samples(paste0(md$out_dir, "/", subdir))
        
        topo_mods1 <- c("mean", mod_vars)
        params <- tibble(mod = rep(topo_mods1, length(topo_vars)),
                         topo = rep(topo_vars, each = length(topo_mods1)))
        
        deltas <- d %>% 
                filter(str_detect(param, "delta")) %>%
                mutate(param = str_remove(param, "delta\\."),
                       param = str_replace_all(param, "\\.", "_")) %>%
                separate(param, c("param", "var"))
        
        n_bases <- (md$s_knots + md$s_degree - 1) ^ 2
        
        deltas <- deltas %>%
                rename(clim_var = var, basis = param) %>% 
                mutate(basis = as.integer(basis),
                       topo_var = topo_vars[ceiling(basis / n_bases)],
                       basis = basis %% n_bases,
                       basis = ifelse(basis == 0, n_bases, basis),
                       clim_var = vars[as.integer(clim_var)])
        
        # spline basis functions
        z <- tensor_splines(ecdf(dmd[[mod_vars[1]]])(g[[mod_vars[1]]]),
                            ecdf(dmd[[mod_vars[2]]])(g[[mod_vars[2]]]),
                            knots = md$s_knots, degree = md$s_degree,
                            xbounds = 0:1, ybounds = 0:1) %>% 
                as.data.frame() %>% as_tibble() %>%
                mutate(id = 1:nrow(.)) %>%
                gather(basis, basis_value, -id) %>%
                mutate(basis = as.integer(str_remove(basis, "V")))
        
        # full posterior
        gg <- deltas %>% 
                left_join(z) %>%
                mutate(delta = value * basis_value) %>%
                group_by(clim_var, topo_var, id, i) %>%
                summarize(delta = sum(delta)) %>%
                left_join(g) %>%
                ungroup()
        
        # rescale deltas to native units
        gg$delta_rscl <- NA
        for(v in topo_vars) gg$delta_rscl[gg$topo_var == v] <- 
                gg$delta[gg$topo_var == v] / dmv[[paste0(v, "_sd")]] 
        for(v in vars) gg$delta_rscl[gg$clim_var == v] <- 
                gg$delta_rscl[gg$clim_var == v] * dmv[[paste0(v, "_sd")]] 
        gg <- gg %>%
                mutate(delta_rscl = ifelse(clim_var == "bio12", log10inc(delta_rscl), delta_rscl),
                       delta_rscl = ifelse(clim_var == "bio12", delta_rscl * 100, delta_rscl),
                       delta_rscl = ifelse(topo_var == "tpi", delta_rscl * 100, delta_rscl),
                       bio1_rscl = bio1 * dmv$bio1_sd + dmv$bio1_mean,
                       bio12_rscl = bio12 * dmv$bio12_sd + dmv$bio12_mean,
                       bio12_rscl = 10^bio12_rscl - 1)
        return(gg)
}

plot_delta_splines <- function(md = slice(metadata, 5)){
        
        
        par <- basename(md$data_file) %>% str_remove("\\.rds")
        
        ## macroclimate grid ##
        g <- readRDS("data/derived/macroclimate_grid.rds") %>%
                select(bio1 = bio1r, bio12 = bio12r, color)
        for(b12 in unique(g$bio12)){
                x <- range(g$bio1[g$bio12 == b12]) 
                g <- bind_rows(g, tibble(bio1 = seq(x[1], x[2], 1), bio12 = b12))
        }
        dmv <- readRDS(md$data_file)$mv
        g <- distinct(g) %>%
                mutate(id = 1:nrow(.)) %>%
                mutate(bio1 = (bio1 - dmv$bio1_mean) / dmv$bio1_sd,
                       bio12 = (log10(bio12+1) - dmv$bio12_mean) / dmv$bio12_sd)
        
        
        ## deltas ##
        gg <- predict_deltas(md, g)
        
        
        
        ## posterior summaries ##
        ggg <- gg %>%
                group_by(clim_var, topo_var, bio1, bio12, bio1_rscl, bio12_rscl) %>%
                summarize(q025 = quantile(delta_rscl, .025),
                          q975 = quantile(delta_rscl, .975),
                          q025_raw = quantile(delta, .025),
                          q975_raw = quantile(delta, .975),
                          sig = sign(q025) == sign(q975),
                          delta = median(delta),
                          delta_rscl = median(delta_rscl),
                          color = color[1]) %>%
                ungroup()
        
        
        
        ## plot ##
        
        
        ref_lines <- ggg %>%
                filter(sig) %>%
                mutate(sign = sign(delta_rscl)) %>%
                group_by(clim_var, topo_var, sign) %>%
                filter(abs(delta) == max(abs(delta))) %>%
                mutate(round = case_when(clim_var == "bio12" ~ round(delta_rscl, 0),
                                         TRUE ~ round(delta_rscl, 1)),
                       clim_units = case_when(clim_var == "bio12" ~ "%",
                                              TRUE ~ "°C"),
                       topo_units = case_when(topo_var == "tpi" ~ "100m",
                                              topo_var == "tpis" ~ "sd",
                                              topo_var == "northness" ~ "N",
                                              topo_var == "eastness" ~ "E",
                                              topo_var == "windward" ~ "WE"),
                       text = paste0(ifelse(sign > 0, "+", ""),
                                     round, clim_units, " / ", topo_units))
        
        
        
        # p <- ggplot() +
        #         facet_grid(clim_var ~ topo_var, scales = "free",
        #                    labeller = label_vars) +
        #         geom_hline(data = ref_lines, 
        #                    aes(yintercept = delta), 
        #                    color = "gray60", size = .25, linetype = "dotted") +
        #         geom_label(data = ref_lines, 
        #                    aes(x = mean(range(ggg$bio1_rscl)), y = delta, label = text), 
        #                    color = "gray60", label.size = NA, size = 2.75) +
        #         geom_hline(yintercept = 0, size = .5, color = "gray40") +
        #         geom_ribbon(data = confidence,
        #                     aes(bio1_rscl, ymin = q025_raw, ymax = q975_raw,
        #                         fill = bio12_rscl, group = bio12_rscl),
        #                     alpha = .5) +
        #         geom_line(data = ggg, aes(bio1_rscl, delta, 
        #                                   color = bio12_rscl, group = bio12_rscl)) +
        #         scale_fill_gradientn(colors = c("orange", "forestgreen", "dodgerblue"),
        #                              trans = "log10", limits = range(ggg$bio12_rscl)) +
        #         scale_color_gradientn(colors = c("orange", "forestgreen", "dodgerblue"),
        #                               trans = "log10", limits = range(ggg$bio12_rscl)) +
        #         style +
        #         theme(panel.grid = element_blank()) +
        #         labs(x = "macro annual temperature (°C)",
        #              color = "macro annual precipitation (mm)",
        #              fill = "macro annual precipitation (mm)",
        #              y = "standardized effect of topographic variable on climate variable")
        # ggsave(paste0("figures/deltas/splines/", par, "_spline_lines_std.pdf"),
        #        p, width = 2 * length(topo_vars), height = 6, units = "in")
        
        
        ref_lines <- ref_lines %>% set_order()
        ggg <- ggg %>% set_order()
        
        
        p <- ggplot() +
                facet_grid(clim_var ~ topo_var, scales = "free",
                           labeller = label_vars) +
                geom_hline(data = ref_lines,
                           aes(yintercept = delta),
                           color = "gray60", size = .25, linetype = "dotted") +
                geom_label(data = ref_lines,
                           aes(x = mean(range(ggg$bio1_rscl)), y = delta, label = text),
                           color = "gray60", label.size = NA, size = 2.75) +
                geom_hline(yintercept = 0, size = .5, color = "gray40") +
                geom_ribbon(data = ggg,
                            aes(bio1_rscl, ymin = q025_raw, ymax = q975_raw,
                                # fill = bio12_rscl,
                                group = bio12_rscl),
                            alpha = .2) +
                geom_line(data = ggg, aes(bio1_rscl, delta,
                                          # color = bio12_rscl,
                                          group = bio12_rscl)) +
                geom_point(data = ggg %>% filter(!is.na(color)) %>%
                                   arrange(desc(bio12_rscl)),
                           aes(bio1_rscl, delta, fill = color),
                           shape = 21, size = 2) +
                scale_fill_identity() +
                style +
                theme(panel.grid = element_blank()) +
                labs(x = "macro annual temperature (°C)",
                     color = "macro annual precipitation (mm)",
                     fill = "macro annual precipitation (mm)",
                     y = "standardized effect of topographic variable on climate variable")
        ggsave(paste0("figures/deltas/splines/", par, "_spline_lines.pdf"),
               p, width = 2 * length(unique(ggg$topo_var)), height = 6, units = "in")
        
        
        
        ## marginal means ##
        
        for(fun in c("mean", "median")){
                fx <- switch(fun, "mean" = mean, "median" = median)
                
                mag_ct <- gg %>%
                        mutate(comp = ifelse(is.na(color), "spline", "point")) %>%
                        group_by(bio1_rscl, bio12_rscl, i, comp) %>%
                        summarize(delta = fx(abs(delta)),
                                  color = color[1]) %>%
                        group_by(bio1_rscl, bio12_rscl, comp) %>%
                        summarize(q025 = quantile(delta, .025),
                                  q975 = quantile(delta, .975),
                                  sig = sign(q025) == sign(q975),
                                  delta = median(delta),
                                  color = color[1]) %>%
                        ungroup() %>%
                        mutate(clim_var = "MEAN", topo_var = "MEAN")
                
                mag_c <- gg %>%
                        mutate(comp = ifelse(is.na(color), "spline", "point")) %>%
                        group_by(bio1_rscl, bio12_rscl, i, comp, clim_var) %>%
                        summarize(delta = fx(abs(delta)),
                                  color = color[1]) %>%
                        group_by(bio1_rscl, bio12_rscl, comp, clim_var) %>%
                        summarize(q025 = quantile(delta, .025),
                                  q975 = quantile(delta, .975),
                                  sig = sign(q025) == sign(q975),
                                  delta = median(delta),
                                  color = color[1]) %>%
                        ungroup() %>%
                        mutate(topo_var = "MEAN")
                
                mag_t <- gg %>%
                        mutate(comp = ifelse(is.na(color), "spline", "point")) %>%
                        group_by(bio1_rscl, bio12_rscl, i, comp, topo_var) %>%
                        summarize(delta = fx(abs(delta)),
                                  color = color[1]) %>%
                        group_by(bio1_rscl, bio12_rscl, comp, topo_var) %>%
                        summarize(q025 = quantile(delta, .025),
                                  q975 = quantile(delta, .975),
                                  sig = sign(q025) == sign(q975),
                                  delta = median(delta),
                                  color = color[1]) %>%
                        ungroup() %>%
                        mutate(clim_var = "MEAN")
                
                mag <- bind_rows(mag_ct, mag_c, mag_t) %>% set_order()
                
                
                ## plot ##
                p <- ggplot() +
                        facet_grid(clim_var ~ topo_var, scales = "free",
                                   labeller = label_vars) +
                        
                        geom_rect(data = mag,
                                  xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
                                  fill = "gray80") +
                        geom_rect(data = mag,
                                  xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0,
                                  fill = "gray50") +
                        
                        geom_hline(data = ref_lines, 
                                   aes(yintercept = delta), 
                                   color = "gray60", size = .25, linetype = "dotted") +
                        geom_label(data = ref_lines, 
                                   aes(x = mean(range(ggg$bio1_rscl)), y = delta, label = text), 
                                   color = "gray60", label.size = NA, size = 2.75) +
                        geom_hline(yintercept = 0, size = .5, color = "gray40") +
                        geom_ribbon(data = ggg,
                                    aes(bio1_rscl, ymin = q025_raw, ymax = q975_raw, group = bio12_rscl),
                                    alpha = .2) +
                        geom_line(data = ggg, 
                                  aes(bio1_rscl, delta, group = bio12_rscl)) +
                        geom_point(data = ggg %>% filter(!is.na(color)) %>%
                                           arrange(desc(bio12_rscl)), 
                                   aes(bio1_rscl, delta, fill = color),
                                   shape = 21, size = 2) +
                        
                        geom_ribbon(data = mag %>% filter(comp == "spline"),
                                    aes(bio1_rscl, ymin = q025, ymax = q975, group = bio12_rscl),
                                    alpha = .2) +
                        geom_line(data = mag %>% filter(comp == "spline"),
                                  aes(bio1_rscl, delta, group = bio12_rscl)) +
                        geom_point(data = mag %>% filter(comp == "point"),
                                   aes(bio1_rscl, delta, fill = color),
                                   shape = 21, size = 2) +
                        
                        scale_fill_identity() +
                        style +
                        theme(panel.grid = element_blank()) +
                        labs(x = "macro annual temperature (°C)",
                             color = "macro annual precipitation (mm)",
                             fill = "macro annual precipitation (mm)",
                             y = "standardized effect of topographic variable on climate variable")
                
                ggsave(paste0("figures/deltas/splines/", par, "_spline_lines_marginal_", fun, ".pdf"),
                       p, width = 2 * length(unique(ggg$topo_var)), height = 6, units = "in")
        }
        
        
        
        # aspect arrows
        ad <- ggg %>%
                select(clim_var, topo_var, bio1_rscl, bio12_rscl, delta) %>%
                mutate(delta = ifelse(clim_var == "bio12", -delta, delta)) %>%
                spread(topo_var, delta) %>%
                mutate(angle = atan2(northness, eastness),
                       magnitude = sqrt(northness^2 + eastness^2)) %>%
                filter(bio1_rscl %in% unique(bio1_rscl)[seq(1, 100, 5)])
        p <- ggplot(mapping = aes(bio1_rscl, bio12_rscl, angle = angle, mag = 1)) +
                geom_point(data = ad, aes(color = clim_var), size = 3, alpha = .6) +
                geom_point(data = ad, color = "white", size = 3) +
                geom_arrow2(data = ad %>% filter(clim_var == "bio12"),
                            color = "darkgreen", fill = "darkgreen", scale = .34, length = 2, theta = 10, alpha = .6) +
                geom_arrow2(data = ad %>% filter(clim_var == "bio5"),
                            color = "darkred", fill = "darkred", scale = .34, length = 2, theta = 10, alpha = .6) +
                geom_arrow2(data = ad %>% filter(clim_var == "bio6"),
                            color = "darkblue", fill = "darkblue", scale = .34, length = 2, theta = 10, alpha = .6) +
                annotate(geom = "rect", xmin = 13, xmax = 24.5, ymin = 2600, ymax = 7500,
                         fill = "white", color = "black", size = .25) +
                geom_point(data = ad, color = "black", fill = "white", shape = 21, size = 4.5) +
                geom_arrow2(data = data.frame(bio1_rscl = 15, bio12_rscl = 3200, angle = pi/2),
                            fill = "black", scale = .34, length = 2, theta = 10) +
                annotate(geom = "point", x = 15, y = 3200, 
                         color = "black", fill = "white", shape = 21, size = 4.5) +
                annotate(geom = "text", x = 15, y = 6000, label = "N",
                         fontface = "bold", size = 5) +
                style +
                coord_cartesian(ylim = c(100, 8100),
                                xlim = range(ad$bio1_rscl)) +
                scale_x_continuous(expand = c(.2, .2)) +
                scale_y_log10(expand = c(0, 0)) +
                scale_color_manual(labels = c("dryness", "heat", "cold"),
                                   values = c("darkgreen", "darkred", "darkblue")) +
                labs(x = "macro annual temperature (°C)",
                     y = "macro annual precipitation (mm)",
                     color = NULL) +
                theme(legend.direction = "vertical",
                      legend.position = c(.86, .83),
                      legend.background = element_blank(),
                      legend.key.size = unit(.02, "in"))
        ggsave(paste0("figures/deltas/splines/", par, "_aspect_arrows.pdf"),
               p, width = 4, height = 4, units = "in")
        
        
        # version showing posterior uncertainty
        ii <- sample(1:200, 100)
        ad <- gg %>%
                select(clim_var, topo_var, bio1_rscl, bio12_rscl, delta, id, i) %>%
                mutate(delta = ifelse(clim_var == "bio12", -delta, delta)) %>%
                spread(topo_var, delta) %>%
                mutate(angle = atan2(northness, eastness),
                       magnitude = sqrt(northness^2 + eastness^2)) %>%
                filter(bio1_rscl %in% unique(bio1_rscl)[seq(1, 100, 5)],
                       i %in% ii) %>%
                mutate(clim_var = recode(clim_var,
                                         "bio5" = "heat",
                                         "bio6" = "cold",
                                         "bio12" = "dryness"))
        p <- ggplot(mapping = aes(bio1_rscl, bio12_rscl, angle = angle, mag = 1)) +
                facet_wrap(~ clim_var, nrow = 1, labeller = label_vars) +
                geom_arrow2(data = ad,
                            color = "gray70", fill = "gray70", scale = .34, length = 2, theta = 3, alpha = .1) +
                geom_point(data = ad, color = "black", fill = "white", shape = 21, size = 4.5) +
                style +
                coord_cartesian(ylim = c(100, 8100),
                                xlim = range(ad$bio1_rscl)) +
                scale_x_continuous(expand = c(.2, .2)) +
                scale_y_log10(expand = c(0, 0)) +
                labs(x = "macro annual temperature (°C)",
                     y = "macro annual precipitation (mm)")
        ggsave(paste0("figures/deltas/splines/", par, "_aspect_arrows_posterior.png"),
               p, width = 12, height = 4, units = "in")
        
        
        
        # ordination
        vol <- ggg %>%
                filter(bio1_rscl %in% unique(bio1_rscl)[seq(1, 100, 5)]) %>%
                split(paste(.$bio1, .$bio12)) %>%
                map_df(function(x){
                        pc <- x %>%
                                select(clim_var, topo_var, delta) %>%
                                spread(topo_var, delta) %>%
                                select(-clim_var) %>%
                                prcomp(center = TRUE, scale. = FALSE)
                        data.frame(bio1 = x$bio1[1],
                                   bio12 = x$bio12[1],
                                   gini = ineq::Gini(pc$sdev),
                                   entropy = ineq::entropy(pc$sdev),
                                   ratio21 = pc$sdev[2] / pc$sdev[1],
                                   p1 = pc$sdev[1] / sum(pc$sdev),
                                   pc2 = pc$sdev[2],
                                   volume = 4/3 * pi * prod(pc$sdev / pc$sdev))
                })
        
        ggplot(vol, aes(bio1, bio12, fill = pc2)) +
                geom_tile() +
                scale_fill_viridis_c(limits = c(0,NA))
        
        
        
        
        
        
        ### strongest effect ###
        
        g2 <- expand_grid(bio1 = seq(min(g$bio1), max(g$bio1), length.out = 20),
                          bio12 = seq(min(g$bio12), max(g$bio12), length.out = 20)) %>%
                mutate(id = 1:nrow(.))
        gg2 <- predict_deltas(md, g2)
        ggg2 <- gg2 %>%
                group_by(clim_var, topo_var, bio1, bio12, bio1_rscl, bio12_rscl) %>%
                summarize(q025 = quantile(delta_rscl, .025),
                          q975 = quantile(delta_rscl, .975),
                          q025_raw = quantile(delta, .025),
                          q975_raw = quantile(delta, .975),
                          sig = sign(q025) == sign(q975),
                          delta = median(delta),
                          delta_rscl = median(delta_rscl)) %>%
                ungroup()
        ggg2 %>%
                group_by(bio1, bio12) %>%
                mutate(delta = abs(delta)) %>%
                filter(delta == max(delta)) %>%
                ggplot(aes(bio1, bio12, fill = paste0(clim_var, "\n", topo_var))) +
                geom_tile() +
                style
        
        
        
        
        return(NULL)
        
        
        p <- ggg %>%
                mutate(delta = ifelse(sig, delta, NA)) %>%
                ggplot(aes(bio1_rscl, bio12_rscl, fill = delta)) +
                facet_grid(clim_var ~ topo_var) +
                geom_tile() +
                scale_fill_gradientn(colors = c("darkred", "orange", "gold", "gray70",
                                                "lightblue", "dodgerblue", "darkblue"),
                                     values = c(0, .4, .46, .5, .54, .6, 1),
                                     limits = max(abs(gg$delta), na.rm = T) * c(-1, 1),
                                     na.value = "gray85") +
                scale_y_log10() +
                style +
                labs(x = "mean annual temperature",
                     y = "total annual precipitation")
        ggsave(paste0("figures/deltas/splines/", par, "_spline_heatmaps.pdf"),
               p, width = 2 * length(topo_vars), height = 6, units = "in")
        
        
        # mag <- gg %>%
        #         group_by(bio1_rscl, bio12_rscl, i) %>%
        #         summarize(c_var = clim_var[abs(delta) == max(abs(delta))],
        #                   t_var = topo_var[abs(delta) == max(abs(delta))],
        #                   var = paste(t_var, c_var))
        # 
        # mag %>%
        #         ggplot(aes(bio1_rscl, bio12_rscl, fill = var)) +
        #         geom_tile() +
        #         scale_y_log10()
        
        
        ## PCA ##
        
        library(patchwork)
        
        x <- gg %>%
                mutate(var = paste0(topo_var, clim_var)) %>%
                group_by(var, bio1_rscl, bio12_rscl) %>%
                summarize(delta = median(delta)) %>%
                mutate(var = str_replace(var, "northness", "N"),
                       var = str_replace(var, "eastness", "E"),
                       var = str_replace(var, "windward", "WW"),
                       var = str_replace(var, "bio", "b")) %>%
                spread(var, delta) %>%
                ungroup()
        pc <- prcomp(select(x, -bio1_rscl, -bio12_rscl))
        x <- x %>%
                mutate(pc1 = pc$x[,1],
                       pc2 = pc$x[,2])
        x$color <- colormap::colors2d(as.data.frame(select(x, bio1_rscl, bio12_rscl)),
                                      xtrans = "rank", ytrans = "rank",
                                      colors = c("green", "gold", "magenta", "blue"))
        rot <- as.data.frame(pc$rotation) %>%
                rownames_to_column("var")
        p_pc <- ggplot(x, aes(pc1, pc2)) +
                # geom_hline(yintercept = 0, color = "gray") +
                # geom_vline(xintercept = 0, color = "gray") +
                geom_path(aes(group = bio1_rscl, color = color)) +
                geom_path(aes(group = bio12_rscl, color = color)) +
                geom_point(aes(color = color)) +
                scale_color_identity() +
                geom_segment(data = rot, 
                             aes(x = 0, xend = PC1/6, y = 0, yend = PC2/6),
                             color = "black") +
                geom_text(data = rot, 
                          aes(x = PC1/6, y = PC2/6, label = var),
                          color = "black") +
                theme_classic() +
                coord_fixed()
        p_clim <- ggplot(x, aes(bio1_rscl, bio12_rscl, color = color)) +
                geom_path(aes(group = bio1_rscl)) +
                geom_path(aes(group = bio12_rscl)) +
                geom_point(aes(color = color)) +
                scale_color_identity() +
                scale_y_log10() +
                theme_classic()
        p <- p_clim + p_pc + plot_layout(nrow = 1, widths = c(1, 2))
        ggsave(paste0("figures/deltas/splines/", par, "_spline_PCA.pdf"),
               p, width = 7, height = 4, units = "in")
        
        
        # library(plotly)
        # gge <- expand(ggg, bio1, bio12) %>% left_join(ggg)
        # x <- gge %>% filter(clim_var == "bio12", topo_var == "northness") %>%
        #         select(bio1, bio12, delta) %>%
        #         spread(bio12, delta) %>%
        #         select(-bio1) %>%
        #         as.matrix()
        # fig <- plot_ly(x = unique(gge$bio12), y = unique(gge$bio1), z = x) %>% 
        #         add_surface()
        # fig
        
        
        # ggg %>%
        #         select(clim_var, topo_var, bio1_rscl, bio12_rscl, delta) %>%
        #         spread(topo_var, delta) %>%
        #         ggplot(aes(eastness, northness)) +
        #         facet_grid(. ~ clim_var) +
        #         geom_vline(xintercept = 0) +
        #         geom_hline(yintercept = 0) +
        #         geom_point() +
        #         coord_fixed() +
        #         style
        
}

metadata %>% slice(c(2, 5)) %>%
        # filter(str_detect(data_file, "param_9")) %>%
        split(1:nrow(.)) %>%
        map(plot_delta_splines)
# stop()






############### niches ##################

# unpack parameters
get_niches <- function(fitdir, vars = c("bio5", "bio6", "bio12")){
        
        load_samples(fitdir) %>%
                filter(str_detect(param, "alpha|mu|sigma|tau|Lomega")) %>%
                separate(param, c("param", "species", "var", "var2"), sep = "\\.") %>% 
                mutate(v = define(var, vars),
                       v2 = define(var2, vars),
                       var = define(var, vars, define = F),
                       var2 = define(var2, vars, define = F),
                       species = as.integer(species))
}


predict_suitability <- function(n, x){
        require(mvtnorm)
        fmu <- n %>% filter(param == "mu") %>% arrange(var) %>% pull(value)
        fsigma <- n %>% filter(param == "sigma") %>% 
                arrange(var, var2) %>% pull(value) %>% 
                matrix(nrow = length(fmu), byrow = T)
        falpha <- n %>% filter(param == "alpha") %>% pull(value)
        pmax <- dmvnorm(fmu, fmu, fsigma)
        dmvnorm(x, fmu, fsigma) / pmax * falpha
}


plot_topo_niche <- function(md = slice(metadata, 5)){
        
        
        sp <- readRDS(md$data_file)$spid
        d <- readRDS(md$data_file)
        mv <- d$mv
        
        
        # for every FIA plot including absences, predict suitability
        # then for each macroclimate bin, calcualte mean northness weighted by prob, 
        # vs mean northness weighted by true occ (i.e. just mean of presences)
        
        delta <- d$md %>%
                mutate(batch = rep(1:10, each = ceiling(nrow(.)/10))[1:nrow(.)]) %>%
                split(.$batch) %>%
                map_df(function(x){
                        x %>% mutate(id = 1:nrow(.)) %>%
                                predict_deltas(md, ., subdir = "pml")
                })
        
        m <- delta %>%
                mutate(id = as.integer(factor(paste(batch, id)))) %>%
                select(-bio1, -i, -batch, -delta_rscl) %>%
                unite(effect, topo_var, clim_var) %>%
                spread(effect, delta) %>%
                mutate(b5 = bio5 + northness * northness_bio5 + 
                               eastness * eastness_bio5 + 
                               windward * windward_bio5 + 
                               tpi * tpi_bio5,
                       b6 = bio6 + northness * northness_bio6 + 
                               eastness * eastness_bio6 + 
                               windward * windward_bio6 + 
                               tpi * tpi_bio6,
                       b12 = bio12 + northness * northness_bio12 + 
                               eastness * eastness_bio12 + 
                               windward * windward_bio12 + 
                               tpi * tpi_bio12)
        
        # load niche params
        niches <- get_niches(paste0(md$out_dir, "/pml"), vars = d$vars) %>% 
                filter(param %in% c("mu", "sigma", "alpha")) %>%
                left_join(rename(sp, gs = species, species = id)) %>%
                select(-species) %>% rename(species = gs)
        
        # occurrence probabilities, predicted and true
        m$pred <- m %>%
                split(.$species) %>%
                map(function(x){
                        niche <- filter(niches, species == x$species[1])
                        predict_suitability(niche, x %>% select(b12, b5, b6))
                }) %>%
                do.call("c", .)
        m$pred_macro <- m %>%
                split(.$species) %>%
                map(function(x){
                        niche <- filter(niches, species == x$species[1])
                        predict_suitability(niche, x %>% select(bio12, bio5, bio6))
                }) %>%
                do.call("c", .)
        m$true <- m$npres / m$n
        
        weighted.mean(m$pred[m$species == "Pseudotsuga menziesii"] / 
                              m$pred_macro[m$species == "Pseudotsuga menziesii"], 
                      m$npres[m$species == "Pseudotsuga menziesii"])
        
        
        # add niche means
        m <- niches %>% filter(param == "mu") %>%
                select(var, value, species) %>%
                mutate(var = paste0(var, "mu")) %>%
                spread(var, value) %>%
                left_join(m, .)
        
        # add niche variances
        m <- niches %>% filter(param == "sigma",
                               var == var2) %>%
                select(var, value, species) %>%
                mutate(var = paste0(var, "sigma"),
                       value = sqrt(value)) %>%
                spread(var, value) %>%
                left_join(m, .)
        
        # macroclimate diff relative to niche mean and sd
        m <- m %>% mutate(bio5d = (bio5 - bio5mu) / bio5sigma,
                          bio6d = (bio6 - bio6mu) / bio6sigma,
                          bio12d = (bio12 - bio12mu) / bio12sigma)
        
        sp_plot2 <- function(sp = "Pseudotsuga menziesii"){
                
                niche <- filter(niches, species == sp)
                mu <- filter(niche, param == "mu")
                sigma <- filter(niche, param == "sigma", var == var2) %>%
                        mutate(value = sqrt(value))
                
                g <- expand_grid(bio5 = seq(max(mu$value[mu$var == "bio5"] - 2 * sigma$value[sigma$var == "bio5"],
                                                min(d$md$bio5)),
                                            min(mu$value[mu$var == "bio5"] + 2 *sigma$value[sigma$var == "bio5"],
                                                max(d$md$bio5)),
                                            length.out = 15),
                                 bio12 = seq(max(mu$value[mu$var == "bio12"] - 2 * sigma$value[sigma$var == "bio12"],
                                                 min(d$md$bio12)),
                                             min(mu$value[mu$var == "bio12"] + 2 *sigma$value[sigma$var == "bio12"],
                                                 max(d$md$bio12)),
                                             length.out = 15),
                                 bio6 = c(-1, 0, 1),
                                 northness = seq(-4, 4, length.out = 11),
                                 eastness = seq(-4, 4, length.out = 11),
                                 tpi = 0,
                                 windward = 0) %>%
                        filter(sqrt(northness^2 + eastness^2) <= 4) %>%
                        mutate(bio1 = (bio5 + bio6) /2)
                
                g <- g %>% 
                        mutate(id = 1:nrow(.)) %>%
                        predict_deltas(md, ., subdir = "pml")
                
                g <- g %>%
                        select(-bio1, -i, -delta_rscl) %>%
                        unite(effect, topo_var, clim_var) %>%
                        spread(effect, delta) %>%
                        mutate(b5 = bio5 + northness * northness_bio5 + 
                                       eastness * eastness_bio5 + 
                                       windward * windward_bio5 + 
                                       tpi * tpi_bio5,
                               b6 = bio6 + northness * northness_bio6 + 
                                       eastness * eastness_bio6 + 
                                       windward * windward_bio6 + 
                                       tpi * tpi_bio6,
                               b12 = bio12 + northness * northness_bio12 + 
                                       eastness * eastness_bio12 + 
                                       windward * windward_bio12 + 
                                       tpi * tpi_bio12)
                
                g$pred <- predict_suitability(niche, g %>% select(b12, b5, b6))
                g$pred_macro <- predict_suitability(niche, g %>% select(bio12, bio5, bio6))
                
                gg <- g %>%
                        group_by(bio5, bio12, bio6) %>%
                        summarize(northness = weighted.mean(northness, pred),
                                  eastness = weighted.mean(eastness, pred),
                                  b5 = weighted.mean(b5, pred),
                                  b12 = weighted.mean(b12, pred),
                                  pred = mean(pred),
                                  pred_macro = mean(pred_macro)) %>%
                        ungroup() %>%
                        mutate(angle = atan2(northness, eastness),
                               pred_macro = pred_macro / max(pred),
                               pred = pred / max(pred),
                               
                               b5 = b5 * mv$bio5_sd + mv$bio5_mean,
                               b12 = 10 ^ (b12 * mv$bio12_sd + mv$bio12_mean) - 1,
                               
                               bio5 = bio5 * mv$bio5_sd + mv$bio5_mean,
                               bio6 = bio6 * mv$bio6_sd + mv$bio6_mean,
                               bio12 = 10 ^ (bio12 * mv$bio12_sd + mv$bio12_mean) - 1,
                               bio6 = paste0("macro cold = ", round(bio6, 0), " °C"),
                               bio6 = factor(bio6, levels = rev(unique(bio6))))
                
                
                p <- ggplot(gg, aes(bio5, bio12, angle = angle, fill = angle)) +
                        facet_grid(bio6 ~ .) +
                        geom_arrow2(data = filter(gg, pred_macro > .001),
                                    mag = 1, scale = .25, length = 1) +
                        geom_contour(aes(z = pred_macro), 
                                     breaks = c(.001, .01, .1, .5, .9),
                                     color = "black") +
                        style +
                        scale_fill_gradientn(colors = c("forestgreen", "dodgerblue", "red", 
                                                        "gold", "limegreen", "forestgreen")) +
                        coord_cartesian() +
                        scale_x_continuous(expand = c(0, 0)) +
                        scale_y_continuous(expand = c(0, 0), trans = "log10") +
                        theme(legend.position = "none") +
                        labs(x = "macro heat (°C)",
                             y = "macro moisture (mm)")
                ggsave(paste0("figures/pinwheel_", sp, ".png"), p,
                       width = 3, height = 7, units = "in")
                
                p <- ggplot(gg, aes(bio5, bio12, angle = angle, fill = angle)) +
                        facet_grid(bio6 ~ .) +
                        # geom_arrow2(data = filter(gg, pred_macro > .001),
                        #             mag = 1, scale = .1, length = 1) +
                        geom_contour(aes(z = pred_macro), 
                                     breaks = c(.001, .01, .1, .5, .9),
                                     color = "black") +
                        geom_segment(aes(bio5, bio12, xend = b5, yend = b12),
                                     arrow = arrow(angle = 15, type = "closed",
                                                   length = unit(.2, "in"))) +
                        style +
                        scale_fill_gradientn(colors = c("forestgreen", "dodgerblue", "red", 
                                                        "gold", "limegreen", "forestgreen"),
                                             limits = c(-pi, pi)) +
                        coord_cartesian() +
                        scale_x_continuous(expand = c(0, 0)) +
                        scale_y_continuous(expand = c(0, 0), trans = "log10") +
                        theme(legend.position = "none") +
                        labs(x = "macro heat (°C)",
                             y = "macro moisture (mm)")
                ggsave(paste0("figures/pinwheel_", sp, "_delta.png"), p,
                       width = 3, height = 7, units = "in")
        }
        
        c("Pseudotsuga menziesii", "Acer rubrum",
          "Pinus contorta", "Pinus ponderosa", "Nyssa sylvatica") %>%
                map(sp_plot2)
        
        
        
        sp_plot <- function(sp = "Pseudotsuga menziesii"){
                
                # summarize
                if(sp == "All species"){
                        means <- m
                        res <- 8
                        nmin <- 1000
                }else{
                        means <- m %>% filter(species == sp)
                        res <- 3
                        nmin <- 10
                }
                
                means <- means %>%
                        mutate(bio5 = plyr::round_any(bio5d, diff(range(bio5d))/res),
                               bio12 = plyr::round_any(bio12d, diff(range(bio12d))/res)) %>%
                        gather(model, prob, pred, true) %>%
                        mutate(weight = n * prob) %>%
                        group_by(model, bio5, bio12) %>%
                        summarize(northness = weighted.mean(northness, weight),
                                  eastness = weighted.mean(eastness, weight),
                                  windward = weighted.mean(windward, weight),
                                  tpi = weighted.mean(tpi, weight),
                                  bio6d = weighted.mean(bio6d, weight),
                                  prob = mean(prob),
                                  n = sum(weight)) %>%
                        ungroup()
                
                pd <- means %>%
                        select(-prob, -bio6d) %>%
                        group_by(bio5, bio12) %>%
                        mutate(n = n[model == "true"][1]) %>%
                        filter(n > 0) %>%
                        ungroup() %>%
                        gather(topo_var, value, northness:tpi) %>%
                        spread(model, value) %>% 
                        filter(n > nmin) %>%
                        mutate(color = colormap::colors2d(
                                select(., bio5, bio12) %>% as.data.frame(),
                                colors = c("forestgreen", "gold", "red", "cornflowerblue")))
                
                pd <- pd %>%
                        mutate(color = colormap::colorwheel2d(
                                select(., bio5, bio12) %>% as.data.frame(),
                                colors = c("gray60", "gold", "forestgreen", "cyan", "blue",
                                           "magenta", "red"),
                                origin = c(0, 0)))
                
                ppp <- pd %>%
                        gather(model, mean, pred, true) %>%
                        spread(topo_var, mean) %>%
                        mutate(angle = atan2(northness, eastness)) %>%
                        select(bio5, bio12, model, angle, color) %>%
                        spread(model, angle)
                
                pad <- .2
                p_angles <- ppp %>%
                        ggplot(mapping = aes(bio5, bio12, fill = color, mag = 1)) +
                        facet_wrap(~ "direction") +
                        geom_arrow2(aes(angle = true),
                                    fill = "black", 
                                    scale = .3, length = 1) +
                        geom_arrow2(aes(angle = pred),
                                    fill = ppp$color,
                                    scale = .3, length = 1) +
                        scale_fill_identity() +
                        style +
                        theme(panel.grid = element_blank()) +
                        theme(axis.title.y = element_blank(),
                              axis.text.y = element_blank(),
                              axis.ticks.y = element_blank(),
                              axis.title.x = element_text(hjust = 1)) +
                        labs(x = "macro heat niche position (sd)        ") +
                        scale_x_continuous(limits = range(pd$bio5) +
                                                   diff(range(pd$bio5)) * c(-pad, pad)) +
                        scale_y_continuous(limits = range(pd$bio12) +
                                                   diff(range(pd$bio12)) * c(-pad, pad))
                
                p_scatter <- pd %>%
                        mutate(clim_var = NA) %>% set_order() %>%
                        arrange(desc(n)) %>%
                        ggplot(aes(true, pred, color = color, fill = color)) +
                        facet_wrap(~ topo_var, labeller = label_vars, scales = "free", 
                                   ncol = 1, strip.position = "right") +
                        geom_vline(xintercept = 0) +
                        geom_hline(yintercept = 0) +
                        geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
                        geom_point(aes(size = n), shape = 21, color = "black") +
                        scale_color_identity() +
                        scale_fill_identity() +
                        style +
                        theme(axis.title.y = element_blank(),
                              axis.text.y = element_blank(),
                              axis.ticks.y = element_blank()) +
                        labs(x = "observed mean\ntopographic position")
                
                p_lines <- pd %>% 
                        mutate(clim_var = NA) %>% set_order() %>%
                        arrange(desc(n)) %>%
                        ggplot(aes(bio5, pred, color = color, fill = color)) +
                        facet_wrap(~ topo_var, labeller = label_vars, scales = "free", 
                                   ncol = 1, strip.position = "right") +
                        geom_hline(yintercept = 0) +
                        geom_line(aes(group = bio12)) +
                        geom_point(aes(size = n), shape = 21, color = "black") +
                        scale_color_identity() +
                        scale_fill_identity() +
                        style +
                        theme(legend.position = "none",
                              strip.text = element_blank()) +
                        labs(y = "predicted mean topographic position",
                             x = "macro heat\nniche position (sd)")
                
                p_clim <- pd %>% 
                        group_by(bio5, bio12) %>%
                        summarize(n = mean(n), color = color[1]) %>%
                        ggplot(aes(bio5, bio12, color = color, fill = color)) +
                        facet_grid(. ~ "color key") +
                        geom_path(aes(group = bio12)) +
                        geom_point(aes(size = n), shape = 21, color = "black") +
                        scale_color_identity() +
                        scale_fill_identity() +
                        style +
                        labs(x = NULL,
                             y = "macro moisture\nniche position (sd)") +
                        scale_x_continuous(limits = range(pd$bio5) +
                                                   diff(range(pd$bio5)) * c(-pad, pad)) +
                        scale_y_continuous(limits = range(pd$bio12) +
                                                   diff(range(pd$bio12)) * c(-pad, pad))
                
                p <- p_lines + p_scatter + p_clim + p_angles +
                        plot_annotation(
                                title = sp,
                                theme = theme(plot.title = element_text(hjust = .5))) +
                        plot_layout(design = "AB
                            CD",
                                    widths = c(1, 1),
                                    heights = c(5, 1)) &
                        theme(legend.position = "none",
                              axis.title = element_text(size = 10)) &
                        scale_size_continuous(trans = "log10")
                p
                
                ggsave(paste0("figures/topo_niche_", sp, ".png"), p,
                       width = 4, height = 8, units = "in")
        }
        
        c("All species", 
          "Pseudotsuga menziesii", "Acer rubrum",
          "Pinus contorta", "Pinus ponderosa", "Nyssa sylvatica")[2] %>%
                map(sp_plot)
        
        
        
        
        sp_plot <- function(sp = "Pseudotsuga menziesii"){
                
                means <- m %>% filter(species == sp)
                res <- 4
                nmin <- 100
                
                
                means <- means %>%
                        mutate(bio5 = plyr::round_any(bio5d, diff(range(bio5d))/res),
                               bio12 = plyr::round_any(bio12d, diff(range(bio12d))/res),
                               bio6 = plyr::round_any(bio6d, diff(range(bio6d))/3)) %>%
                        gather(model, prob, pred, true) %>%
                        mutate(weight = n * prob) %>%
                        group_by(model, bio5, bio12, bio6) %>%
                        summarize(northness = weighted.mean(northness, weight),
                                  eastness = weighted.mean(eastness, weight),
                                  windward = weighted.mean(windward, weight),
                                  tpi = weighted.mean(tpi, weight),
                                  # bio6d = weighted.mean(bio6d, weight),
                                  prob = mean(prob),
                                  n = sum(weight)) %>%
                        ungroup()
                
                pd <- means %>%
                        select(-prob) %>%
                        group_by(bio5, bio12, bio6) %>%
                        mutate(n = n[model == "true"][1]) %>%
                        filter(n > 0) %>%
                        ungroup() %>%
                        gather(topo_var, value, northness:tpi) %>%
                        spread(model, value) %>% 
                        filter(n > nmin)
                
                ppp <- pd %>%
                        gather(model, mean, pred, true) %>%
                        spread(topo_var, mean) %>%
                        mutate(angle = atan2(northness, eastness)) %>%
                        select(bio5, bio12, bio6, model, angle) %>%
                        spread(model, angle)
                
                
                pad <- .2
                p_angles <- ppp %>%
                        ggplot(mapping = aes(bio5, bio12, mag = 1)) +
                        facet_wrap(~ bio6) +
                        geom_arrow2(aes(angle = true),
                                    fill = "black", 
                                    scale = .3, length = 1) +
                        geom_arrow2(aes(angle = pred, fill = pred),
                                    scale = .3, length = 1) +
                        style +
                        scale_fill_gradientn(colors = c("forestgreen", "dodgerblue", "red", 
                                                        "gold", "forestgreen"),
                                             limits = c(-pi, pi)) +
                        theme(panel.grid = element_blank()) +
                        theme(axis.title.y = element_blank(),
                              axis.text.y = element_blank(),
                              axis.ticks.y = element_blank(),
                              axis.title.x = element_text(hjust = 1)) +
                        labs(x = "macro heat niche position (sd)        ") +
                        scale_x_continuous(limits = range(pd$bio5) +
                                                   diff(range(pd$bio5)) * c(-pad, pad)) +
                        scale_y_continuous(limits = range(pd$bio12) +
                                                   diff(range(pd$bio12)) * c(-pad, pad))
                p_angles
                
                
                
        }
        
        
}









####################################

plot_niches <- function(md = slice(metadata, 1)){
        
        d <- readRDS(md$data_file)
        
        # n <- get_niches(paste0(md$out_dir, "/hmc"))
        n <- get_niches(paste0("data/derived/stan/mu_con_par_1", "/hmc"))
        
        s <- n %>%
                filter(param == "mu") %>%
                group_by(species, var) %>%
                summarize(value = mean(value)) %>%
                rename(id = species) %>%
                left_join(d$spid)
        
        s %>%
                spread(var, value) %>%
                ecoclim::pairsData(xy_vars = unique(s$var),
                                   z_vars = "species") %>%
                mutate(group = case_when(str_detect(species, "Quercus") ~ "Quercus",
                                         str_detect(species, "Pinus") ~ "Pinus",
                                         TRUE ~ "other")) %>%
                ggplot(aes(x_value, y_value, color = group)) +
                facet_grid(y_var ~ x_var) +
                geom_point() +
                scale_color_manual(values = c("black", "red", "dodgerblue")) +
                theme_minimal() +
                labs(color = NULL)
        
        
        s <- n %>%
                filter((param == "sigma" & var == var2)
                       | param == "mu") %>%
                group_by(species, param, var) %>%
                summarize(value = median(value)) %>%
                rename(id = species) %>%
                left_join(d$spid)
        s %>%
                mutate(group = case_when(str_detect(species, "Quercus") ~ "Quercus",
                                         str_detect(species, "Pinus") ~ "Pinus",
                                         TRUE ~ "other")) %>%
                spread(param, value) %>%
                ggplot(aes(mu, sqrt(sigma), color = group, group = group)) +
                facet_wrap(~ var) +
                geom_point() +
                stat_ellipse(level = .75) +
                scale_y_log10() +
                scale_color_manual(values = c("black", "red", "dodgerblue")) +
                theme_minimal() +
                labs(color = NULL)
        
        
        
        n %>%
                filter((param == "sigma" & var == var2)
                       | param == "tau") %>%
                select(-var2, -v2) %>%
                spread(param, value) %>%
                sample_n(1000) %>%
                ggplot(aes(sqrt(sigma), tau)) +
                geom_point()
        
        
        
        predict_suitability <- function(sp, iter, samples, x){
                require(mvtnorm)
                
                f <- filter(samples, species == sp, i == iter)
                
                fmu <- f %>% filter(param == "mu") %>% arrange(var) %>% pull(value)
                fsigma <- f %>% filter(param == "sigma") %>% 
                        arrange(var, var2) %>% pull(value) %>% 
                        matrix(nrow = length(fmu), byrow = T)
                falpha <- f %>% filter(param == "alpha") %>% arrange(var) %>% pull(value)
                
                pmax <- dmvnorm(fmu, fmu, fsigma)
                
                x %>%
                        mutate(prob = dmvnorm(., fmu, fsigma) / pmax * falpha,
                               species = sp, i = iter)
        }
        
        
        g <- expand_grid(bio5 = seq(-5, 5, .4),
                         bio6 = seq(-5, 5, .4),
                         bio12 = seq(-5, 5, .4))
        
        # n0 <- n # original
        # n1 <- n # tau switched
        # n2 <- n # tau switched, Lomega constrained
        
        # high_tau <- n0 %>%
        #         filter(param == "tau", var %in% c("1", "3")) %>%
        #         group_by(species, var) %>%
        #         summarize(value = median(value)) %>% ungroup() %>%
        #         filter(value > quantile(value, .98))
        
        # f0 <- n0 %>%
        #         select(sp = species, iter = i) %>%
        #         filter(sp %in% high_tau$species) %>%
        #         distinct() %>%
        #         group_by(sp) %>%
        #         sample_n(10) %>%
        #         ungroup() %>%
        #         pmap_df(predict_suitability, samples = n0, x = g)
        
        n %>%
                filter(param == "sigma",
                       var != 2, var2 != 2) %>%
                group_by(i, species) %>%
                summarize(r = value[var == 1 & var2 == 3] /
                                  sqrt(prod(value[var == var2]))) %>% 
                pull(r) %>% hist()
        
        r1 <- n1 %>%
                filter(param == "sigma",
                       var != 2, var2 != 2) %>%
                group_by(i, species) %>%
                summarize(r = value[var == 1 & var2 == 3] /
                                  sqrt(prod(value[var == var2]))) %>%
                ungroup() %>%
                filter(abs(r) > .95)
        
        f <- n %>% 
                select(sp = species, iter = i) %>% 
                # filter(sp %in% high_tau$species) %>%
                filter(paste(iter, sp) %in% paste(r1$i, r1$species)) %>%
                # filter(iter == 1, sp == 83) %>%
                distinct() %>%
                group_by(sp) %>%
                sample_n(1) %>% 
                ungroup() %>%
                pmap_df(predict_suitability, samples = n, x = g)
        
        # f %>% group_by(bio5, bio6, bio12) %>% 
        #         summarize(prob = mean(prob)) %>% ungroup() %>%
        #         filter(bio6 == bio6[prob == max(prob)]) %>%
        #         ggplot(aes(bio5, bio12, fill = prob)) +
        #         geom_raster()
        
        colors <- colormap::distant_colors(length(unique(f$species)))
        
        f %>% group_by(species) %>%
                filter(bio6 == bio6[prob == max(prob)]) %>%
                mutate(prob = prob / max(prob)) %>%
                ggplot(aes(bio5, bio12, z = prob, 
                           color = factor(species), group = paste(species, i))) +
                geom_contour(breaks = .25, alpha = 1) +
                scale_color_manual(values = colors) +
                xlim(-5, 5) + ylim(-5, 5) +
                theme_minimal() +
                theme(legend.position = "none")
        
        bind_rows(mutate(n0, version = 0),
                  mutate(n1, version = 1),
                  mutate(n, version = 2)) %>%
                filter(param == "tau") %>%
                group_by(version) %>%
                summarize(value = quantile(value, .5))
        
        
        
        
        ### niche rotations ###
        
        stop("be sure that variable numbers are still b5, 6, 12")
        
        sig <- n %>%
                filter(param == "sigma",
                       var == 1, var2 == 3)
        sig %>%
                filter(abs(value) < 10) %>%
                ggplot(aes(value)) +
                geom_vline(xintercept = 0) +
                geom_histogram(boundary = 0) +
                geom_vline(xintercept = mean(sig$value), color = "red")
        
        
        # correlations 
        
        r <- n %>%
                filter(param == "sigma",
                       var != 2, var2 != 2) %>%
                group_by(i, species) %>%
                summarize(r = value[var == 1 & var2 == 3] /
                                  sqrt(prod(value[var == var2])))
        mur <- n %>% 
                filter(param == "mu", var != 2) %>%
                group_by(species) %>%
                filter(abs(mean(value)) < 2) %>% # rm extreme species
                ungroup() %>%
                select(i, species, var, value) %>%
                mutate(var = paste0("v", var)) %>%
                spread(var, value) %>%
                left_join(r)
        fit <- lm(r ~ v1 + v3 + v1*v3, data = mur)
        
        expand_grid(v1 = seq(-3, 3, .1), v3 = seq(-3, 3, .1)) %>%
                mutate(r = predict(fit, .)) %>%
                ggplot(aes(v1, v3, fill = r, z = r)) +
                geom_raster() +
                geom_contour(color = "black", 
                             breaks = seq(-1, 1, .1)) +
                scale_fill_gradient2(midpoint = 0)
        
        
        
        r <- n %>%
                filter(param == "sigma",
                       var != 3, var2 != 3) %>%
                group_by(i, species) %>%
                summarize(r = value[var == 1 & var2 == 2] /
                                  sqrt(prod(value[var == var2])))
        mur <- n %>% 
                filter(param == "mu", var != 3) %>%
                group_by(species) %>%
                filter(abs(mean(value)) < 3) %>% # rm extreme species
                ungroup() %>%
                select(i, species, var, value) %>%
                mutate(var = paste0("v", var)) %>%
                spread(var, value) %>%
                left_join(r)
        fit <- lm(r ~ v1 + v2 + v1*v2, data = mur)
        
        expand_grid(v1 = seq(-3, 3, .1), v2 = seq(-3, 3, .1)) %>%
                mutate(r = predict(fit, .)) %>%
                ggplot(aes(v1, v2, fill = r, z = r)) +
                geom_raster() +
                geom_contour(color = "black", 
                             breaks = seq(-1, 1, .1)) +
                scale_fill_gradient2(midpoint = 0)
        
        
        ### how many optima occur outside obs range? ###
        
        bbox <- d$md %>%
                filter(npres > 0) %>%
                select(species, bio5, bio6, bio12) %>%
                gather(var, value, -species) %>%
                group_by(species, var) %>%
                summarize(min = min(value), max = max(value)) %>%
                ungroup() %>%
                mutate(var = recode(var, "bio5" = "1", "bio6" = "2", "bio12" = "3"))
        
        mu <- n %>% 
                filter(param == "mu") %>%
                rename(id = species) %>%
                left_join(d$spid) %>%
                left_join(bbox) %>%
                mutate(inbounds = value > min & value < max,
                       dist = case_when(inbounds ~ 0,
                                        value < min ~ min - value,
                                        value > max ~ value - max))
        
        mean(mu$inbounds)
        hist(mu$dist[mu$dist < 5])
        
        
        
}



niche_data <- function(md){
        
        d <- readRDS(md$data_file)
        
        n <- get_niches(paste0(md$out_dir, "/hmc"))
        
        s <- n %>%
                group_by(param, species, var) %>%
                summarize(value = median(value)) %>%
                ungroup() %>%
                mutate(model = basename(md$out_dir))
        s
}

n <- metadata %>% #slice(1:3) %>% 
        split(1:nrow(.)) %>% 
        map_df(niche_data)

n %>%
        filter(param %in% c("mu", "tau"),
               species %in% unique(species)) %>%
        ggplot(aes(species, value, color = model)) +
        facet_grid(param + var ~ ., scales = "free") +
        geom_point()




