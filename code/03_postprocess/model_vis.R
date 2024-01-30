
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
        mutate(out_dir = str_replace(out_dir, "param", "pbs_d2_k8d3_le5"),
               out_dir = str_replace(out_dir, "le5_9", "le5_bulk_9"),
               s_knots = 8, s_degree = 3) %>%
        filter(file.exists(paste0(out_dir, "/pml/fit.rds")))



# load stan samples from MCMC or MLE 
load_samples <- function(fitdir){
        f <- list.files(fitdir, pattern = "csv", full.names = T)
        f <- f[!grepl("fit_summary", f)]
        f %>%
                map_df(read_csv, comment = "#") %>%
                mutate(i = 1:nrow(.)) %>%
                gather(param, value, -i)
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
                # mutate(id = 1:nrow(.)) %>%
                mutate(id = g$id) %>%
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


plot_trig <- function(){
        
        sa <- expand_grid(slope = seq(0, 90, 5),
                          aspect = seq(0, 180, 15)) %>%
                filter(! aspect %in% c(15, 165)) %>%
                mutate(northness = cos(aspect  / 180 * pi) * sin(slope / 180 * pi)) 
        
        p <- sa %>%
                ggplot(aes(slope, northness, group = aspect)) +
                geom_line(color = "darkred") +
                geom_text(data = filter(sa, slope == max(slope)) %>%
                                  mutate(aspect = paste0(aspect, "°"),
                                         aspect = ifelse(aspect == "0°", "0° (N)", aspect),
                                         aspect = ifelse(aspect == "180°", "180° (S)", aspect)),
                          aes(label = aspect),
                          hjust = 0, nudge_x = 1, color = "darkred") +
                annotate(geom = "text", x = 100, y = 0, label = "aspect", angle = -90, color = "darkred") +
                scale_x_continuous(breaks = seq(0, 90, 15),
                                   labels = function(x) paste0(x, "°"),
                                   expand = c(0, 0)) +
                coord_cartesian(clip = "off",
                                xlim = c(0, 90)) +
                style +
                theme(legend.position = "none",
                      plot.margin = unit(c(.01, .11, .01, .01), "npc"))
        
        sa <- expand_grid(slope = seq(0, 90, 15),
                          aspect = seq(0, 360, 1)) %>%
                filter(slope != 75) %>%
                mutate(northness = cos(aspect / 180 * pi) * sin(slope / 180 * pi),
                       eastness = sin(aspect / 180 * pi) * sin(slope / 180 * pi),
                       x = sin(aspect / 180 * pi) * sin(slope / 180 * pi),
                       y = cos(aspect / 180 * pi) * sin(slope / 180 * pi))
        
        s2 <- sqrt(2)/2
        cmps <- tibble(x = c(0, s2, 1, s2, 0, -s2, -1, -s2), 
                       y = c(1, s2, 0, -s2, -1, -s2, 0, s2),
                       expand = ifelse(x==0 | y==0, 1.1, 1.125),
                       label = c("N", "NE", "E", "SE", "S", "SW", "W", "NW"))
        
        p <- ggplot() +
                geom_segment(data = cmps, aes(x = 0, y = 0, xend = x, yend = y),
                             color = "darkred", alpha = .5) +
                geom_path(data = sa, aes(eastness, northness, group = slope)) +
                geom_label(data = cmps, aes(x*expand, y*expand, label = label),
                           color = "darkred") +
                geom_label(data = filter(sa, aspect == 30),
                           aes(x, y, 
                               label = paste0(slope, "°"))) +
                geom_text(data = filter(sa, aspect == 30, slope == 90),
                          aes(x*1.15, y*1.15), fontface = "bold",
                          label = "slope") +
                geom_text(data = filter(sa, aspect == 10, slope == 90),
                          aes(x*1.12, y*1.12),
                          fontface = "bold", color = "darkred",
                          label = "aspect") +
                scale_color_gradient(breaks = seq(0, 90, 15)) +
                scale_fill_gradientn(colors = c("black", "gray20", "gray50", "white"),
                                     values = c(0, .33, .66, 1)) +
                coord_fixed() +
                theme_minimal() +
                theme(legend.position = "none",
                      axis.title = element_text(face = "bold"),
                      axis.text = element_text(face = "bold")) +
                labs(x = "eastness", 
                     y = "northness")
        
        ggsave("figures/manuscript/northness_eastness.pdf", p, width = 6, height = 6, units = "in")
}
#plot_trig()


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
        
        ## plot splines ##
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
        
        ref_lines <- ref_lines %>% set_order()
        ggg <- ggg %>% set_order()
        
        p <- ggplot() +
                facet_grid(topo_var ~ clim_var, scales = "free",
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
        ggsave(paste0("figures/manuscript/", par, "_spline_lines.pdf"),
               p, width = 2 * length(unique(ggg$topo_var)), height = 6, units = "in")
        
        # annotated version (fig. 4d)
        p <- ggplot() +
                facet_grid(topo_var ~ clim_var, scales = "free",
                           labeller = label_vars) +
                geom_hline(data = ref_lines,
                           aes(yintercept = delta),
                           color = "gray70", size = .5) +
                geom_hline(data = ref_lines,
                           aes(yintercept = delta),
                           color = "gray30", size = .5, linetype = "dashed") +
                geom_label(data = ref_lines %>% filter(delta > 0),
                           aes(x = mean(range(ggg$bio1_rscl)), y = delta, label = text),
                           color = "white", label.size = NA, size = 2.5,
                           fill = "gray30", vjust = .3) +
                geom_label(data = ref_lines %>% filter(delta < 0),
                           aes(x = mean(range(ggg$bio1_rscl)), y = delta, label = text),
                           color = "white", label.size = NA, size = 2.5,
                           fill = "gray30", vjust = .7) +
                geom_hline(yintercept = 0, size = .5, color = "gray30") +
                geom_ribbon(data = ggg,
                            aes(bio1_rscl, ymin = q025_raw, ymax = q975_raw,
                                group = bio12_rscl),
                            alpha = .2) +
                geom_line(data = ggg, aes(bio1_rscl, delta,
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
                     y = "standardized effect of topographic variable on bioclimate variable")
        ggsave(paste0("figures/manuscript/", par, "_spline_lines_annotated.pdf"),
               p, width = 5, height = 9, units = "in")
        
        
        ## linear effect lines for example macroclimates (fig. 4c) ##
        z <- ggg %>% 
                filter(!is.na(color),
                       clim_var == "bio12",
                       topo_var == "northness") %>% 
                filter(bio12 == min(bio12) & bio1 == max(bio1) |
                               bio1 == min(bio1) & bio12 == min(bio12[bio1 == min(bio1)]) |
                               bio1 == max(bio1) & bio12 == max(bio12[bio1 == max(bio1)])) %>%
                arrange(desc(bio12_rscl)) %>%
                expand_grid(northness = seq(-1.5, 1.5, .1)) %>%
                mutate(clim = northness * delta_rscl,
                       clim_q025 = northness * q025,
                       clim_q975 = northness * q975)
        p <- z %>%
                ggplot(aes(northness, clim, 
                           ymin = clim_q025, ymax = clim_q975,
                           group = paste(bio1, bio12),
                           color = color)) +
                geom_vline(xintercept = 0, size = .5, color = "gray30") +
                geom_hline(yintercept = 0, size = .5, color = "gray30") +
                geom_ribbon(alpha = .2, color = NA) +
                geom_line(size = 1.5, color = "black") +
                geom_line(size = 1) +
                scale_color_identity() +
                style +
                scale_y_continuous(labels = function(x) paste0(ifelse(sign(x) == 1, "+",
                                                                      ifelse(sign(x) == 0, "", "-")), 
                                                               abs(x), "%")) +
                scale_x_continuous(expand = c(0, 0),
                                   breaks = c(-1, -.5, 0, .5, 1),
                                   labels = c("-1", "-0.5\n(S 30° or\nSE/SW 45°)", "0\n(level or\ndue E/W)", "0.5\n(N 30° or\nNE/NW 45°)", "1")) +
                coord_cartesian(xlim = c(-1, 1),
                                ylim = c(-40, 40)) +
                labs(x = "northness",
                     y = "bioclimate moisture anomaly\nfrom landscape macroclimate")
        ggsave(paste0("figures/manuscript/", par, "_topoclim_eg.pdf"),
               p, width = 4, height = 3.2, units = "in")
        
        
        ## aspect arrows ##
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
        ggsave(paste0("figures/manuscript/", par, "_aspect_arrows.pdf"),
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
        ggsave(paste0("figures/manuscript/", par, "_aspect_arrows_posterior.png"),
               p, width = 12, height = 4, units = "in")
}

metadata %>% 
        split(1:nrow(.)) %>%
        map(plot_delta_splines)




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
        dmvnorm(x[,sort(colnames(x))], fmu, fsigma) / pmax * falpha
}

# relationships between topography, macroclimate, and species occurrence
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
                mutate(b5 = bio5 + 
                               northness * northness_bio5 + 
                               eastness * eastness_bio5 + 
                               windward * windward_bio5 + 
                               tpi * tpi_bio5,
                       b6 = bio6 + 
                               northness * northness_bio6 + 
                               eastness * eastness_bio6 + 
                               windward * windward_bio6 + 
                               tpi * tpi_bio6,
                       b12 = bio12 + 
                               northness * northness_bio12 + 
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
        
        bell_curves <- function(){
                pal <- c("blue", "dodgerblue", "gray70", "orange", "red")        
                weighted.sd <- function(x, wt) sqrt(Hmisc::wtd.var(x, wt))
                msp <- m
                
                # all species
                p_a_d <- msp %>%
                        group_by(species) %>%
                        mutate(w = npres^2 / n,
                               bio5d = (bio5 - weighted.mean(bio5, w)) ,
                               bio6d = (bio6 - weighted.mean(bio6, w)) ,
                               bio12d = (bio12 - weighted.mean(bio12, w)) ) %>%
                        select(species,
                               n, npres, true, pred,
                               northness:windward,
                               bio5d:bio12d) %>%
                        gather(clim_var, clim, bio5d:bio12d) %>%
                        gather(topo_var, topo, northness:windward) %>%
                        mutate(clim = plyr::round_any(clim, .2),
                               topo = plyr::round_any(topo, 2)) %>%
                        group_by(clim_var, topo_var, clim, topo) %>%
                        summarize(p = sum(npres) / sum(n)) %>%
                        group_by(clim_var, topo, topo_var) %>%
                        mutate(p = (p + lag(p) + lead(p) + lag(p, 2) + lead(p, 2)) / 5) %>%
                        filter(abs(topo) <= 4) %>%
                        mutate(topo = factor(topo,
                                             levels = c(-4, -2, 0, 2, 4),
                                             labels = c("< -3", "-3 to -1", "-1 to 1", "1 to 3", "> 3"))) %>%
                        mutate(clim_var = str_remove(clim_var, "d")) %>%
                        mutate(clim = case_when(clim_var == "bio5" ~ clim * mv$bio5_sd,
                                                clim_var == "bio6" ~ clim * mv$bio6_sd,
                                                clim_var == "bio12" ~ clim * mv$bio12_sd)) %>%
                        filter(clim_var == "bio12" & abs(clim) <= .6 |
                                       clim_var == "bio5" & abs(clim) <= 8 |
                                       clim_var == "bio6" & abs(clim) <= 12) %>%
                        group_by(clim_var, topo, topo_var) %>%
                        mutate(p = p / max(p, na.rm = T)) %>%
                        set_order()
                
                padl <- p_a_d %>%
                        group_by(clim_var, topo_var, topo) %>%
                        summarize(clim = weighted.mean(clim, p))
                
                p_all <- p_a_d %>%
                        ggplot(aes(clim, p, color = topo, fill = topo, group = topo)) +
                        facet_grid(topo_var ~ clim_var, labeller = label_vars, scales = "free") +
                        geom_point(size = .35) +
                        geom_area(position = "identity", alpha = .1) +
                        scale_color_manual(values = pal) +
                        scale_fill_manual(values = pal) +
                        scale_x_continuous(expand = c(0, 0)) +
                        scale_y_continuous(expand = c(0, 0), limits = c(0, 1.05), breaks = c(0, .5, 1)) +
                        style +
                        theme(legend.position = "right",
                              panel.grid = element_blank()) +
                        guides(color = guide_legend(reverse = T),
                               fill= guide_legend(reverse = T)) +
                        labs(y = "relative occurrence frequency",
                             x = expression(paste(Delta, " moisture (log10 mm)       ", Delta, " high temperature (°C)      ", Delta, " low temperature (°C)")),
                             color = "topography\nz-score\n(stdev)",
                             fill = "topography\nz-score\n(stdev)")
                
                spp <- c("Quercus garryana", "Pinus contorta", "Ulmus rubra")
                p_spp <- function(sp){
                        msp %>%
                                filter(species == sp) %>%
                                group_by(species) %>%
                                mutate(bio5d = bio5,
                                       bio6d = bio6,
                                       bio12d = bio12) %>%
                                select(species,
                                       n, npres, true, pred,
                                       northness:windward,
                                       bio5d:bio12d) %>%
                                gather(clim_var, clim, bio5d:bio12d) %>%
                                gather(topo_var, topo, northness:windward) %>%
                                mutate(clim = plyr::round_any(clim, .2),
                                       topo = plyr::round_any(topo, 2)) %>%
                                group_by(species, clim_var, topo_var, clim, topo) %>%
                                summarize(p = sum(npres) / sum(n)) %>%
                                group_by(species, clim_var, topo, topo_var) %>%
                                arrange(clim) %>%
                                mutate(p = (p + lag(p) + lead(p) + lag(p, 2) + lead(p, 2)) / 5) %>%
                                filter(abs(topo) <= 4) %>%
                                mutate(p = p / max(p, na.rm = T)) %>%
                                mutate(topo = factor(topo,
                                                     levels = c(-4, -2, 0, 2, 4),
                                                     labels = c("< -3", "-3 to -1", "-1 to 1", "1 to 3", "> 3"))) %>%
                                mutate(clim_var = str_remove(clim_var, "d")) %>%
                                mutate(clim = case_when(clim_var == "bio5" ~ clim * mv$bio5_sd + mv$bio5_mean,
                                                        clim_var == "bio6" ~ clim * mv$bio6_sd + mv$bio6_mean,
                                                        clim_var == "bio12" ~ 10^(clim * mv$bio12_sd + mv$bio12_mean))) %>%
                                filter(species == spp[1] & between(clim, 350, 2500) |
                                               species == spp[2] & between(clim, 350, 1800) |
                                               species == spp[3] & between(clim, 700, 1800)) %>%
                                filter(clim_var == "bio12",
                                       topo_var == "northness") %>%
                                set_order() %>%
                                mutate(species = factor(species, level = spp)) %>%
                                ggplot(aes(clim, p, color = topo, fill = topo, group = topo)) +
                                facet_grid("northness ~ moisture" ~ species) +
                                geom_point(size = .35) +
                                geom_area(position = "identity", alpha = .1) +
                                scale_color_manual(values = pal) +
                                scale_fill_manual(values = pal) +
                                scale_x_continuous(expand = c(0, 0), trans = "log10", breaks = c(500, 700, 1000, 1400, 2000)) +
                                scale_y_continuous(expand = c(0, 0), limits = c(0, 1.05), breaks = c(0, .5, 1)) +
                                style +
                                theme(strip.text.x = element_text(face = "italic"),
                                      legend.position = "none") +
                                guides(color = guide_legend(reverse = T),
                                       fill= guide_legend(reverse = T)) +
                                labs(y = "relative occurrence frequency",
                                     x = "macroclimate moisture (mm)",
                                     color = "northness\n(stdev)",
                                     fill = "northness\n(stdev)")
                }
                
                p <- (p_spp(spp[1]) + theme(axis.title.x = element_blank(),
                                            axis.title.y = element_blank())) + 
                        (p_spp(spp[2]) + theme(axis.title.x = element_blank())) + 
                        (p_spp(spp[3]) + theme(axis.title.y = element_blank())) + 
                        (p_all + theme(axis.title.y = element_blank())) + 
                        plot_layout(design = "AD
                                    BD
                                    CD", 
                                    widths = c(1, 3))
                
                ggsave(paste0("figures/manuscript/fig3.pdf"), p,
                       width = 9.5, height = 7, units = "in")
                
                
                
                #### true vs predicted bell curves ####
                
                pp <- msp %>%
                        group_by(species) %>%
                        mutate(w = npres^2 / n,
                               bio5d = (bio5 - weighted.mean(bio5, w)) ,
                               bio6d = (bio6 - weighted.mean(bio6, w)) ,
                               bio12d = (bio12 - weighted.mean(bio12, w)) ) %>%
                        select(species,
                               n, npres, true, pred,
                               northness:windward,
                               bio5d:bio12d) %>%
                        gather(clim_var, clim, bio5d:bio12d) %>%
                        gather(topo_var, topo, northness:windward) %>%
                        mutate(clim = plyr::round_any(clim, .2),
                               topo = plyr::round_any(topo, 2)) %>%
                        group_by(clim_var, topo_var, clim, topo) %>%
                        summarize(true = weighted.mean(true, n),
                                  pred = weighted.mean(pred, n)) %>%
                        gather(model, p, true, pred) %>%
                        group_by(clim_var, topo, topo_var, model) %>%
                        mutate(p = (p + lag(p) + lead(p) + lag(p, 2) + lead(p, 2)) / 5) %>%
                        filter(abs(topo) <= 4)
                
                ppp <- pp %>%
                        mutate(topo = factor(topo,
                                             levels = c(-4, -2, 0, 2, 4),
                                             labels = c("< -3", "-3 to -1", "-1 to 1", "1 to 3", "> 3"))) %>%
                        mutate(clim_var = str_remove(clim_var, "d")) %>%
                        mutate(clim = case_when(clim_var == "bio5" ~ clim * mv$bio5_sd,
                                                clim_var == "bio6" ~ clim * mv$bio6_sd,
                                                clim_var == "bio12" ~ clim * mv$bio12_sd)) %>%
                        filter(clim_var == "bio12" & abs(clim) <= .6 |
                                       clim_var == "bio5" & abs(clim) <= 8 |
                                       clim_var == "bio6" & abs(clim) <= 12) %>%
                        group_by(clim_var, topo, topo_var, model) %>%
                        mutate(p = p / max(p, na.rm = T)) %>%
                        set_order() %>%
                        mutate(model = factor(model, 
                                              levels = c("true", "pred"),
                                              labels = c("observed", "predicted")))
                
                p1 <- ppp %>%
                        ggplot(aes(clim, p, color = topo, fill = topo, group = topo,
                                   linetype = model, shape = model)) +
                        facet_grid(topo_var+model ~ clim_var, 
                                   labeller = label_vars2,
                                   scales = "free") +
                        geom_point(size = .35) +
                        geom_area(position = "identity", alpha = .1, size = .25) +
                        scale_color_manual(values = pal) +
                        scale_fill_manual(values = pal) +
                        scale_shape_manual(values = c(1, 3)) +
                        scale_x_continuous(expand = c(0, 0)) +
                        scale_y_continuous(expand = c(0, 0), limits = c(0, 1.05), breaks = c(0, .5, 1)) +
                        style +
                        guides(color = guide_legend(reverse = T, order = 2),
                               fill = guide_legend(reverse = T, order = 2),
                               shape = guide_legend(override.aes = list(size = 1),
                                                    order = 1),
                               linetype = guide_legend(order = 1)) +
                        labs(y = "relative occurrence frequency",
                             x = expression(paste(Delta, " moisture (log10 mm)       ", Delta, " high temperature (°C)      ", Delta, " low temperature (°C)")),
                             color = "topographic\nvalue (stdev)  ",
                             fill = "topographic\nvalue (stdev)  ",
                             shape = "realm  ", linetype = "realm  ")
                
                p2 <- ppp %>%
                        spread(model, p) %>%
                        ggplot(aes(observed, predicted, group = clim)) +
                        facet_grid(topo_var ~ clim_var, 
                                   labeller = label_vars) +
                        geom_abline(size = .5, color = "gray") +
                        geom_point(size = .25, color = "gray40") +
                        geom_smooth(method = lm, se = F, 
                                    size = .4, color = "forestgreen") +
                        geom_smooth(method = lm, se = F, aes(group = 1),
                                    size = .5, color = "purple") +
                        scale_x_continuous(breaks = c(0, .5, 1),
                                           labels = c(0, "0.5", 1)) +
                        scale_y_continuous(breaks = c(0, .5, 1),
                                           labels = c(0, "0.5", 1)) +
                        style +
                        labs(x = "observed relative occurrence frequency",
                             y = "predicted relative occurrence probability")
                
                p <- p1 + p2 + 
                        plot_layout(nrow = 1, widths = c(2, 1.5)) +
                        plot_annotation(tag_levels = "a") &
                        theme(legend.box = "vertical",
                              legend.direction = "horizontal")
                ggsave(paste0("figures/manuscript/topo_occ_curves_obs_pred.pdf"), p,
                       width = 11, height = 9, units = "in")
                
        }
        bell_curves()
        
}
metadata %>% slice(2) %>% plot_topo_niche()

# species niche ellipsoids
plot_niches <- function(md){
        
        d <- readRDS(md$data_file)
        n <- get_niches(paste0(md$out_dir, "/pml"))
        
        s <- n %>%
                filter(param == "mu") %>%
                group_by(species, var) %>%
                summarize(value = mean(value)) %>%
                rename(id = species) %>%
                left_join(d$spid)
        
        predict_suitability <- function(sp, iter, samples, x){
                require(mvtnorm)
                
                f <- filter(samples, species == sp, i == iter)
                
                fmu <- f %>% filter(param == "mu") %>% arrange(var) %>% pull(value)
                fsigma <- f %>% filter(param == "sigma") %>% 
                        arrange(var, var2) %>% pull(value) %>% 
                        matrix(nrow = length(fmu), byrow = T)
                falpha <- f %>% filter(param == "alpha") %>% arrange(var) %>% pull(value)
                
                pmax <- dmvnorm(fmu, fmu, fsigma)
                
                x[,sort(colnames(x))] %>%
                        mutate(prob = dmvnorm(., fmu, fsigma) / pmax * falpha,
                               species = sp, i = iter)
        }
        
        sq <- seq(-3, 3, .1)
        g <- expand_grid(bio5 = sq,
                         bio6 = sq,
                         bio12 = sq)
        
        f <- n %>% 
                select(sp = species, iter = i) %>% 
                distinct() %>%
                group_by(sp) %>%
                sample_n(1) %>% 
                ungroup() %>%
                pmap_df(predict_suitability, samples = n, x = g) %>% 
                mutate(bio5 = bio5 * d$mv$bio5_sd + d$mv$bio5_mean,
                       bio6 = bio6 * d$mv$bio6_sd + d$mv$bio6_mean,
                       bio12 = bio12 * d$mv$bio12_sd + d$mv$bio12_mean,
                       bio12 = 10^bio12 - 1)
        
        
        cpal <- c("darkorchid4", "blue", "dodgerblue", "gray70", "orange", "red", "darkred", "darkred")
        fopt <- f %>%
                group_by(species) %>%
                filter(bio6 == bio6[prob == max(prob)])
        cval <- seq(0, 1, length.out = length(cpal))
        
        p_all <- f %>%
                group_by(species) %>%
                filter(bio6 == bio6[prob == max(prob)]) %>%
                mutate(prob = prob / max(prob)) %>%
                ungroup() %>%
                sample_n(nrow(.)) %>%
                ggplot(aes(bio5, bio12, z = prob, 
                           color = bio6, 
                           fill = bio6, 
                           group = species)) +
                facet_wrap(~ "All species") +
                stat_contour(geom = "polygon", alpha = .1, breaks = .5) +
                scale_color_gradientn(colors = cpal, values = cval, limits = range(fopt$bio6)) +
                scale_fill_gradientn(colors = cpal, values = cval, limits = range(fopt$bio6)) +
                scale_x_continuous(expand = c(0, 0)) +
                scale_y_continuous(expand = c(0, 0), trans = "log10") +
                guides(color = guide_colorbar(barwidth = 15)) +
                style +
                theme(axis.title.y = element_blank()) +
                labs(x = "micro high temperature (°C)",
                     y = "micro moisture (mm)",
                     color = "micro low temperature (°C)",
                     fill = "micro low temperature (°C)")
        
        
        sq <- seq(-3, 3, .02)
        g <- expand_grid(bio5 = sq,
                         bio6 = seq(-3, 3, .1),
                         bio12 = sq)
        
        spp <- c("Quercus garryana", "Pinus contorta", "Ulmus rubra")
        
        fs <- n %>% 
                select(sp = species, iter = i) %>% 
                distinct() %>%
                filter(sp %in% d$spid$id[d$spid$species %in% spp]) %>%
                group_by(sp) %>%
                sample_n(1) %>% 
                ungroup() %>%
                pmap_df(predict_suitability, samples = n, x = g) %>% 
                mutate(bio5 = bio5 * d$mv$bio5_sd + d$mv$bio5_mean,
                       bio6 = bio6 * d$mv$bio6_sd + d$mv$bio6_mean,
                       bio12 = bio12 * d$mv$bio12_sd + d$mv$bio12_mean,
                       bio12 = 10^bio12 - 1)
        
        p_spp <- fs %>% 
                rename(id = species) %>%
                left_join(d$spid) %>%
                filter(species %in% spp) %>%
                mutate(species = factor(species, levels = spp)) %>%
                group_by(species) %>%
                mutate(prob = prob / max(prob)) %>%
                ungroup() %>%
                sample_n(nrow(.)) %>%
                ggplot(aes(bio5, bio12, z = prob, 
                           color = bio6, 
                           fill = bio6, 
                           group = bio6)) +
                facet_wrap(~species, ncol = 1, scales = "free") +
                stat_contour(geom = "polygon", alpha = .05, breaks = .5) +
                scale_color_gradientn(colors = cpal, values = cval, limits = range(fopt$bio6)) +
                scale_fill_gradientn(colors = cpal, values = cval, limits = range(fopt$bio6)) +
                scale_y_continuous(trans = "log10") +
                style +
                labs(x = "micro high temperature (°C)",
                     y = "micro moisture (mm)",
                     color = "micro low temperature (°C)",
                     fill = "micro low temperature (°C)")
        
        p <- p_spp + theme(legend.position = "none") + 
                p_all +
                plot_layout(nrow = 1,
                            widths = c(1, 3),
                            guides = "collect") +
                plot_annotation(theme = theme(legend.position = "bottom"))
        
        ggsave(paste0("figures/manuscript/niches.pdf"), p,
               width = 9, height = 7, units = "in")
}
metadata %>% slice(2) %>% plot_niches()
