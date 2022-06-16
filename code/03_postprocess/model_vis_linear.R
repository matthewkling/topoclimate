


# unpack delta parameters
get_deltas <- function(fitdir, vars, topo_vars, mod_vars, 
                       translate = T, spline = T){
        
        d <- load_samples(fitdir)
        
        topo_mods1 <- c("mean", mod_vars)
        params <- tibble(mod = rep(topo_mods1, length(topo_vars)),
                         topo = rep(topo_vars, each = length(topo_mods1)))
        
        
        deltas <- d %>% 
                filter(str_detect(param, "delta")) %>%
                mutate(param = str_remove(param, "delta\\."),
                       param = str_replace_all(param, "\\.", "_")) %>%
                separate(param, c("param", "var"))
        
        if(translate) deltas <- deltas %>%
                mutate(var = define(var, vars),
                       topo = params$topo[as.integer(param)],
                       mod = params$mod[as.integer(param)],
                       mod = factor(mod, levels = topo_mods1))
        
        if(spline) deltas <- deltas %>%
                rename(clim_var = var, basis = param) %>% 
                mutate(basis = as.integer(basis),
                       topo_var = topo_vars[ceiling(basis / max(basis))],
                       basis = basis %% max(basis),
                       basis = ifelse(basis == 0, max(basis), basis),
                       clim_var = vars[as.integer(clim_var)])
        
        deltas
}

plot_deltas <- function(md = slice(metadata, 11)){
        
        ## unpack model metadata ##
        
        vars <- md$vars[[1]]
        mod_vars <- md$mod_vars[[1]]
        topo_vars <- md$topo_vars[[1]]
        subplots <- md$subplots[[1]]
        d <- readRDS(md$data_file)
        
        par <- basename(md$data_file) %>% str_remove("\\.rds")
        tag <- paste0(#"climate: ", paste(vars, collapse = ", "), "\n",
                #"topo: ", paste(topo_vars, collapse = ", "), "\n",
                "subplots: ", paste(subplots, collapse = ", "))
        
        
        ## load parameter estimates ##
        
        post <- get_deltas(paste0(md$out_dir, "/hmc"), 
                           vars, topo_vars, mod_vars)
        # mle <- get_deltas(paste0(md$out_dir, "/pml"), topo_vars, mod_vars)
        
        
        
        ## grid points plots ##
        
        # grid points covering forest climate space
        if("wspeed" %in% mod_vars){
                g <- expand_grid(bio1 = seq(-2, 2, 2),
                                 bio12 = seq(-2, 2, 2),
                                 wspeed = seq(-2, 2, 2))
        }else{
                g <- d$md %>%
                        filter(npres > 0) %>%
                        select(all_of(mod_vars)) %>%
                        mutate_all(round) %>%
                        distinct()
        }
        
        dd <- post %>%
                select(-param) %>%
                mutate(mod = paste0(mod, "_effect")) %>%
                spread(mod, value) %>%
                expand_grid(g) %>%
                mutate(effect = mean_effect +
                               bio1*bio1_effect + bio12*bio12_effect)
        if("wspeed" %in% mod_vars) dd <- mutate(dd, effect = effect + wspeed * wspeed_effect)
        
        stop("have coefs been misinterpreted? do they double-count topography since it's been premultiplied into the predictors?")
        
        s <- dd %>%
                group_by_at(c(mod_vars, "var", "topo")) %>%
                summarize(effect = median(effect)) %>% ungroup()
        s$hex <- colormap::colors2d(select(s, bio1, bio12) %>% as.data.frame(),
                                    c("forestgreen", "gold", "red", "dodgerblue"),
                                    xtrans = "rank", ytrans = "rank")
        
        if("wspeed" %in% mod_vars){
                p <- s %>%
                        ggplot(aes(bio1, effect, fill = hex, alpha = wspeed)) +
                        facet_grid(var ~ topo) +
                        geom_hline(yintercept = 0) +
                        geom_path(aes(group = paste(wspeed, bio12)), color = "black", size = .25) +
                        geom_path(aes(group = paste(wspeed, bio1)), color = "black", size = .25) +
                        geom_point(shape = 21)
                
        }else{
                p <- s %>%
                        ggplot(aes(bio1, effect, fill = hex)) +
                        facet_grid(var ~ topo) +
                        geom_hline(yintercept = 0) +
                        geom_path(aes(group = bio12), color = "black", size = .25) +
                        geom_path(aes(group = bio1), color = "black", size = .25) +
                        geom_point(shape = 21)
        }
        p <- p + scale_fill_identity() +
                scale_alpha_continuous(range = c(.4, 1), breaks = sq) +
                scale_x_continuous(breaks = sq) +
                style +
                theme(legend.position = "right",
                      panel.grid = element_blank()) +
                labs(y = "standardized effect",
                     x = "temperature",
                     color = "precipitation",
                     title = tag)
        ggsave(paste0("figures/deltas/", par, "_points_median.png"),
               p, width = 2 * length(topo_vars), height = 6, units = "in")
        ggsave(paste0("figures/deltas/", par, "_points_median.pdf"),
               p, width = 2 * length(topo_vars), height = 6, units = "in")
        
        
        # heatmaps
        g <- d$md %>%
                filter(npres > 0) %>%
                select(all_of(mod_vars)) %>%
                mutate_all(function(x) plyr::round_any(x, .25)) %>%
                distinct()
        dd <- post %>%
                select(-param) %>%
                mutate(mod = paste0(mod, "_effect")) %>%
                spread(mod, value) %>%
                expand_grid(g) %>%
                mutate(effect = mean_effect +
                               bio1*bio1_effect + bio12*bio12_effect)
        s <- dd %>%
                group_by_at(c(mod_vars, "var", "topo")) %>%
                summarize(sig = sign(quantile(effect, .025)) == sign(quantile(effect, .975)),
                          effect = median(effect)) %>% 
                ungroup()
        p <- s %>%
                ggplot(aes(bio1, bio12, 
                           z = effect,
                           fill = ifelse(sig, effect, NA))) +
                facet_grid(var ~ topo) +
                geom_tile() +
                # scale_fill_gradientn(colors = c("orange", "red", "black", "dodgerblue", "cyan"),
                #                      limits = max(abs(dd$effect)) * c(-1, 1)) +
                scale_fill_gradientn(colors = c("darkred", "orange", "gray70", "dodgerblue", "darkblue"),
                                     na.value = "gray90",
                                     limits = max(abs(dd$effect)) * c(-1, 1)) +
                theme_minimal() +
                theme(panel.grid = element_blank(),
                      strip.text = element_text(color = "white"),
                      strip.background = element_rect(fill = "black")) +
                labs(fill = "effect",
                     x = "temperature",
                     y = "precipitation",
                     title = tag)
        ggsave(paste0("figures/deltas/", par, "_heatmap.png"),
               p, width = 2 * length(topo_vars), height = 6, units = "in")
        
        
        return("mmmmkay")
        
        #############
        
        sq <- seq(-2, 2, 2)
        dd <- post %>%
                select(-param) %>%
                mutate(mod = paste0(mod, "_effect")) %>%
                spread(mod, value) %>%
                expand_grid(bio1 = seq(min(sq), max(sq), length.out = 20), 
                            bio12 = sq) %>%
                mutate(effect = mean_effect +
                               bio1*bio1_effect + 
                               bio12*bio12_effect) %>%
                arrange(sample(1:nrow(.)))
        
        p <- dd %>%
                ggplot(aes(bio1, effect, color = bio12, 
                           group = paste(i, bio12))) +
                facet_grid(var ~ topo) +
                geom_hline(yintercept = 0) +
                geom_line(alpha = .025) +
                scale_color_gradientn(colors = c("red", "black", "dodgerblue")) +
                scale_alpha_continuous(range = c(.4, 1), breaks = sq) +
                scale_x_continuous(breaks = sq) +
                style +
                theme(legend.position = "right") +
                labs(y = "standardized effect",
                     x = "temperature",
                     color = "precipitation",
                     title = tag)
        # ggsave(paste0("figures/deltas/", par, "_samples.png"),
        #        p, width = 6, height = 6, units = "in")
        
        p <- dd %>%
                group_by(bio1, bio12, var, topo) %>%
                summarize(effect = median(effect)) %>%
                ggplot(aes(bio1, effect, color = bio12, 
                           group = bio12)) +
                facet_grid(var ~ topo) +
                geom_hline(yintercept = 0) +
                geom_line() +
                scale_color_gradientn(colors = c("red", "black", "dodgerblue")) +
                scale_alpha_continuous(range = c(.4, 1), breaks = sq) +
                scale_x_continuous(breaks = sq) +
                style +
                theme(legend.position = "right") +
                labs(y = "standardized effect",
                     x = "temperature",
                     color = "precipitation",
                     title = tag)
        # ggsave(paste0("figures/deltas/", par, "_median.png"),
        #        p, width = 6, height = 6, units = "in")
        
        ##################
        
        
}

# metadata %>% split(1:nrow(.)) %>% map(plot_deltas)


#### splines, backup. ########


library(tidyverse)


source("code/utils.R")


# metadata on pre-processed model inputs
metadata <- readRDS("data/derived/binned/metadata.rds") %>%
        rename(data_file = file_path) %>%
        mutate(mod_vars = map(avars, function(x) c("bio12", x)),
               out_dir = paste0("data/derived/stan/", 
                                str_remove(basename(data_file), "\\.rds"))) %>%
        filter(dir.exists(out_dir))



# load stan samples from MCMC or MLE 
load_samples <- function(fitdir){
        list.files(fitdir, pattern = "csv", full.names = T) %>%
                map_df(read_csv, comment = "#") %>%
                mutate(i = 1:nrow(.)) %>%
                gather(param, value, -i)
}

define <- function(var, vars){
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
}

get_spline_deltas <- function(fitdir, vars, topo_vars, mod_vars,
                              g, dmd, dmv){
        
        d <- load_samples(fitdir)
        
        topo_mods1 <- c("mean", mod_vars)
        params <- tibble(mod = rep(topo_mods1, length(topo_vars)),
                         topo = rep(topo_vars, each = length(topo_mods1)))
        
        
        deltas <- d %>% 
                filter(str_detect(param, "delta")) %>%
                mutate(param = str_remove(param, "delta\\."),
                       param = str_replace_all(param, "\\.", "_")) %>%
                separate(param, c("param", "var"))
        
        deltas <- deltas %>%
                rename(clim_var = var, basis = param) %>% 
                mutate(basis = as.integer(basis),
                       topo_var = topo_vars[ceiling(basis / max(basis))],
                       basis = basis %% max(basis),
                       basis = ifelse(basis == 0, max(basis), basis),
                       clim_var = vars[as.integer(clim_var)])
        
        
        # spline basis functions
        z <- tensor_splines(ecdf(dmd[[mod_vars[1]]])(g[[mod_vars[1]]]),
                            ecdf(dmd[[mod_vars[2]]])(g[[mod_vars[2]]]),
                            knots = md$s_knots, degree = md$s_degree) %>% 
                as.data.frame() %>% as_tibble() %>%
                mutate(id = 1:nrow(.)) %>%
                gather(basis, basis_value, -g) %>%
                mutate(basis = as.integer(str_remove(basis, "V")))
        
        # full posterior
        gg <- deltas %>% 
                left_join(z) %>%
                mutate(delta = value * basis_value) %>%
                group_by(clim_var, topo_var, g, i) %>%
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
                       # delta_rscl = ifelse(topo_var == "tpis", delta_rscl * 4, delta_rscl),
                       bio1_rscl = bio1 * dmv$bio1_sd + dmv$bio1_mean,
                       bio12_rscl = bio12 * dmv$bio12_sd + dmv$bio12_mean,
                       bio12_rscl = 10^bio12_rscl - 1)
        return(gg)
        
}

spline_metadata <- metadata %>%
        # select(data_file, out_dir, vars, topo_vars, mod_vars) %>%
        mutate(out_dir = str_replace(out_dir, "param", "pbs_d2_k8d3_le5"),
               s_knots = 8, s_degree = 3) %>%
        filter(file.exists(paste0(out_dir, "/pml/fit.rds")))


plot_delta_splines <- function(md = slice(spline_metadata, 3)){
        
        
        ## unpack model metadata ##
        
        vars <- md$vars[[1]]
        mod_vars <- md$mod_vars[[1]]
        topo_vars <- md$topo_vars[[1]]
        subplots <- md$subplots[[1]]
        d <- readRDS(md$data_file)
        
        dmd <- d$md %>% filter(species %in% unique(species))
        dmv <- d$mv
        spid <- d$spid
        dmd$sp_id <- as.integer(factor(dmd$species))
        dmd <- left_join(dmd, select(spid, species))
        
        par <- basename(md$data_file) %>% str_remove("\\.rds")
        tag <- paste0(#"climate: ", paste(vars, collapse = ", "), "\n",
                #"topo: ", paste(topo_vars, collapse = ", "), "\n",
                "subplots: ", paste(subplots, collapse = ", "))
        
        
        ## load parameter estimates ##
        
        out_dir2 <- paste0(md$out_dir, "/hmc")
        if(!file.exists(paste0(out_dir2, "/fit.rds"))){
                out_dir2 <- paste0(md$out_dir, "/pml")
        }
        
        
        
        # mle <- get_deltas(paste0(md$out_dir, "/pml"), topo_vars, mod_vars)
        
        
        
        ## construct spline bases for mod_var grid ##
        
        # g <- dmd %>% 
        #         filter(npres > 0) %>%
        #         select(bio1, bio12) %>%
        #         mutate(bio1 = plyr::round_any(bio1, .1),
        #                bio12 = plyr::round_any(bio12, .5)) %>%
        #         distinct()
        # for(b12 in unique(g$bio12)){
        #         x <- range(g$bio1[g$bio12 == b12]) 
        #         g <- bind_rows(g, tibble(bio1 = seq(x[1], x[2], .1), bio12 = b12))
        # }
        # g <- distinct(g) %>%
        #         mutate(g = 1:nrow(.))
        
        # prediction grid
        g <- readRDS("data/derived/macroclimate_grid.rds") %>%
                select(bio1 = bio1r, bio12 = bio12r, color)
        for(b12 in unique(g$bio12)){
                x <- range(g$bio1[g$bio12 == b12]) 
                g <- bind_rows(g, tibble(bio1 = seq(x[1], x[2], 1), bio12 = b12))
        }
        g <- distinct(g) %>%
                mutate(id = 1:nrow(.)) %>%
                mutate(bio1 = (bio1 - dmv$bio1_mean) / dmv$bio1_sd,
                       bio12 = (log10(bio12+1) - dmv$bio12_mean) / dmv$bio12_sd)
        
        
        gg <- get_spline_deltas(out_dir2, 
                                vars, topo_vars, mod_vars, 
                                g, dmd, dmv)
        browser()
        
        
        ## fitted delta values ##
        
        
        
        
        # posterior summaries
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
        
        confidence <- ggg# %>%
        #filter(bio12 %in% c(-2, 0, 2))
        
        
        p <- ggplot() +
                facet_grid(clim_var ~ topo_var, scales = "free") +
                geom_hline(yintercept = 0, size = .5, color = "gray40") +
                geom_ribbon(data = confidence,
                            aes(bio1_rscl, ymin = q025, ymax = q975,
                                fill = bio12_rscl, group = bio12_rscl),
                            alpha = .5) +
                geom_line(data = ggg, aes(bio1_rscl, delta_rscl, 
                                          color = bio12_rscl, group = bio12_rscl)) +
                scale_fill_gradientn(colors = c("orange", "forestgreen", "dodgerblue"),
                                     trans = "log10", limits = range(ggg$bio12_rscl)) +
                scale_color_gradientn(colors = c("orange", "forestgreen", "dodgerblue"),
                                      trans = "log10", limits = range(ggg$bio12_rscl)) +
                style +
                labs(x = "macro annual temperature (°C)",
                     color = "macro annual precipitation (mm)",
                     fill = "macro annual precipitation (mm)",
                     y = "effect of topographic variable on climate variable")
        ggsave(paste0("figures/deltas/splines/", par, "_spline_lines.pdf"),
               p, width = 2 * length(topo_vars), height = 6, units = "in")
        
        
        
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
        
        
        label_vars <- function(x){
                data.frame(z = recode(as.character(x[,1]),
                                      "windward" = "windward exposure",
                                      "tpi" = "prominence",
                                      "tpis" = "prominence",
                                      "bio12" = "precipitation",
                                      "bio5" = "summer max temp",
                                      "bio6" = "winter min temp"))
        }
        
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
                geom_ribbon(data = confidence,
                            aes(bio1_rscl, ymin = q025_raw, ymax = q975_raw,
                                fill = bio12_rscl, group = bio12_rscl),
                            alpha = .5) +
                geom_line(data = ggg, aes(bio1_rscl, delta, 
                                          color = bio12_rscl, group = bio12_rscl)) +
                scale_fill_gradientn(colors = c("orange", "forestgreen", "dodgerblue"),
                                     trans = "log10", limits = range(ggg$bio12_rscl)) +
                scale_color_gradientn(colors = c("orange", "forestgreen", "dodgerblue"),
                                      trans = "log10", limits = range(ggg$bio12_rscl)) +
                style +
                theme(panel.grid = element_blank()) +
                labs(x = "macro annual temperature (°C)",
                     color = "macro annual precipitation (mm)",
                     fill = "macro annual precipitation (mm)",
                     y = "standardized effect of topographic variable on climate variable")
        ggsave(paste0("figures/deltas/splines/", par, "_spline_lines_std.pdf"),
               p, width = 2 * length(topo_vars), height = 6, units = "in")
        
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
                geom_ribbon(data = confidence,
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
        ggsave(paste0("figures/deltas/splines/", par, "_spline_lines_colors.pdf"),
               p, width = 2 * length(topo_vars), height = 6, units = "in")
        
        
        
        
        
        mag_ct <- gg %>%
                mutate(comp = ifelse(is.na(color), "spline", "point")) %>%
                group_by(bio1_rscl, bio12_rscl, i, comp) %>%
                summarize(delta = mean(abs(delta)),
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
                summarize(delta = mean(abs(delta)),
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
                summarize(delta = mean(abs(delta)),
                          color = color[1]) %>%
                group_by(bio1_rscl, bio12_rscl, comp, topo_var) %>%
                summarize(q025 = quantile(delta, .025),
                          q975 = quantile(delta, .975),
                          sig = sign(q025) == sign(q975),
                          delta = median(delta),
                          color = color[1]) %>%
                ungroup() %>%
                mutate(clim_var = "MEAN")
        
        
        set_order <- function(x) x %>%
                mutate(clim_var = factor(clim_var, 
                                         levels = c("bio12", "bio5", "bio6", "MEAN")),
                       topo_var = factor(topo_var, 
                                         levels = c("northness", "eastness", "windward", "tpi", "MEAN")))
        
        
        mag <- bind_rows(mag_ct, mag_c, mag_t) %>% set_order()
        ref_lines <- ref_lines %>% set_order()
        confidence <- confidence %>% set_order()
        ggg <- ggg %>% set_order()
        
        # mag %>%
        #         filter(comp == "spline") %>%
        #         ggplot() +
        #         facet_grid(clim_var ~ topo_var) +
        #         geom_hline(yintercept = 0, size = .5, color = "gray40") +
        #         geom_ribbon(aes(bio1_rscl, ymin = q025, ymax = q975,
        #                         group = bio12_rscl),
        #                     alpha = .2) +
        # geom_line(aes(bio1_rscl, delta,
        #               group = bio12_rscl)) +
        # geom_point(data = mag %>% filter(comp == "point"),
        #            aes(bio1_rscl, delta, fill = color),
        #            shape = 21, size = 3) +
        #         scale_fill_identity() +
        #         style
        
        
        
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
                geom_ribbon(data = confidence,
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
        ggsave(paste0("figures/deltas/splines/", par, "_spline_lines_margins.pdf"),
               p, width = 2 * length(topo_vars), height = 6, units = "in")
        
        
        
        
        
        
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

spline_metadata %>% #slice(3:4) %>% 
        filter(str_detect(data_file, "param_9")) %>%
        split(1:nrow(.)) %>%
        map(plot_delta_splines)
stop()





############### niches ##################

# unpack parameters
get_niches <- function(fitdir, vars = c("bio5", "bio6", "bio12")){
        
        load_samples(fitdir) %>%
                filter(str_detect(param, "alpha|mu|sigma|tau|Lomega")) %>%
                separate(param, c("param", "species", "var", "var2"), sep = "\\.") %>% 
                mutate(v = define(var, vars),
                       v2 = define(var2, vars),
                       species = as.integer(species))
}

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
        split(1:nrow(.)) %>% map_df(niche_data)

n %>%
        filter(param %in% c("mu", "tau"),
               species %in% unique(species)) %>%
        ggplot(aes(species, value, color = model)) +
        facet_grid(param + var ~ ., scales = "free") +
        geom_point()




