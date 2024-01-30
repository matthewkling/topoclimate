
# Simplest version of model: a single-species analysis with stationary topoclimate effects

library(cmdstanr)
library(tidyverse)

source("code/utils.R")


minimal_model_mcmc <- function(sp = "Juniperus osteosperma",
                      data_file = "data/derived/binned/param_9.rds",
                      mod_file = "code/02_model/minimal_model.stan",
                      vars = c("bio5", "bio6", "bio12"),
                      mod_vars = c("bio12", "bio1"),
                      topo_vars = c("northness", "eastness", "tpi", "windward")
){
        d <- readRDS(data_file)
        md <- d$md %>% filter(species == sp)
        mv <- d$mv
        
        # macroclimate and topo matrices
        m <- md %>% select(all_of(vars)) %>% as.matrix()
        z <- md %>% select(all_of(topo_vars)) %>% as.matrix()
        
        # package data for stan
        dat <- list(
                K = ncol(m),
                N = nrow(md),
                D = ncol(z),
                z = z,
                m = m,
                nn = md$n,
                y = md$npres
        )
        
        # fit model
        model <- cmdstan_model(mod_file)
        init <- function() list(tau = rep(1, dat$K),
                                alpha = .25,
                                delta = matrix(0, dat$D, dat$K),
                                Lomega = matrix(0, dat$K, dat$K))
        fit <- model$sample(data = dat, init = init,
                            chains = 5, parallel_chains = 5)
        
        # format posterior samples
        fit$draws() %>% 
                as.data.frame() %>% as_tibble() %>%
                mutate(species = sp,
                       iter = 1:nrow(.)) %>%
                gather(var, value, -species, -iter) %>%
                separate(var, c("chain", "var"), "\\.") %>%
                filter(str_detect(var, "delta")) %>%
                mutate(var = str_remove(var, "delta\\["),
                       var = str_remove(var, "]")) %>%
                separate(var, c("topo", "clim"), remove = F) %>%
                mutate(topo_var = topo_vars[as.integer(topo)],
                       clim_var = vars[as.integer(clim)])
}


minimal_model_map <- function(sp = "Juniperus osteosperma",
                         topo_vars = c("northness", "eastness", "tpi", "windward")){
        
        ## setup ##
        # message(sp, fold)
        temp_dir <- "data/temp"
        data_file <- "data/derived/binned/param_9.rds"
        mod_file <- "code/02_model/minimal_model.stan"
        clim_vars <- c("bio5", "bio6", "bio12")
        mod_vars <- c("bio12", "bio1")
        d <- readRDS(data_file)
        md <- d$md %>% filter(species == sp)
        m <- md %>% select(all_of(clim_vars)) %>% as.matrix()
        z <- md %>% select(all_of(topo_vars)) %>% as.matrix()
        
        md$id <- 1:nrow(m)
        
        train <- rep(T, nrow(m)) 
        
        ##  step 1: fit bioclimate model ##
        dat <- list(
                K = ncol(m),
                N = nrow(m[train, ]),
                D = ncol(z),
                z = z[train, ],
                m = m[train, ],
                nn = md$n[train],
                y = md$npres[train]
        )
        init <- function() list(mu = rep(0, dat$K),
                                tau = rep(1, dat$K),
                                alpha = .25,
                                delta = matrix(0, dat$D, dat$K),
                                Lomega = matrix(0, dat$K, dat$K))
        model <- cmdstan_model(mod_file)
        bioclim_fit <- model$optimize(data = dat, 
                                      init = list(init()),
                                      output_dir = temp_dir,
                                      refresh = 0)
        
        bioclim_fit %>% 
                get_samples() %>% 
                get_deltas(topo_vars, clim_vars) %>%
                mutate(species = sp,
                       npres = sum(md$npres),
                       bio12 = weighted.mean(md$bio12, md$npres),
                       bio1 = weighted.mean((md$bio5 + md$bio6)/2, md$npres))
}




# the three example species featured in paper (MCMC)
mcmc4 <- c("Quercus garryana", "Pinus contorta", "Ulmus rubra") %>%
        map_df(minimal_model_mcmc,
               topo_vars = c("northness", "eastness", "tpi", "windward")) %>%
        mutate(topo_vars = "all")
mcmc3 <- c("Quercus garryana", "Pinus contorta", "Ulmus rubra") %>%
        map_df(minimal_model_mcmc,
               topo_vars = c("northness", "eastness", "tpi")) %>%
        mutate(topo_vars = "no wind")

# all species (MAP)
plan(multisession, workers = 8)
spp <- readRDS("data/derived/binned/param_9.rds")$md %>% 
        pull(species) %>% unique()
modes4 <- future_map_dfr(spp, 
                         topo_effects,
                         topo_vars = c("northness", "eastness", "tpi", "windward")) %>%
        mutate(topo_vars = "all",
               species = "all")
modes3 <- future_map_dfr(spp, 
                         topo_effects,
                         topo_vars = c("northness", "eastness", "tpi")) %>%
        mutate(topo_vars = "no wind",
               species = "all")



p <- bind_rows(modes3, modes4) %>%
        rename(value = topo_effect) %>%
        bind_rows(mcmc3, mcmc4) %>%
        set_order() %>%
        mutate(weight = ifelse(is.na(npres), 1, npres),
               topo_vars = ifelse(topo_vars == "all", "all four", "no wind. exp.")) %>%
        ggplot(aes(y = value, x = species, fill = species, color = species, weight = weight,
                   alpha = topo_vars, group = paste(species, topo_vars))) +
        facet_grid(topo_var~clim_var, labeller = label_vars) +
        geom_hline(yintercept = 0, linewidth = .5) +
        geom_boxplot(outlier.shape = NA, linewidth = .25) +
        coord_cartesian(ylim = c(-.16, .18)) +
        scale_fill_manual(values = c("black", "darkred", "darkblue", "darkgreen", "orange")) +
        scale_color_manual(values = c("black", "darkred", "darkblue", "darkgreen", "orange")) +
        scale_alpha_manual(values = c(.7, .1)) +
        guides(alpha = guide_legend(override.aes = list(fill = "black"))) +
        theme_bw() +
        theme(strip.background = element_rect(fill = "black"),
              strip.text = element_text(color = "white"),
              axis.text.x = element_text(hjust = 1, angle = 45, color = c("black", "darkred", "darkblue", "darkgreen"))) +
        labs(y = "standardized effect of topographic variable on bioclimate variable",
             x = NULL,
             alpha = "topographic\nvariables\nin model")
ggsave("figures/manuscript/minimal_species_posterior.pdf", 
       p, width = 6, height = 7, units = "in")



weighted.quantile <- function(x, w, q){
        tibble(x = x, w = w) %>%
                arrange(x) %>%
                mutate(p = cumsum(w) / sum(w),
                       d = abs(p - q)) %>%
                filter(d == min(d)) %>% 
                pull(x)
}

p <- modes4 %>%
        set_order() %>%
        mutate(ppt = case_when(bio12 < weighted.quantile(bio12, npres, .333) ~ "low",
                               bio12 > weighted.quantile(bio12, npres, .666) ~ "high",
                               TRUE ~ "medium"),
               ppt = factor(ppt, levels = c("low", "medium", "high") %>% rev())) %>%
        ggplot(aes(bio1, topo_effect, weight = npres, color = ppt, fill = ppt)) +
        facet_grid(topo_var~clim_var, labeller = label_vars, scales = "free") +
        geom_hline(yintercept = 0, linewidth = .5) +
        geom_point(size = .25) +
        geom_smooth(method = lm, alpha = .3, linewidth = .5) +
        coord_cartesian(ylim = c(-.3, .3)) +
        scale_color_manual(values = c("orangered", "gold4", "turquoise4") %>% rev()) +
        scale_fill_manual(values = c("orangered", "gold4", "turquoise4") %>% rev()) +
        theme_bw() +
        theme(strip.background = element_rect(fill = "black"),
              strip.text = element_text(color = "white")) +
        labs(y = "standardized effect of topographic variable on bioclimate variable",
             x = "species average mean annual temperature",
             color = "species\naverage\nprecipitation",
             fill = "species\naverage\nprecipitation")
ggsave("figures/manuscript/minimal_species_mods.pdf", 
       p, width = 6, height = 7, units = "in")




