
# this script tests the statistical model by simulating a dataset
# and evaluating the fitted versus true topoclimate coefficients

# 1) simulate data with known niches and topo effects
# 2) summarize using actual binning function
# 3) fit stan model
# 4) test model ability to recover known deltas

library(tidyverse)
library(furrr)
library(parallel)
library(conflicted)
library(rethinking)

conflict_prefer("map", "purrr")
conflict_prefer("filter", "dplyr")
conflict_prefer("lag", "dplyr")
conflict_prefer("extract", "rstan")
conflict_prefer("Position", "ggplot2")
conflict_prefer("rstudent", "stats")


# functions from elsewhere in the project
source("code/02_model/model_functions.r")
source("code/03_bin_functions.r")
source("code/03_postprocess/utils.R")


# 1: simulate input data ######################################################

# variables
clim_vars <- c("bio5", "bio6", "bio12")
topo_vars <- c("northness", "eastness", "tpi")
mod_vars <- c("bio1", "bio12")

set.seed(123)


## landscape topography and macroclimate ############

n <- 1000000
plots <- data.frame(bio5 = rnorm(n),
                    bio6 = rnorm(n),
                    bio12 = rnorm(n),
                    northness = rnorm(n),
                    eastness = rnorm(n),
                    tpi = rnorm(n)) %>%
        mutate(bio1 = (bio5 + bio6) / 2,
               plot_id = 1:nrow(.)) %>%
        as_tibble()


## topoclimate ###########

# random "true" coefficients for each modifier * topo * climate variable combination
deltas <- expand_grid(clim_var = clim_vars,
                      topo_var = topo_vars,
                      mod_var = c("intercept", mod_vars)) %>%
        mutate(effect = rnorm(nrow(.), 0, .1))

# calculate plot topoclimate from topography, macroclimate, and these deltas
# bio5micro = bio5 + northness_heat_effect * northness + eastness_heat_effect * eastness + tpi_heat_effect * tpi
# northness_heat_effect = delta1 + delta2 * bio1 + delta3 * bio12
topoclim <- plots %>%
        mutate(bio12x = bio12,
               intercept = 1) %>%
        gather(clim_var, clim_value, all_of(clim_vars)) %>%
        rename(bio12 = bio12x) %>%
        gather(topo_var, topo_value, all_of(topo_vars)) %>%
        gather(mod_var, mod_value, all_of(c("intercept", mod_vars))) %>%
        full_join(deltas) %>%
        
        mutate(value = mod_value * effect) %>%
        group_by(plot_id, clim_var, clim_value, topo_var, topo_value) %>%
        summarize(topo_effect = sum(value)) %>%
        
        mutate(value = topo_value * topo_effect) %>%
        group_by(plot_id, clim_var, clim_value) %>%
        summarize(delta = sum(value)) %>%
        mutate(micro_value = clim_value + delta) %>%
        select(plot_id, clim_var, micro_value) %>%
        spread(clim_var, micro_value) %>%
        ungroup()


## species niches ########

# simulate random parameters for a multidimensional Gaussian niche
make_species <- function(id, clim_vars){
        
        k <- length(clim_vars)
        
        # niche means
        mu <- rnorm(k, 0, 2) %>% setNames(clim_vars)
        
        # niche variance-covariance
        tau <- rgamma(k, shape = 2, scale = 1)
        omega <- rlkjcorr(1, k, 1.5)
        diag_tau <- matrix(0, k, k)
        diag(diag_tau) <- tau
        sigma <- diag_tau %*% omega %*% diag_tau
        rownames(sigma) <- colnames(sigma) <- clim_vars
        
        # prevalence at niche optimum
        alpha <- runif(1, .1, 1)
        
        list(id = id, mu = mu, sigma = sigma, alpha = alpha)
}

nspp <- 20
species <- letters[1:nspp] %>% map(make_species, clim_vars = clim_vars)


## species occurrences across landscape ######

# use species niche and plot microclimate to simulate occurrences across plots
populate_landscape <- function(sp, ls){
        
        # unpack niche parameters
        mu <- sp$mu
        sigma <- sp$sigma
        alpha <- sp$alpha
        invsigma <- solve(sigma)
        k <- length(mu)
        
        # microclimate data matrix
        x <- ls %>% select(all_of(names(mu))) %>% as.matrix()
        n <- nrow(x)
        
        # occurrence probability on each plot (Gaussian probability surface)
        a <- 1 / (2 * pi) ^ (k / 2)
        c <- x - matrix(mu, n, k, byrow = T)
        p <- rep(0, n)
        for(i in 1:n){
                ci <- matrix(c[i,], 1)
                p[i] <- ci %*% invsigma %*% t(ci)
        }
        p <- a * 1/sqrt(determinant(sigma, logarithm = F)$modulus) * exp(-.5 * p)
        ci <- matrix(rep(0, k), 1)
        pmu <- ci %*% invsigma %*% t(ci)
        pmu <- a * 1/sqrt(determinant(sigma, logarithm = F)$modulus) * exp(-.5 * pmu)
        p <- p / pmu[1,1] * alpha
        
        # sample random presences per occurrence probabilities
        ls$p <- p
        ls$occ <- rbernoulli(nrow(ls), ls$p)
        ls$species <- sp$id
        return(ls %>% filter(occ) %>% select(plot_id, species))
}


occ <- species %>% map_df(populate_landscape, ls = topoclim)

plots <- left_join(plots, occ)



## binned summaries ################

d <- plots %>%
        mutate(subplot_id = plot_id, subplot = 1,
               lon = runif(nrow(.)), lat = runif(nrow(.)))

b <- bin_data(spp = na.omit(unique(d$species)),
              nbins = 4, 
              vars = clim_vars,
              topo_vars = topo_vars,
              avars = "bio1",
              ncores = 5,
              file_path = NULL)



# 2: fit model ################################################################

fitfile <- "data/derived/stan/simulation.rds"
fit <- b %>% 
        fit_model(modfile = "code/02_model/model_optim.stan",
                  outfile = fitfile,
                  vars = clim_vars,
                  topo_vars = topo_vars,
                  mod_vars = mod_vars,
                  engine = "pml", 
                  iter = 10000,
                  return = "object")


# 3: compare true and fitted effects ##########################################

# d <- b
md <- b$md
mv <- b$mv
dc <- b$dc

lookup <- function(var){
        list(bio5 = "heat", 
             bio6 = "cold", 
             bio12 = "moisture",
             bio1 = "tmean")[[var]]
}

topo_mods1 <- c("int", mod_vars)
params <- tibble(mod = rep(topo_mods1, length(topo_vars)),
                 topo = rep(topo_vars, each = length(topo_mods1)))


get_deltas <- function(fitfile){
        fit <- readRDS(fitfile)
        f <- as.data.frame(fit$par) %>%
                rownames_to_column("param") %>% as_tibble()
        names(f)[2] <- "value"
        
        deltas <- f %>% 
                filter(str_detect(param, "delta")) %>%
                mutate(param = str_remove(param, "delta\\["),
                       param = str_remove(param, "\\]"),
                       param = str_replace(param, ",", "_")) %>%
                separate(param, c("param", "var")) %>%
                mutate(var = recode(var, 
                                    "1" = lookup(clim_vars[1]), 
                                    "2" = lookup(clim_vars[2]), 
                                    "3" = lookup(clim_vars[3])),
                       topo = params$topo[as.integer(param)],
                       mod = params$mod[as.integer(param)])#,
                       # mod = factor(mod, levels = params$mod[1:(length(topo_vars)+1)]))
        deltas
}

fitted_deltas <- get_deltas(fitfile) %>%
        rename(clim_var = var,
               topo_var = topo,
               mod_var = mod,
               fitted = value) %>%
        select(-param) %>%
        mutate(clim_var = recode(clim_var,
                                 "heat" = "bio5",
                                 "cold" = "bio6",
                                 "moisture" = "bio12"),
               mod_var = ifelse(mod_var == "int", "intercept", mod_var))

p <- left_join(deltas, fitted_deltas) %>%
        ggplot(aes(effect, fitted)) +
        geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
        geom_point(color = "darkred") +
        geom_smooth(method = lm, se = F, color = "darkred", size = .25) +
        theme_minimal() +
        coord_fixed() +
        scale_x_continuous(breaks = seq(-.1, .1, .1)) +
        scale_y_continuous(breaks = seq(-.1, .1, .1)) +
        labs(x = "simulated true effect",
             y = "fitted model estimate")
ggsave("figures/simulation/true_fitted_scatter.png", 
       p, width = 4, height = 4, units = "in")


