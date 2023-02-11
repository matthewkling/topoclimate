
# this script tests the statistical model by simulating a dataset
# and evaluating the fitted versus true topoclimate coefficients

# 1) simulate data with known niches and topo effects
# 3) fit stan model
# 4) test model ability to recover known params

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
source("code/01_preprocess/03_bin_functions.R")
source("code/02_model/model_functions.r")
source("code/utils.R")

set.seed(123)


# 1: simulate input data ######################################################

# variables
clim_vars <- c("bio5", "bio6", "bio12")
topo_vars <- c("northness", "eastness", "windward", "tpi")
mod_vars <- c("bio1", "bio12")

# load empirical data, to use as a structural template
d <- readRDS("data/derived/binned/param_9.rds")$md %>% select(-wspeed)


## microclimate ###########

### generate basis functions of modifier variables ####
z <- make_splines(x = ecdf(d[[mod_vars[1]]])(d[[mod_vars[1]]]),
                  y = ecdf(d[[mod_vars[2]]])(d[[mod_vars[2]]]),
                  knots = 8, degree = 3)
adj <- tensor_adj(z, d = 2)

adj <- map(1:length(topo_vars), function(i) adj + (i - 1) * ncol(z)) %>%
        do.call("rbind", .)
z <- topo_vars %>%
        map(function(v) z * matrix(d[[v]], nrow(z), ncol(z))) %>%
        do.call("cbind", .)



### simulate random "true" known coefficients ####

lambda <- 5 # smoothing parameter
set.seed(1)

# random starting values for deltas
rdelta <- matrix(rnorm(ncol(z) * length(clim_vars), 0, 0.1), 
                 ncol(z), length(clim_vars))
adj2 = adj[,3]
adj1 = adj[,2]
adj0 = adj[,1]
diffs <- (rdelta[adj0,] - (2 * rdelta[adj1,]) + rdelta[adj2,]) ^ 2

# iteratively adjust delta values to match qualitative goals
fx <- function(x){
        delta <- matrix(x, ncol(z), length(clim_vars))
        diffs <- (delta[adj0,] - (2 * delta[adj1,]) + delta[adj2,]) ^ 2
        dnorm(mean(x), 0, .01, log = T) + 
                dgamma(sd(x), 2, 100, log = T) +
                dnorm(mean(diffs), 2, 1e6, log = T) + 
                dgamma(sd(diffs), 2, 1e6, log = T)
}
opt <- optim(as.vector(rdelta), 
             fx,
             method = "BFGS",
             control = list(fnscale = -1,
                            trace = T,
                            maxit = 1000))
delta <- matrix(opt$par, ncol(z), length(clim_vars)) / 10


### calculate microclimate for every data bin ####

m <- d %>% select(all_of(clim_vars)) %>% as.matrix() # macroclimate
topoclimate <- m + z %*% delta # microclimate
topoclimate <- as.data.frame(topoclimate) %>% as_tibble()
names(topoclimate) <- paste0(names(topoclimate), "m")
d <- bind_cols(d, topoclimate)


## populate plots with species ####

populate <- function(s){
        
        message(s)
        ds <- filter(d, species == s)
        
        # approximate niche params
        mu <- ds %>%
                summarize(bio5 = weighted.mean(bio5m, npres),
                          bio6 = weighted.mean(bio6m, npres),
                          bio12 = weighted.mean(bio12m, npres)) %>%
                unlist()
        
        sigma <- ds %>% select(bio5m, bio6m, bio12m) %>% cov.wt(ds$npres^2 / ds$n)
        sigma <- sigma$cov
        
        alpha <- ds %>%
                mutate(delta = sqrt((bio5 - mu["bio5"])^2 + (bio6 - mu["bio6"])^2 + (bio12 - mu["bio12"])^2)) %>%
                filter(delta < .5) %>%
                summarize(p = sum(npres)/sum(n)) %>%
                pull(p)
        alpha <- pmin(alpha * 3, 1)
        
        invsigma <- solve(sigma)
        k <- length(mu)
        
        # microclimate data matrix
        x <- ds %>% select(bio5m, bio6m, bio12m) %>% as.matrix()
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
        ds$p <- p
        ds$occ <- rbinom(nrow(ds), ds$n, ds$p)
        return(ds)
}

b <- map_df(unique(na.omit(d$species)), populate) %>%
        rename(npres_true = npres,
               npres = occ) %>%
        na.omit()



# 2: fit model ################################################################

fit_model <- function(md,
                      mod_file = "code/02_model/model.stan",
                      out_dir = "data/derived/simulation",
                      s_knots = 8,
                      s_degree = 3,
                      pml_iters = 50000,
                      lambda = 5){
        
        require(tidyverse)
        require(cmdstanr)
        
        # file admin
        pml_dir <- paste0(out_dir, "/pml")
        pml_file <- paste0(pml_dir, "/fit.rds")
        if(!dir.exists(out_dir)) dir.create(out_dir)
        if(!dir.exists(pml_dir)) dir.create(pml_dir)
        
        # macroclimate matrix
        vars <- clim_vars
        m <- md %>% select(all_of(vars)) %>% as.matrix()
        
        # basis splines
        z <- make_splines(ecdf(md[[mod_vars[1]]])(md[[mod_vars[1]]]),
                          ecdf(md[[mod_vars[2]]])(md[[mod_vars[2]]]),
                          knots = s_knots, degree = s_degree)
        adj <- tensor_adj(z, d = 2)
        
        # topo matrix
        adj <- map(1:length(topo_vars), function(i) adj + (i - 1) * ncol(z)) %>%
                do.call("rbind", .)
        z <- topo_vars %>%
                map(function(v) z * matrix(md[[v]], nrow(z), ncol(z))) %>%
                do.call("cbind", .)
        
        # data ranges
        md$sp_id <- as.integer(factor(md$species))
        m_min <- md %>% filter(npres > 0) %>% group_by(sp_id) %>% 
                summarize_at(all_of(vars), min) %>% select(-sp_id) %>% as.matrix()
        m_max <- md %>% filter(npres > 0) %>% group_by(sp_id) %>% 
                summarize_at(all_of(vars), max) %>% select(-sp_id) %>% as.matrix()
        
        # model data
        dat <- list(
                K = ncol(m),
                N = nrow(md),
                S = length(unique(md$species)),
                ss = md$sp_id,
                D = ncol(z),
                z = z,
                
                Nadj = nrow(adj),
                adj2 = adj[,3],
                adj1 = adj[,2],
                adj0 = adj[,1],
                lambda = lambda,
                
                y = md$npres,
                m = m,
                nn = md$n,
                
                s0 = md %>% mutate(i = 1:nrow(.)) %>% group_by(species) %>%
                        summarize(start = min(i)) %>% pull(start),
                s1 = md %>% mutate(i = 1:nrow(.)) %>% group_by(species) %>%
                        summarize(end = max(i)) %>% pull(end),
                
                m_span = m_max - m_min,
                m_min = m_min 
        )
        
        # compile model
        model <- cmdstan_model(mod_file)
        
        # fit with ML
        file.remove(list.files(pml_dir, full.names = T))
        fit <- model$optimize(
                data = dat,
                init = list(list(tau = matrix(1, dat$S, dat$K),
                                 alpha = rep(.5, dat$S),
                                 delta = matrix(0, dat$D, dat$K),
                                 Lomega = array(0, c(dat$S, dat$K, dat$K)))),
                sig_figs = 18,
                iter = pml_iters,
                seed = 123
        )
        
        # save
        fit$save_output_files(dir = pml_dir)
        saveRDS(fit, pml_file)
}

fit_model(b)



# 3: compare true and fitted effects ##########################################

## true effects ####

g <- expand_grid(bio1 = seq(-2, 2, .1),
                 bio12 = seq(-2, 2, 1))
zz <- make_splines(x = ecdf(d[[mod_vars[1]]])(g[[mod_vars[1]]]),
                   y = ecdf(d[[mod_vars[2]]])(g[[mod_vars[2]]]),
                   knots = 8, degree = 3)
gzz <- zz %>%
        as.data.frame() %>% as_tibble() %>%
        setNames(paste0("b", 1:100)) %>%
        bind_cols(g, .) %>%
        gather(basis, basis_value, -bio1, -bio12)

true <- delta %>%
        as.data.frame() %>% as_tibble() %>%
        setNames(clim_vars) %>%
        mutate(basis = paste0("b", rep(1:100, length(topo_vars))),
               topo_var = rep(topo_vars, each = 100)) %>%
        gather(clim_var, delta, bio5:bio12)

gzt <- full_join(gzz, true) %>%
        group_by(bio1, bio12, topo_var, clim_var) %>%
        summarize(effect = sum(basis_value * delta))



## fitted effects ####

fit <- as_cmdstan_fit("data/derived/simulation/pml/p_spline_d2-202301061343-1-229f68.csv")

f <- as.data.frame(fit$mle()) %>%
        setNames("value") %>%
        rownames_to_column("param") %>%
        filter(grepl("delta", param)) %>%
        as_tibble() %>%
        mutate(param = str_remove(param, "delta\\["),
               param = str_remove(param, "\\]")) %>%
        mutate(clim_var = rep(clim_vars, each = nrow(delta))) %>%
        mutate(basis = rep(paste0("b", rep(1:100, length(topo_vars))), length(clim_vars)),
               topo_var = rep(rep(topo_vars, each = 100), length(clim_vars))) %>%
        select(-param) %>%
        rename(delta = value)

gzf <- full_join(gzz, f) %>%
        group_by(bio1, bio12, topo_var, clim_var) %>%
        summarize(effect = sum(basis_value * delta))

tf <- bind_rows(gzt %>% mutate(model = "simulated true effect"),
                gzf %>% mutate(model = "fitted effect")) %>%
        mutate(clim_var = factor(clim_var,
                                 levels = c("bio12", "bio5", "bio6"),
                                 labels = c("moisture", "high temperature", "low temperature")),
               topo_var = factor(topo_var,
                                 levels = c("northness", "eastness", "windward", "tpi"),
                                 labels = c("northness", "eastness", "windward exposure", "elevational position")),
               model = factor(model,
                              levels = c("simulated true effect", "fitted effect")))


## plots ####

p <- tf %>%
        ggplot(aes(bio1, effect, color = model, alpha = bio12, 
                   group = paste(model, bio12))) + 
        facet_grid(topo_var ~ clim_var) +
        geom_hline(yintercept = 0, color = "gray") +
        geom_line() +
        theme_bw() +
        theme(panel.grid = element_blank(),
              strip.background = element_rect(fill = "black"),
              strip.text = element_text(color = "white")) +
        scale_alpha_continuous(range = c(.3, 1)) +
        scale_color_manual(values = c("black", "red")) +
        scale_x_continuous(expand = c(0,0), breaks = -1:1) +
        labs(x = "macro annual temperature (std)",
             y = "standardized effect of terrain variable on climate variable",
             alpha = "macro annual\nprecipitation (std)")
ggsave("figures/manuscript/simulation_true_fitted_curves.pdf", 
       p, width = 8, height = 6, units = "in")

p <- full_join(f %>% rename(fitted = delta),
          true %>% rename(true = delta)) %>%
        ggplot(aes(true, fitted)) +
        geom_point(color = "dodgerblue", size = .5) +
        geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
        geom_smooth(method = lm, se = F, color = "red", size = .25) +
        theme_minimal() +
        coord_fixed() +
        scale_x_continuous(breaks = seq(-.2, .2, .1)) +
        scale_y_continuous(breaks = seq(-.2, .2, .1)) +
        labs(x = "simulated true effect",
             y = "fitted model estimate")
ggsave("figures/manuscript/simulation_true_fitted_scatter.pdf", 
       p, width = 4, height = 4, units = "in")


# calculate fit statistics
full_join(f %>% rename(fitted = delta),
          true %>% rename(true = delta)) %>%
        summarize(r2 = cor(fitted, true),
                  rmse = sqrt(mean((fitted - true)^2)),
                  rmse_std = rmse / sd(true),
                  nrmse = rmse / diff(range(true)))
