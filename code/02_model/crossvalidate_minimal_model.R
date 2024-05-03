
library(cmdstanr)
library(tidyverse)
library(furrr)
source("code/utils.R")

get_samples <- function(fit, seed = NULL){
        samples <- fit$draws() %>% 
                as.data.frame() %>% as_tibble() %>%
                mutate(iter = 1:nrow(.)) %>%
                gather(var, value, -iter)
        return(samples)
}

get_deltas <- function(samples, topo_vars, clim_vars){
        samples %>%
                filter(str_detect(var, "delta")) %>%
                mutate(var = str_remove(var, "delta\\["),
                       var = str_remove(var, "]")) %>%
                separate(var, c("topo", "clim"), remove = F) %>%
                mutate(topo_var = topo_vars[as.integer(topo)],
                       clim_var = clim_vars[as.integer(clim)]) %>%
                select(topo_var, clim_var, topo_effect = value)
}

predict_microclimate <- function(d, samples, topo_vars, clim_vars){
        deltas <- get_deltas(samples, topo_vars, clim_vars)
        d %>%
                select(n, npres, all_of(clim_vars), all_of(topo_vars)) %>%
                mutate(datum = 1:nrow(.)) %>%
                gather(clim_var, macro, all_of(clim_vars)) %>%
                gather(topo_var, topo_value, all_of(topo_vars)) %>%
                left_join(deltas, by = join_by(clim_var, topo_var)) %>%
                arrange(datum, clim_var) %>%
                group_by(datum, clim_var, macro) %>%
                summarize(delta = sum(topo_value * topo_effect), .groups = "drop") %>%
                mutate(micro = macro + delta) %>%
                select(datum, clim_var, micro) %>%
                spread(clim_var, micro) %>%
                ungroup() %>%
                select(all_of(clim_vars)) %>%
                as.matrix()
}

get_niches <- function(samples, clim_vars){
        samples %>%
                filter(str_detect(var, "alpha|mu|sigma")) %>%
                mutate(var = str_remove(var, "]"),
                       var = str_replace(var, "\\[", ",")) %>%
                separate(var, c("param", "var", "var2"), sep = ",") %>% 
                mutate(var = clim_vars[as.integer(var)],
                       var2 = clim_vars[as.integer(var2)])
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

predict_occ <- function(d, fit, clim_vars){
        fit %>%
                get_samples() %>%
                get_niches(clim_vars) %>%
                predict_suitability(d)
}

compare_predictors <- function(sp = "Juniperus osteosperma", 
                               nfolds = 5, 
                               fold = 1){ # 0 performs fitting and evaluation on full data set, other values perform cross-validation
        
        ## setup ##
        message(sp, fold)
        temp_dir <- "data/temp"
        data_file <- "data/derived/binned/param_9.rds"
        mod_file <- "code/02_model/minimal_model.stan"
        clim_vars <- c("bio5", "bio6", "bio12")
        mod_vars <- c("bio12", "bio1")
        topo_vars <- c("northness", "eastness", "tpi", "windward")
        d <- readRDS(data_file)
        md <- d$md %>% filter(species == sp)
        m <- md %>% select(all_of(clim_vars)) %>% as.matrix()
        z <- md %>% select(all_of(topo_vars)) %>% as.matrix()
        
        md$id <- 1:nrow(m)
        
        set.seed(123)
        partition <- sample(1:nfolds, nrow(m), replace = T)
        train <- partition != fold
        test <- partition == fold
        
        # to get likelihood of training dataset without cross-validation (i.e. to compute vanilla delta AIC)
        if(fold == 0) test <- train <- rep(T, nrow(m)) 
        
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
        
        ##  step 2: make bioclimate predictions ##
        bioclim_samples <- get_samples(bioclim_fit)
        micro_train <- md %>%
                filter(train) %>%
                predict_microclimate(bioclim_samples, topo_vars, clim_vars)
        
        ## step 3: fit competing niche models, based on macroclimate vs bioclimate ##
        dat_micro <- dat_macro <- list(
                K = ncol(m),
                N = nrow(m[train, ]),
                m = m[train, ],
                nn = md$n[train],
                y = md$npres[train]
        )
        dat_micro$m <- micro_train
        
        mod_file <- "code/02_model/niche_model.stan"
        model <- cmdstan_model(mod_file)
        macro_fit <- model$optimize(data = dat_macro, 
                                    init = list(init()),
                                    output_dir = temp_dir,
                                    refresh = 0)
        micro_fit <- model$optimize(data = dat_micro, 
                                    init = list(init()),
                                    output_dir = temp_dir,
                                    refresh = 0)
        
        ## step 3: make predictions on independent data for competing models ##
        test_data <- md %>% filter(test)
        test_data$occ_macro <- test_data %>% select(all_of(clim_vars)) %>% 
                predict_occ(macro_fit, clim_vars)
        test_data_micro <- md %>% filter(test) %>% 
                predict_microclimate(bioclim_samples, topo_vars, clim_vars)
        colnames(test_data_micro) <- paste0("micro_", colnames(test_data_micro))
        test_data$occ_micro <- test_data_micro %>% predict_occ(micro_fit, clim_vars)
        
        ## step 4: evaluate predicted vs true occurrences ##
        ll <- test_data %>%
                gather(predictors, pred, occ_macro, occ_micro) %>%
                mutate(species = sp,
                       predictors = str_remove(predictors, "occ_"),
                       domain = ifelse(fold == 0, "within-sample", "out-of-sample"),
                       loglike = dbinom(npres, n, pred, log = T)) %>%
                select(species, predictors, domain, id, pred, n, npres, loglike)
        return(ll)
}



# fit and evaluate models
nfolds <- 10
plan(multisession, workers = 8)
spp <- readRDS("data/derived/binned/param_9.rds")$md %>% 
        pull(species) %>% unique()
eval <- expand_grid(sp = spp, nfolds = nfolds, fold = 0:nfolds) %>% 
        future_pmap_dfr(compare_predictors,
                        .options = furrr_options(seed = 123))
write_csv(eval, "data/derived/likelihoods_minimal_model.csv")


