
# test whether bioclimate outperforms macroclimate
# in predicting occurrences of withheld test species


library(cmdstanr)
library(tidyverse)
library(furrr)

source("code/02_model/model_functions.R")
source("code/utils.R")


## model fitting ####

fit_pspline_model_xval <- function(data_file = metadata[1,]$data_file,
                                   mod_file = "code/02_model/model.stan",
                                   out_dir = metadata[1,]$out_dir,
                                   vars = c("bio5", "bio6", "bio12"),
                                   topo_vars = c("northness", "eastness", "tpi"),
                                   mod_vars = c("bio1", "bio12"),
                                   s_knots = 2,
                                   s_degree = 2,
                                   lambda = 5,
                                   pml_iters = 50000,
                                   nfolds = 5, fold = 1
){
        
        require(tidyverse)
        require(cmdstanr)
        
        start <- Sys.time()
        
        if(!dir.exists(out_dir)) dir.create(out_dir)
        out_file <- paste0(out_dir, "/fit_fold", fold, "of", nfolds, ".rds")
        if(file.exists(out_file)) return("skipping")
        
        d <- readRDS(data_file)
        md <- d$md
        
        set.seed(123)
        sp <- unique(md$species)
        folds <- sample(1:nfolds, length(sp), replace = T)
        training_sp <- sp[folds != fold]
        if(fold == 0) training_sp <- sp
        md <- filter(md, species %in% training_sp)
        
        mv <- d$mv
        # spid <- d$spid
        md$sp_id <- as.integer(factor(md$species))
        # md <- left_join(md, select(spid, species))
        
        # macroclimate matrix
        m <- md %>% select(all_of(vars)) %>% as.matrix()
        
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
        
        
        m_min <- md %>% filter(npres > 0) %>% group_by(sp_id) %>% 
                summarize_at(all_of(vars), min) %>% select(-sp_id) %>% as.matrix()
        m_max <- md %>% filter(npres > 0) %>% group_by(sp_id) %>% 
                summarize_at(all_of(vars), max) %>% select(-sp_id) %>% as.matrix()
        
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
        
        model <- cmdstan_model(mod_file)
        
        
        # run ML optimize
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
        
        fit$save_object(out_file)
        fit$save_output_files(dir = out_dir)
        
        return(Sys.time() - start)
}




## bioclimate predictions ####

predict_deltas <- function(md = slice(metadata, 1), # model metadata 
                           g, # macroclimate data for which to predict deltas
                           fit){
        vars <- md$vars[[1]]
        mod_vars <- md$mod_vars[[1]]
        topo_vars <- md$topo_vars[[1]]
        dmd <- readRDS(md$data_file)$md
        
        d <- fit$draws() %>%
                as.data.frame() %>% as.tibble() %>%
                gather(param, value, -lp__)
        
        topo_mods1 <- c("mean", mod_vars)
        params <- tibble(mod = rep(topo_mods1, length(topo_vars)),
                         topo = rep(topo_vars, each = length(topo_mods1)))
        
        deltas <- d %>% 
                filter(str_detect(param, "delta")) %>%
                mutate(param = str_remove(param, "delta\\["),
                       param = str_remove(param, "]"),
                       param = str_replace_all(param, "\\,", "_")) %>%
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
                group_by(clim_var, topo_var, id) %>%
                summarize(delta = sum(delta)) %>%
                left_join(g) %>%
                ungroup()
        return(gg)
}


predict_bioclimate <- function(md = metadata,
                               nfolds = 5, fold = 1){
        
        # fitted model
        fit <- paste0(md$out_dir, "fit_fold", fold, "of", nfolds, ".rds") %>%
                readRDS()
        
        # test data to make predictions for
        d <- readRDS(md$data_file)
        dmd <- d$md %>% mutate(id = 1:nrow(.))
        set.seed(123)
        sp <- unique(dmd$species)
        folds <- sample(1:nfolds, length(sp), replace = T)
        test_sp <- sp[folds == fold]
        if(fold == 0) test_sp <- sp
        g <- filter(dmd, species %in% test_sp)
        
        
        # calcualte deltas and microclimate
        delta <- predict_deltas(md, g, fit)
        m <- delta %>%
                select(-bio1) %>%
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
        
        m %>%
                mutate(model = "full multi-species model",
                       domain = ifelse(fold == 0, "within-sample", "out-of-sample"),
                       id = g$id) %>%
                select(model, domain, species, id, 
                       n, npres,
                       bio5, bio6, bio12,
                       b5, b6, b12)
}


## occurrence predictions ####

get_samples <- function(fit, seed = NULL){
        if(!is.null(seed)) set.seed(seed)
        samples <- fit$draws() %>% 
                as.data.frame() %>% as_tibble() %>%
                mutate(iter = 1:nrow(.)) %>%
                gather(var, value, -iter)
        if(any(grepl("\\.", samples$var))){ # i.e. if data is not ML
                samples <- samples %>%
                        separate(var, c("chain", "var"), "\\.") %>%
                        mutate(chain = as.integer(chain)) %>%
                        filter(chain == sample(chain, 1), 
                               iter == sample(iter, 1))
        }
        return(samples)
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
        samples <- get_samples(fit)
        niche <- get_niches(samples, clim_vars)
        pred <- predict_suitability(niche, d)
}

compare_predictors <- function(sp, d, md = slice(metadata, 1),
                               nfolds = 5, fold = 1){
        
        temp_dir <- "data/temp"
        vars <- md$vars[[1]]
        mod_vars <- md$mod_vars[[1]]
        topo_vars <- md$topo_vars[[1]]
        clim_vars <- md$vars[[1]]
        
        if(fold == 0) dmd <- d %>% filter(species == sp, domain == "within-sample")
        if(fold > 0) dmd <- d %>% filter(species == sp, domain == "out-of-sample")
        M <- dmd %>% select(all_of(clim_vars)) %>% as.matrix()
        m <- dmd %>% select(all_of(str_remove(clim_vars, "io"))) %>% as.matrix()
        
        # set.seed(123)
        partition <- sample(1:nfolds, nrow(m), replace = T)
        train <- partition != fold
        test <- partition == fold
        if(fold == 0) test <- train <- rep(T, nrow(m)) 
        
        # training data for candidate models
        dat_macro <- list(
                K = ncol(m),
                N = nrow(m[train, ]),
                nn = dmd$n[train],
                y = dmd$npres[train]
        )
        dat_micro <- dat_macro
        dat_macro$m <- M[train, ]
        dat_micro$m <- m[train, ]
        
        # fit candidate models
        model <- cmdstan_model("code/02_model/niche_model.stan")
        init <- function() list(mu = rep(0, dat_micro$K),
                                tau = rep(1, dat_micro$K),
                                alpha = .25,
                                Lomega = matrix(0, dat_micro$K, dat_micro$K))
        macro_fit <- model$optimize(data = dat_macro, 
                                    init = list(init()),
                                    output_dir = temp_dir,
                                    refresh = 0)
        micro_fit <- model$optimize(data = dat_micro, 
                                    init = list(init()),
                                    output_dir = temp_dir,
                                    refresh = 0)
        
        # species occurrence predictions for candidate models
        test_data <- dmd %>% filter(test)
        test_data$occ_macro <- test_data %>% select(all_of(clim_vars)) %>% predict_occ(macro_fit, clim_vars)
        test_data$occ_micro <- test_data %>% select(str_remove(clim_vars, "io")) %>% predict_occ(micro_fit, clim_vars)
        
        # log likelihood
        ll <- test_data %>%
                gather(predictors, pred, occ_macro, occ_micro) %>%
                mutate(species = sp,
                       predictors = str_remove(predictors, "occ_"),
                       domain = ifelse(fold == 0, "within-sample", "out-of-sample"),
                       loglike = dbinom(npres, n, pred, log = T)) %>%
                select(species, predictors, domain, id, pred, n, npres, loglike)
        return(ll)
}


## execution ########

# metadata on pre-processed model inputs
metadata <- readRDS("data/derived/binned/metadata.rds") %>%
        filter(file.exists(file_path)) %>%
        rename(data_file = file_path) %>%
        mutate(mod_vars = map(avars, function(x) c("bio12", x)),
               out_dir = paste0("data/derived/stan/", 
                                str_remove(basename(data_file), "\\.rds"))) %>%
        slice(c(9)) %>%
        select(data_file, out_dir, vars, topo_vars, mod_vars) %>%
        mutate(mod_file = "code/02_model/model.stan",
               out_dir = "data/derived/stan_xval/full_topo4/",
               s_knots = 8, s_degree = 3, lambda = 5)

# fit bioclimate models
nfolds <- 10
metadata %>%
        expand_grid(nfolds = nfolds, fold = 0:nfolds) %>%
        pmap(fit_pspline_model_xval)

# make bioclimate predictions for holdout species
bc <- expand_grid(nfolds = nfolds, fold = 0:nfolds) %>%
        pmap_dfr(predict_bioclimate, 
                 md = slice(metadata, 1))

# evaluate bioclimate vs macroclimate predictors
nfolds_sp <- 10
plan(multisession, workers = 8)
eval <- expand_grid(sp = unique(bc$species),
                    nfolds = nfolds_sp, fold = 0:nfolds_sp) %>%
        future_pmap_dfr(possibly(compare_predictors), 
                        md = slice(metadata, 1), 
                        d = bc,
                        .options = furrr_options(seed = 123))
write_csv(eval, "data/derived/likelihoods_full_model.csv")


