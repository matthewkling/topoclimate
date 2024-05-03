

# compare BISHOP to SDM with macroclimate and topo predictors but no bioclimate interaction

# this script also uses functions defined in crossvalidate_minimal_model.R


auc <- function(nsamp, npres, pred){
        n = as.double(sum(nsamp)) # as.double to prevent integer underflow
        n1 = as.double(sum(npres))
        n0 = as.double(n - n1)
        
        d <- tibble(nsamp = nsamp,
                    npres = npres,
                    pred = pred) %>%
                arrange(pred) %>%
                mutate(rank = cumsum(nsamp) - (nsamp/2),
                       rank0 = rank * (nsamp - npres)) %>%
                summarize(U = n0*n1 + (n0*(n0+1))/2 - sum(rank0),
                          AUC = U/(n0*n1))
        return(d$AUC)
}


weighted.median <- function(x, w){
        tibble(x = x, w = w) %>%
                arrange(x) %>%
                mutate(q = cumsum(w) / sum(w)) %>%
                filter(q == max(q[q < .5]) | q == min(q[q > .5])) %>%
                summarize(m = x[1] + (.5 - q[1])/diff(q) * diff(x)) %>%
                pull(m)
}


get_niches_topo <- function(samples, clim_vars, topo_vars){
        samples %>%
                filter(str_detect(var, "alpha|mu|sigma|beta")) %>%
                mutate(var = str_remove(var, "]"),
                       var = str_replace(var, "\\[", ",")) %>%
                separate(var, c("param", "var", "var2"), sep = ",") %>% 
                mutate(var = ifelse(param == "beta", topo_vars[as.integer(var)], clim_vars[as.integer(var)]),
                       var2 = ifelse(param == "beta", topo_vars[as.integer(var2)], clim_vars[as.integer(var2)]))
}

paraboloid <- function(x, mu, sigma){
        xc <- apply(x, 1, function(z) z - mu)
        p <- t(xc) %*% solve(sigma)
        rowSums(p * t(xc))
}


predict_suitability_topo <- function(n, x, z){
        fmu <- n %>% filter(param == "mu") %>% arrange(var) %>% pull(value)
        fsigma <- n %>% filter(param == "sigma") %>% 
                arrange(var, var2) %>% pull(value) %>% 
                matrix(nrow = length(fmu), byrow = T)
        falpha <- n %>% filter(param == "alpha") %>% pull(value)
        fbeta <- n %>% filter(param == "beta") %>% arrange(var) %>% pull(value) %>% matrix(ncol = 1)
        lp <- falpha - 
                paraboloid(x[,sort(colnames(x))], fmu, fsigma) +
                as.matrix(z[,sort(colnames(z))]) %*% fbeta
        inv_logit(lp)
}

predict_occ_topo <- function(d, z, fit, clim_vars, topo_vars){
        fit %>%
                get_samples() %>%
                get_niches_topo(clim_vars, topo_vars) %>%
                predict_suitability_topo(x = d, z = z)
}


compare_predictors_topo <- function(sp = "Juniperus osteosperma", 
                                     nfolds = 5, 
                                     fold = 1, # 0 performs fitting and evaluation on full data set, other values perform cross-validation
                                     apply_delta = T,
                                     apply_beta = T,
                                     apply_epsilon = T,
                                     mod_file = "code/02_model/minimal_model_topo.stan",
                                     e_scale = .1
){ 
        
        ## setup ##
        temp_dir <- "data/temp"
        data_file <- "data/derived/binned/param_9.rds"
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
        test <- !train
        
        # to get likelihood of training dataset without cross-validation (i.e. to compute vanilla delta AIC)
        if(fold == 0) test <- train <- rep(T, nrow(m)) 
        
        model <- cmdstan_model(mod_file)
        dat <- list(
                K = ncol(m),
                N = nrow(m[train, ]),
                D = ncol(z),
                z = z[train, ],
                m = m[train, ],
                nn = md$n[train],
                y = md$npres[train],
                apply_beta = as.integer(apply_beta),
                apply_delta = as.integer(apply_delta),
                apply_epsilon = as.integer(apply_epsilon),
                e_scale = e_scale
        )
        
        init <- function() list(mu = rep(0, dat$K),
                                tau = rep(1, dat$K),
                                alpha = 0,
                                beta = rep(0, dat$D),
                                delta = matrix(0, dat$D, dat$K),
                                Lomega = matrix(0, dat$K, dat$K),
                                # epsilon = matrix(0, dat$N, dat$K),
                                # e_scale = .1,
                                sigma_m = .1,
                                m_true = dat$m
        )
        
        fit <- function(data, model){
                model$optimize(data = data, 
                               init = list(init()),
                               output_dir = temp_dir,
                               algorithm = "lbfgs", tol_obj = 1e-15, iter = 10000,
                               refresh = 0, show_messages = F, show_exceptions = F)
        }
        
        pred <- function(fit, td){
                samples <- get_samples(fit)
                test_data <- td %>% predict_microclimate(samples, topo_vars, clim_vars)
                colnames(test_data) <- paste0("micro_", colnames(test_data))
                micro <- td %>% predict_microclimate(samples, topo_vars, clim_vars)
                topo <- td %>% select(all_of(topo_vars))
                predict_occ_topo(micro, topo, fit, clim_vars, topo_vars) %>% as.vector()
        }
        
        # fit and predict 
        td <- md %>% filter(test)
        td$pred <- dat %>% fit(model) %>% pred(td)
        
        # pointwise likelihoods
        ll <- td %>%
                mutate(domain = ifelse(fold == 0, "within-sample", "out-of-sample"),
                       model = mod_file,
                       incl_bio = apply_delta,
                       incl_topo = apply_beta,
                       incl_err = apply_epsilon,
                       e_scale = e_scale,
                       loglike = dbinom(npres, n, pred, log = T)) %>%
                select(species, id, pred, n, npres, domain:loglike)
        return(ll)
}


#### production run #####

plan(multisession, workers = 4)
eval <- expand_grid(sp = spp,
                    nfolds = nfolds, fold = 1:nfolds,
                    apply_beta = T:F,
                    apply_delta = T:F,
                    apply_epsilon = F,
                    mod_file = "code/02_model/minimal_model_topo.stan",
                    e_scale = 1) %>%
        future_pmap_dfr(possibly(compare_predictors_topo),
                        .options = furrr_options(seed = 123))
write_csv(eval, "data/derived/likelihoods_macro_topo.csv")















