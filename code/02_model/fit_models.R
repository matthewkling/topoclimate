
library(cmdstanr)
library(tidyverse)
library(callr)


# metadata on pre-processed model inputs
metadata <- readRDS("data/derived/binned/metadata.rds") %>%
        filter(file.exists(file_path)) %>%
        rename(data_file = file_path) %>%
        mutate(mod_vars = map(avars, function(x) c("bio12", x)),
               out_dir = paste0("data/derived/stan/", 
                                str_remove(basename(data_file), "\\.rds")))



# meta <- slice(metadata, 4)

fit_stan_model <- function(data_file = meta$data_file,
                           mod_file = "code/02_model/model.stan",
                           out_dir = meta$out_dir,
                           pml_redo = F, hmc_redo = F, # whether to refit model if outfile already exists
                           vars = c("bio5", "bio6", "bio12"),
                           topo_vars = c("northness", "eastness", "tpi"),
                           mod_vars = c("bio1", "bio12"),
                           pml_iters = 50000,
                           hmc_warmups = 100, 
                           hmc_samples = 100, hmc_chains = 5
){
        
        require(tidyverse)
        require(cmdstanr)
        
        ## file admin ##############
        
        start <- Sys.time()
        
        pml_dir <- paste0(out_dir, "/pml")
        hmc_dir <- paste0(out_dir, "/hmc")
        pml_file <- paste0(pml_dir, "/fit.rds")
        hmc_file <- paste0(hmc_dir, "/fit.rds")
        
        if(!dir.exists(out_dir)) dir.create(out_dir)
        if(!dir.exists(pml_dir)) dir.create(pml_dir)
        if(!dir.exists(hmc_dir)) dir.create(hmc_dir)
        
        
        
        ## data formatting ###########
        
        d <- readRDS(data_file)
        
        md <- d$md %>% filter(species %in% unique(species))
        mv <- d$mv
        spid <- d$spid
        md$sp_id <- as.integer(factor(md$species))
        md <- left_join(md, select(spid, species))
        
        # macroclimate matrix
        m <- md %>% select(all_of(vars)) %>% as.matrix()
        
        # topo matrix
        z <- cbind(rep(1, nrow(md)), md %>% select(all_of(mod_vars)) %>% as.matrix())
        z <- topo_vars %>%
                map(function(v) z * matrix(md[[v]], nrow(z), ncol(z))) %>%
                do.call("cbind", .)
        
        dat <- list(
                K = ncol(m),
                N = nrow(md),
                S = length(unique(md$species)),
                ss = md$sp_id,
                D = ncol(z),
                z = z,
                y = md$npres,
                m = m,
                nn = md$n,
                
                si = 1:length(unique(md$species)),
                sp1 = md %>% mutate(i = 1:nrow(.)) %>% group_by(species) %>%
                        summarize(start = min(i)) %>% pull(start),
                spn = md %>% count(species) %>% pull(n)
        )
        
        
        ## ML optim ############
        
        if((!file.exists(pml_file)) | pml_redo){
                
                model <- cmdstan_model(mod_file)
                
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
                # ~ 40 minutes
                
                fit$save_output_files(dir = pml_dir)
                saveRDS(fit, pml_file)
        }
        
        
        ## HMC sampling ###############
        
        if((!file.exists(hmc_file)) | hmc_redo){
                
                # start sampling at MLE
                fit <- readRDS(pml_file)
                init <- function() list(tau = matrix(fit$mle("tau") , dat$S, dat$K),
                                        alpha = fit$mle("alpha") %>% setNames(NULL),
                                        delta = matrix(fit$mle("delta"), dat$D, dat$K),
                                        Lomega = array(fit$mle("Lomega") , c(dat$S, dat$K, dat$K)))
                
                
                fit <- model1$sample(
                        data = dat,
                        init = init,
                        sig_figs = 18,
                        iter_warmup = hmc_warmups, iter_sampling = hmc_samples, 
                        refresh = 10,
                        chains = hmc_chains, parallel_chains = hmc_chains
                )
                # ~ 7.5 hours
                
                fit$save_output_files(dir = hmc_dir)
                saveRDS(fit, hmc_file)
                
        }
        
        
        return(start - Sys.time())
        
}

# run models as background processes
jobs <- metadata %>%
        slice(1:4) %>%
        select(data_file, out_dir, vars, topo_vars, mod_vars) %>%
        split(1:nrow(.)) %>%
        lapply(function(x) r_bg(function(par, fun){
                require(purrr)
                pmap(par, fun)}, 
                args = list(par = x, fun = fit_stan_model)))

jobs

