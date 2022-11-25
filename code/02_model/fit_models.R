
library(cmdstanr)
library(tidyverse)


# metadata on pre-processed model inputs
metadata <- readRDS("data/derived/binned/metadata.rds") %>%
        filter(file.exists(file_path)) %>%
        rename(data_file = file_path) %>%
        mutate(mod_vars = map(avars, function(x) c("bio12", x)),
               out_dir = paste0("data/derived/stan/", 
                                str_remove(basename(data_file), "\\.rds")))



fit_stan_model <- function(data_file = metadata[1,]$data_file,
                           mod_file = "code/02_model/model.stan",
                           out_dir = metadata[1,]$out_dir,
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
        
        
        ## ML optim ############
        
        if((!file.exists(pml_file)) | pml_redo){
                
                file.remove(list.files(pml_dir, full.names = T))
                
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
                # ~ 40 minutes
                
                fit$save_output_files(dir = pml_dir)
                saveRDS(fit, pml_file)
        }
        
        
        ## HMC sampling ###############
        
        if((!file.exists(hmc_file)) | hmc_redo){
                
                file.remove(list.files(hmc_dir, full.names = T))
                
                # use MLE as initial values
                fit <- readRDS(pml_file)
                init <- function() list(tau = matrix(fit$mle("tau") , dat$S, dat$K),
                                        alpha = fit$mle("alpha") %>% setNames(NULL),
                                        delta = matrix(fit$mle("delta"), dat$D, dat$K),
                                        Lomega = array(fit$mle("Lomega") , c(dat$S, dat$K, dat$K)))
                
                # run MCMC
                fit <- model$sample(
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
        
        return(Sys.time() - start)
}



# run models as parallel background jobs

bg_job <- function(job, wd, model_function = fit_stan_model){
        job <<- job # export to global env so it's available to job.R
        model_function <<- model_function
        rstudioapi::jobRunScript("code/02_model/job.R",
                                 name = paste(basename(job$out_dir), collapse = ", "),
                                 workingDir = wd,
                                 importEnv = T)
}

# n_proc <- 3 # number of simultaneous processes
# metadata %>%
#         #slice(c(1)) %>%
#         select(data_file, out_dir, vars, topo_vars, mod_vars) %>%
#         mutate(pml_redo = T, hmc_redo = T) %>%
#         split(rep(1:n_proc, nrow(.))[1:nrow(.)]) %>%
#         map(bg_job, wd = getwd())


# metadata %>% 
#         slice(c(1)) %>%
#         select(data_file, out_dir, vars, topo_vars, mod_vars) %>%
#         mutate(pml_redo = T, hmc_redo = T, 
#                mod_file = "code/02_model/mu_constrain.stan",
#                out_dir = str_replace(out_dir, "param", "mu_con_par")) %>%
#         split(rep(1:n_proc, nrow(.))[1:nrow(.)]) %>%
#         map(bg_job, wd = getwd())






# build basis splines
make_splines <- function(x, y, xbounds = NULL, ybounds = NULL,
                         knots = 2, degree = 3){
        k = seq(0, 1, length.out = knots)
        k = k[2:(length(k)-1)]
        if(knots == 2) k <- NULL
        
        library(splines)
        if(is.null(xbounds)) xbounds <- range(x)
        if(is.null(ybounds)) xbounds <- range(y)
        B1 <- bs(x, knots = k, Boundary.knots = xbounds, degree=degree, intercept = TRUE)
        B2 <- bs(y, knots = k, Boundary.knots = xbounds, degree=degree, intercept = TRUE)
        B <- rep(1, nrow(B1))
        for(i in 1:ncol(B1)) for(j in 1:ncol(B2)) B <- cbind(B, B1[,i] * B2[,j])
        B <- B[,2:ncol(B)]
        B
}

fit_spline_model <- function(data_file = metadata[1,]$data_file,
                           mod_file = "code/02_model/model.stan",
                           out_dir = metadata[1,]$out_dir,
                           pml_redo = F, hmc_redo = F, # whether to refit model if outfile already exists
                           vars = c("bio5", "bio6", "bio12"),
                           topo_vars = c("northness", "eastness", "tpi"),
                           mod_vars = c("bio1", "bio12"),
                           s_knots = 2,
                           s_degree = 2,
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
        
        
        s <- make_splines(ecdf(md[[mod_vars[1]]])(md[[mod_vars[1]]]),
                          ecdf(md[[mod_vars[2]]])(md[[mod_vars[2]]]),
                          knots = s_knots, degree = s_degree)
        
        # x <- data.frame(x = ecdf(md[[mod_vars[1]]])(md[[mod_vars[1]]]),
        #                 y = ecdf(md[[mod_vars[2]]])(md[[mod_vars[2]]]),
        #                 s = s[,2]) %>%
        #         filter((1:nrow(.)) %% 100 == 0)
        # ggplot(x, aes(x, y, color = s)) + geom_point()
        
        # topo matrix
        # z <- cbind(rep(1, nrow(md)), md %>% select(all_of(mod_vars)) %>% as.matrix())
        z <- s
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
        
        
        ## ML optim ############
        
        if((!file.exists(pml_file)) | pml_redo){
                
                file.remove(list.files(pml_dir, full.names = T))
                
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
                # ~ 40 minutes
                
                fit$save_output_files(dir = pml_dir)
                saveRDS(fit, pml_file)
        }
        
        
        ## HMC sampling ###############
        
        if((!file.exists(hmc_file)) | hmc_redo){
                
                file.remove(list.files(hmc_dir, full.names = T))
                
                # use MLE as initial values
                fit <- readRDS(pml_file)
                init <- function() list(tau = matrix(fit$mle("tau") , dat$S, dat$K),
                                        alpha = fit$mle("alpha") %>% setNames(NULL),
                                        delta = matrix(fit$mle("delta"), dat$D, dat$K),
                                        Lomega = array(fit$mle("Lomega") , c(dat$S, dat$K, dat$K)))
                
                # run MCMC
                fit <- model$sample(
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
        
        return(Sys.time() - start)
}



# n_proc <- 3 # number of simultaneous processes
# mm <- metadata %>%
#         slice(c(3, 7, 9)) %>%
#         select(data_file, out_dir, vars, topo_vars, mod_vars) %>%
#         mutate(pml_redo = T, hmc_redo = T,
#                out_dir = str_replace(out_dir, "param", "splines_k3d3"),
#                s_knots = 3, s_degree = 3) %>%
#         split(rep(1:n_proc, nrow(.))[1:nrow(.)]) %>%
#         map(bg_job, wd = getwd(), model_function = fit_spline_model)




# indices of adjacent bases
tensor_adj <- function(s, d = 1){
        b <- matrix(1:ncol(s), sqrt(ncol(s)))
        f <- function(x, d){
                n <- length(x)
                y <- lapply(1:(d+1), function(i) x[i:(n+i-1)])
                y <- do.call("paste", y)
                y[!grepl("NA", y)]
        }
        bx <- as.vector(apply(b, 1, f, d = d))
        by <- as.vector(apply(b, 2, f, d = d))
        tibble(x = c(bx, by)) %>% 
                separate(x, paste0("i", 1:(d+1)), convert = T) %>% 
                as.matrix()
}

fit_pspline_model <- function(data_file = metadata[1,]$data_file,
                             mod_file = "code/02_model/p_spline.stan",
                             out_dir = metadata[1,]$out_dir,
                             pml_redo = F, hmc_redo = F, # whether to refit model if outfile already exists
                             vars = c("bio5", "bio6", "bio12"),
                             topo_vars = c("northness", "eastness", "tpi"),
                             mod_vars = c("bio1", "bio12"),
                             s_knots = 2,
                             s_degree = 2,
                             lambda = 5,
                             pml_iters = 50000,
                             hmc_warmups = 200, 
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
        
        # z <- cbind(rep(1, nrow(md)), md %>% select(all_of(mod_vars)) %>% as.matrix())
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
        
        
        ## ML optim ############
        
        if((!file.exists(pml_file)) | pml_redo){
                
                file.remove(list.files(pml_dir, full.names = T))
                
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
                # ~ 40 minutes
                
                fit$save_output_files(dir = pml_dir)
                saveRDS(fit, pml_file)
        }
        
        ## HMC sampling ###############
        
        if((!file.exists(hmc_file)) | hmc_redo){
                
                file.remove(list.files(hmc_dir, full.names = T))
                
                # use MLE as initial values
                fit <- readRDS(pml_file)
                init <- function() list(tau = matrix(fit$mle("tau") , dat$S, dat$K),
                                        alpha = fit$mle("alpha") %>% setNames(NULL),
                                        delta = matrix(fit$mle("delta"), dat$D, dat$K),
                                        Lomega = array(fit$mle("Lomega") , c(dat$S, dat$K, dat$K)))
                
                # run MCMC
                fit <- model$sample(
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
        
        return(Sys.time() - start)
}

# n_proc <- 3
# metadata %>%
#         #slice(3) %>%
#         select(data_file, out_dir, vars, topo_vars, mod_vars) %>%
#         mutate(pml_redo = T, hmc_redo = T,
#                out_dir = str_replace(out_dir, "param", "psplines_k8d3_diffbeta"),
#                s_knots = 8, s_degree = 3, lambda = 1000) %>%
#         # pmap(fit_pspline_model) %>%
#         split(rep(1:n_proc, nrow(.))[1:nrow(.)]) %>%
#         map(bg_job, wd = getwd(), model_function = fit_pspline_model)


# ## version used on aws for recent ms draft
# n_proc <- 1
# metadata %>%
#         slice(c(9, 11)) %>%
#         select(data_file, out_dir, vars, topo_vars, mod_vars) %>%
#         mutate(pml_redo = F, hmc_redo = T, 
#                mod_file = "code/02_model/p_spline_d2.stan",
#                out_dir = str_replace(out_dir, "param", "pbs_d2_k8d3_le5"),
#                s_knots = 8, s_degree = 3, lambda = 5) %>%
#         # pmap(fit_pspline_model) %>%
#         split(rep(1:n_proc, nrow(.))[1:nrow(.)]) %>%
#         map(bg_job, wd = getwd(), model_function = fit_pspline_model)


