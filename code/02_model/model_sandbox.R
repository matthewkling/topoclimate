


library(tidyverse)
library(cmdstanr)
library(tictoc)


# metadata on preprocessed model inputs
metadata <- readRDS("data/derived/binned/metadata.rds") %>%
        filter(file.exists(file_path)) %>%
        rename(data_file = file_path) %>%
        mutate(mod_vars = map(avars, function(x) c("bio12", x)),
               out_dir = paste0("data/derived/stan/", 
                                str_remove(basename(data_file), "\\.rds")),
               pml_file = "fit_pml.rds",
               stan_file = "fit_stan.rds")

meta <- slice(metadata, 4)




data_file = meta$data_file
vars = c("bio5", "bio6", "bio12")
topo_vars = c("northness", "eastness", "tpi")
mod_vars = c("bio1", "bio12")

##


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

# model and inputs
# model <- stan_model(mod_file)
dat <- list(
        K = ncol(m),
        N = nrow(md),
        S = length(unique(md$species)),
        ss = md$sp_id,
        D = ncol(z),
        z = z,
        y = md$npres,
        m = m,
        nn = md$n
)


### base model ####

model0 <- cmdstan_model("code/02_model/model_stan.stan")
tic()
fit <- model0$optimize(
        data = dat,
        init = list(list(tau = matrix(1, dat$S, dat$K),
                         alpha = rep(.5, dat$S),
                         delta = matrix(0, dat$D, dat$K),
                         Lomega = array(0, c(dat$S, dat$K, dat$K)))),
        sig_figs = 18,
        iter = 500
)
toc()
# 724 sec



# ### threaded #####
# 
# model <- cmdstan_model("code/02_model/reduce_model.stan", 
#                        cpp_options = list(stan_threads = TRUE))
# tic()
# fit <- model$optimize(
#         data = dat,
#         init = list(list(tau = matrix(1, dat$S, dat$K),
#                          alpha = rep(.5, dat$S),
#                          delta = matrix(0, dat$D, dat$K),
#                          Lomega = array(0, c(dat$S, dat$K, dat$K)))),
#         sig_figs = 18,
#         threads = 7,
#         iter = 500
# )
# toc()
# # 1306 s




### threaded, by species #####

modelS <- cmdstan_model("code/02_model/reduce_sp_model.stan", 
                        cpp_options = list(stan_threads = TRUE))
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

tic()
fit <- modelS$optimize(
        data = dat,
        init = list(list(tau = matrix(1, dat$S, dat$K),
                         alpha = rep(.5, dat$S),
                         delta = matrix(0, dat$D, dat$K),
                         Lomega = array(0, c(dat$S, dat$K, dat$K)))),
        sig_figs = 18,
        threads = 7,
        iter = 500
)
toc()
# 1206 s











### dot product model, not threaded ####

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
model0 <- cmdstan_model("code/02_model/dotprod_model.stan")
tic()
fit <- model0$optimize(
        data = dat,
        init = list(list(tau = matrix(1, dat$S, dat$K),
                         alpha = rep(.5, dat$S),
                         delta = matrix(0, dat$D, dat$K),
                         Lomega = array(0, c(dat$S, dat$K, dat$K)))),
        sig_figs = 18,
        iter = 500
)
toc()
# 232 sec



### dot product model, threaded ####

model0 <- cmdstan_model("code/02_model/dotprod_threaded_model.stan", 
                        cpp_options = list(stan_threads = TRUE))
tic()
fit <- model0$optimize(
        data = dat,
        init = list(list(tau = matrix(1, dat$S, dat$K),
                         alpha = rep(.5, dat$S),
                         delta = matrix(0, dat$D, dat$K),
                         Lomega = array(0, c(dat$S, dat$K, dat$K)))),
        sig_figs = 18,
        threads = 3,
        iter = 500
)
toc()
#  sec







#######################################################

fit_mod <- function(data_file, mod_file, out_dir, out_file = NULL, 
                    vars = c("bio5", "bio6", "bio12"),
                    topo_vars = c("northness", "eastness", "tpi", "windward"),
                    mod_vars = c("bio1", "bio12", "wspeed"),
                    engine = "stan", # "stan" or "pml"
                    iter = 2000,
                    return = "path", # "path" or "object"
                    
                    # used only for stan
                    optim_file = NULL, # point estimates will be used as inits
                    warmup = 500, nchains = 5, treemax = 10, adelta = .8){
        
        # outdir <- paste0(dirname(outfile), "/", str_remove(basename(outfile), "\\.rds"))
        # outfile <- paste0(outdir, "/fit.rds")
        if(!dir.exists(out_dir)) dir.create(out_dir)
        
        # input data
        if(class(data_file) == "character") d <- readRDS(data_file)
        if(class(data_file) == "list") d <- data_file
        
        md <- d$md# %>% filter(species %in% unique(species)[1:5])
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
        
        # model and inputs
        # model <- stan_model(mod_file)
        model <- cmdstan_model(mod_file)
        dat <- list(
                K = ncol(m),
                N = nrow(md),
                S = length(unique(md$species)),
                ss = md$sp_id,
                D = ncol(z),
                z = z,
                y = md$npres,
                m = m,
                nn = md$n
        )
        
        if(engine == "pml"){
                fit <- model$optimize(
                        data = dat,
                        init = list(list(tau = matrix(1, dat$S, dat$K),
                                         alpha = rep(.5, dat$S),
                                         delta = matrix(0, dat$D, dat$K),
                                         Lomega = array(0, c(dat$S, dat$K, dat$K)))),
                        sig_figs = 18,
                        iter = iter
                )
        }
        
        if(engine == "stan"){
                
                #### stan #####
                ofit <- readRDS(paste0(out_dir, "/", optim_file))
                init <- function() list(tau = matrix(ofit$mle("tau") , dat$S, dat$K),
                                        alpha = ofit$mle("alpha") %>% setNames(NULL),
                                        delta = matrix(ofit$mle("delta"), dat$D, dat$K),
                                        Lomega = array(ofit$mle("Lomega") , c(dat$S, dat$K, dat$K)))
                browser()
                # working
                fit <- model$sample(
                        data = dat, 
                        init = init,
                        refresh = 1,
                        sig_figs = 18,
                        iter_warmup = warmup, iter_sampling = iter, 
                        chains = 1, #nchains, parallel_chains = nchains, 
                        max_treedepth = treemax, adapt_delta = adelta
                )
                fit$profiles()
                # x <- fit$draws(variables = "delta",
                #                format = "df")
                
                
                
                #### variational #####
                init <- function() list(tau = matrix(1, dat$S, dat$K),
                                        alpha = rep(.1, dat$S),
                                        delta = matrix(0, dat$D, dat$K),
                                        Lomega = array(0, c(dat$S, dat$K, dat$K)))
                vfit <- model$variational(
                        data = dat, 
                        init = init,
                        sig_figs = 18
                )
                
                
        }
        
        
        # # save result
        fit$save_output_files(dir = out_dir)
        saveRDS(fit, paste0(out_dir, "/", out_file))
        # if(is.null(outfile)) outfile <- paste0(outdir, basename(data_file), basename(mod_file))
        # saveRDS(fit, outfile)
        if(return == "path") return(out_file)
        if(return == "object") return(fit)
}



pml <- metadata %>%
        slice(4) %>%
        select(data_file, out_dir, out_file = pml_file,
               vars, topo_vars, mod_vars) %>%
        pmap(fit_mod,
             mod_file = "code/02_model/profiling_model.stan",
             engine = "pml",
             iter = 50000)
stop()

stn <- metadata %>%
        slice(4) %>%
        select(data_file, out_dir, out_file = stan_file, optim_file = pml_file,
               vars, topo_vars, mod_vars) %>%
        pmap(fit_mod,
             mod_file = "code/02_model/profiling_model.stan",
             engine = "stan",
             iter = 10, warmup = 20)
