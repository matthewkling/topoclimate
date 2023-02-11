
library(cmdstanr)
library(tidyverse)

source("code/02_model/model_functions.R")

# metadata on pre-processed model inputs
metadata <- readRDS("data/derived/binned/metadata.rds") %>%
        filter(file.exists(file_path)) %>%
        rename(data_file = file_path) %>%
        mutate(mod_vars = map(avars, function(x) c("bio12", x)),
               out_dir = paste0("data/derived/stan/", 
                                str_remove(basename(data_file), "\\.rds"))) %>%
        slice(c(3, 9))

# model fitting
n_proc <- 1
metadata %>%
        select(data_file, out_dir, vars, topo_vars, mod_vars) %>%
        mutate(pml_redo = F, hmc_redo = F,
               mod_file = "code/02_model/model.stan",
               out_dir = str_replace(out_dir, "param", "pbs_d2_k8d3_le5_bulk"),
               s_knots = 8, s_degree = 3, lambda = 5,
               hmc_warmups = 500, hmc_samples = 1000, hmc_chains = 5) %>%
        pmap(fit_pspline_model) 

