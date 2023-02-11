
library(tidyverse)
library(posterior)

# fit <- readRDS("data/derived/stan/pbs_d2_k8d3_le5_bulk_9/hmc/fit.rds")
fit <- readRDS("data/derived/stan/pbs_d2_k8d3_le5_9/hmc/fit.rds")

sd <- fit$sampler_diagnostics()
cd <- fit$cmdstan_diagnose()
s <- fit$summary()

quantile(s$rhat[grepl("delta", s$variable)], c(.25, .50, .75))
quantile(na.omit(s$rhat[!grepl("lp|delta", s$variable)]), c(.25, .50, .75))

quantile(s$ess_bulk[grepl("delta", s$variable)], c(.25, .50, .75))
quantile(na.omit(s$ess_bulk[!grepl("lp|delta", s$variable)]), c(.25, .50, .75))

quantile(s$ess_tail[grepl("delta", s$variable)], c(.25, .50, .75))
quantile(na.omit(s$ess_tail[!grepl("lp|delta", s$variable)]), c(.25, .50, .75))

write_csv(s, "data/derived/stan/pbs_d2_k8d3_le5_bulk_9/hmc/fit_summary.csv")
