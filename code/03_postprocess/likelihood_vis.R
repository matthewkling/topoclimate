
library(tidyverse)

evals <- bind_rows(read_csv("data/derived/likelihoods_minimal_model.csv") %>% 
                           mutate(model = "minimal single-species model"),
                   read_csv("data/derived/likelihoods_full_model.csv") %>% 
                           mutate(model = "full multi-species model"))

delta <- evals %>%
        group_by(species, model, domain, predictors) %>%
        summarize(mean_loglike = weighted.mean(loglike, n),
                  loglike = sum(loglike),
                  n = sum(n),
                  npres = sum(npres)) %>%
        group_by(model, species, domain, n, npres) %>%
        arrange(predictors) %>%
        summarize(delta = diff(loglike),
                  delta_mean = diff(mean_loglike),
                  .groups = "drop") %>%
        rename(species_sum = delta, species_mean = delta_mean) %>%
        gather(summary, delta, species_sum, species_mean) %>%
        mutate(summary = str_replace(summary, "_", " "),
               preferred = case_when(delta > 0 ~ "bioclimate",
                                     delta < 0 ~ "macroclimate",
                                     TRUE ~ "neither"),
               nspecies = 1,
               nplots = n,
               npresences = npres,
               model = str_remove(model, " model")) %>%
        gather(weighting, w, npresences, nspecies) %>%
        mutate(weighting = str_replace(weighting, "n", "# ")) %>%
        rename(prediction = domain,
               count = weighting)

# summary stats
stats <- delta %>%
        group_by(model, prediction, summary, count) %>%
        summarize(prop_bio = weighted.mean(delta > 0, w),
                  wilcox_p = wilcox.test(delta * w)$p.value,
                  xmax = max(delta)) %>%
        group_by(model, summary) %>%
        mutate(xmax = max(xmax))

p <- ggplot() +
        facet_grid(prediction + count ~ model + summary, scales = "free",
                   labeller = label_both) +
        geom_histogram(data = delta,
                       aes(delta, fill = preferred, weight = w),
                       boundary = 0, position = "stack", bins = 50) +
        geom_text(data = stats, 
                  aes(xmax, 0, label = paste0(round(prop_bio * 100, 1), "% bioclim;  \np = ",
                                              signif(wilcox_p, 2), "  ")), 
                  hjust = 1, vjust = -1, lineheight = .8, size = 3, color = "dodgerblue") +
        scale_fill_manual(values = c("dodgerblue", "darkred")) +
        theme_bw() +
        theme(strip.text = element_text(color = "white"),
              strip.background = element_rect(fill = "black"),
              legend.position = "top") +
        labs(x = "delta log likelihood (bioclimate minus macroclimate)",
             y = "count (# species or # plots)",
             fill = "higher-performing predictor set  ")
ggsave("figures/manuscript/validation_histograms.pdf",
       p, width = 9, height = 9, units = "in")



# analysis 2: comparison to macro + topo #########################

eval <- read_csv("data/derived/likelihoods_macro_topo.csv")

pd <- eval %>%
        mutate(model = case_when(incl_bio == 0 & incl_topo == 0 ~ "macroclim",
                                 incl_bio == 0 & incl_topo == 1 ~ "macroclim_topo",
                                 incl_bio == 1 & incl_topo == 0 ~ "bioclim",
                                 incl_bio == 1 & incl_topo == 1 ~ "bioclim_topo")) %>%
        group_by(model, domain, species) %>%
        summarize(AUC = auc(n, npres, pred),
                  loglike = weighted.mean(loglike, n),
                  npres = sum(npres)) %>%
        gather(metric, value, AUC, loglike) %>%
        group_by(species, domain, metric) %>%
        mutate(bioclim_topo = value[model == "bioclim_topo"] - value,
               bioclim = value[model == "bioclim"] - value,
               macroclim = value[model == "macroclim"] - value,
               macroclim_topo = value[model == "macroclim_topo"] - value) %>%
        gather(dmodel, dvalue, bioclim_topo:macroclim_topo) %>% 
        filter(domain == "out-of-sample") %>%
        mutate(model = factor(model, levels = c("bioclim_topo", "bioclim", "macroclim_topo", "macroclim"),
                              labels = c("bioclim + topo", "bioclim", "macroclim + topo", "macroclim")),
               dmodel = factor(dmodel, levels = c("bioclim_topo", "bioclim", "macroclim_topo", "macroclim"),
                               labels = c("bioclim + topo", "bioclim", "macroclim + topo", "macroclim"))) %>%
        group_by(metric, dmodel, model) %>%
        mutate(sign = mean(dvalue) > 0,
               p = paste0(round(mean(dvalue > 0), 2) * 100, "%"),
               p = ifelse(model == dmodel, NA, p),
               p = ifelse(species != species[1], NA, p),
               y = quantile(dvalue, .75)) %>%
        filter(dmodel == "bioclim",
               model == "macroclim + topo")

stats <- pd %>%
        mutate(delta = dvalue) %>%
        group_by(metric) %>%
        summarize(prop_bio = mean(delta > 0),
                  mean = mean(delta),
                  wilcox_p = wilcox.test(delta)$p.value,
                  xmax = max(delta))

pd <- pd %>% mutate(sign = factor(sign(dvalue), levels = c(-1, 1), labels = c("macroclim + topo", "BISHOP"))) 

p <- ggplot() +
        facet_wrap(~metric, scales = "free_x") +
        geom_histogram(data = pd, aes(dvalue, fill = sign, group = sign), boundary = 0) +
        geom_text(data = stats, 
                  aes(xmax, 0, label = paste0("mean = ", signif(mean, 2), ";  \n",
                                              round(prop_bio * 100, 1), "% positive;  \n",
                                              "p = ", signif(wilcox_p, 2), "  ")), 
                  hjust = 1, vjust = -1, lineheight = .8, size = 3, color = "dodgerblue") +
        scale_fill_manual(values = c("dodgerblue", "darkred")[2:1]) +
        theme_bw() +
        theme(strip.text = element_text(color = "white"),
              strip.background = element_rect(fill = "black"),
              legend.position = "top") +
        labs(x = "difference in out-of-sample AUC or pointwise log likelihood (BISHOP minus standard SDM)",
             y = "number of species",
             fill = "higher-performing model  ")

ggsave("figures/manuscript/validation_histograms_topo.pdf",
       p, width = 9, height = 5, units = "in")



