
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
        labs(#x = "\u0394 log likelihood (bioclimate minus macroclimate)",
             x = "delta log likelihood (bioclimate minus macroclimate)",
             y = "count (# species or # plots)",
             fill = "higher-performing predictor set  ")
ggsave("figures/manuscript/validation_histograms.pdf",
       p, width = 9, height = 9, units = "in")
