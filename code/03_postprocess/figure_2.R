
library(tidyverse)
library(patchwork)

pal <- c("blue", "dodgerblue", "gray60", "orange", "red")        

d <- expand_grid(macroclimate = seq(-3, 4, .1),
                 n = c(-.66, -.33, 0, .33, .66)) %>%
        mutate(microclimate = macroclimate + n * (1.5 + .3 * macroclimate),
               occurrence = dnorm(microclimate))

p1 <- ggplot(d, aes(macroclimate, microclimate, color = factor(n))) +
        geom_line() +
        # annotate(geom = "text", x = -2.5, hjust = 0, size = 3,
        #          y = c(2.5, 3, 3.5),
        #          color = c("blue", "gray50", "red"), 
        #          label = c("steep north-facing", "level", "steep south-facing")) +
        coord_fixed() +
        scale_x_continuous(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0)) +
        coord_cartesian(ylim = range(d$macroclimate)) +
        scale_color_manual(values = pal) +
        theme_classic() +
        theme(legend.position = "none")

p2 <- ggplot(d %>% filter(n == 0), 
             aes(microclimate, ymin = 0, ymax = occurrence)) +
        geom_ribbon(alpha = .25, fill = "gray", color = "gray") +
        scale_x_continuous(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0), limits = c(0, max(d$occurrence) * 1.1)) +
        scale_fill_manual(values = pal) +
        scale_color_manual(values = pal) +
        theme_classic() +
        theme(legend.position = "none") +
        labs(y = "occurrence")

p3 <- ggplot(d, aes(macroclimate, ymin = 0, ymax = occurrence, fill = factor(n), color = factor(n))) +
        geom_ribbon(alpha = .25) +
        scale_x_continuous(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0), limits = c(0, max(d$occurrence) * 1.1)) +
        scale_fill_manual(values = pal) +
        scale_color_manual(values = pal) +
        theme_classic() +
        theme(legend.position = "none") +
        labs(y = "occurrence")

p <- p1 + p2 + p3 &
        theme_classic() +
        theme(legend.position = "none",
              axis.text = element_blank(),
              axis.ticks = element_blank())
p

ggsave("figures/manuscript/fig2.pdf", p,
       width = 6, height = 2.1, units = "in")
