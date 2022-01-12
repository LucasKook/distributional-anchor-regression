# vis scenario iv1
# Lucas Kook
# 05.11.2020

# Dependencies ------------------------------------------------------------

library(tidyverse)
library(ggpubr)
library(patchwork)
theme_set(theme_pubr())

bpath <- "results"
xilabs <- c("-infinity", round(log10(c(1, 3, 7, 1000, 10000)), 2))

# Read results ------------------------------------------------------------

cfx <- read.csv(file.path(bpath, "scenario-iv1-cfx.csv")) %>% 
  gather(key = "model.xi", value = "coef") %>% 
  separate(model.xi, into = c("model", "xi"), sep = "\\.") %>% 
  mutate(xi = ordered(xi, levels = c(0, 1, 3, 7, 1000, 10000),
                      labels = xilabs))

nll <- read.csv(file.path(bpath, "scenario-iv1-nll.csv")) %>%
  gather(key = "model.xi", value = "logLik") %>% 
  separate(model.xi, into = c("model", "xi"), sep = "\\.") %>% 
  mutate(xi = ordered(xi, levels = c(0, 1, 3, 7, 1000, 10000),
                      labels = xilabs))

# Vis ---------------------------------------------------------------------

probs <- seq(0, 1, length.out = 1e3)
p1 <- nll %>% 
  mutate(p = rep(probs, length(xilabs))) %>% 
  ggplot(aes(x = p, y = logLik, group = xi, color = xi)) +
  geom_line(show.legend = TRUE) +
  scale_color_viridis_d(name = parse(text = "log[10]~xi"), 
                        labels = parse(text = xilabs)) +
  theme(legend.position = c(0.35, 0.75)) +
  labs(x = parse(text = "alpha"),
       y = parse(text = "alpha*'-'*quantile~of~-logLik[i]")) +
  guides(color = guide_legend(ncol = 2)) +
  labs(tag = "A")

p2 <- cfx %>% 
  ggplot(aes(x = xi, y = -coef, color = xi)) +
  geom_boxplot(show.legend = FALSE) +
  geom_hline(aes(yintercept = 0.3), lty = 2, show.legend = FALSE) +
  labs(x = parse(text = "log[10]~xi"), y = parse(text = "hat(beta)")) +
  scale_x_discrete(labels = parse(text = xilabs)) +
  labs(tag = "B")

p1 + p2

ggsave("figures/vis-iv1.pdf", height = 3, width = 7)
