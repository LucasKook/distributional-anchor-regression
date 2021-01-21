# Visualize scenario la
# Lucas Kook
# 05.11.2020

# Dependencies ------------------------------------------------------------

library(tidyverse)
library(ggpubr)
library(patchwork)
theme_set(theme_pubr())

bpath <- "../experiments/"

# Read results ------------------------------------------------------------

res_lm <- read.csv(file.path(bpath, "scenario-la-lin.csv")) %>% 
  mutate(quantile = (0:10)/10) %>% 
  gather(key = "model", value = "mean_ape", mean_ape_anchor:mean_ape_plain)
res_Lm <- read.csv("results/la/scenario-la-Lm.csv") %>%
  gather(key = "model.run", value = "logLik") %>% 
  separate(model.run, into = c("model", "run"), sep = "\\.")

# Vis ---------------------------------------------------------------------

p1 <- ggplot(res_lm, aes(x = quantile, y = mean_ape, color = model)) +
  geom_point(show.legend = TRUE) +
  labs(y = parse(text = "alpha*'-'*quantile~of~APE"), x = parse(text = "alpha")) +
  scale_color_discrete(name = element_blank(), labels = c("anchor", "plain")) +
  theme(legend.position = c(0.2, 0.9)) +
  labs(tag = "A")

probs <- seq(0, 1, length.out = 1e3)

p2 <- res_Lm %>% 
  group_by(model, run) %>% 
  summarise(qfun = quantile(-logLik, probs = probs)) %>%
  mutate(p = probs) %>% 
  group_by(model, p) %>% 
  summarise(qfun = mean(qfun)) %>% 
  ggplot(aes(x = p, y = qfun, group = model, color = model)) +
  geom_line(show.legend = FALSE) +
  labs(x = parse(text = "alpha"), 
       y = parse(text = "alpha*'-'*quantile~of~-logLik[i]")) +
  labs(tag = "B")

p1 + p2

ggsave("figures/vis-la.pdf", height = 3, width = 7)
