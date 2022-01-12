# Visualize scenario la
# Lucas Kook
# 05.11.2020

# Dependencies ------------------------------------------------------------

library(tidyverse)
library(ggpubr)
library(patchwork)
theme_set(theme_pubr())

# Read results ------------------------------------------------------------

res_lm <- read.csv("results/la/scenario-la-lin.csv") %>% 
  gather(key = "model", value = "ape", ape_anchor:ape_plain) %>% 
  separate(model, into = c("jnk", "model"), sep = "_")
res_Lm <- read.csv("results/la/scenario-la-Lm.csv") %>%
  gather(key = "model.run", value = "logLik") %>% 
  separate(model.run, into = c("model", "run"), sep = "\\.")

probs <- seq(0, 1, length.out = 1e3)

# Vis ---------------------------------------------------------------------

lmdat <- res_lm %>% 
  group_by(model, run) %>% 
  summarise(qfun = quantile(ape, probs = probs)) %>%
  mutate(p = probs)

slmdat <- lmdat %>% 
  group_by(model, p) %>% 
  summarise(median = median(qfun))

p1 <- lmdat %>%  
  ggplot(aes(x = p, y = qfun, group = paste(model, run), color = model)) +
  geom_line(show.legend = FALSE, alpha = 0.05) +
  geom_line(aes(x = p, y = median, group = model, color = model), 
            data = slmdat, show.legend = TRUE, inherit.aes = FALSE) +
  labs(x = parse(text = "alpha"), y = parse(text = "alpha*'-'*quantile~of~APE[i]")) +
  theme(legend.position = c(0.2, 0.8), legend.key.size = ggplot2::unit(0.8, "line")) +
  labs(tag = "a", color = element_blank())

Lmdat <- res_Lm %>% 
  group_by(model, run) %>% 
  summarise(qfun = quantile(-logLik, probs = probs)) %>%
  mutate(p = probs)

sLmdat <- Lmdat %>% 
  group_by(model, p) %>% 
  summarise(median = median(qfun))

p3 <- Lmdat %>%  
  ggplot(aes(x = p, y = qfun, group = paste(model, run), color = model)) +
  geom_line(show.legend = FALSE, alpha = 0.05) +
  geom_line(aes(x = p, y = median, group = model, color = model), 
            data = sLmdat, show.legend = FALSE, inherit.aes = FALSE) +
  labs(x = parse(text = "alpha"), y = parse(text = "alpha*'-'*quantile~of~NLL[i]")) +
  labs(tag = "b") +
  ylim(0, 200)

p1 + p3

ggsave("vis-la.jpg", height = 3, width = 7, dpi = 1200)
