# Visualize scenario nla
# Lucas Kook
# 05.11.2020

# Dependencies ------------------------------------------------------------

library(tidyverse)
library(ggpubr)
library(patchwork)
theme_set(theme_pubr())

bpath <- "../experiments/"

# Read results ------------------------------------------------------------

res_Lm <- read.csv(file.path(bpath, "scenario-nla-Lm.csv")) %>% 
  gather(key = "model.run", value = "logLik") %>% 
  separate(model.run, into = c("model", "run"), sep = "\\.")

# Vis ---------------------------------------------------------------------

p2 <- res_Lm %>% 
  group_by(model, run) %>% 
  summarise(nll = - mean(logLik)) %>% 
  ggplot(aes(x = model, y = nll, fill = model, color = model)) +
  geom_dotplot(show.legend = TRUE, binaxis = "y", stackdir = "center", 
               binwidth = 0.4, dotsize = 2) +
  labs(x = element_blank(), y = "NLL") +
  theme(legend.position = c(0.2, 0.8), legend.key.size = unit(0.8, "line")) +
  scale_color_discrete(name = element_blank()) +
  scale_fill_discrete(name = element_blank()) +
  labs(tag = "A")

probs <- seq(0, 1, length.out = 1e3)

p3 <- res_Lm %>% 
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

p2 + p3

ggsave("figures/vis-nla.pdf", height = 3, width = 7)
