# Visualize scenario nla
# Lucas Kook
# 05.11.2020

# Dependencies ------------------------------------------------------------

library(tidyverse)
library(ggpubr)
library(patchwork)
theme_set(theme_pubr())

# Read results ------------------------------------------------------------

res_Lm <- read.csv("scenario-nla.csv") %>%
  gather(key = "model.run", value = "logLik") %>% 
  separate(model.run, into = c("model", "run"), sep = "\\.")

# Vis ---------------------------------------------------------------------

probs <- seq(0, 1, length.out = 1e3)

pdat <- res_Lm %>% 
  group_by(model, run) %>% 
  summarise(qfun = quantile(-logLik, probs = probs)) %>%
  mutate(p = probs)

spdat <- pdat %>% 
  group_by(model, p) %>% 
  summarise(median = median(qfun))

mmax <- max(pdat$qfun)

p2 <- res_Lm %>% 
  group_by(model, run) %>% 
  summarise(nll = - mean(logLik)) %>% 
  ggplot(aes(x = model, y = nll, fill = model, color = model)) +
  geom_dotplot(show.legend = TRUE, binaxis = "y", stackdir = "center", 
               binwidth = 0.4, dotsize = 2) +
  labs(x = element_blank(), y = "mean NLL") +
  theme(legend.position = c(0.2, 0.8), 
        legend.key.size = ggplot2::unit(0.8, "line")) +
  scale_color_discrete(name = element_blank()) +
  scale_fill_discrete(name = element_blank()) +
  ylim(0, mmax) +
  labs(tag = "a")

p3 <- pdat %>% 
  ggplot(aes(x = p, y = qfun, group = paste(model, run), color = model)) +
  geom_line(show.legend = FALSE, alpha = 0.05) +
  geom_line(aes(x = p, y = median, group = model, color = model),
            data = spdat, show.legend = FALSE, inherit.aes = FALSE) +
  labs(x = parse(text = "alpha"),
       y = parse(text = "alpha*'-'*quantile~of~NLL[i]")) +
  ylim(0, mmax) +
  labs(tag = "b")

p2 + p3

ggsave("figures/vis-nla.jpg", height = 3, width = 7, dpi = 1200)
