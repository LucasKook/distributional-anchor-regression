# Vis boosting vs c-probit
# Lucas Kook
# Oct 2021

# Deps --------------------------------------------------------------------

library(tidyverse)
library(ggpubr)
library(patchwork)
theme_set(theme_pubr())

# Read --------------------------------------------------------------------

out <- read_csv("results/nla/scenario-nla-boosting.csv")

out1 <- out %>% 
  gather("tmp", "val", -1) %>% 
  separate("tmp", into = c("model", "loss"), sep = "_") %>% 
  separate("model", into = c("type", "model"), sep = "\\.") %>% 
  mutate(model = factor(case_when(
    model == "cprobit" ~ "c*'-'*probit",
    model == "boosting" ~ "L[2]~anchor~boosting",
    model == "lmrf" ~ "L[2]~anchor~boosting"
  )),
  loss = factor(case_when(
    loss == "ape" ~ "APE",
    loss == "logLik" ~ "NLL"
  )))

# Vis ---------------------------------------------------------------------

probs <- seq(0, 1, length.out = 1e3)
pdat <- out1 %>% 
  filter(loss == "APE") %>% 
  group_by(run, model, type) %>% 
  summarise(qfun = quantile(val, probs = probs)) %>% 
  mutate(p = probs)

spdat <- pdat %>% 
  group_by(model, type, p) %>% 
  summarise(median = median(qfun))

ggplot(pdat, aes(x = p, y = qfun, color = type, group = paste(run, model, type))) +
  geom_line(alpha = 0.05, show.legend = FALSE) +
  geom_line(aes(x = p, y = median, color = type), inherit.aes = FALSE, data = spdat) +
  facet_wrap(~ model, labeller = label_parsed) +
  labs(x = expression(alpha), y = expression(alpha*'-'*quantile~of~APE[i])) +
  scale_color_discrete(name = element_blank()) +
  theme(legend.position = c(0.1, 0.8))

ggsave("figures/vis-nla-boosting.jpg", height = 3, width = 7, dpi = 1200)
