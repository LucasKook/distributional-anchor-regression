# Vis scenario iv1 lm
# Lucas Kook
# Sep 21

# Deps --------------------------------------------------------------------

library(tidyverse)
library(ggpubr)
library(patchwork)
theme_set(theme_pubr())

# Read --------------------------------------------------------------------

# Mis-specified lm/Lm results
out1 <- read_csv("results/nla/scenario-nla-lm.csv")
out2 <- read_csv("results/nla/scenario-nla-Lma.csv")

# Anchor cprobit results
res_cprobit <- read.csv("results/nla/scenario-nla.csv") %>%
  gather(key = "model.run", value = "logLik") %>% 
  separate(model.run, into = c("model", "run"), sep = "\\.")

probs <- seq(0, 1, length.out = 1e3)

cprobit_pdat <- res_cprobit %>% 
  group_by(model, run) %>% 
  summarise(qfun = quantile(-logLik, probs = probs)) %>%
  mutate(p = probs)

cprobit_spdat <- cprobit_pdat %>% 
  group_by(model, p) %>% 
  summarise(median = median(qfun))

mmax <- max(cprobit_pdat$qfun)

# Vis ---------------------------------------------------------------------

pdat <- out1 %>% group_by(run) %>% 
  summarise(anchor = quantile(ape_anchor, probs), 
            plain = quantile(ape_plain, probs)) %>% 
  mutate(p = probs) %>% 
  gather("model", "quantile", anchor:plain)

spdat <- pdat %>% group_by(model, p) %>% summarise(median = median(quantile))

p1 <- ggplot(pdat, aes(x = p, y = quantile, group = paste(run, model), color = model)) +
  geom_line(alpha = 0.1, show.legend = FALSE) +
  geom_line(aes(x = p, y = median, color = model), inherit.aes = FALSE, data = spdat, lwd = 1,
            show.legend = FALSE) +
  labs(x = expression(alpha), y = expression(alpha*'-'*quantile~of~APE[i])) +
  scale_color_discrete(name = element_blank())

pdat2 <- out2 %>% group_by(run) %>% 
  summarise(anchor = quantile(-anchor, probs), plain = quantile(-plain, probs)) %>% 
  mutate(p = probs) %>% 
  gather("model", "quantile", anchor:plain)

spdat2 <- pdat2 %>% group_by(model, p) %>% summarise(median = median(quantile))

p2 <- ggplot(pdat2, aes(x = p, y = quantile, group = paste(run, model), color = model)) +
  geom_line(alpha = 0.1, show.legend = FALSE) +
  geom_line(aes(x = p, y = median, color = model, linetype = "Lm"), inherit.aes = FALSE, data = spdat2, 
            lwd = 1, show.legend = TRUE) +
  theme(legend.position = c(0.4, 0.9), legend.box = "horizontal", 
        legend.text = element_text(size = 9)) +
  geom_line(aes(x = p, y = median, group = model, color = model, linetype = "c-probit"),
            data = cprobit_spdat, show.legend = TRUE, inherit.aes = FALSE,
            lwd = 1) +
  scale_linetype_manual(values = c("Lm" = 1, "c-probit" = 4)) +
  labs(x = expression(alpha), y = expression(alpha*'-'*quantile~of~NLL[i]),
       color = element_blank(), linetype = element_blank())

((p1 + labs(tag = "a")) + (p2 + labs(tag = "b")))

ggsave("figures/vis-nla-lm.jpg", height = 3, width = 7, dpi = 1200)
