# Vis scenario iv1 lm
# Lucas Kook
# Sep 21

# Deps --------------------------------------------------------------------

library(tidyverse)
library(ggpubr)
library(patchwork)
theme_set(theme_pubr())

xilabs <- c("-infinity", round(log10(c(1, 3, 7, 1000, 10000)), 2))
probs <- seq(0, 1, length.out = 1e3)

# Read --------------------------------------------------------------------

out1 <- read_csv("results/iv1/scenario-iv1-lm.csv")
out2 <- read_csv("results/iv1/scenario-iv1-Lma.csv")
nll <- read.csv("results/iv1/scenario-iv1-nll.csv") %>%
  gather(key = "model.xi", value = "logLik") %>% 
  separate(model.xi, into = c("model", "xi"), sep = "\\.") %>% 
  mutate(p = rep(probs, length(xilabs))) %>% 
  mutate(xi = ordered(xi, levels = c(0, 1, 3, 7, 1000, 10000),
                      labels = xilabs))

p0 <- nll %>% 
  ggplot(aes(x = p, y = logLik, group = xi, color = xi)) +
  geom_line(show.legend = TRUE) +
  scale_color_viridis_d(name = parse(text = "log[10]~xi"), labels = parse(text = xilabs)) +
  labs(x = parse(text = "alpha"),
       y = parse(text = "alpha*'-'*quantile~of~NLL[i]")) +
  theme(legend.position = c(0.35, 0.75)) +
  guides(color = guide_legend(ncol = 2)) +
  labs(tag = "a")

# Vis ---------------------------------------------------------------------

pdat <- out1 %>% group_by(run) %>% 
  summarise(anchor = quantile(ape_anchor, probs), plain = quantile(ape_plain, probs)) %>% 
  mutate(p = probs) %>% 
  gather("model", "quantile", anchor:plain)

spdat <- pdat %>% group_by(model, p) %>% summarise(median = median(quantile))

p1 <- ggplot(pdat, aes(x = p, y = quantile, group = paste(run, model), color = model)) +
  geom_line(alpha = 0.1) +
  geom_line(aes(x = p, y = median, color = model), inherit.aes = FALSE, data = spdat, lwd = 1) +
  labs(x = expression(alpha), y = expression(alpha*'-'*quantile~of~APE[i])) +
  theme(legend.position = c(0.2, 0.9)) +
  scale_color_discrete(name = element_blank())

pdat2 <- out2 %>% group_by(run) %>% 
  summarise(anchor = quantile(-anchor, probs), plain = quantile(-plain, probs)) %>% 
  mutate(p = probs) %>% 
  gather("model", "quantile", anchor:plain)

spdat2 <- pdat2 %>% group_by(model, p) %>% summarise(median = median(quantile))

pp2 <- ggplot(pdat2, aes(x = p, y = quantile, group = paste(run, model), color = model)) +
  geom_line(alpha = 0.1, show.legend = FALSE) +
  geom_line(aes(x = p, y = median, color = model), inherit.aes = FALSE, data = spdat2, 
            show.legend = FALSE) +
  labs(x = expression(alpha), y = expression(alpha*'-'*quantile~of~NLL[i]),
       color = element_blank()) +
  geom_line(aes(x = p, y = logLik, group = xi, linetype = xi), data = nll,
            inherit.aes = FALSE, show.legend = FALSE, alpha = 0.75) +
  guides(linetype = guide_legend(ncol = 3), color = "none")

sub <- ggplot(spdat2, aes(x = p, y = median, color = model)) +
  geom_line(show.legend = FALSE) +
  labs(x = element_blank(), y = element_blank(),
       color = element_blank(), linetype = element_blank()) +
  scale_linetype_discrete(name = parse(text = "c*'-'*probit:~log[10]~xi"), 
                          labels = parse(text = xilabs)) +
  geom_line(aes(x = p, y = logLik, group = xi, linetype = xi), data = nll,
            inherit.aes = FALSE, show.legend = TRUE, alpha = 0.75) +
  xlim(0.9, 1) +
  theme(text = element_text(size = 9), legend.position = c(0.42, 0.77)) +
  guides(linetype = guide_legend(ncol = 3), color = "none")

p2 <- pp2 + annotation_custom(ggplotGrob(sub),
                       xmin = 0, xmax = 0.9, ymin = 7, ymax = 25)

(p1 + labs(tag = "a")) + (p2 + labs(tag = "b"))

ggsave("figures/vis-iv1-lm.jpg", height = 3.5, width = 8, dpi = 1200)
