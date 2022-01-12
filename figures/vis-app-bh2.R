# vis application BostonHousing2
# Lucas Kook
# 06.11.2020

# Dependencies ------------------------------------------------------------

library(tidyverse)
library(ggpubr)
library(patchwork)
theme_set(theme_pubr())
data("BostonHousing2", package = "mlbench")

bpath <- "results/app-bh2"
mlabs <- c("Lm", "c*'-'*probit", "c*'-'*logiti", "c*'-'*logit~exact", 
           "c*'-'*logit~censored")

# Read --------------------------------------------------------------------

cfs <- data.frame(path = list.files(bpath, pattern = "cfx")) %>% 
  separate(path, into = c("app", "bh", "mod", "xi", "cfx"), sep = "-", 
           remove = FALSE) %>% 
  separate(xi, into = c("parm", "xi"), sep = "xi") %>% 
  select(-app, -bh, -parm, -cfx) %>% 
  mutate(dat = lapply(paste0(bpath, path), read.csv)) %>% 
  unnest(dat) %>% 
  gather("predictor", "estimate", crim:lstat) %>% 
  mutate(mod = factor(mod, levels = c("Lm", "BoxCox", "Colr", "Colre", "Colrr"), 
                      labels = mlabs),
         xi = factor(xi, levels = c(0, 10^(0:4)), labels = c("-infinity", 0:4)),
         estimate = ifelse(str_detect(mod, "logit"), -estimate, estimate))

p1 <- cfs %>% 
  filter(model == "anchor", predictor %in% c("indus", "lstat", "nox", "rm"),
         mod != mlabs[3]) %>% 
  ggplot(aes(x = predictor, y = estimate, color = xi)) +
  geom_hline(aes(yintercept = 0), col = "gray65", lty = 2) +
  geom_boxplot(outlier.size = rel(0.5)) +
  facet_grid(~ mod, labeller = label_parsed) + 
  scale_color_viridis_d(name = parse(text = "log[10]~xi"),
                        labels = parse(text = levels(cfs$xi))) +
  theme(text = element_text(size = 12)) +
  labs(y = expression(hat(beta)[CV]), x = element_blank()) + 
  guides(color = guide_legend(nrow = 1))

lls <- data.frame(path = list.files(bpath, pattern = "logLik")) %>% 
  separate(path, into = c("app", "bh", "mod", "xi", "cfx"), sep = "-",
           remove = FALSE) %>% 
  separate(xi, into = c("parm", "xi"), sep = "xi") %>% 
  select(-app, -bh, -parm, -cfx) %>% 
  mutate(dat = lapply(paste0(bpath, path), function(x) {
    tmp <- read.csv(x)
    tmp$env <- levels(BostonHousing2$town)
    tmp
  })) %>% 
  unnest(dat) %>% 
  gather("model", "logLik", plain:anchor) %>%
  mutate(mod = factor(mod, levels = c("Lm", "BoxCox", "Colr", "Colre", "Colrr"), 
                      labels = mlabs),
         xi = factor(xi, levels = c(0, 10^(0:4)), labels = paste0(c("-infinity", 0:4))))

outliers <- lls %>% 
  filter(model == "anchor", mod != mlabs[3]) %>% 
  group_by(mod, xi) %>% 
  filter(env %in% c("Boston Beacon Hill", "Boston Back Bay", "Boston North End"))
  
p2 <- lls %>% 
  filter(model == "anchor", mod != mlabs[3]) %>%
  ggplot(aes(x = xi, y = -logLik, color = xi)) +
  geom_boxplot(outlier.size = rel(0.5), show.legend = TRUE) +
  facet_grid(~ mod, labeller = label_parsed) +
  labs(x = expression(log[10]~xi), y = expression(NLL[CV])) +
  theme(text = element_text(size = 12)) +
  scale_x_discrete(labels = parse(text = levels(lls$xi))) +
  scale_color_viridis_d(name = parse(text = "log[10]~xi"),
                        labels = parse(text = levels(cfs$xi))) +
  geom_path(aes(group = env), data = outliers, lwd = 0.5) +
  ggrepel::geom_text_repel(aes(label = env),
                           data = outliers %>% filter(xi == "-infinity"),
                           nudge_y = 1.5, nudge_x = 2, show.legend = FALSE,
                           cex = 2.8, color = "black") +
  guides(color = guide_legend(nrow = 1))
        

ggarrange(p2 + labs(tag = "a"), p1 + labs(tag = "b"), ncol = 1, align = "v",
          common.legend = TRUE, legend = "top")

ggsave("figures/vis-app-bh2-xi.eps", height = 5, width = 9)
