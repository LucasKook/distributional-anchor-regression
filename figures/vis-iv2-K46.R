# Evaluate scenario iv2 K=4,6
# Lucas Kook
# Oct 2021

# Dependencies ------------------------------------------------------------

library(tidyverse)
library(ggpubr)
library(patchwork)
theme_set(theme_pubr())

bpath <- "results/iv2-K46" 

# Read coefs --------------------------------------------------------------

cfx_files <- 
  data.frame(files = list.files(bpath, pattern = "cfx"), stringsAsFactors = FALSE)

lvn <- paste0("n", ns <- c(300, 500, 1000))
lbn <- paste0("n==", ns)
lvncl <- paste0("ncl", ncls <- c(4, 6, 10))
lbncl <- paste0("K==", ncls)
lvs <- paste0("shift", ss <- c(0, 1, 1.8, 3))
lbs <- paste0("do(A==", ss, ")")
lvba <- paste0("ba", bas <- c(-1, 0.5, 1, 2))
lbba <- paste0("M[X]==", bas)
lvxi <- paste0("xi", xis <- c(0, 10^(0:4)))
lbxi <- c("-infinity", 0:4)

dat <- cfx_files %>% 
  mutate(files_sans_ext = tools::file_path_sans_ext(files)) %>% 
  separate(files_sans_ext, into = c("type", "model", "sim", "parms"), sep = "_") %>% 
  separate(parms, into = c("n", "ncl", "bx", "ba", "bh", "sd", "shift", "xi"), sep = "-(?![0-9])") %>%
  mutate(n = factor(n, levels = lvn, labels = lbn),
         ncl = factor(ncl, levels = lvncl, labels = lbncl),
         shift = factor(shift, levels = lvs, labels = lbs),
         ba = factor(ba, levels = lvba, labels = lbba),
         xi = factor(xi, levels = lvxi, labels = lbxi)) %>% 
  mutate(dat = lapply(paste(bpath, files, sep = "/"), read_csv)) %>% 
  unnest(cols = dat)

p1 <- dat %>% 
  filter(n == "n==1000", ba != "M[X]==-1", model == "anchor") %>% 
  ggplot(aes(x = xi, y = x, color = xi)) +
  geom_boxplot(show.legend = FALSE) +
  facet_grid(ncl ~ ba, labeller = label_parsed) +
  geom_hline(aes(yintercept = 0.5), lty = 2, show.legend = FALSE) +
  labs(x = expression(log[10]~xi), y = expression(hat(beta))) +
  theme(legend.position = c(0.1, 0.9), text = element_text(size = 9)) +
  scale_color_viridis_d() +
  scale_x_discrete(labels = parse(text = levels(droplevels(dat$xi))))

# Read logLiks ------------------------------------------------------------

ll_files <- 
  data.frame(files = list.files(bpath, pattern = "ll"), stringsAsFactors = FALSE)

dat_ll <- ll_files %>% 
  mutate(files_sans_ext = tools::file_path_sans_ext(files)) %>% 
  separate(files_sans_ext, into = c("type", "model", "sim", "parms"), sep = "_") %>% 
  separate(parms, into = c("n", "ncl", "bx", "ba", "bh", "sd", "shift", "xi"), sep = "-(?![0-9])") %>%
  mutate(n = factor(n, levels = lvn, labels = lbn),
         ncl = factor(ncl, levels = lvncl, labels = lbncl),
         shift = factor(shift, levels = lvs, labels = lbs),
         ba = factor(ba, levels = lvba, labels = lbba),
         xi = factor(xi, levels = lvxi, labels = lbxi)) %>% 
  mutate(dat = lapply(paste(bpath, files, sep = "/"), read_csv)) %>% 
  unnest(cols = dat)

probs <- seq(0, 1, length.out = 1e3)
sum_dat <- dat_ll %>%
  gather("run", "logLik", starts_with("V")) %>%
  group_by(model, n, ncl, bx, ba, bh, sd, shift, xi) %>%
  summarise(qn = quantile(logLik, p = probs)) %>%
  mutate(p = probs)

sum_dat2 <- dat_ll %>% 
  gather("run", "logLik", starts_with("V")) %>% 
  group_by(model, n, ncl, bx, ba, bh, sd, shift, xi, run) %>% 
  summarise(NLL = mean(logLik))

p3 <- sum_dat2 %>% filter(shift != "do(A==0)", ba != "M[X]==-1", model == "anchor") %>% 
  ggplot(aes(x = xi, y = NLL, color = xi)) +
  geom_boxplot(show.legend = FALSE, outlier.size = 1) +
  facet_grid(ncl ~ ba, labeller = label_parsed) +
  scale_color_viridis_d(name = expression(log[10]~xi), 
                        labels = parse(text = levels(droplevels(sum_dat2$xi)))) +
  labs(x = expression(log[10]~xi), y = "NLL") +
  theme(text = element_text(size = 9)) +
  scale_x_discrete(labels = parse(text = levels(droplevels(sum_dat2$xi)))) +
  guides(color = guide_legend(direction = "horizontal", nrow = 1))

ggarrange(p3 + labs(tag = "a"), p1 + labs(tag = "b"), nrow = 2, 
          heights = c(2, 2), common.legend = FALSE, align = "v")

ggsave("figures/vis-iv2-K46.eps", height = 7.5, width = 8.5)
