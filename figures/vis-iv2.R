# Evaluate scenario iv2
# Lucas Kook
# 30.10.2020

# Dependencies ------------------------------------------------------------

library(tidyverse)
library(ggpubr)
library(patchwork)
theme_set(theme_pubr())

# Read coefs --------------------------------------------------------------

bpath <- "../experiments/"

cfx_files <- 
  data.frame(files = list.files(bpath, pattern = "cfx"), stringsAsFactors = FALSE)

lvn <- paste0("n", ns <- c(300, 500, 1000))
lbn <- paste0("n==", ns)
lvncl <- paste0("ncl", ncls <- c(4, 6, 10))
lbncl <- paste0("K==", ncls)
lvs <- paste0("shift", ss <- c(0, 1, 1.8, 3))
lbs <- paste0("do(A==", ss, ")")
lvba <- paste0("ba", bas <- c(-1, 0.5, 1, 2))
lbba <- paste0("zeta[2]==", bas)
lvxi <- paste0("xi", xis <- c(0, 10^(0:4)))
lbxi <- c("-infinity", 0:4)

dat <- cfx_files %>% 
  mutate(files_sans_ext = tools::file_path_sans_ext(files)) %>% 
  separate(files_sans_ext, into = c("type", "model", "sim", "parms"), sep = "_") %>% 
  separate(parms, into = c("n", "ncl", "bx", "ba", "bh", "sd", "shift", "xi"), 
           sep = "-(?![0-9])") %>%
  mutate(n = factor(n, levels = lvn, labels = lbn),
         ncl = factor(ncl, levels = lvncl, labels = lbncl),
         shift = factor(shift, levels = lvs, labels = lbs),
         ba = factor(ba, levels = lvba, labels = lbba),
         xi = factor(xi, levels = lvxi, labels = lbxi)) %>% 
  mutate(dat = lapply(paste(bpath, files, sep = "/"), read_csv)) %>% 
  unnest(cols = dat)

p1 <- dat %>% 
  filter(ncl == "K==10", n == "n==1000", ba != "zeta[2]==-1", model == "anchor") %>% 
  mutate(shift = "do(A%in%'{'*0*','*1*','*1.8*','*3*'}')") %>% 
  ggplot(aes(x = xi, y = x, color = xi)) +
  geom_boxplot(show.legend = FALSE) +
  facet_grid(shift ~ ba, labeller = label_parsed) +
  geom_hline(aes(yintercept = 0.5), lty = 2, show.legend = FALSE) +
  labs(x = expression(log[10]~xi), y = expression(hat(beta))) +
  theme(legend.position = c(0.1, 0.9), text = element_text(size = 9)) +
  scale_color_viridis_d() +
  scale_x_discrete(labels = parse(text = c("-infinity", 0:4)))

# Read logLiks ------------------------------------------------------------

ll_files <- 
  data.frame(files = list.files(bpath, pattern = "ll"), stringsAsFactors = FALSE)

dat_ll <- ll_files %>% 
  mutate(files_sans_ext = tools::file_path_sans_ext(files)) %>% 
  separate(files_sans_ext, into = c("type", "model", "sim", "parms"), sep = "_") %>% 
  separate(parms, into = c("n", "ncl", "bx", "ba", "bh", "sd", "shift", "xi"), 
           sep = "-(?![0-9])") %>%
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
  group_by(model, n, ncl, bx, ba, bh, sd, shift, xi, run) %>% 
  summarise(NLL = mean(logLik))

p3 <- sum_dat %>% filter(shift != "do(A==0)", ba != "zeta[2]==-1", 
                         model == "anchor") %>% 
  ggplot(aes(x = xi, y = NLL, color = xi)) +
  geom_boxplot(show.legend = FALSE, outlier.size = 1) +
  facet_grid(shift ~ ba, labeller = label_parsed) +
  scale_color_viridis_d(name = expression(log[10]~xi), 
                        labels = parse(text = levels(sum_dat$xi))) +
  labs(x = expression(log[10]~xi), y = "NLL") +
  theme(text = element_text(size = 9)) +
  scale_x_discrete(labels = parse(text = c("-infinity", 0:4))) +
  guides(color = guide_legend(direction = "horizontal", nrow = 1))

ggarrange(p3 + labs(tag = "A"), p1 + labs(tag = "B"), nrow = 2, 
          heights = c(2, 1), common.legend = FALSE, align = "v")
ggsave("figures/vis-iv2.pdf", height = 5.5, width = 6.5)
