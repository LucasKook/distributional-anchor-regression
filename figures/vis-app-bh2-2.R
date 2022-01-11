# Vis application BH2 predicted densities
# Lucas Kook
# 09.11.2020

# Dependencies ------------------------------------------------------------

library(anchor)
library(tram)
library(survival)
library(CVXR)
library(tidyverse)
library(ggpubr)
theme_set(theme_pubr())

# Data --------------------------------------------------------------------

data("BostonHousing2", package = "mlbench")
nep <- 8e3

# Funs --------------------------------------------------------------------

Lm_anchor <- function(m0, xi, data) {
  A <- data$A
  dat <- data.frame(y = data$Y, x = data$X)
  trdat <- tramnet:::.get_tram_data(m0)
  xe <- data$X[trdat$exact$which, , drop = FALSE]
  nth <- trdat$npar
  nb <- ncol(data$X)
  theta <- Variable(nth)
  beta <- Variable(nb)
  Pia <- A %*% solve(t(A) %*% A) %*% t(A)
  resids <- (theta[1] + theta[2] * trdat$exact$ay[, 2] - xe %*% beta)
  v <- Pia %*% resids
  z <- trdat$exact$ay %*% theta - xe %*% beta
  const <- -log(sqrt(2 * pi)) * trdat$nobs
  ll <- const - sum(z^2/2) + sum_entries(log(trdat$exact$aypr %*% theta))
  obj <- - ll + xi * sum_entries(power(v, 2))
  const <- list(trdat$const$ui %*% theta >= trdat$const$ci)
  prob <- Problem(Minimize(obj), constraints = const)
  res <- solve(prob)
  class(res) <- "Lma"
  return(res)
}

logLik.Lma <- function(obj, newdata = NULL, ...) {
  if (is.null(newdata))
    return(-obj$value)
  ret <- dnorm(obj[[1]][1] + obj[[1]][2] * c(newdata$Y) -
                 c(newdata$X %*% obj[[2]]), log = TRUE) + log(obj[[1]][2])
  return(ret)
}

logLik.pLm <- function(obj, newdata = NULL, ...) {
  if (is.null(newdata)) {
    class(obj) <- class(obj)[-1L]
    return(logLik(obj))
  } else {
    cfx <- coef(obj, with_baseline = TRUE)
    cfb <- cfx[1:2]
    cfs <- cfx[-1:-2]
    ret <- dnorm(cfb[1] + cfb[2] * c(newdata$Y) -
                   c(newdata$X %*% cfs), log = TRUE) + log(cfb[2])
    return(ret)
  }
}

BoxCox_anchor <- function(m0, xi, data) {
  A <- data$A
  dat <- data.frame(y = data$Y, x = data$X)
  trdat <- tramnet:::.get_tram_data(m0)
  xe <- data$X[trdat$exact$which, , drop = FALSE]
  nth <- trdat$npar
  nb <- ncol(data$X)
  theta <- Variable(nth)
  beta <- Variable(nb)
  Pia <- A %*% solve(t(A) %*% A) %*% t(A)
  resids <- (trdat$exact$ay %*% theta - xe %*% beta)
  v <- Pia %*% resids
  z <- trdat$exact$ay %*% theta - xe %*% beta
  const <- -log(sqrt(2 * pi)) * trdat$nobs
  ll <- const - sum(z^2/2) + sum_entries(log(trdat$exact$aypr %*% theta))
  obj <- - ll + xi * sum_entries(power(v, 2))
  const <- list(trdat$const$ui %*% theta >= trdat$const$ci)
  prob <- Problem(Minimize(obj), constraints = const)
  res <- solve(prob)
  class(res) <- "BC"
  return(res)
}

logLik.BC <- function(obj, newdata = NULL, m0) {
  if (is.null(newdata))
    return(-obj$value)
  yb <- model.matrix(m0$model$bases$response, data = data.frame(cmedv = newdata$Y))
  ybp <- model.matrix(m0$model$bases$response, data = data.frame(cmedv = newdata$Y),
                      deriv = c("cmedv" = 1L))
  ret <- dnorm(c(yb %*% obj[[1]]) -
                 c(newdata$X %*% obj[[2]]), log = TRUE) + log(c(ybp %*% obj[[1]]))
  return(ret)
}

logLik.pBC <- function(obj, newdata = NULL, ...) {
  if (is.null(newdata)) {
    class(obj) <- class(obj)[-1L]
    return(logLik(obj))
  } else {
    cfx <- coef(obj, with_baseline = TRUE)
    nt <- length(cfx) - length(coef(obj))
    cfb <- cfx[1:nt]
    cfs <- cfx[-1:-nt]
    yb <- model.matrix(obj$model$bases$response,
                       data = data.frame(cmedv = newdata$Y, newdata$X))
    ybp <- model.matrix(obj$model$bases$response,
                        data = data.frame(cmedv = newdata$Y, newdata$X),
                        deriv = c("cmedv" = 1L))
    ret <- dnorm(c(yb %*% cfb) - c(newdata$X %*% cfs), log = TRUE) +
      log(c(ybp %*% cfb))
    return(ret)
  }
}

to_list <- function(fml, data) {
  mf <- model.frame(fml, data)
  Y <- model.response(mf)
  X <- anchor:::.rm_int(model.matrix(fml, data))
  return(list(Y = matrix(Y, ncol = 1), X = X))
}

.theta_to_gamma <- function(gammas) {
  thetas <- gammas
  thetas[2:length(thetas)] <- log(diff(gammas))
  thetas[is.nan(thetas)] <- - 1e20
  return(thetas)
}

warm_start <- function(m, mt, ord) {
  tmp <- get_weights(m$mod_baseline)
  cfx <- coef(mt, with_baseline = TRUE)
  tmp[[1]][] <- .theta_to_gamma(cfx[1:(ord + 1L)])
  set_weights(m$mod_baseline, tmp)
  tmp2 <- get_weights(m$mod_shift)
  tmp2[[1]][] <- coef(mt)
  set_weights(m$mod_shift, tmp2)
  return(invisible(m))
}

pred_dens <- function(m, cfx, newdata, type = "density", K = 1e4) {
  coef(m) <- cfx
  predict(m, type = type, newdata = newdata, K = K)
}

# Parms -------------------------------------------------------------------

txi <- 10

# Lm ----------------------------------------------------------------------

prs <- c("crim", "zn", "indus", "nox", "rm", "age", "lstat")
fml <- as.formula(paste("cmedv ~", paste(prs, collapse = "+")))
mLm <- as.mlt(Lm(fml, data = BostonHousing2, extrapolate = TRUE, 
                 support = supp <- range(BostonHousing2$cmedv), add = c(-5, 0)))
m0Lm <- Lm(cmedv ~ 1, data = BostonHousing2, extrapolate = TRUE, 
           support = supp, add = add <- c(-5, 0))
dat <- to_list(fml, BostonHousing2)
dat$A <- model.matrix(~ town, data = BostonHousing2)
aLm <- Lm_anchor(m0Lm, xi = txi, data = dat)
cfxLm <- c(aLm[[1]], aLm[[2]])
nd <- BostonHousing2[BostonHousing2$town == "Boston Beacon Hill", -6]
densLm <- pred_dens(mLm, cfxLm, nd)

mmLm <- as.mlt(Lm(fml, data = BostonHousing2, extrapolate = TRUE, 
                 support = supp <- range(BostonHousing2$cmedv), add = c(-5, 0)))
mdat <- to_list(fml, BostonHousing2)
mdat$A <- model.matrix(~ town, data = BostonHousing2)
maLm <- Lm_anchor(m0Lm, xi = 0, data = dat)
mcfxLm <- c(maLm[[1]], maLm[[2]])
mdensLm <- pred_dens(mmLm, mcfxLm, nd)

# BoxCox ------------------------------------------------------------------

mBC <- as.mlt(BoxCox(fml, data = BostonHousing2, extrapolate = TRUE, 
                     support = supp, add = add))
m0BC <- BoxCox(cmedv ~ 1, data = BostonHousing2, extrapolate = TRUE, 
               support = supp, add = add)
aBC <- BoxCox_anchor(m0BC, xi = txi, data = dat)
cfxBC <- c(aBC[[1]], aBC[[2]])
densBC <- pred_dens(mBC, cfxBC, nd)

mmBC <- as.mlt(BoxCox(fml, data = BostonHousing2, extrapolate = TRUE, 
                      support = supp, add = add))
maBC <- BoxCox_anchor(m0BC, xi = 0, data = dat)
mcfxBC <- c(maBC[[1]], maBC[[2]])
mdensBC <- pred_dens(mmBC, mcfxBC, nd)

# Exact (no) censoring ----------------------------------------------------

mCoc <- Colr(fml, data = BostonHousing2, 
             add = add,
             extrapolate = TRUE, support = supp)
m0Coc <- Colr(cmedv ~ 1, data = BostonHousing2, 
              add = add,
              extrapolate = TRUE, support = supp)
trdat <- tramnet:::.get_tram_data(m0Coc)
dat$y <- trdat$exact$ay
dat$yp <- trdat$exact$aypr
dat$cens <- rep(0, nrow(BostonHousing2))

mbl <- mod_baseline(8L)
msh <- mod_shift(ncol(dat$X))
me <- ontram(mod_bl = mbl, mod_sh = msh, x_dim = ncol(dat$X), y_dim = 7L,
            method = "logit", n_batches = 1, epochs = nep)
class(me) <- c("colrnn", class(me))
warm_start(me, mCoc, 6L)
mhe <- fit_colr_r(me, x_train = dat$X, y_train = dat$y, yp_train = dat$yp,
                 cens_train = dat$cens, history = FALSE, xi = txi, a_train = dat$A)
cfxCoc <- coef(mhe, with_baseline = TRUE)
densCoc <- pred_dens(as.mlt(mCoc), cfxCoc, nd)

mmbl <- mod_baseline(8L)
mmsh <- mod_shift(ncol(dat$X))
mme <- ontram(mod_bl = mmbl, mod_sh = mmsh, x_dim = ncol(dat$X), y_dim = 7L,
            method = "logit", n_batches = 1, epochs = nep)
class(me) <- c("colrnn", class(me))
warm_start(mme, mCoc, 6L)
mmhe <- fit_colr_r(mme, x_train = dat$X, y_train = dat$y, yp_train = dat$yp,
                 cens_train = dat$cens, history = FALSE, xi = 0, a_train = dat$A)
mcfxCoc <- coef(mmhe, with_baseline = TRUE)
mdensCoc <- pred_dens(as.mlt(mCoc), mcfxCoc, nd)

# Right censoring ---------------------------------------------------------

BostonHousing2$cmedvv <- with(BostonHousing2, Surv(cmedv, cmedv != 50))
fml3 <- as.formula(paste("cmedvv ~", paste(prs, collapse = "+")))
mCocr <- Colr(fml3, data = BostonHousing2, 
             add = add,
             extrapolate = TRUE, support = supp)
m0Cocr <- Colr(cmedv ~ 1, data = BostonHousing2, 
              add = add,
              extrapolate = TRUE, support = supp)
trdatr <- tramnet:::.get_tram_data(m0Cocr)
dat$y <- trdatr$exact$ay
dat$yp <- trdatr$exact$aypr
dat$cens <- as.numeric(BostonHousing2$cmedv == 50)

mbl <- mod_baseline(8L)
msh <- mod_shift(ncol(dat$X))
mr <- ontram(mod_bl = mbl, mod_sh = msh, x_dim = ncol(dat$X), y_dim = 7L,
            method = "logit", n_batches = 1, epochs = nep)
class(mr) <- c("colrnn", class(mr))
warm_start(mr, mCocr, 6L)
mhr <- fit_colr_r(mr, x_train = dat$X, y_train = dat$y, yp_train = dat$yp,
                 cens_train = dat$cens, history = FALSE, xi = txi, a_train = dat$A)
cfxCocr <- coef(mhr, with_baseline = TRUE)
densCocr <- pred_dens(as.mlt(mCocr), cfxCocr, nd)

mmbl <- mod_baseline(8L)
mmsh <- mod_shift(ncol(dat$X))
mmr <- ontram(mod_bl = mmbl, mod_sh = mmsh, x_dim = ncol(dat$X), y_dim = 7L,
            method = "logit", n_batches = 1, epochs = nep)
class(mr) <- c("colrnn", class(mr))
warm_start(mmr, mCocr, 6L)
mmhr <- fit_colr_r(mmr, x_train = dat$X, y_train = dat$y, yp_train = dat$yp,
                 cens_train = dat$cens, history = FALSE, xi = 0, a_train = dat$A)
mcfxCocr <- coef(mmhr, with_baseline = TRUE)
mdensCocr <- pred_dens(as.mlt(mCocr), mcfxCocr, nd)

# Vis ---------------------------------------------------------------------

dLm <- data.frame(cmedv = as.numeric(row.names(densLm)), densLm) %>% 
  gather(key = "obs", value = "dens", X1:X3) %>% 
  mutate(model = "Lm", xi = 10)
mdLm <- data.frame(cmedv = as.numeric(row.names(mdensLm)), mdensLm) %>% 
  gather(key = "obs", value = "dens", X1:X3) %>% 
  mutate(model = "Lm", xi = 0)
dBC <- data.frame(cmedv = as.numeric(row.names(densBC)), densBC) %>% 
  gather(key = "obs", value = "dens", X1:X3) %>% 
  mutate(model = "c-probit", xi = 10)
mdBC <- data.frame(cmedv = as.numeric(row.names(mdensBC)), mdensBC) %>% 
  gather(key = "obs", value = "dens", X1:X3) %>% 
  mutate(model = "c-probit", xi = 0)
dCoc <- data.frame(cmedv = as.numeric(row.names(densCoc)), densCoc) %>%
  gather(key = "obs", value = "dens", X1:X3) %>%
  mutate(model = "c-logit (exact)", xi = 10)
mdCoc <- data.frame(cmedv = as.numeric(row.names(mdensCoc)), mdensCoc) %>%
  gather(key = "obs", value = "dens", X1:X3) %>%
  mutate(model = "c-logit (exact)", xi = 0)
dCocr <- data.frame(cmedv = as.numeric(row.names(densCocr)), densCocr) %>% 
  gather(key = "obs", value = "dens", X1:X3) %>% 
  mutate(model = "c-logit (censored)", xi = 10)
mdCocr <- data.frame(cmedv = as.numeric(row.names(mdensCocr)), mdensCocr) %>% 
  gather(key = "obs", value = "dens", X1:X3) %>% 
  mutate(model = "c-logit (censored)", xi = 0)

pdat <- full_join(dLm, dBC) %>% full_join(., dCoc) %>% 
  full_join(., dCocr) %>% full_join(., mdLm) %>% 
  full_join(., mdBC) %>% full_join(., mdCocr) %>% full_join(., mdCoc) %>% 
  mutate(model = factor(model, levels = c("Lm", "c-probit", "c-logit (exact)", "c-logit (censored)"),
                        labels = c("Lm", "c-probit", "c-logit (exact)", "c-logit (censored)")),
         obs = factor(obs, levels = c("X1", "X2", "X3"), labels = c("Loc 1", "Loc 2", "Loc 3")))
pdat2 <- BostonHousing2[BostonHousing2$town == "Boston Beacon Hill", ]

ggplot(pdat, aes(x = cmedv, y = dens, color = ordered(xi), linetype = obs)) +
  geom_line() +
  labs(y = parse(text = "hat(f)[Y*'|'*x](y*'|'*x)")) +
  scale_linetype_discrete(name = element_blank()) +
  scale_color_viridis_d(name = parse(text = "xi"), option = "C", end = 0.8) +
  scale_x_continuous(breaks = seq(0, 50, 10)) +
  facet_wrap(~ model, ncol = 2) +
  geom_vline(xintercept = 50, data = NULL, lty = 2, color = "gray60")

ggsave("figures/vis-app-bh2-dens-revised.eps", width = 7, height = 6)
