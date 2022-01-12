# Simulation scenario nla but with misspecified la
# 04.11.2020
# Lucas Kook

# Dependencies ------------------------------------------------------------

library(anchor)
library(tram)
library(CVXR)
library(tidyverse)
library(ggpubr)
library(patchwork)
theme_set(theme_pubr())

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

logLik.Lma <- function(obj, newdata = NULL) {
  if (is.null(newdata))
    return(-obj$value)
  ret <- dnorm(obj[[1]][1] + obj[[1]][2] * c(newdata$Y) -
                 c(newdata$X %*% obj[[2]]), log = TRUE) + log(obj[[1]][2])
  return(ret)
}

logLik.plain <- function(obj, newdata = NULL) {
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

from <- function(todistr, p, dd, cfx, ord, support, add) {
  yvar <- numeric_var("y", support = support, add = add)
  yB <- Bernstein_basis(yvar, order = ord, ui = "increasing", log_first = FALSE)
  st <- as.basis(as.formula(
    paste("~", paste0("X", seq_len(p), collapse = " + "), "- 1") ), data = dd)
  m <- ctm(response = yB, shifting = st, todistr = todistr)
  coef(m) <- cfx
  return(m)
}

qFUN <- function(x, dist) {
  ret <- switch(dist, "Normal" = qnorm(x), "Logistic" = qlogis(x),
                "MinExtrVal" = log(-log(1 - x)),
                "MaxExtrVal" = -log(-log(x)))
  return(ret)
}

gen_dat_iv1 <- function(n = 1000, sdX = 0.75, sdH = 0.75, ba = 0.3,
                        bh = 0.6, ord = 6, bx = cfs <- 0.3, shift = FALSE,
                        perturb = 12, r = 1, q = 1, gammas = 1, shift_strength = 1) {
  epsX <- rnorm(n, sd = sdX)
  H <- rnorm(n, sd = sdH)
 
  if (!shift) {
    A <- ifelse(runif(n) <= 0.5, -1, 1)
  } else {
    A <- perturb
  }
  
  X <- ba * A + bh * H + epsX
  
  pY <- function(x) pchisq(x, df = 3)
  qY <- function(p) qchisq(p, df = 3)
  
  supp <- qY(c(1e-12, 1 - 1e-12))
  add <- c(-1, 1) * 5
  distr <- "Normal"
  
  xx <- seq(from = supp[1], to = supp[2], length.out = ord + 1)
  bl_cfx <- qFUN(pY(xx), dist = distr)
  cfs <- c(bx, # X
           bh) # H
  
  dat <- data.frame(X1 = X, X2 = H)
  
  simm <- from(distr, 2, dat, c(bl_cfx, cfs), ord, supp, add)
  Y <- trtf:::.R2vec(simulate(simm, newdata = dat))
  ret <- list(Y = matrix(Y, ncol = 1),
              X = matrix(X, ncol = 1),
              H = matrix(H, ncol = 1),
              A = matrix(A, ncol = 1))
  return(ret)
}

# Simulation --------------------------------------------------------------

nsim <- 100
res_lin <- anchor_simulation(anchor_model = anchor_regression,
                             plain_model = plain_linear, gammas = rnorm(2),
                             generate_data = gen_dat_iv1,
                             gamma = 15, seed = 7, nsim = nsim)

# Write -------------------------------------------------------------------

out1 <- bind_rows(res_lin, .id = "run")
write.csv(out1, "results/iv1/scenario-iv1-lm.csv", quote = FALSE, row.names = FALSE)

# Lm ----------------------------------------------------------------------

vlla <- list()

set.seed(7)
for (iter in seq_len(nsim)) {
  message(paste0("Iteration ", iter))
  train_dat <- gen_dat_iv1(shift = FALSE)
  test_dat <- gen_dat_iv1(n = 2000, shift = TRUE)
  A <- train_dat$A
  dat <- data.frame(y = train_dat$Y, x = train_dat$X)
  m0 <- Lm(y ~ 1, data = dat)
  m_plain <- Lm(y ~ ., data = dat)
  class(m_plain) <- c("plain", class(m_plain))
  m_anchor <- Lm_anchor(m0 = m0, xi = 7, data = train_dat)
  vlla[[iter]] <- data.frame(anchor = logLik(m_anchor, newdata = test_dat),
                             plain = logLik(m_plain, newdata = test_dat))
}

# Write -------------------------------------------------------------------

out2 <- bind_rows(vlla, .id = "run")
write.csv(out2, "results/iv1/scenario-iv1-Lma.csv", quote = FALSE, row.names = FALSE)
