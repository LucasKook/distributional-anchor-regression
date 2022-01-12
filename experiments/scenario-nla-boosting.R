# scenario iv1
# Lucas Kook
# 05.11.2020

# Dependencies ------------------------------------------------------------

library(anchor)
library(tram)
library(trtf)
library(CVXR)
library(tidyverse)
library(ggpubr)
library(patchwork)
theme_set(theme_pubr())

# Funs --------------------------------------------------------------------

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

BoxCox_anchor <- function(m0, xi, data) {
  A <- data$A
  dat <- data.frame(y = data$Y, x = data$X)
  trdat <- anchor:::.get_tram_data(m0)
  xe <- data$X[trdat$exact$which, , drop = FALSE]
  nth <- trdat$npar
  nobs <- NROW(xe)
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
  prob <- Problem(Minimize(obj/nobs), constraints = const)
  res <- solve(prob)
  class(res) <- "Lma"
  return(res)
}

logLik.Lma <- function(obj, m0, newdata = NULL) {
  if (is.null(newdata))
    return(-obj$value)
  yb <- model.matrix(m0$model$bases$response, data = data.frame(y = newdata$Y))
  ybp <- model.matrix(m0$model$bases$response, data = data.frame(y = newdata$Y),
                      deriv = c("y" = 1L))
  ret <- dnorm(c(yb %*% obj[[1]]) -
                 c(newdata$X %*% obj[[2]]), log = TRUE) + log(c(ybp %*% obj[[1]]))
  return(ret)
}

logLik.plain <- function(obj, newdata = NULL) {
  if (is.null(newdata)) {
    class(obj) <- class(obj)[-1L]
    return(logLik(obj))
  } else {
    cfx <- coef(obj, with_baseline = TRUE)
    nt <- length(cfx) - length(coef(obj))
    cfb <- cfx[1:nt]
    cfs <- cfx[-1:-nt]
    yb <- model.matrix(obj$model$bases$response,
                       data = data.frame(y = newdata$Y, x = newdata$X))
    ybp <- model.matrix(obj$model$bases$response,
                        data = data.frame(y = newdata$Y, x = newdata$X),
                        deriv = c("y" = 1L))
    ret <- dnorm(c(yb %*% cfb) - c(newdata$X %*% cfs), log = TRUE) +
      log(c(ybp %*% cfb))
    return(ret)
  }
}

# Sim ---------------------------------------------------------------------

nsim <- 100
ntr <- 300
nte <- 2000
xis <- 6
out <- list()
cfxp <- cfxa <- list()

set.seed(241068)
for (iter in seq_len(nsim)) {
  message(paste0("Iteration ", iter))
  train_dat <- generate_data_m2(n = ntr)
  test_dat <- generate_data_m2(n = nte, shift = TRUE)
  A <- train_dat$A
  dat <- data.frame(y = train_dat$Y, x = train_dat$X)
  nd <- data.frame(y = test_dat$Y, x = test_dat$X)
  m0 <- BoxCox(y ~ 1, data = dat)
  m_plain <- BoxCox(y ~ ., data = dat, supp = range(dat$y, nd$y))
  class(m_plain) <- c("plain", class(m_plain))
  vllp <- logLik(m_plain, newdata = test_dat)
  ppr <- as.vector(trtf:::.R2vec(predict(as.mlt(m_plain), newdata = nd,
                                         type = "quantile", prob = 0.5)))
  ppreds <- abs(ppr - nd$y)
  cfxp[[iter]] <- coef(m_plain)
  pboost <- plain_lmrf(train_dat)
  bpr <- as.vector(evaluate(pboost, newdata = test_dat$X, target = test_dat$Y)$ape)
  for (xi in xis) {
    xidx <- which(xi == xis)
    m_anchor <- BoxCox_anchor(m0 = m0, xi = xi, data = train_dat)
    vlla <- logLik(m_anchor, m0 = m0, newdata = test_dat)
    cfxa[[iter]] <- m_anchor[[2]]
    mpp <- as.mlt(m_plain)
    cfa <- c(m_anchor[[1]], m_anchor[[2]])
    names(cfa) <- names(coef(mpp))
    coef(mpp) <- cfa
    ppa <- as.vector(trtf:::.R2vec(predict(mpp, newdata = nd, 
                                           type = "quantile", prob = 0.5)))
    apreds <- abs(ppa - nd$y)
    aboost <- anchor_boosting(lmrf_base_learner, train_dat, mstop = 50, 
                              gamma = 2 * xi + 1, nu = 0.01)
    abpr <- as.vector(evaluate(aboost, newdata = test_dat$X, 
                               target = test_dat$Y)$ape)
  }
  out[[iter]] <- data.frame(run = iter, plain.cprobit_logLik = -vllp, 
                            anchor.cprobit_logLik = -vlla, 
                            plain.cprobit_ape = ppreds, anchor.cprobit_ape = apreds,
                            plain.lmrf_ape = bpr, anchor.boosting_ape = abpr)
}

# Write -------------------------------------------------------------------

out1 <- dplyr::bind_rows(out)
write.csv(out1, "results/nla/scenario-nla-boosting.csv", row.names = FALSE, quote = FALSE)
