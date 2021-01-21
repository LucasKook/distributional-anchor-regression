# Application: BostonHousing2
# Lucas Kook
# 06.11.2020

# Dependencies ------------------------------------------------------------

library(anchor)
library(tram)
library(survival)
library(CVXR)

# Data --------------------------------------------------------------------

data("BostonHousing2", package = "mlbench")

# Funs --------------------------------------------------------------------

Lm_anchor <- function(m0, xi, data) {
  A <- data$A
  dat <- data.frame(y = data$Y, x = data$X)
  trdat <- anchor:::.get_tram_data(m0)
  neg <- ifelse(trdat$negative, -1, 1)
  xe <- data$X[trdat$exact$which, , drop = FALSE]
  nth <- trdat$npar
  nb <- ncol(data$X)
  theta <- Variable(nth)
  beta <- Variable(nb)
  Pia <- A %*% solve(t(A) %*% A) %*% t(A)
  if (xi != 0) {
    resids <- (theta[1] + theta[2] * trdat$exact$ay[, 2] - xe %*% beta)
    v <- Pia %*% resids
    pen <- sum_entries(power(v, 2))
  } else {
    pen <- 0
  }
  z <- trdat$exact$ay %*% theta + neg * xe %*% beta
  const <- -log(sqrt(2 * pi)) * trdat$nobs
  ll <- const - sum(power(z, 2)/2) + sum_entries(log(trdat$exact$aypr %*% theta))
  obj <- - ll + xi * pen
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
  trdat <- anchor:::.get_tram_data(m0)
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

loeo_cv <- function(formula = fml, anchor = "town", data = BostonHousing2,
                    modp = Lm, moda = Lm_anchor, cls = "pLm", xi = 1) {
  nenv <- length(levels(data[[anchor]]))
  nparm <- length(strsplit(as.character(formula)[3], split = "\\+")[[1]])
  cfx_plain <- cfx_anchor <- matrix(nrow = nenv, ncol = nparm)
  ll_plain <- ll_anchor <- list()
  
  for (env in levels(data[[anchor]])) {
    run <- which(levels(data[[anchor]]) == env)
    message("run ", run, "/", nenv)
    idx <- which(data[[anchor]] == env)
    train <- data[-idx, ]
    train[anchor] <- droplevels(train[[anchor]])
    test <- to_list(formula, data[idx, ])
    A <- model.matrix(as.formula(paste("~ ", anchor)), data = train)
    X <- model.matrix(formula, data = train)[, -1L]
    mp <- modp(formula, data = train)
    class(mp) <- c(cls, class(mp))
    cfx_plain[run, ] <- coef(mp)
    ll_plain[[run]] <- logLik(mp, newdata = test)
    m0 <- modp(update(formula, . ~ 1), data = train)
    ltrain <- to_list(formula, data = train)
    ltrain$A <- A
    ma <- moda(m0, xi = xi, data = ltrain)
    cfx_anchor[run, ] <- ma[[2]]
    ll_anchor[[run]] <- logLik(ma, newdata = test, m0 = m0)
  }
 
  return(list(cfx_plain = cfx_plain, cfx_anchor = cfx_anchor,
              ll_plain = ll_plain, ll_anchor = ll_anchor))
}

write_loeo <- function(res, nm = "app-bh2", cnms, mod) {
  cfx <- data.frame(rbind(res$cfx_plain, res$cfx_anchor))
  cfx$model = rep(c("plain", "anchor"), each = nrow(cfx) / 2)
  colnames(cfx) <- c(cnms, "model")
  write.csv(cfx, paste0(nm, "-", mod, "-cfx.csv"), quote = FALSE, row.names = FALSE)
  ll <- data.frame(plain = do.call("c", lapply(res$ll_plain, mean)),
                   anchor = do.call("c", lapply(res$ll_anchor, mean)))
  write.csv(ll, paste0(nm, "-", mod, "-logLik.csv"), quote = FALSE, row.names = FALSE)
  return(list(cfx = cfx, ll = ll))
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

loeo_cv_nn_r <- function(formula, anchor = "town", data = BostonHousing2,
                         xi = 1, nep = 16e3, nb = 1, 
			 cens = as.numeric(BostonHousing2$cmedv == 50),
                         m0 = Colr(cmedv ~ 1, data = BostonHousing2)) {
  nenv <- length(levels(data[[anchor]]))
  nparm <- length(strsplit(as.character(formula)[3], split = "\\+")[[1]])
  cfx_plain <- cfx_anchor <- matrix(nrow = nenv, ncol = nparm)
  ll_plain <- ll_anchor <- list()
  trdat <- anchor:::.get_tram_data(m0)
  ord <- 6L
  
  for (env in levels(data[[anchor]])) {
    run <- which(levels(data[[anchor]]) == env)
    message("run ", run, "/", nenv)
    idx <- which(data[[anchor]] == env)
    train <- data[-idx, , drop = FALSE]
    train[anchor] <- droplevels(train[[anchor]])
    tr <- to_list(update(formula, cmedv ~ .), train)
    tr$y <- trdat$exact$ay[-idx, , drop = FALSE]
    tr$yp <- trdat$exact$aypr[-idx, , drop = FALSE]
    tr$cens <- cens[-idx]
    A <- model.matrix(as.formula(paste("~ ", anchor)), data = train)
    te <- to_list(update(formula, cmedv ~ .), data[idx, ])
    te$y <- trdat$exact$ay[idx, , drop = FALSE]
    te$yp <- trdat$exact$aypr[idx, , drop = FALSE]
    te$cens <- cens[idx]
    mCoc0 <- Colr(formula, train)
    mbl <- mod_baseline(ord + 2L)
    msh <- mod_shift(ncol(tr$X))
    m <- ontram(mod_bl = mbl, mod_sh = msh, x_dim = ncol(tr$X), y_dim = ord + 1L,
                method = "logit", n_batches = nb, epochs = nep)
    class(m) <- c("colrnn_r", class(m))
    warm_start(m, mCoc0, ord)
    mh <- fit_colr_r(m, x_train = tr$X, y_train = tr$y, yp_train = tr$yp, cens_train = tr$cens,
                   history = histo <- FALSE, xi = xi, a_train = A)
    mbl1 <- mod_baseline(ord + 2L)
    msh1 <- mod_shift(ncol(tr$X))
    mc <- ontram(mod_bl = mbl1, mod_sh = msh1, x_dim = ncol(tr$X), y_dim = ord + 1L,
                 method = "logit", n_batches = nb, epochs = nep)
    class(mc) <- c("colrnn_r", class(mc))
    warm_start(mc, mCoc0, ord)
    cfx_plain[run, ] <- coef(mc)
    ll_plain[[run]] <- logLik(mc, x = te$X, y = te$y, yp = te$yp, cens = te$cens, indiv = TRUE)
    cfx_anchor[run, ] <- coef(m)
    ll_anchor[[run]] <- logLik(m, x = te$X, y = te$y, yp = te$yp, cens = te$cens, indiv = TRUE)
  }
  return(list(cfx_plain = cfx_plain, cfx_anchor = cfx_anchor,
              ll_plain = ll_plain, ll_anchor = ll_anchor))
}

# Params ------------------------------------------------------------------

xis <- c(0, 10^(0:4))

# Lm ----------------------------------------------------------------------

prs <- c("crim", "zn", "indus", "nox", "rm", "age", "lstat")
fml <- as.formula(paste("cmedv ~", paste(prs, collapse = "+")))
mLm <- Lm(fml, data = BostonHousing2)

for (txi in xis) {
  rLm <- loeo_cv(modp = Lm, moda = Lm_anchor, cls = "pLm", xi = txi)
  dLm <- write_loeo(rLm, cnms = prs, mod = paste0("Lm-xi", txi))
}

# c-probit ----------------------------------------------------------------

mBC <- BoxCox(fml, data = BostonHousing2)
for (txi in xis) {
  rBC <- loeo_cv(modp = BoxCox, moda = BoxCox_anchor, cls = "pBC", xi = txi)
  dBC <- write_loeo(rBC, cnms = prs, mod = paste0("BoxCox-xi", txi))
}

# c-logit -----------------------------------------------------------------

BostonHousing2$cmedvv <- with(BostonHousing2, Surv(cmedv, as.numeric(cmedv != 50)))

fml3 <- as.formula(paste("cmedvv ~", paste(prs, collapse = "+")))
mCoc <- Colr(fml3, data = BostonHousing2)
for (txi in xis) {
  rCocr <- loeo_cv_nn_r(formula = fml3, anchor = "town", data = BostonHousing2, 
                     xi = txi)
  dCocr <- write_loeo(rCocr, cnms = prs, mod = paste0("Colrr-xi", txi))
}
