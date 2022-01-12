# Simulation scenario la
# 04.11.2020
# Lucas Kook

# Dependencies ------------------------------------------------------------

library(anchor)
library(tram)
library(CVXR)

# Funs --------------------------------------------------------------------

Lm_anchor <- function(m0, xi, data) {
  A <- data$A
  dat <- data.frame(y = data$Y, x = data$X)
  trdat <- anchor:::.get_tram_data(m0)
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

# Simulation --------------------------------------------------------------

nsim <- 100
res_lin <- anchor_simulation(anchor_model = anchor_regression,
                             plain_model = plain_linear, gammas = rnorm(2),
                             generate_data = generate_data_lin,
                             gamma = 13, shift_strength = sqrt(10), sdx = 1,
                             seed = 7, nsim = nsim)

# Write -------------------------------------------------------------------

out1 <- dplyr::bind_rows(res_lin, .id = "run")
write.csv(out1, "results/la/scenario-la-lin.csv", quote = FALSE, row.names = FALSE)

# Lm ----------------------------------------------------------------------

nsim <- 100
vllp <- vlla <- matrix(nrow = n <- 2000, ncol = nsim)

set.seed(7)
for (iter in seq_len(nsim)) {
  message(paste0("Iteration ", iter))
  train_dat <- generate_data_lin(sdx = 1, gammas = rnorm(2), shift = FALSE)
  test_dat <- generate_data_lin(n = n, shift = TRUE, 
                                shift_strength = sqrt(10), sdx = 1)
  A <- train_dat$A
  dat <- data.frame(y = train_dat$Y, x = train_dat$X)
  m0 <- Lm(y ~ 1, data = dat)
  m_plain <- Lm(y ~ ., data = dat)
  class(m_plain) <- c("plain", class(m_plain))
  vllp[, iter] <- logLik(m_plain, newdata = test_dat)
  m_anchor <- Lm_anchor(m0 = m0, xi = 6, data = train_dat)
  vlla[, iter] <- logLik(m_anchor, newdata = test_dat)
}

# Write -------------------------------------------------------------------

out2 <- data.frame(anchor = vlla, plain = vllp)
write.csv(out2, "results/la/scenario-la-Lm.csv", quote = FALSE, row.names = FALSE)
