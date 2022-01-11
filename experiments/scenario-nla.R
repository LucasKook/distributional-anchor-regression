# Scenario nla
# Lucas Kook
# 05.11.2020

# Dependencies ------------------------------------------------------------

library(anchor)
library(tram)
library(CVXR)

# Functions ---------------------------------------------------------------

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
    cfb <- cfx[1:7]
    cfs <- cfx[-1:-7]
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
vllp <- vlla <- matrix(nrow = n <- 2000, ncol = nsim)

set.seed(1)
for (iter in seq_len(nsim)) {
  message(paste0("Iteration ", iter))
  train_dat <- generate_data_m2(shift = FALSE)
  test_dat <- generate_data_m2(n = 2000, shift = TRUE)
  A <- train_dat$A
  dat <- data.frame(y = train_dat$Y, x = train_dat$X)
  m0 <- BoxCox(y ~ 1, data = dat)
  m_plain <- BoxCox(y ~ ., data = dat)
  class(m_plain) <- c("plain", class(m_plain))
  vllp[, iter] <- logLik(m_plain, newdata = test_dat)
  m_anchor <- BoxCox_anchor(m0 = m0, xi = 6, data = train_dat)
  vlla[, iter] <- logLik(m_anchor, m0 = m0, newdata = test_dat)
}

# Write -------------------------------------------------------------------

out <- data.frame(anchor = vlla, plain = vllp)
write.csv(out, "scenario-nla.csv", quote = FALSE, row.names = FALSE)
