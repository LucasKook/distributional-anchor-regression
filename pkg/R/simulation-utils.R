#' @export
evaluate <- function(mod, newdata, target) {
  UseMethod("evaluate")
}

#' @method evaluate ontram
#' @export
evaluate.ontram <- function(mod, newdata, target) {
  preds <- predict(mod, x = newdata, y = target)
  lli <- -log(c(preds$pdf * target)[as.logical(c(target))])
  return(lli)
}

#' @method evaluate default
#' @export
evaluate.default <- function(mod, newdata, target) {
  predicted <- predict(mod, newdata = data.frame(x = newdata))
  ape <- abs(target - predicted)
  mse <- mean((target - predicted)^2)
  ret <- list(ape = ape, mse = mse)
  return(ret)
}

#' @method evaluate keras.engine.sequential.Sequential
#' @export
evaluate.keras.engine.sequential.Sequential <- function(mod, newdata, target) {
  predicted <- mod(newdata)$numpy()
  ape <- abs(target - predicted)
  mse <- mean((target - predicted)^2)
  ret <- list(ape = ape, mse = mse)
  return(ret)
}

# Simple linear anchor regression scenario from Buehlmann (2018)
#' @export
generate_data_lin <- function(n = 300, r = 2, q = 1, p = 10,
                              sdy = 0.25, sdx = 1, gammas = rnorm(r),
                              shift = FALSE, shift_strength = 1) {
  epsY <- rnorm(n, sd = sdy)
  epsXj <- matrix(rnorm(n * p, sd = sdx), nrow = n, ncol = p)
  H <- matrix(rnorm(n * q), nrow = n, ncol = q)
  A <- matrix(rnorm(n * r), nrow = n, ncol = r)
  if (shift)
    A <- shift_strength * A

  X <- c(A %*% gammas + H) + epsXj
  Y <- 3 * X[, 2] + 3 * X[, 3] + H - 2 * A[, 1] + epsY

  ret <- list(Y = Y,
              X = X,
              A = A,
              H = H)
  return(ret)
}

# Scenario M2 from Buehlmann (2018)
#' @importFrom mvtnorm rmvnorm
#' @export
generate_data_m2 <- function(n = 300, r = 2, q = 1, p = 10, sdy = 0.25, sdx = 0.5,
                             gammas = rep(1, r), shift = FALSE, shift_strength = NULL) {
  epsY <- rnorm(n, sd = sdy)
  epsXj <- rmvnorm(n = n, sigma = sdx^2 * diag(p))
  H <- matrix(rnorm(n * q), nrow = n, ncol = q)
  if (shift) {
    mu <- rnorm(n = r, mean = 10, sd = 2)
    A <- rmvnorm(n = n, mean = mu, sigma = diag(r))
  } else {
    A <- rmvnorm(n = n, sigma = diag(r))
  }

  X <- c(A %*% gammas + 2 * H) + epsXj

  f <- function(X) {
    X[, 2] + X[, 3] + (X[, 2] > 0) + (X[, 2] <= -0.5) * (X[, 3] <= 1)
  }

  Y <- f(X) + 3 * H - 2 * A[, 1] + epsY

  ret <- list(Y = Y,
              X = X,
              A = A,
              H = H)
  return(ret)
}

#' ordinal IV scenario
#' @export
gen_dat_ordinal <- function(n = 1000, sd = 1, ba = 1, bx = 1, bh = 1,
                            perturb = FALSE, shift = 1.8, nclasses = 4) {
  epsX <- rnorm(n, sd = sd)
  H <- rnorm(n, sd = sd)
  if (!perturb) {
    A <- ifelse(runif(n) <= 0.5, -ba, ba)
  } else {
    A <- rep(shift, n)
  }
  X <- A + bh * H + epsX

  ord <- nclasses - 1L
  bl_cfx <- qlogis(seq(0, 1, length.out = ord + 2))[-c(1, ord + 2)]

  lY <- matrix(bl_cfx, ncol = ord, nrow = n, byrow = TRUE) - bx * X - bh * H
  pY <- t(apply(cbind(0, plogis(lY), 1), 1, diff))

  Y <- t(apply(pY, 1, function(p) rmultinom(n = 1, size = 1, prob = p)))

  return(list(X = matrix(X, ncol = 1), Y = Y, A = matrix(A, ncol = 1), H = matrix(H, ncol = 1)))
}

#' @export
anchor_simulation <- function(nsim = 100, alps = seq(0, 1, length.out = 11),
                              n = 300, nout = 2000, r = 2, q = 1, p = 10,
                              sdy = 0.25, sdx = 0.5,
                              gammas = rep(1, r), shift_strength = NULL,
                              generate_data = generate_data_m2,
                              anchor_model, plain_model, seed = NULL, ...) {
  ape_anchor <- list()
  if (!is.null(seed))
    set.seed(seed)
  for (iter in seq_len(nsim)) {
    message(paste0("Iteration ", iter))
    train_dat <- generate_data(n = n, r = r, q = q, p = p,
                               gammas = gammas, shift = FALSE)
    test_dat <- generate_data(n = nout, r = r, q = q, p = p,
                              shift = TRUE, gammas = gammas, shift_strength = shift_strength)
    mod_anchor <- anchor_model(data = train_dat, ...)
    mod_plain <- plain_model(data = train_dat, ...)
    met_anchor <- evaluate(mod = mod_anchor, newdata = test_dat$X,
                           target = test_dat$Y)
    met_plain <- evaluate(mod = mod_plain, newdata = test_dat$X,
                          target = test_dat$Y)
    ape_anchor[[iter]] <- data.frame(ape_anchor = met_anchor$ape,
                                     ape_plain = met_plain$ape)
  }
  return(ape_anchor)
}
