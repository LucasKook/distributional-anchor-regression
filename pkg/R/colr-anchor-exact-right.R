#' anchor logLik colr exact right
#' @examples
#' library(tram)
#' library(survival)
#' data("BostonHousing2", package = "mlbench")
#' BostonHousing2$cmedvv <- with(BostonHousing2, Surv(cmedv, cmedv != 50))
#' m <- Colr(cmedv ~ 1, data = BostonHousing2, order = 4)
#' mn <- Colr(cmedvv ~ nox, data = BostonHousing2, order = 4)
#' trdat <- tramnet:::.get_tram_data(m)
#' ntheta <- trdat$npar
#' y <- trdat$exact$ay
#' yp <- trdat$exact$aypr
#' x <- model.matrix(~ 0 + nox, data = BostonHousing2)
#' mbl <- mod_baseline(ntheta + 1L)
#' msh <- mod_shift(1L)
#' theta_to_gamma <- function(gammas) {
#'   thetas <- gammas
#'   thetas[2:length(thetas)] <- log(diff(gammas))
#'   thetas[is.nan(thetas)] <- - 1e20
#'   return(thetas)
#' }
#' gms <- theta_to_gamma(coef(mn, with_baseline = TRUE)[-(ntheta + 1)])
#' tmp <- get_weights(mbl)
#' tmp[[1]][] <- 0
#' tmp[[1]][] <- gms
#' set_weights(mbl, tmp)
#' tmp2 <- get_weights(msh)
#' tmp2[[1]][] <- coef(mn)
#' set_weights(msh, tmp2)
#' mo <- ontram(mbl, msh, y_dim = ntheta, x_dim = 1L, n_batches = 1, epochs = 1e3, lr = 0)
#' gammas <- mbl(matrix(1.0, ncol = 1L, nrow = 1L))
#' betas <- k_squeeze(msh(x), axis = 0L)
#' cens <- as.numeric(BostonHousing2$cmedv == 50)
#' colr_loss_r(gammas, betas, NULL, k_constant(y), k_constant(yp), k_constant(cens), prm = 0, xi = 0)
#' - logLik(mn) / nrow(BostonHousing2)
#' coef(mn, w = T); coef(mo, w = T)
#' colr_gradient_r(k_constant(x), k_constant(y), k_constant(yp), k_constant(cens), mo, verbose = TRUE)
#' logLik.colrnn_r(mo, x = x, y = y, yp = yp, cens = cens)
#' @export
colr_loss_r <- function(gammas, beta = NULL, eta = NULL, y, yp, cens, prm,
                            distr = k_sigmoid, xi, indiv = FALSE) {
  thetas <- gamma_to_theta_colr(gammas)
  ay <- tf$linalg$matvec(y, tf$squeeze(thetas, axis = 0L))
  aypr <- tf$linalg$matvec(yp, tf$squeeze(thetas, axis = 0L))
  if (is.null(beta))
    beta <- k_zeros_like(ay)
  if (is.null(eta))
    eta <- k_zeros_like(ay)
  if (length(dim(beta)) == 1L || length(dim(eta)) == 1L) {
    if (length(dim(ay)) > 1L) {
      ay <- k_squeeze(ay, 2L)
      aypr <- k_squeeze(aypr, 2L)
    }
  }
  if (length(dim(ay)) != length(dim(beta)))
    beta <- tf$squeeze(beta, 1L)
  pl <- distr(ay + beta + eta)
  li <- (1 - cens) * (pl * (1 - pl) * aypr) + cens * (1 - pl)
  if (indiv)
    return(tf$math$log(li + 1e-16))
  if (xi != 0) {
    resids <- (1 - cens) * (1 - exp(ay + beta)) / (1 + exp(ay + beta)) +
      cens * (pl * (1 - pl)) / (1 - pl)
    pen <- tf$reduce_mean(tf$square(tf$linalg$matvec(prm, resids)))
  } else {
    pen <- 0
  }
  nll <- - tf$reduce_mean(tf$math$log(li + 1e-16)) + xi * pen
  return(nll)
}

#' anchor fit colr
#' @examples
#' library(tram)
#' library(survival)
#' data("BostonHousing2", package = "mlbench")
#' BostonHousing2$cmedvv <- with(BostonHousing2, Surv(cmedv, cmedv != 50))
#' m <- Colr(cmedv ~ 1, data = BostonHousing2, order = 4, log_first = TRUE)
#' mn <- Colr(cmedvv ~ rm, data = BostonHousing2, order = 4, log_first = TRUE)
#' trdat <- tramnet:::.get_tram_data(m)
#' ntheta <- trdat$npar
#' y <- trdat$exact$ay
#' yp <- trdat$exact$aypr
#' x <- model.matrix(~ 0 + rm, data = BostonHousing2)
#' mbl <- mod_baseline(ntheta + 1L)
#' msh <- mod_shift(1L)
#' theta_to_gamma <- function(gammas) {
#'   thetas <- gammas
#'   thetas[2:length(thetas)] <- log(diff(gammas))
#'   thetas[is.nan(thetas)] <- - 1e20
#'   return(thetas)
#' }
#' gms <- theta_to_gamma(coef(mn, with_baseline = TRUE)[-(ntheta+1)])
#' tmp <- get_weights(mbl)
#' tmp[[1]][] <- 0
#' tmp[[1]][] <- gms
#' set_weights(mbl, tmp)
#' tmp2 <- get_weights(msh)
#' tmp2[[1]][] <- coef(mn)
#' set_weights(msh, tmp2)
#' mo <- ontram(mbl, msh, y_dim = ntheta, x_dim = 1L, n_batches = 1, epochs = 1e3, lr = 1e-3)
#' A <- model.matrix(~ town, data = BostonHousing2)
#' cens <- as.numeric(BostonHousing2$cmedv == 50)
#' fit_colr_r(mo, x_train = x, y_train = y, yp_train = yp, cens_train = cens, a_train = A, xi = 0)
#' coef(mn, w = T); coef(mo, w = T)
#' @export
fit_colr_r <- function(model, history = FALSE, x_train = NULL,
                     y_train, yp_train, cens_train, img_train = NULL,
                     save_model = FALSE, x_test = NULL, y_test = NULL,
                     img_test = NULL, xi, a_train, verbose = 1) {
  stopifnot(nrow(x_train) == nrow(y_train))
  stopifnot(ncol(y_train) == model$y_dim)
  apply_gradient_tf <- tf_function(colr_gradient_r)
  n <- nrow(y_train)
  start_time <- Sys.time()
  if (verbose != 0) {
    message("Training ordinal transformation model neural network.")
    message(paste0("Training samples: ", nrow(y_train)))
    message(paste0("Batches: ", bs <- model$n_batches))
    message(paste0("Batch size: ", ceiling(n/bs)))
    message(paste0("Epochs: ", nep <- model$epoch))
    pb <- txtProgressBar(min = 1, max = nep, style = 3)
  }
  prm <- a_train %*% solve(t(a_train) %*% a_train) %*% t(a_train)
  if (history) {
    model_history <- list(train_loss = c(), test_loss = c())
    class(model_history) <- "ontram_history"
  }
  for (epo in seq_len(nep)) {
    if (verbose == 2) {
      message(paste0("Training epoch: ", epo))
    } else if (verbose == 1) {
      setTxtProgressBar(pb, epo)
    }
    batch_idx <- sample(rep(seq_len(bs), ceiling(n/bs)), n)
    for (bat in seq_len(bs)) {
      idx <- which(batch_idx == bat)
      y_batch<- tf$constant(.batch_subset(y_train, idx, dim(y_train)),
                                   dtype = "float32")
      yp_batch<- tf$constant(.batch_subset(yp_train, idx, dim(yp_train)),
                                   dtype = "float32")
      cens_batch <- tf$constant(cens_train[idx], dtype = "float32")
      if (!is.null(x_train)) {
        x_batch <- tf$constant(ontram:::.batch_subset(x_train, idx, dim(x_train)),
                               dtype = "float32")
      } else {
        x_batch <- NULL
      }
      if (!is.null(img_train)) {
        img_batch <- lapply(img_train, function(x) tf$constant(ontram:::.batch_subset(x, idx, dim(x)),
                                                               dtype = "float32"))
      } else {
        img_batch <- NULL
      }
      if (!is.null(a_train)) {
        batch_prm <- k_constant(prm[,idx][idx,], dtype = "float32")
      } else {
        batch_prm <- 0
      }
      apply_gradient_tf(x_batch, y_batch, yp_batch, cens_batch, model, img_batch,
                        response_varying = model$response_varying,
                        xi = xi, prm = batch_prm)
    }
    if (history) {
      train_loss <- predict(model, x = x_train, y = y_train, cens_train = cens_train, im = img_train)$negLogLik
      test_loss <- predict(model, x = x_test, y = y_test, cens_train = cens_train, im = img_test)$negLogLik
      model_history$train_loss <- append(model_history$train_loss, train_loss)
      model_history$test_loss <- append(model_history$test_loss, test_loss)
    }
  }
  end_time <- Sys.time()
  if (verbose != 0)
    message(paste0("Training took ", end_time - start_time))
  if (history)
    return(model_history)
  return(invisible(model))
}

#' anchor colr gradient
#' @export
colr_gradient_r <- function(x_train, y_train, yp_train, cens_train, model,
                          img_train = NULL, verbose = FALSE,
                          response_varying = FALSE, xi = 0, prm = 0) {
  with(tf$GradientTape() %as% tape, {
    if (response_varying) {
      fwd_gamma <- model$mod_baseline(img_train)
    } else {
      fwd_gamma <- model$mod_baseline(matrix(1))
      fwd_gamma <- k_reshape(fwd_gamma, c(1L, model$y_dim))
    }
    fwd_beta <- NULL
    if (!is.null(x_train) & !is.null(model$mod_shift)) {
      fwd_beta <- model$mod_shift(x_train)
    }
    fwd_eta <- NULL
    if (!is.null(img_train) & !is.null(model$mod_image)) {
      fwd_eta <- model$mod_image(img_train)
    }
    nll <- colr_loss_r(gammas = fwd_gamma, beta = fwd_beta, eta = fwd_eta,
                           y = y_train, yp = yp_train, cens = cens_train,
                           distr = model$distr, xi = xi, prm = prm)
  })
  train_parms <- list(model$mod_baseline$trainable_variables)
  if (!is.null(model$mod_shift))
    train_parms <- append(train_parms, list(model$mod_shift$trainable_variables))
  if (!is.null(model$mod_image))
    train_parms <- append(train_parms, list(model$mod_image$trainable_variables))
  nabla <- tape$gradient(nll, train_parms)
  if (verbose)
    print(nabla)
  nparm <- length(train_parms)
  for (i in seq_len(nparm)) {
    model$optimizer$apply_gradients(
      purrr::transpose(list(nabla[[i]], train_parms[[i]]))
    )
  }
}

#' Function for predicting cdf, pdf and class and compute logLik
#' @method logLik colrnn_r
#' @export
logLik.colrnn_r <- function(model, x = NULL, y, yp, cens, im = NULL, indiv = FALSE) {
  y <- tf$constant(y, dtype = "float32")
  yp <- tf$constant(yp, dtype = "float32")
  gammas <- model$mod_baseline(matrix(1))
  if (!is.null(x)) {
    x <- tf$constant(x, dtype = "float32")
    betas <- k_squeeze(model$mod_shift(x), 2L)
  } else {
    betas <- tf$zeros_like(thetas)
  }
  if (!is.null(im)) {
    im <- tf$constant(im, dtype = "float32")
    etas <- model$mod_image(im)
  } else {
    etas <- NULL
  }
  nll <- colr_loss_r(gammas, betas, etas, y, yp, cens, prm = 0, xi = 0, indiv = indiv)$numpy()
  return(nll)
}

gamma_to_theta_colr <- function(gammas) {
  grows <- as.integer(nrow(gammas))
  gcols <- as.integer(ncol(gammas))
  # if (gcols == 1L)
  # return(gammas)
  theta1 <- k_reshape(gammas[, 1L], shape = c(grows, 1L))
  thetas <- k_exp(k_reshape(gammas[, 2L:gcols], shape = c(grows, gcols - 1L)))
  ret <- k_concatenate(c(theta1, theta1 + tf$math$cumsum(thetas, axis = 1L)), axis = 0L)
  return(ret)
}

# not very elegant, look for different solution
.batch_subset <- function(obj, idx, dim) {
  ndim <- length(dim)
  if (ndim == 2L) {
    ret <- obj[idx, , drop = FALSE]
  } else if (ndim == 3L) {
    ret <- obj[idx, , , drop = FALSE]
  } else if (ndim == 4L) {
    ret <- obj[idx, , , , drop = FALSE]
  }
  return(ret)
}
