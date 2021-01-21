#' anchor logLik ontram
#' @examples
#' y <- k_constant(matrix(c(0, 0, 1, 0, 0, 0,
#'                          0, 1, 0, 0, 0, 0), nrow = 2L, ncol = 6L, byrow = TRUE))
#' gammas <- k_constant(matrix(c(-2, 0, 0, 0, 0), nrow = 1L, ncol = 5L, byrow = TRUE))
#' prm <- k_constant(matrix(c(0.46, 0.04, 0.04, 0.1), nrow = 2L, ncol = 2L, byrow = TRUE))
#' compute_logLik_a(gammas, y_train = y, prm = prm, xi = 10)
#' @export
compute_logLik_a <- function(gammas, beta = NULL, eta = NULL, y_train, prm,
                            distr = k_sigmoid, xi) {
  thetas <- gamma_to_theta(gammas)
  yu <- tf$pad(y_train, matrix(c(0L, 1L, 0L, 0L), ncol = 2))
  yl <- tf$pad(y_train, matrix(c(0L, 0L, 0L, 1L), ncol = 2))
  intr_upper <- tf$linalg$matmul(yu, tf$transpose(thetas))
  intr_lower <- tf$linalg$matmul(yl, tf$transpose(thetas))
  if (all(dim(intr_upper) == c(nrow(yu), nrow(yl)))) {
    intr_upper <- tf$linalg$diag_part(intr_upper)
    intr_lower <- tf$linalg$diag_part(intr_lower)
  }
  if (is.null(beta))
    beta <- k_zeros_like(intr_upper)
  if (is.null(eta))
    eta <- k_zeros_like(intr_upper)
  if (length(dim(beta)) == 1L || length(dim(eta)) == 1L) {
    if (length(dim(intr_upper)) > 1L) {
      intr_upper <- k_squeeze(intr_upper, 2L)
      intr_lower <- k_squeeze(intr_lower, 2L)
    }
  }
  pu <- distr(intr_upper - beta - eta)
  pl <- distr(intr_lower - beta - eta)
  li <- pu - pl
  resids <- (pl * (1 - pl) - pu * (1 - pu)) / (pu - pl)
  pen <- tf$reduce_mean(tf$square(tf$matmul(prm, resids)))
  nll <- - tf$reduce_mean(tf$math$log(li + 1e-16)) + xi * pen
  return(nll)
}

#' anchor fit ontram
#' @export
fit_ontram_a <- function(model, history = FALSE, x_train = NULL,
                        y_train, img_train = NULL, save_model = FALSE,
                        x_test = NULL, y_test = NULL, img_test = NULL,
                        xi, a_train, verbose = 1) {
  stopifnot(nrow(x_train) == nrow(y_train))
  stopifnot(ncol(y_train) == model$y_dim)
  apply_gradient_tf <- tf_function(apply_gradient_a)
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
      y_batch <- tf$constant(ontram:::.batch_subset(y_train, idx, dim(y_train)),
                             dtype = "float32")
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
      apply_gradient_tf(x_batch, y_batch, model, img_batch,
                        response_varying = model$response_varying,
                        xi = xi, prm = batch_prm)
    }
    if (history) {
      train_loss <- predict(model, x = x_train, y = y_train, im = img_train)$negLogLik
      test_loss <- predict(model, x = x_test, y = y_test, im = img_test)$negLogLik
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

#' anchor ontram gradient
#' @export
apply_gradient_a <- function(x_train, y_train, model, img_train = NULL,
                            verbose = FALSE, response_varying = FALSE,
                            xi = 0, prm = 0) {
  with(tf$GradientTape() %as% tape, {
    if (response_varying) {
      fwd_gamma <- model$mod_baseline(img_train)
    } else {
      fwd_gamma <- model$mod_baseline(matrix(1))
      fwd_gamma <- k_reshape(fwd_gamma, c(1L, model$y_dim - 1L))
    }
    fwd_beta <- NULL
    if (!is.null(x_train) & !is.null(model$mod_shift)) {
      fwd_beta <- model$mod_shift(x_train)
    }
    fwd_eta <- NULL
    if (!is.null(img_train) & !is.null(model$mod_image)) {
      fwd_eta <- model$mod_image(img_train)
    }
    nll <- compute_logLik_a(gammas = fwd_gamma, beta = fwd_beta, eta = fwd_eta,
                           y_train = y_train, distr = model$distr, xi = xi,
                           prm = prm)
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
