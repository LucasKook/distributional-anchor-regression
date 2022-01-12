#' @export
anchor_gradient <- function(Wg, target, predicted) {
  Wg2 <- Wg %*% Wg
  n <- nrow(target)
  ret <- Wg2 %*% (target - predicted) / n
  return(ret)
}

#' @export
anchor_boosting_step <- function(fprev, fupdate, nu = 0.1) {
  ret <- fprev + nu * fupdate
  return(ret)
}

#' @export
rf_base_learner <- function(pseudo_resp, X) {
  mdat <- data.frame(y = pseudo_resp,
                     x = X)
  m <- randomForest(y ~ ., data = mdat)
  class(m) <- c("rf_bl", class(m))
  return(m)
}

#' @method predict rf_bl
#' @export
predict.rf_bl <- function(obj, newdata = NULL, ...) {
  if (!is.null(newdata)) {
    ret <- predict(obj, newdata = newdata, ...)
  } else {
    ret <- predict(obj, ...)
  }
  return(ret)
}

#' @importFrom randomForest randomForest
#' @export
lmrf_base_learner <- function(pseudo_resp, X) {
  mdat <- data.frame(y = pseudo_resp,
                     x = X)
  mlm <- lm(y ~ ., data = mdat)
  mdat$res <- residuals(mlm)
  mdat <- mdat[,-1L]
  mrf <- randomForest(res ~ ., data = mdat)
  ret <- list(lm = mlm, rf = mrf)
  class(ret) <- "lmrf_bl"
  return(ret)
}

#' @method predict lmrf_bl
#' @export
predict.lmrf_bl <- function(obj, newdata = NULL, ...) {
  mlm <- obj[[1]]
  mrf <- obj[[2]]
  if (!is.null(newdata)) {
    ret <- predict(mlm, newdata = newdata, ...) +
      predict(mrf, newdata = newdata, ...)
  } else {
    ret <- predict(mlm, ...) + predict(mrf, ...)
  }
  return(ret)
}

#' @export
plain_lmrf <- function(data, gamma = NULL, base_learner = NULL, mstop = NULL) {
  ret <- lmrf_base_learner(data$Y, data$X)
  return(ret)
}

#' @export
anchor_boosting <- function(base_learner, data, mstop = 2, gamma = 7, nu = 0.1) {
  n <- nrow(data$X)
  fout <- matrix(0, nrow = n, ncol = mstop)
  bls <- list()
  Wg <- get_Wg(A = data$A, gamma = 7)
  target <- data$Y
  pb <- txtProgressBar(min = 0, max = mstop, width = 60)
  cat("Anchor Boosting:")
  for (m in seq_len(mstop)) {
    setTxtProgressBar(pb = pb, value = m)
    predicted <- fout[, max(c(1, m - 1)), drop = FALSE]
    nabla <- anchor_gradient(Wg = Wg, target = target, predicted = predicted)
    bl <- base_learner(pseudo_resp = nabla, X = data$X)
    fupdate <- as.matrix(predict(bl, newdata = data.frame(x = data$X)),
                         nrow = n, ncol = 1)
    fout[, m] <- anchor_boosting_step(fprev = predicted, fupdate = fupdate, nu = nu)
    bls[[m]] <- bl
  }
  close(pb)
  ret <- list(base_learners = bls, fitted = fout, nu = nu)
  class(ret) <- "anchor_boost"
  return(ret)
}

#' @method predict anchor_boost
#' @export
predict.anchor_boost <- function(obj, newdata, ...) {
  preds <- lapply(obj$base_learners, predict, newdata = newdata)
  aggr <- do.call("cbind", preds)
  ret <- apply(aggr, 1, sum) * obj$nu
  return(ret)
}
