#' @export
get_Wg <- function(A, gamma) {
  PiA <- A %*% solve(t(A) %*% A) %*% t(A)
  Wg <- diag(nrow(A)) - (1 - sqrt(gamma)) * PiA
  return(Wg)
}

#' @export
anchor_trafo <- function(data, gamma = 7) {
  A <- data$A
  Y <- data$Y
  X <- data$X
  Wg <- get_Wg(A, gamma)
  tY <- Wg %*% Y
  tX <- Wg %*% X
  ret <- list(tY = tY, tX = tX)
  return(ret)
}

#' @export
anchor_regression <- function(data, gamma = 7) {
  trafo <- anchor_trafo(data = data, gamma = gamma)
  tY <- trafo$tY
  tX <- trafo$tX
  dat <- data.frame(y = tY, x = tX)
  ret <- lm(y ~ 0 + ., data = dat)
  class(ret) <- c("lin_ar", class(ret))
  return(ret)
}

#' @export
plain_linear <- function(data, gamma = NULL) {
  dat <- data.frame(y = data$Y, x = data$X)
  ret <- lm(y ~ 0 + ., data = dat)
  return(ret)
}

#' @method predict lin_ar
#' @export
predict.lin_ar <- function(obj, newdata = NULL, ...) {
  class(obj) <- class(obj)[-1L]
  if (!is.null(newdata)) {
    ret <- predict(obj, newdata = newdata, ...)
  } else {
    ret <- predict(obj, ...)
  }
  return(ret)
}
