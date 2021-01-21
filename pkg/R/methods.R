#' residuals ontram
#' @method residuals ontram
#' @export
residuals.ontram <- function(object, x = NULL, y, im = NULL) {
  cdf <- t(predict(object, x = x, y = y, im = im)$cdf)
  yu <- t(cbind(0, y))
  yl <- t(cbind(y, 0))
  pu <- c(cdf)[as.logical(c(yu))]
  pl <- c(cdf)[as.logical(c(yl))]
  resids <- (pl * (1 - pl) - pu * (1 - pu)) / (pu - pl)
  return(resids)
}

.get_tram_data <- function (mod) {
  dat <- mget(c("iY", "eY", "offset"), envir = environment(mod$parm),
              ifnotfound = list(NULL))
  out <- list()
  np <- length(mod$par)
  out$npar <- np
  out$nobs <- NROW(dat$iY$Yleft) + NROW(dat$eY$Y)
  if (is.null(dat$iY)) {
    out$censl <- list(ay = matrix(0, ncol = np, nrow = 0),
                      which = NULL)
    out$censr <- list(ay = matrix(0, ncol = np, nrow = 0),
                      which = NULL)
    out$censi <- list(ayl = matrix(0, ncol = np, nrow = 0),
                      ayr = matrix(0, ncol = np, nrow = 0), which = NULL)
  }
  else {
    idxr <- which(is.finite(dat$iY$Yleft[, 1]) & !is.finite(dat$iY$Yright[,
                                                                          1]))
    idxl <- which(!is.finite(dat$iY$Yleft[, 1]) & is.finite(dat$iY$Yright[,
                                                                          1]))
    idxi <- which(is.finite(dat$iY$Yleft[, 1]) & is.finite(dat$iY$Yright[,
                                                                         1]))
    out$censl <- list(ay = dat$iY$Yright[idxl, , drop = FALSE],
                      which = dat$iY$which[idxl])
    out$censr <- list(ay = dat$iY$Yleft[idxr, , drop = FALSE],
                      which = dat$iY$which[idxr])
    out$censi <- list(ayl = dat$iY$Yleft[idxi, , drop = FALSE],
                      ayr = dat$iY$Yright[idxi, , drop = FALSE], which = dat$iY$which[idxi])
  }
  if (is.null(dat$eY)) {
    out$exact <- list(ay = matrix(0, ncol = np, nrow = 0),
                      aypr = matrix(0, ncol = np, nrow = 0), which = NULL)
  }
  else {
    out$exact <- list(ay = dat$eY$Y, aypr = dat$eY$Yprime,
                      which = dat$eY$which)
  }
  out$offset <- dat$offset
  out$weights <- mod$weights
  out$negative <- mod$negative
  out$errdistr <- mod$todistr$name
  if (!is.null(dat$eY)) {
    out$constr <- attr(dat$eY$Y, "constraint")
  }
  else {
    out$constr <- attr(dat$iY$Yleft, "constraint")
  }
  idxi <- which(is.finite(dat$iY$trunc$left[, 1]) & !is.finite(dat$iY$trunc$right[,
                                                                                  1]))
  idxe <- which(is.finite(dat$eY$trunc$left[, 1]) & !is.finite(dat$eY$trunc$right[,
                                                                                  1]))
  out$truncl <- list(ay = rbind(dat$iY$trunc$left[idxi, , drop = FALSE],
                                dat$eY$trunc$left[idxe, , drop = FALSE]), which = c(dat$iY$which[idxi],
                                                                                    dat$eY$which[idxe]))
  idxi <- which(!is.finite(dat$iY$trunc$left[, 1]) & is.finite(dat$iY$trunc$right[,
                                                                                  1]))
  idxe <- which(!is.finite(dat$eY$trunc$left[, 1]) & is.finite(dat$eY$trunc$right[,
                                                                                  1]))
  out$truncr <- list(ay = rbind(dat$iY$trunc$right[idxi, ,
                                                   drop = FALSE], dat$eY$trunc$right[idxe, , drop = FALSE]),
                     which = c(dat$iY$which[idxi], dat$eY$which[idxe]))
  idxi <- which(is.finite(dat$iY$trunc$left[, 1]) & is.finite(dat$iY$trunc$right[,
                                                                                 1]))
  idxe <- which(is.finite(dat$eY$trunc$left[, 1]) & is.finite(dat$eY$trunc$right[,
                                                                                 1]))
  out$trunci <- list(ayl = rbind(dat$iY$trunc$left[idxi, ,
                                                   drop = FALSE], dat$eY$trunc$left[idxe, , drop = FALSE]),
                     ayr = rbind(dat$iY$trunc$right[idxi, , drop = FALSE],
                                 dat$eY$trunc$right[idxe, , drop = FALSE]), which = c(dat$iY$which[idxi],
                                                                                      dat$eY$which[idxe]))
  return(out)
}

.rm_int <- function(x) {
  if (all(x[, 1] == 1))
    return(x[, -1L, drop = FALSE])
  return(x)
}
