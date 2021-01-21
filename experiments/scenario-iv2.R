# Simulation anchor ontrams
# Lucas Kook
# 29.10.2020

set.seed(1)

# Dependencies ------------------------------------------------------------

library(anchor)
library(tram)

# Helpers -----------------------------------------------------------------

warm_start <- function(m, mt, ncl) {
  tmp <- get_weights(m$mod_baseline)
  tmp[[1]][] <- ontram:::.to_gamma(coef(mt, with_baseline = TRUE)[-(ncl)])
  set_weights(m$mod_baseline, tmp)
  tmp2 <- get_weights(m$mod_shift)
  tmp2[[1]][] <- coef(mt)
  set_weights(m$mod_shift, tmp2)
  return(invisible(m))
}

mplt <- function(q_plain, q_anchor, alps) {
  matplot(alps, q_plain, col = "orange", pch = 20, las = 1, xlab = expression(alpha),
       ylab = expression(alpha~"quantile of logLik"),
       ylim = c(min(c(q_anchor, q_plain)), max(c(q_anchor, q_plain))))
  matpoints(alps, q_anchor, col = "cornflowerblue", pch = 20)
  legend("topleft", pch = 20, col = c("cornflowerblue", "orange"),
         legend = c("Anchor", "Plain"), bty = "n")
}

make_file_name <- function(bn, parnames, parvals, ext, include_date = FALSE) {
  fn <- paste0(bn, "_", paste0(parnames, parvals, collapse = "-"))
  dte <- Sys.Date()
  if (include_date)
    fn <- paste0(fn, "_", dte)
  ret <- paste0(fn, ".", ext)
  return(ret)
}

my_write <- function(obj, fn) {
  write.csv(x = obj, file = fn, quote = FALSE, row.names = FALSE)
}

# Params ------------------------------------------------------------------

setup <- data.frame(
  expand.grid(
    n = 1e3,
    ncl = 10,
    bx = 0.5,
    ba = c(-1, 0.5, 1, 2),
    bh = 1.5,
    sdxyh = 0.75^2,
    shift = c(0, 1, 1.8, 3),
    xi = c(0, 1, 10, 100, 1000, 10000)
  )
)

bn <- "oann-sim"
ext <- "csv"

# Simulation --------------------------------------------------------------
  
nsim <- 200
alps <- (0:10)/10

for (scenario in 1:nrow(setup)) {
  scen <- setup[scenario, ]
  cat(round(do.call("c", scen), 2))
  ll_anchor <- ll_plain <- matrix(NA, nrow = scen$n, ncol = nsim)
  cfx_anchor <- cfx_plain <- vector(mode = "numeric", length = nsim)
  fn <- make_file_name(bn = bn, parnames = colnames(scen), parvals = scen[1, ],
                       ext = ext, include_date = FALSE)
  
  for (run in 1:nsim) {
    # Generate data -----------------------------------------------------------
    
    tr <- gen_dat_ordinal(n = scen$n, nclasses = scen$ncl,
                          bx = scen$bx, ba = scen$ba, bh = scen$bh, 
                          sd = scen$sdxyh)
    te <- gen_dat_ordinal(n = 2 * scen$n, nclasses = scen$ncl, bx = scen$bx, 
                          ba = scen$ba, bh = scen$bh, sd = scen$sdxyh,
                          perturb = TRUE, shift = scen$shift)
    dat <- data.frame(y = ordered(apply(tr$Y, 1, which.max)), x = c(tr$X))
    mt <- Polr(y ~ x, data = dat)
    
    # Fit anchor ontram -------------------------------------------------------
    
    m <- ontram_polr(x_dim = ncol(tr$X), y_dim = ncol(tr$Y),
                     method = "logit", n_batches = 4, epochs = nep <- 200)
    warm_start(m, mt, scen$ncl)
    mh <- fit_ontram_a(m, x_train = tr$X, y_train = tr$Y, history = histo <- FALSE,
                       x_test = te$X, y_test = te$Y, a_train = tr$A, xi = scen$xi)
    mc <- ontram_polr(x_dim = ncol(tr$X), y_dim = ncol(tr$Y),
                      method = "logit", n_batches = 4, epochs = nep)
    warm_start(mc, mt, scen$ncl)
    mch <- fit_ontram2(mc, x_train = tr$X, y_train = tr$Y, history = histo,
                       x_test = te$X, y_test = te$Y)
    
    # Write out ---------------------------------------------------------------
    
    ll_anchor[, run] <- evaluate(m, te$X, te$Y)
    ll_plain[, run] <- evaluate(mc, te$X, te$Y)
    
    cfx_anchor[run] <- coef(m)
    cfx_plain[run] <- coef(mc)
  }
  
  my_write(ll_anchor, paste0("ll_anchor", "_", fn))
  my_write(ll_plain, paste0("ll_plain", "_", fn))
  my_write(cfx_anchor, paste0("cfx_anchor", "_", fn))
  my_write(cfx_plain, paste0("cfx_plain", "_", fn))
  
  qq_anchor <- apply(ll_anchor, 2, quantile, probs = alps, na.rm = TRUE)
  qq_plain <- apply(ll_plain, 2, quantile, probs = alps, na.rm = TRUE)
}
