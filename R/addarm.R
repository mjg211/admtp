addarm <- function(ratio = 1, narms = 2, sd = 10, diff = 3, power = 0.9,
                   nstage1 = 100, alpha = 0.05) {

  n          <- ((ratio + 1)/ratio*(sd)^2*(stats::qnorm(1 - alpha) +
                                             stats::qnorm(power))^2) / (diff^2)
  message("For a single two-arm parallel group trial we need ", ceiling(n),
          " on the experimental arm\n")
  message("We need ", ratio*ceiling(n), " on the control arm\n")
  message("Assuming a one-sided test, CRD = ", diff, ", SD = ", sd,
          ", power = ", power, " and allocation ratio of ", ratio, "\n")
  if (nstage1 > ratio*n) {
    warning("Overlap bigger than required sample size for one arm\n")
  }
  newcorr    <- 0
  corr       <- (ratio*n/(ratio + 1))*(ratio*n - nstage1)/(ratio*n)^2
  while (abs(newcorr - corr) > 0.0001) {
     sig     <- rbind(c(1, corr),
                      c(corr, 1))
     mid     <- mvtnorm::qmvnorm(1 - alpha, sigma = sig)$quantile
     n       <-
       ((ratio + 1)/(ratio)*(sd)^2*(mid + stats::qnorm(power))^2)/(diff)^2
     newcorr <- (ratio*n/(ratio + 1))*(ratio*n - nstage1)/(ratio*n)^2
     corr    <- newcorr
  }
  finalen    <- ceiling(n)
  message("Sample size per experimental arm: ", finalen, "\n")
  finalcn    <- ratio*ceiling(n)
  message("Sample size per control arm: ", finalcn, "\n")
  message("Therefore the total sample size is ", ceiling(finalen),
          " per exp. arm, plus ", finalcn, " on the concurrent control arm for",
          " treatment 2, plus ", nstage1, " controls\n")
  message("The full sample size is: ", 2*finalen + finalcn + nstage1)
}
