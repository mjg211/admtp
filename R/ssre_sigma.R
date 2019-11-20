ssre_sigma <- function(n1 = 50, alpha = 0.025, beta = 0.2, delta = 3.5,
                       tau = 0, sigma = 10, estimator = "one_sample", nB = 2,
                       inflation = F, Nmin = 1, Nmax = 1e10,
                       replicates = 10000) {
  nmin             <- Nmin
  nmax             <- nmax
  sigma_hat        <- reject <- numeric(replicates)
  two_n1           <- 2*n1
  N                <- rep(two_n1, replicates)
  if (estimator == "one_sample") {
    df             <- two_n1 - 1
  }
  else if (estimator == "block") {
    B              <- two_n1/nB
    df             <- two_n1 - nB
    TB             <- numeric(B)
    half_nB        <- 0.5*nB
  }
  else if (estimator == "unblinded") {
    df             <- two_n1 - 2
  }
  sqrt_inv_df      <- sqrt(1/df)
  n_factor         <- 4*((qnorm(1 - alpha) + qnorm(1 - beta))/delta)^2
  if (!inflation) {
    inflation_fac  <- 1
  }
  else {
    inflation_fac  <- (qt(1 - alpha, df) + qt(1 - beta, df))^2/
      (qnorm(1 - alpha) + qnorm(1 - beta))^2
  }
  for (i in 1:replicates) {
    x01            <- rnorm(n1, sd = sigma)
    x11            <- rnorm(n1, tau, sigma)
    if (estimator == "one_sample") {
      sigma_hat[i] <- sqrt_inv_df*sqrt(sum((c(x01, x11) -
                                              mean(c(x01, x11)))^2))
    }
    else if (estimator == "block") {
      for (b in 1:B) {
        range      <- (1 + half_nB*(b - 1)):(half_nB*b)
        TB[b]      <- sum(c(x01[range], x11[range]))
      }
      sigma_hat[i] <- sqrt_inv_df*sqrt(sum((TB - mean(TB))^2))
    }
    else if (estimator == "unblinded") {
      sigma_hat[i] <- sqrt_inv_df*sqrt(sum((x01 - mean(x01))^2) +
                                         sum((x11 - mean(x11))^2))
    }
    n              <- ceiling(inflation_fac*n_factor*sigma_hat[i]^2)
    if (estimator != "block") {
      if (n < nmin) {
        n2         <- ceiling(0.5*(nmin - two_n1))
      }
      else {
        if (nmax != Inf) {
          if (n > nmax) {
            n2     <- ceiling(0.5*(nmax - two_n1))
          }
          else {
            n2     <- ceiling(0.5*(n - two_n1))
          }
        }
        else {
          n2       <- ceiling(0.5*(n - two_n1))
        }
      }
    }
    else {
      while (n%%nB != 0) {
        n          <- n + 1
      }
      if (n < nmin) {
        n2         <- ceiling(0.5*(nmin - two_n1))
        while ((2*n2)%%nB != 0) {
          n2       <- n2 + 1
        }
      }
      else {
        if (nmax != Inf) {
          if (n > nmax) {
            n2     <- floor(0.5*(nmax - two_n1))
          }
          else {
            n2     <- ceiling(0.5*(n - two_n1))
          }
          while ((2*n2)%%nB != 0) {
            n2     <- n2 + 1
          }
        }
        else {
          n2       <- ceiling(0.5*(n - two_n1))
          while (n2%%nB != 0) {
            n2     <- n2 + 1
          }
        }
      }
    }
    if (n2 > 0) {
      N[i]         <- N[i] + 2*n2
      x0           <- c(x01, rnorm(n2, sd = sigma))
      x1           <- c(x11, rnorm(n2, tau, sigma))
      S02          <- var(x0)
      S12          <- var(x1)
      Spool        <- sqrt((0.5*N[i] - 1)*(S02 + S12)/(N[i] - 1))
      reject[i]    <- (sqrt(1/N[i])*(sum(x1) - sum(x0))/Spool >=
                         qt(1 - alpha, N[i] - 1))
    }
    else {
      S02          <- var(x01)
      S12          <- var(x11)
      Spool        <- sqrt((n1 - 1)*(S02 + S12)/(N[i] - 1))
      reject[i]    <- (sqrt(1/N[i])*(sum(x11) - sum(x01))/Spool >=
                         qt(1 - alpha, N[i] - 1))
    }
  }
  list(empirical_power = mean(reject), average_N_hat = mean(N),
       average_sigma_hat = mean(sigma_hat), sigma_hat = sigma_hat, N_hat = N)
}