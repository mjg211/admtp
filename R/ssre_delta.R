ssre_delta <- function(z1, n1 = 50, f1 = 1, e1 = 2.76, alpha = 0.025,
                       beta = 0.2, Nmax = 1e10) {
  nmax            <- Nmax - n1
  e2              <- uniroot(ssre_delta_find_e2, interval = c(0.2, 4),
                             f1 = f1, e1 = e1, alpha = alpha, beta = beta,
                             n1 = n1, nmax = nmax)$root
  n2              <- numeric(length(z1))
  for (i in 1:length(n2)) {
    if (any(z1[i] <= f1, z1[i] >= e1)) {
      n2[i]       <- 0
    } else {
      if (z1[i] > (e2 + qnorm(1 - beta))*sqrt(n1/(n1 + nmax))) {
        Z_beta_z1 <- qnorm(1 - beta)
      } else if (z1[i] < (e2 + qnorm(1 - beta))*sqrt(n1/(n1 + nmax))) {
        Z_beta_z1 <- z1[i]*sqrt((nmax + n1)/n1) - e2
      }
      n2[i]       <- ((e2 + Z_beta_z1)^2/(z1[i]^2) - 1)*n1
    }
  }
  n2
}