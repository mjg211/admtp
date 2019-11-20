ssre_delta_find_e2 <- function(e2, f1, e1, alpha, beta, n1, nmax) {
  integral <- integrate(ssre_delta_density,
                        lower = f1,
                        upper = min(e1, e2 + qnorm(1 - beta)), e2 = e2,
                        beta = beta, n1 = n1, nmax = nmax)$value
  integral - (1 - pnorm(f1) - alpha)
}