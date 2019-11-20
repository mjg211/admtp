ssre_delta_density <- function(u, e2, beta, n1, nmax) {
  Z_beta_u <- sapply(sqrt((n1 + nmax)/n1)*u - e2,
                     function(x) min(qnorm(1 - beta), x))
  A1       <- e2*(e2 + Z_beta_u) - u^2
  A2       <- sqrt((e2 + Z_beta_u)^2 - u^2)
  pnorm(A1/A2)*dnorm(u)
}