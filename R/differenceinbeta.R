differenceinbeta <- function(a1, b1, a2, b2) {
  # integrates over sample space of first beta distribution and qbeta of the
  # second distribution:
  stats::integrate(difference_in_beta_given_x,
                   lower = stats::qbeta(1e-10, a1, b1),
                   upper = stats::qbeta(1 - 1e-10, a1, b1),
                   a1    = a1,
                   b1    = b1,
                   a2    = a2,
                   b2    = b2)$value

}