#' ...
#'
#' \code{mtp_power}...
#'
#' @param mu .
#' @param sigma .
#' @param alpha Defaults to \code{0.05}.
#' @param correction Defaults to \code{"none"}.
#' @param replicates Defaults to \code{1e5}.
#' @return ...
#' @seealso ...
#' @examples
#' # A basic example
#' basic <- mtp_power()
#' @export
mtp_power <- function(mu, sigma, alpha = 0.05, correction = "none",
                      replicates = 1e5) {

  ##### Check inputs ###########################################################



  ##### Perform main computations ##############################################

  K                <- length(mu)
  trial_sim        <- stats::pnorm(-mvtnorm::rmvnorm(replicates, mu, sigma))
  if (correction == "none") {
    return(Rfast::rowsums(trial_sim <= alpha))
  } else if (correction == "bonferroni") {
    return(Rfast::rowsums(trial_sim <= alpha/K))
  } else if (correction == "sidak") {
    return(Rfast::rowsums(trial_sim <= 1 - (1-alpha)^(1/K)))
  } else if (correction == "dunnett") {
    return(
      Rfast::rowsums(
        trial_sim <= stats::pnorm(mvtnorm::qmvnorm(1 - alpha,
                                                   sigma = sigma)$quantile,
                                  lower.tail = F)
      )
    )
  } else if (correction == "holm") {
    trial_sim_sort <- t(apply(trial_sim, 1, sort))
    return(
      apply(trial_sim_sort, 1,
            function(x) { min(c(which(x > alpha/(K:1)), K + 1)) - 1 })
    )
  }
}