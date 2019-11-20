#' Calculate sample size for a fixed-sample two-arm randomised controlled trial
#'
#' \code{des_fixed} calculates the sample size required by a fixed-sample
#' (single-stage) two-arm randomised controlled trial that tests for
#' superiority, assuming a normally distributed outcome variable. It supports
#' unequal allocation to the two arms along with one- and two-sided null
#' hypotheses.
#'
#' @param alpha A \code{\link{numeric}} giving the desired type-I error rate. It
#' must belong to (0,1). Defaults to \code{0.05}.
#' @param beta A \code{\link{numeric}} giving the desired type-II error rate. It
#' must belong to (0,1). Defaults to \code{0.2}.
#' @param delta A \code{\link{numeric}} giving the treatment effect to power
#' for. It must belong to (0,Inf). Defaults to \code{0.5}.
#' @param sigma A \code{\link{numeric}} giving the assumed value of the standard
#' deviation of the responses. It must belong to (0,Inf). Defaults to \code{1}.
#' @param ratio A \code{\link{numeric}} giving the allocation ratio to the
#' experimental arm relative to the control arm. It must belong to (0,Inf).
#' Defaults to \code{1}.
#' @param two_sided A \code{\link{logical}} variable specifying whether a one-
#' or two-sided null hypothesis should be assumed. Defaults to \code{F}.
#' @return A \code{\link{list}} containing the following elements:
#' \itemize{
#' \item A \code{\link{numeric}} in the slot \code{$n} giving the sample size
#' required in the control arm.
#' \item A \code{\link{numeric}} in the slot \code{$N} giving the total required
#' sample size.
#' \item A \code{\link{numeric}} in the slot \code{$P_HG} giving the probability
#' of rejecting the null hypothesis when there is no treatment effect.
#' \item A \code{\link{numeric}} in the slot \code{$P_LFC} giving the
#' probability of rejecting the null hypothesis when the treatment effect is
#' \code{delta}.
#' \item A \code{\link{numeric}} in the slot \code{$u} giving the critical
#' rejection threshold for the test statistic.
#' }
#' @examples
#' # For the default parameters
#' default <- des_fixed()
#' @export
des_fixed <- function(alpha = 0.05, beta = 0.2, delta = 0.5, sigma = 1,
                      ratio = 1, two_sided = F) {

  ##### Check inputs ###########################################################

  if (any(!is.numeric(alpha), length(alpha) != 1, alpha <= 0, alpha >= 1)) {
    stop("alpha must be a single numeric that belongs to (0, 1)")
  }
  if (any(!is.numeric(beta), length(beta) != 1, beta <= 0, beta >= 1)) {
    stop("beta must be a single numeric that belongs to (0, 1)")
  }
  if (any(!is.numeric(delta), length(delta) != 1, delta <= 0,
          is.infinite(delta))) {
    stop("delta must be a single numeric that belongs to (0, Inf)")
  }
  if (any(!is.numeric(sigma), length(sigma) != 1, sigma <= 0,
          is.infinite(sigma))) {
    stop("sigma must be a single numeric that belongs to (0, Inf)")
  }
  if (any(!is.numeric(ratio), length(ratio) != 1, ratio <= 0,
          is.infinite(ratio))) {
    stop("ratio must be a single numeric that belongs to (0, Inf)")
  }
  if (!is.logical(two_sided)) {
    stop("two_sided must be logical")
  }

  ##### Perform main computations ##############################################

  u       <- stats::qnorm(1 - alpha/ifelse(two_sided, 2, 1))
  n       <- ceiling((((u + stats::qnorm(1 - beta))*sigma*sqrt(1 + 1/ratio))/
                        delta)^2)
  while ((ratio*n)%%1 != 0) {
    n     <- n + 1
  }
  P_HG    <- alpha
  if (two_sided) {
    P_LFC <- stats::pnorm(-u, delta*sqrt(n/(1 + 1/ratio))/sigma) +
      stats::pnorm(u, delta*sqrt(n/(1 + 1/ratio))/sigma, lower.tail = F)
  } else {
    P_LFC <- stats::pnorm(u, delta*sqrt(n/(1 + 1/ratio))/sigma, lower.tail = F)
  }

  ##### Outputting #############################################################

  list(n     = n,
       N     = n*(1 + ratio),
       P_HG  = P_HG,
       P_LFC = P_LFC,
       u     = u)

}