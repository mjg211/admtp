#' Error rates for normally distributed test statistics
#'
#' \code{error_rates} finds the maximum per-hypothesis type I error-rate (PHER),
#' family-wise error-rate (FWER), and false-discovery rate (FDR) for normally
#' distributed test statistics with mean \code{delta}, variance 1, and a
#' specified level of correlation between each pair of test statistics (see
#' \code{correlation}).
#'
#' @param delta A \code{\link{numeric}} \code{\link{vector}} giving the
#' standardised means of the test statistics. It is assumed that the null
#' hypotheses are that each component of \code{delta} is less than or equal to 0
#' (so zero or negative \code{delta} values correspond to true nulls and
#' positive delta values correspond to false nulls). Its elements must be
#' finite. Defaults to \code{numeric(5)}.
#' @param correlation A \code{\link{numeric}} giving the correlation between the
#' test statistics. It must belong to [0,1]. Defaults to \code{0}.
#' @param critical_value A \code{\link{numeric}} giving the critical rejection
#' threshold such that when a test statistic takes a value above
#' \code{critical_value}, its corresponding null hypothesis is rejected.
#' Defaults to \code{1.96}.
#' @return A \code{\link{list}} containing the following elements:
#' \itemize{
#' \item A \code{\link{numeric}} in the slot \code{$fdr} giving the estimated
#' FDR.
#' \item A \code{\link{numeric}} in the slot \code{$fwer} giving the FWER.
#' \item A \code{\link{numeric}} in the slot \code{$max_pher} giving the maximum
#' PHER.
#' }
#' @seealso \code{\link{power_of_test_corrections}},
#' \code{\link{test_corrections}}.
#' @examples
#' # For the default parameters
#' default <- error_rates()
#' @export
error_rates <- function(delta = numeric(5), correlation = 0,
                        critical_value = 1.96) {

  ##### Check inputs ###########################################################

  if (any(!is.numeric(delta), is.infinite(delta))) {
    stop("delta must a numeric vector, whose elements belong to (-Inf,Inf)")
  }
  if (any(!is.numeric(correlation), length(correlation) != 1, correlation < 0,
          correlation > 1)) {
    stop("correlation must be a single numeric that belongs to [0,1]")
  }
  if (any(!is.numeric(critical_value), length(critical_value) != 1)) {
    stop("critical_value must be a single numeric")
  }

  ##### Perform main computations ##############################################

  K                 <- length(delta)
  cov               <- matrix(correlation, K, K) + diag(1 - correlation, K)
  pher              <- (1 - stats::pnorm(critical_value,
                                         delta, 1))*as.integer(delta <= 0)
  fwer              <- 0
  which_null        <- (delta <= 0)
  if (sum(which_null) > 0) {
    fwer            <- 1 - mvtnorm::pmvnorm(upper = rep(critical_value,
                                                        sum(which_null)),
                                            mean  = delta[which_null],
                                            sigma = cov[which_null,
                                                        which_null])[1]
  }
  test_statistics   <- mvtnorm::rmvnorm(100000, delta, cov)
  false_discoveries <- apply(test_statistics, 1,
                             function(x) { length(which(x > critical_value &
                                                          which_null)) })
  discoveries       <- apply(test_statistics, 1,
                             function(x) { length(which(x > critical_value)) })
  fdr               <- mean(replace(false_discoveries/discoveries,
                                    discoveries == 0, 0))

  ##### Outputting #############################################################

  list(fdr      = fdr,
       fwer     = fwer,
       max_pher = max(pher))

}