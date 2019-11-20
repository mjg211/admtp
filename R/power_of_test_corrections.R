#' Calculate power of multiple testing procedures
#'
#' \code{power_of_test_corrections} uses simulation to provide the expected
#' number of correct and incorrect rejections for each of three multiple testing
#' procedures: Benjamini-Hochberg, Bonferroni, and Holm.
#'
#' Since it uses simulation, it may take a few seconds to run. If it takes too
#' long, then change \code{replicates} to a smaller number.
#'
#' As in \code{\link{error_rates}}, \code{delta} is the vector of test statistic
#' means. However, in this case the correlation parameter (see
#' \code{correlation}) is only between the first and second test statistic
#' (to allow us to investigate negative values).
#'
#' @param delta A \code{\link{numeric}} \code{\link{vector}} giving the
#' standardised means of the test statistics. It is assumed that the null
#' hypotheses are that each component of \code{delta} is less than or equal to 0
#' (so zero or negative \code{delta} values correspond to true nulls and
#' positive delta values correspond to false nulls). Defaults to
#' \code{numeric(5)}.
#' @param correlation A \code{\link{numeric}} giving the correlation between
#' the first and second test statistics. Defaults to \code{0}.
#' @param target_alpha A \code{\link{numeric}} giving the target familywise
#' error-rate. It must belong to [-1,1]. Defaults to \code{0.05}.
#' @param replicates A \code{\link{numeric}} giving the number of replicate
#' simulations to use. It must belong to {2,3,...}. Defaults to \code{1e4}.
#' @return A \code{\link{list}} containing the following elements:
#' \itemize{
#' \item A \code{\link{numeric}} \code{\link{vector}} in the slot
#' \code{$average_correct_rejections} giving the average number of correct
#' rejections for each of the three supported methods.
#' \item A \code{\link{numeric}} \code{\link{vector}} in the slot
#' \code{$average_incorrect_rejections} giving the average number of incorrect
#' rejections for each of the three supported methods.
#' }
#' @seealso \code{\link{error_rates}}, \code{\link{test_corrections}}.
#' @examples
#' # For the default parameters
#' default <- power_of_test_corrections()
#' @export
power_of_test_corrections <- function(delta = numeric(5), correlation = 0,
                                      target_alpha = 0.05, replicates = 1e4) {

  ##### Check inputs ###########################################################

  if (any(!is.numeric(delta), is.infinite(delta))) {
    stop("delta must a numeric vector, whose elements belong to (-Inf,Inf)")
  }
  if (any(!is.numeric(correlation), length(correlation) != 1, correlation < -1,
          correlation > 1)) {
    stop("correlation must be a single numeric that belongs to [-1,1]")
  }
  if (any(!is.numeric(target_alpha), length(target_alpha) != 1,
          target_alpha < 0, target_alpha > 1)) {
    stop("target_alpha must be a single numeric that belongs to (0,1)")
  }
  if (any(!is.numeric(replicates), length(replicates) != 1,
          replicates%%1 != 0, replicates < 2)) {
    stop("replicates must be a single numeric that belongs to {2,3,...}")
  }

  ##### Perform main computations ##############################################

  K                            <- length(delta)
  sigma                        <- diag(1, K)
  sigma[1, 2]                  <- sigma[2, 1] <- correlation
  test_statistics              <- mvtnorm::rmvnorm(replicates, delta, sigma)
  pvalues                      <- 1 - stats::pnorm(test_statistics)
  benjamini_hochberg           <- bonferroni <- holm <- matrix(F, replicates, K)
  for (i in 1:replicates) {
    test_i                     <- test_corrections(pvalues[i, ], target_alpha)
    benjamini_hochberg[i, ]    <- test_i$benjamini_hochberg
    bonferroni[i, ]            <- test_i$bonferroni
    holm[i, ]                  <- test_i$holm
  }
  delta_neg                    <- (delta <= 0)
  delta_pos                    <- (delta > 0)
  benjamini_hochberg_correct   <- as.double(benjamini_hochberg%*%delta_pos)
  bonferroni_correct           <- as.double(bonferroni%*%delta_pos)
  holm_correct                 <- as.double(holm%*%delta_pos)
  benjamini_hochberg_incorrect <- as.double(benjamini_hochberg%*%delta_neg)
  bonferroni_incorrect         <- as.double(bonferroni%*%delta_neg)
  holm_incorrect               <- as.double(holm%*%delta_neg)

  ##### Outputting #############################################################

  average_correct_rejections        <- c(mean(benjamini_hochberg_correct),
                                         mean(bonferroni_correct),
                                         mean(holm_correct))
  average_incorrect_rejections      <- c(mean(benjamini_hochberg_incorrect),
                                         mean(bonferroni_incorrect),
                                         mean(holm_incorrect))
  names(average_correct_rejections) <- names(average_incorrect_rejections) <-
    c("Benjamini-Hochberg", "Bonferroni", "Holm")
  list(average_correct_rejections   = average_correct_rejections,
       average_incorrect_rejections = average_incorrect_rejections)

}