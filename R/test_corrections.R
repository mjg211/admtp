#' Apply multiple testing procedures to p-values
#'
#' \code{test_corrections} allows you to apply Benjamini-Hochberg, Bonferroni,
#' and Holm, to a set of p-values. For each method, it returns a vector where
#' elements take value \code{\link[base]{TRUE}} if the corresponding null
#' hypothesis is rejected, and \code{\link[base]{FALSE}} otherwise.
#'
#' Note that \code{\link[stats]{p.adjust}} implements more corrections.
#'
#' @param pvalues A \code{\link{numeric}} \code{\link{vector}} giving the
#' p-values for each hypothesis. Its elements must all belong to (0,1].
#' @param target_alpha A \code{\link{numeric}} giving the target familywise
#' error-rate. It must belong to (0,1). Defaults to \code{0.05}.
#' @return A \code{\link{list}} containing the following elements:
#' \itemize{
#' \item A \code{\link{numeric}} \code{\link{vector}} in the slot
#' \code{$benjamini_hochberg} giving the rejection results when applying
#' Benjamini-Hochberg.
#' \item A \code{\link{numeric}} \code{\link{vector}} in the slot
#' \code{$bonferroni} giving the rejection results when applying
#' Bonferroni.
#' \item A \code{\link{numeric}} \code{\link{vector}} in the slot
#' \code{$holm} giving the rejection results when applying Holm.
#' }
#' @seealso \code{\link{error_rates}}, \code{\link{power_of_test_corrections}}.
#' @examples
#' # A basic example with three hypotheses
#' three <- test_corrections(c(0.01, 0.02, 0.03))
#' @export
test_corrections <- function(pvalues, target_alpha = 0.05) {

  ##### Check inputs ###########################################################

  if (any(!is.numeric(pvalues), pvalues <= 0, pvalues > 1)) {
    stop("pvalues must a numeric vector, whose elements all belong to (0,1]")
  }
  if (any(!is.numeric(target_alpha), length(target_alpha) != 1,
          target_alpha < 0, target_alpha > 1)) {
    stop("target_alpha must be a single numeric that belongs to (0,1)")
  }

  ##### Perform main computations ##############################################

  K                                      <- length(pvalues)
  bonferroni                             <- (pvalues <= target_alpha/K)
  pvalues_ordered                        <- pvalues[order(pvalues)]
  holm_rejections_ordered                <- (pvalues_ordered <=
                                               target_alpha/(K + 1 - (1:K)))
  min_entry                              <-
    suppressWarnings(min(which(holm_rejections_ordered == F)))
  if (!is.infinite(min_entry)) {
    holm_rejections_ordered[min_entry:K] <- F
  }
  holm                                   <- logical(K)
  holm[order(pvalues)]                   <- holm_rejections_ordered
  bh_rejections_ordered                  <- (pvalues_ordered <=
                                               target_alpha*(1:K)/K)
  max_entry                              <-
    suppressWarnings(max(which(bh_rejections_ordered == T)))
  if(!is.infinite(max_entry)) {
    bh_rejections_ordered[1:max_entry]   <- T
  }
  benjamini_hochberg                     <- logical(K)
  benjamini_hochberg[order(pvalues)]     <- bh_rejections_ordered

  ##### Outputting #############################################################

  list(benjamini_hochberg = benjamini_hochberg,
       bonferroni         = bonferroni,
       holm               = holm)

}