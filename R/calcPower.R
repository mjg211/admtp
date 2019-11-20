#' ...
#'
#' \code{calcPower}...
#'
#' @param graph .
#' @param mean .
#' @param f Defaults to \code{NULL}.
#' @param corr_sim .
#' @param alpha Defaults to \code{0.025}.
#' @param n_sim Defaults to \code{1e4}.
#' @return ...
#' @seealso \code{\link{g_mtp}}, ...
#' @examples
#' # A basic example
#' basic <- calcPower(...)
#' @export
calcPower <- function(graph, mean, f = NULL, corr_sim, alpha = 0.025,
                      n_sim = 1e4) {

  ##### Check inputs ###########################################################



  ##### Perform main computations ##############################################

  sims           <- mvtnorm::rmvnorm(n_sim, mean, corr_sim)
  pvals          <- stats::pnorm(sims, lower.tail = F)
  out            <- matrix(0, n_sim, ncol(graph$G))
  for (i in 1:n_sim) {
    out[i, ]     <- gMCP(graph, pvalue = pvals[i,], alpha)$rejected
  }
  pow            <- Rfast::colmeans(out)
  avg_pow        <- sum(out)/nrow(out)
  at_least_1     <- mean(Rfast::rowsums(out) > 0)
  all_pow        <- mean(Rfast::rowsums(out) == ncol(out))
  result         <- list(exp_rejections = avg_pow,
                         local_power    = pow,
                         pow_at_least_1 = at_least_1,
                         reject_all     = all_pow)
  if(is.list(f)) {
    fn           <- names(f)
    result[[fn]] <- sum(apply(out, 1, f[[fn]]))/nrow(out)
  }

  ##### Outputting #############################################################

  result

}