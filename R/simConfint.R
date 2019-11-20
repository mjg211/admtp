#' ...
#'
#' \code{simConfInt}...
#'
#' @param graph .
#' @param pvalues .
#' @param confint .
#' @param alternative Defaults to \code{c("less", "greater")}.
#' @param estimates .
#' @param df .
#' @param alpha Defaults to \code{0.05}.
#' @param mu Defaults to \code{0}.
#' @return ...
#' @seealso ...
#' @examples
#' # A basic example
#' basic <- simConfint()
#' @export

simConfint <- function(graph, pvalues, confint, alternative=c("less", "greater"),
                      estimates, df, alpha=0.05, mu=0) {

  ##### Check input variables ##################################################



  ##### Perform main computations ##############################################

  result      <- gMCP(graph, pvalues, alpha=alpha, fweights = TRUE)
  if (all(result$rejected==1)) {
    alpha     <- graph$weights*alpha
  } else {
    alpha     <- result$weights*alpha
  }

  if (confint == "t") {
    dist      <- function(x) { qt(p = x, df = df) }
  } else if (confint=="normal") {
    dist      <- qnorm
  } else {
    stop('Parameter confint has to be a function or "t" or "normal"')
  }
  if (alternative == "greater") {
    stderr    <- abs(estimates/dist(1 - pvalues))
    lb        <- estimates + dist(alpha)*stderr
    lb        <- ifelse(result$rejected, max(0, lb), lb)
    ub        <- rep(Inf, length(lb))
  } else if (alternative == "less") {
    stderr    <- abs(estimates/dist(pvalues))
    ub        <- estimates + dist(1 - alpha)*stderr
    ub        <- ifelse(result$rejected, min(0, ub), ub)
    lb        <- rep(-Inf, length(ub))
  } else {
    stop("Specify alternative as \"less\" or \"greater\".")
  }
  m           <- matrix(c(lb, estimates, ub), ncol=3)
  colnames(m) <- c("lower bound", "estimate", "upper bound")
  return(m)

}