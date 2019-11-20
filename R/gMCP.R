#' ...
#'
#' \code{gMCP}...
#'
#' @param graph .
#' @param pvalues .
#' @param alpha Defaults to \code{0.05}.
#' @param fweights Defaults to \code{F}.
#' @return ...
#' @seealso \code{\link{calcPower}}, ...
#' @examples
#' # A basic example
#' basic <- gMCP(...)
#' @export
gMCP <- function(graph, pvalues, alpha = 0.05, fweights = F) {

  ##### Check inputs ###########################################################



  ##### Perform main computations ##############################################

  G                     <- graph$G
  n                     <- ncol(G)
  h                     <- numeric(n)
  names(h)              <- paste0("H", 1:n)
  a                     <- alpha*graph$weights
  crit                  <- 0
  while (crit == 0) {
    test                <- (pvalues <= a)
    if (any(test)) {
      rej               <- which.max(test)
      h[rej]            <- 1
      Gtemp             <- matrix(0, n, n)
      for (i in 1:n) {
        a[i]            <- a[i] + a[rej]*G[rej, i]
        if (G[i, rej]*G[rej, i] < 1) {
          for (j in 1:n) {
            Gtemp[i, j] <- (G[i, j] + G[i, rej]*G[rej, j])/
                             (1 - G[i, rej]*G[rej, i])
          }
        }
        Gtemp[i, i]     <- 0
      }
      G                 <- Gtemp
      G[rej, ]          <- G[, rej] <- 0
      a[rej]            <- 0
    } else {
      crit              <- 1
    }
  }

  ##### Outputting #############################################################

  if (fweights) {
    list(pvalues = pvalues, alpha = alpha, rejected = h, weights = a/alpha)
  } else {
    list(pvalues = pvalues, alpha = alpha, rejected = h)
  }

}