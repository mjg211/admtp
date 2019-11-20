# Support functions to find probability of one beta distribution being greater
# than another P(X > Y) where X ~ Beta(a2, b2), Y ~ Beta(a1, b1)
differenceinbeta_givenx <- function(x, a1, b1, a2, b2) {

  dbeta(x, a1, b1)*(1 - pbeta(x, a2, b2))

}