# Code to find posterior distribution for multi-arm trial with binary outcome

# findposterior finds the posterior probability of each experimental arm being
# better than control:
# Arguments:
# K - number of experimental treatments
# nperarm - sample size per arm at the interim
# responses.control - number of responses observed out of controls
# responses.exp - vector of length K with number of responses observed out of
# experimental treatments
# prior.a - value of a parameter in prior Beta(a,b) distribution
# prior.b - value of b parameter in prior Beta(a,b) distribution

findposterior <- function(K = 3, n_per_arm = rep(20, K + 1),
                          responses_control = 5, responses_exp = c(3, 10, 5),
                          prior_a = 0.25, prior_b = 0.75) {

  ##### Check inputs ###########################################################



  ##### Perform main computations ##############################################

  # Posterior distribution is Beta(prior_a + responses, prior_b + n - responses)
  # Find posterior for each experimental
  posterior_exp_a                   <- prior_a + responses_exp
  posterior_exp_b                   <- prior_b + n_per_arm[-1] - responses_exp

  # Find posterior for control
  posterior_control_a               <- prior_a + responses_control
  posterior_control_b               <- prior_b + n_per_arm[1] -
    responses_control

  # Find posterior probability that each experimental is better than control by
  # calling differenceinbeta function
  posterior_exp_better_than_control <-
    sapply(1:K, function(x) { difference_in_beta(posterior_control_a,
                                                 posterior_control_b,
                                                 posterior_exp_a[x],
                                                 posterior_exp_b[x]) })

  ##### Outputting #############################################################

  return(posterior_exp_better_than_control)

}