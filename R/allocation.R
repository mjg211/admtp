# allocation probabilities takes the posteriors found in findposterior and
# converts them to allocation probabilities according to the formulae in the
# lecture
# arguments:
# posterior_exp_better_than_control - output from find_posterior
# gamma - weighting factor (should be greater than 0)
# control_allocation - the desired allocation to control (should be between 0
# and 1)

allocation <- function(posterior_exp_better_than_control, gamma = 0.5,
                       control_allocation =
                         1/(length(posterior_exp_better_than_control) + 1)) {

  ##### Check inputs ###########################################################



  ##### Perform main computations ##############################################

  # Control allocation is set to controlallocation, the remaining arms have
  # allocation set to the posterior raised to the power of gamma

  allocation_exp <- (posterior_exp_better_than_control^gamma)/
    sum(posterior_exp_better_than_control^gamma)
  # Change so that sum of experimental allocation + control allocation = 1
  allocation_exp <- (1 - control_allocation)*allocation_exp

  ##### Outputting #############################################################

  c(control_allocation, allocation_exp)

}