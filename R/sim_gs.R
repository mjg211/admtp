sim_gs <- function(l, u, n = 72, tau = 0, sigma = 10, ratio = 1,
                   sigma_z = sigma, t_test = F, adjusted = F, alpha = 0.05,
                   replicates = 100000) {

  ##### Check input variables ##################################################

  if (length(l) != length(u)) {
    stop("The lower and upper stopping boundaries (l and u) must be vectors of",
         " the same length.")
  }
  J       <- length(l)
  if (any(l[1:(J - 1)] >= u[1:(J - 1)])) {
    stop("Each lower interim stopping boundary (in l) must be strictly less ",
         "then the corresponding upper interim stopping boundary (in u).")
  }
  if (l[J] != u[J]) {
    stop("The final lower stopping boundary (in l) must be equal to the final ",
         "upper stopping boundary (in u).")
  }
  if (any(n < 1, n%%1 != 0)) {
    stop("The stage-wise group size in the control arm (n) must be an integer ",
         "greater than or equal to 1.")
  }
  if (sigma <= 0) {
    stop("The true value of the standard deviation of the responses (sigma) ",
         "must be a real strictly greater than 0.")
  }
  if (ratio <= 0) {
    stop("The allocation ratio between the experimental and control arms ",
         "(ratio) must be a real strictly greater than 0.")
  }
  if (sigma_z <= 0) {
    stop("The assumed value of the standard deviation of the responses in the ",
         "z test-statistics (sigma_z) must be a real strictly greater than 0.")
  }
  if ((n*ratio)%%1 != 0) {
    stop("The product of the allocation ratio between the experimental and ",
         "control arms (ratio) and the stage-wise group size in the control ",
         "arm (n) must be an integer greater than or equal to 1.")
  }
  if (!is.logical(t_test)) {
    stop("t_test must be logical.")
  }
  if (!is.logical(adjusted)) {
    stop("adjusted must be logical.")
  }
  if (all(t_test, adjusted)) {
    stop("Specifying t_test = T and adjusted = T is not currently supported.")
  }
  if (any(alpha <= 0, alpha >= 1)) {
    stop("The significance level to use in confidence interval construction ",
         "(alpha) must be a real strictly between 0 and 1.")
  }
  if (any(replicates < 1, replicates%%1 != 0)) {
    stop("The number of replicate simulations to conduct (replicates) must be ",
         "an integer greater than or equal to 1.")
  }
  if (all(alpha != 0.05, adjusted)) {
    warning("alpha has been changed from default but this will have no effect ",
            "given the choice of adjusted.")
  }
  if (all(sigma != sigma_z, t_test)) {
    warning("sigma_z has been changed from default but this will have no ",
            "effect given the choice of t_test.")
  }

  ##### Initialise internal functions ##########################################

  covariance <- function(sqrt_I) {
    J               <- length(sqrt_I)
    Sigma           <- diag(J)
    for (j1 in 2:J) {
      j2            <- 1:(j1 - 1)
      Sigma[j1, j2] <- Sigma[j2, j1] <- sqrt_I[j2]/sqrt_I[j1]
    }
    Sigma
  }

  information <- function(n, J, sigma, ratio) {
    (1:J)*n*ratio/(sigma^2*(1 + ratio))
  }

  p_value <- function(theta, J, l, u, sqrt_I, Lambda) {
    theta_sqrt_I <- theta*sqrt_I
    P            <- pnorm(u[1], theta_sqrt_I[1], lower.tail = F)
    if (J > 1) {
      for (j in 2:J) {
        P        <- P + mvtnorm::pmvnorm(c(l[1:(j - 1)], u[j]),
                                         c(u[1:(j - 1)], Inf),
                                         theta_sqrt_I[1:j],
                                         sigma = Lambda[1:j, 1:j])
      }
    }
    P
  }

  p_value_root <- function(theta, J, l, u, sqrt_I, Lambda, h) {
    theta_sqrt_I <- theta*sqrt_I
    P            <- pnorm(u[1], theta_sqrt_I[1], lower.tail = F)
    if (J > 1) {
      for (j in 2:J) {
        P        <- P + mvtnorm::pmvnorm(c(l[1:(j - 1)], u[j]),
                                         c(u[1:(j - 1)], Inf),
                                         theta_sqrt_I[1:j],
                                         sigma = Lambda[1:j, 1:j])
      }
    }
    P - h
  }

  sim_gs_internal <- function(tau, n, l, u, sigma, ratio, sigma_z, t_test,
                              adjusted, alpha, replicates, J) {
    seq_J                <- 1:J
    Lambda               <- covariance(sqrt(seq_J))
    n_vec                <- c(0, seq_J*n)
    rn_vec               <- ratio*n_vec
    ss_vec               <- n_vec[seq_J + 1] + rn_vec[seq_J + 1]
    rej                  <- cov_naive <- pval_naive <- 0
    est_naive            <- sample_size <- numeric(replicates)
    if (adjusted) {
      cov_adj            <- pval_adj <- 0
      est_adj            <- est_naive
    } else {
      bias_adj           <- cov_adj <- est_adj <- pval_adj <- rmse_adj <- NA
    }
    Z_j                  <- numeric(J)
    x_0                  <- numeric(rn_vec[J + 1])
    x_1                  <- numeric(n_vec[J + 1])
    if (!t_test) {
      sqrt_I             <- sqrt(information(n, J, sigma_z, ratio))
    } else {
      sqrt_I             <- Z_j
      df                 <- ss_vec - 2
    }
    for (i in 1:replicates) {
      sum_x_0            <- sum_x_1 <- 0
      for (j in 1:J) {
        range_0          <- (1 + n_vec[j]):n_vec[j + 1]
        range_1          <- (1 + rn_vec[j]):rn_vec[j + 1]
        x_0[range_0]     <- rnorm(n, sd = sigma)
        x_1[range_1]     <- rnorm(rn_vec[2], tau, sigma)
        sum_x_0          <- sum_x_0 + sum(x_0[range_0])
        sum_x_1          <- sum_x_1 + sum(x_1[range_1])
        if (t_test) {
          sigma_hat_j2   <-
            (sum((x_0[1:n_vec[j + 1]] - sum_x_0/n_vec[j + 1])^2) +
               sum((x_1[1:rn_vec[j + 1]] - sum_x_1/rn_vec[j + 1])^2))/df[j]
          sqrt_I[j]      <- sqrt(n_vec[j + 1]*ratio/(sigma_hat_j2*(ratio + 1)))
        }
        Z_j[j]           <-
          (sum_x_1/rn_vec[j + 1] - sum_x_0/n_vec[j + 1])*sqrt_I[j]
        if (any(Z_j[j] > u[j], Z_j[j] <= l[j])) {
          sample_size[i] <- ss_vec[j]
          if (Z_j[j] > u[j]) {
            rej          <- rej + 1
          }
          est_naive[i]   <- Z_j[j]/sqrt_I[j]
          pval_naive     <- pval_naive + (1 - pnorm(Z_j[j])[1])
          lci_naive      <- est_naive[i] - qnorm(1 - alpha)/sqrt_I[j]
          cov_naive      <- cov_naive + (lci_naive <= tau)
          if (adjusted) {
            if (j == 1) {
              pval_adj   <- pval_adj + (1 - pnorm(Z_j[j])[1])
              est_adj[i] <- uniroot(f = p_value_root,
                                    interval = c(-10000, 10000),
                                    J = j, l = l[1], u = Z_j[j],
                                    sqrt_I = sqrt_I[1:j],
                                    Lambda = Lambda[1:j, 1:j],
                                    h = 0.5)$root
              lci_adj    <- uniroot(f = p_value_root,
                                    interval = c(-10000, 10000),
                                    J = j, l = l[1], u = Z_j[j],
                                    sqrt_I = sqrt_I[1:j],
                                    Lambda = Lambda[1:j, 1:j],
                                    h = alpha)$root
            } else {
              pval_adj   <- pval_adj + p_value(0, j, l[1:j],
                                               c(u[1:(j - 1)], Z_j[j]),
                                               sqrt_I[1:j], Lambda[1:j, 1:j])
              est_adj[i] <- uniroot(f = p_value_root,
                                    interval = c(-10000, 10000),
                                    J = j, l = l[1:j],
                                    u = c(u[1:(j - 1)], Z_j[j]),
                                    sqrt_I = sqrt_I[1:j],
                                    Lambda = Lambda[1:j, 1:j],
                                    h = 0.5)$root
              lci_adj    <- uniroot(f = p_value_root,
                                    interval = c(-10000, 10000),
                                    J = j, l = l[1:j],
                                    u = c(u[1:(j - 1)], Z_j[j]),
                                    sqrt_I = sqrt_I[1:j],
                                    Lambda = Lambda[1:j, 1:j],
                                    h = alpha)$root
            }
            cov_adj      <- cov_adj + (lci_adj <= tau)
          }
          break
        }
      }
      if (all(adjusted, i%%ceiling(replicates/10) == 0)) {
        message(ceiling(i/10), "% of simulations completed.")
      }
    }
    P                    <- rej/replicates
    ESS                  <- mean(sample_size)
    SDSS                 <- sqrt(sum((sample_size - ESS)^2)/(replicates - 1))
    sort_sample_size     <- sort(sample_size)
    MSS                  <-
      0.5*(sort_sample_size[ceiling(0.5*replicates)] +
             sort_sample_size[ceiling(0.5*replicates + 1)])
    est_naive            <- mean(est_naive)
    bias_naive           <- est_naive - tau
    rmse_naive           <- sqrt(sum((est_naive - tau)^2)/replicates)
    cov_naive            <- cov_naive/replicates
    pval_naive           <- pval_naive/replicates
    if (adjusted) {
      est_adj            <- sum(est_adj)/replicates
      bias_adj           <- est_adj - tau
      rmse_adj           <- sqrt(sum((est_adj - tau)^2)/replicates)
      cov_adj            <- cov_adj/replicates
      pval_adj           <- pval_adj/replicates
    }
    c(P, ESS, SDSS, MSS, est_naive, est_adj, bias_naive, bias_adj,
      rmse_naive, rmse_adj, pval_naive, pval_adj, cov_naive, cov_adj)
  }

  ##### Perform main computations ##############################################

  opchar <- sim_gs_internal(tau, n, l, u, sigma, ratio, sigma_z, t_test,
                            adjusted, alpha, replicates, J)

  ##### Output #################################################################

  return(list(P = opchar[1], ESS = opchar[2], SDSS = opchar[3], MSS = opchar[4],
              est_naive = opchar[5], est_adj = opchar[6],
              bias_naive = opchar[7], bias_adj = opchar[8],
              rmse_naive = opchar[9], rmse_adj = opchar[10],
              pval_naive = opchar[11], pval_adj = opchar[12],
              cov_naive = opchar[13], cov_adj = opchar[14]))
}