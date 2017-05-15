#' @import ExtremeDeconvolution MCMCpack mvtnorm

#' @title Estimate prior
#'
#' @description TODO: fill in
#'
#' @param X n x 2 matrix of estimated effects
#' @param S n x 2 matrix of std. errors
#' @param K number of mixture components to use in EM
#'
#' @return prior_res list from output of ExtremeDeconvolution
#' @export
estimate_prior <- function(X, S, K){
  
  X <- X
  C <- (t(X) %*% X) / (nrow(X) - 1) # sample covariance mat
  mu_fix <- matrix(0, nrow=K, ncol=2) # fix mean at 0
  
  # fix covariance of first component to spike at 0
  sigma_init <- list()
  sigma_init[[1]] <- matrix(0, nrow=2, ncol=2)
  
  # list of booleans for which component has fixed covariance
  sigma_fix <- list()
  sigma_fix[[1]] <- TRUE
  
  # randomly intialize rest of covariance matrices
  for (k in 2:K){
    sigma_init[[k]] <- MCMCpack::riwish(3, C) # inv-wishart distributed with sample covariance mat
    sigma_fix[[k]] <- FALSE
  }
  
  # intialize mixture proportions
  pi_init <- rep(1 / K, K)
  
  # run the EM
  prior_res <- ExtremeDeconvolution::extreme_deconvolution(X, S^2, pi_init,
                                                           mu_fix, sigma_init,
                                                           fixmean=TRUE,
                                                           fixcovar=sigma_fix)
  
  # convert to array and rename
  prior_res$U0 <- array(as.numeric(unlist(prior_res$xcovar)), dim=c(2, 2, K))
  prior_res$pi_hat <- prior_res$xamp
  
  return(prior_res)
}

#' @title Sample prior
#'
#' @description TODO: fill in
#'
#' @param N number of samples
#' @param prior_res list output from extreme deconvolution
#'
#' @return R N x 3 matrix of samples and latent component as columns
#' @export
sample_prior <- function(N, prior_res){
  pi_hat <- prior_res$pi_hat
  U0 <- prior_res$U0
  K <- length(pi_hat)
  R <- matrix(NA, nrow=N, ncol=3)
  for(i in 1:N){
    for(k in 1:K){
      k <- which(rmultinom(1, 1, pi_hat) == 1)
      R[i, 1:2] <- mvtnorm::rmvnorm(n=1, mean=c(0, 0), sigma=U0[,,k])
      R[i, 3] <- k
    }
  }
  return(R)
}

# plot_fitted_prior <- function(){
#   
#   
# }