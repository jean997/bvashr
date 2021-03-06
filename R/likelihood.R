#' @import mvtnorm matrixStats

#' @title Get an observation's likelihood under causal model
#' 
#' @description computes the likelihood for a single observation
#'              for the bvash causal model
#'
#' @param y_i vector containing betahats for two traits
#' @param lambda_12 double causal effect of trait 2 --> trait 1
#' @param lambda_21 double causal effect of trait 1 --> trait 2
#' @param U array prior diagonal 2 x 2 covariance matricies
#' @param rho double correlation of estimated effects
#' @param s_1 double std. error for trait 1
#' @param s_2 double std. error for trait 2
#'
#' @return lik_i vector likelihood under each component for a single observation
#' @export
get_causal_comp_lik_i <- function(y_i, lambda_12, lambda_21, U, rho, s_1, s_2){

    # number of mixture components
    K <- dim(U)[3]
    
    # loadings matrix
    Lambda <- matrix(c(1, lambda_21, lambda_12, 1), nrow=2, byrow=TRUE)
  
    # error matrix
    cov_i <- rho * s_1 * s_2
    V_i <- matrix(c(s_1^2, cov_i, cov_i, s_2^2), nrow=2, byrow=TRUE)
    
    # compute component likelihoods
    lik_i <- rep(NA, K)
    for(k in 1:K){
        Sigma_ik <- Lambda %*% U[,,k] %*% t(Lambda) + V_i
        lik_i[k] <- dmvnorm(y_i, sigma=Sigma_ik)
    }
  
    return(lik_i)

}

#' @title Get likelihood matrix under causal model
#'
#' @description computes likelihood matrix under the bvash causal model
#'
#' @param Y matrix of betahats n x 2 
#' @param lambda_12 double causal effect of trait 2 --> trait 1
#' @param lambda_21 double causal effect of trait 1 --> trait 2
#' @param U array prior diagonal 2 x 2 covariance matricies
#' @param rho double correlation of the estimated effects
#' @param S matrix of std. errors n x 2
#'
#' @return L matrix likelihood matrix n x K
#' @export
get_causal_lik_mat <- function(Y, lambda_12, lambda_21, U, rho, S){
  
    # number of observations
    n <- nrow(Y)
    stopifnot(ncol(Y) == 2)
    stopifnot(all(dim(S) == c(n, 2)))

    # compute tranpose of likelihood matrix (K x n)
    Lt <- sapply(1:n, FUN=function(i){
        get_causal_comp_lik_i(Y[i, ], lambda_12, lambda_21, U, rho, S[i, 1], S[i, 2])
    })
    
    return(t(Lt))

}

#' @title Get an observation's log likelihood under causal model
#'
#' @description computes the log-likelihood for a single observation
#'              for the bvash causal model
#'
#' @param y_i vector containing betahats for two traits
#' @param lambda_12 double causal effect of trait 2 --> trait 1
#' @param lambda_21 double causal effect of trait 1 --> trait 2
#' @param U array prior diagonal 2 x 2 covariance matricies
#' @param w vector prior K mixture proportions (weights)
#' @param rho double correlation of the estimated effects
#' @param s_1 double std. error for trait 1
#' @param s_2 double std. error for trait 2
#'
#' @return log_lik_i double likelihood a single observation
#' @export
get_causal_log_lik_i <- function(y_i, lambda_12, lambda_21, U, w, rho, s_1, s_2){
  
    # number of components
    K <- dim(U)[3]
    stopifnot(length(w) == K)
  
    # loadings matrix
    Lambda <- matrix(c(1, lambda_21, lambda_12, 1), nrow=2, byrow=TRUE)
  
    # error matrix
    cov_i <- rho * s_1 * s_2
    V_i <- matrix(c(s_1^2, rho * s_1 * s_2, rho * s_1 * s_2, s_2^2), nrow=2, byrow=TRUE)
  
    log_lik_i <- rep(NA, K)
    for(k in 1:K){
        Sigma_ik <- Lambda %*% U[,,k] %*% t(Lambda) + V_i
        log_lik_i[k] <- log(w[k]) * dmvnorm(y_i, sigma=Sigma_ik, log=TRUE)
    }
  
    return(matrixStats::logSumExp(log_lik_i))
  
}

#' @title Get log-likelihood under causal model 
#'
#' @description computes log-likelihood of all the data 
#'
#' @param Y matrix of betahats n x 2 
#' @param lambda_12 double causal effect of trait 2 --> trait 1
#' @param lambda_21 double causal effect of trait 1 --> trait 2
#' @param U array prior diagonal 2 x 2 covariance matricies
#' @param w vector prior K mixture proportions (weights)
#' @param rho double correlation of the estimated effects
#' @param S matrix of std. errors n x 2
#'
#' @return log_lik double log likelihood 
#' @export
get_causal_log_lik <- function(Y, lambda_12, lambda_21, U, w, rho, S){
  
    # number of observations
    n <- nrow(Y)
    stopifnot(ncol(Y) == 2)
    stopifnot(all(dim(S) == c(n, 2)))
  
    # compute log-liklihood using all observations
    log_lik <- sapply(1:n, FUN=function(i){
        get_causal_log_lik_i(Y[i, ], lambda_12, lambda_21, U, w, rho, S[i, 1], S[i, 2])
    })

    log_lik <- sum(log_lik)
    
    return(log_lik)
    
}
