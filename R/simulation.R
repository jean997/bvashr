#' @import mvtnorm matrixStats

#' @title Gets prior params for simulation 
#'
#' @description gets the prior (w, mu, U) parameters associated with a paticular 
#'              simulation scenario
#'
#' @param pi0 double proportion of nulls
#' @param scenario char of scenario 
#'
#' @return prior_params list of prior params
#' @export
get_prior_params <- function(pi0, scenario){

    pi1 <- 1 - pi0

    if(scenario == "spiky"){
    
        # mixture proportions and means
        w <- c(pi0, pi1 * .4, pi1 * .2, pi1 * .2, pi1 * .2)
        mu <- array(0, dim=c(2, 1, 5))
    
        # covariance matricies
        U <- array(dim=c(2, 2, 5))
        U[,,1] <- matrix(0, nrow=2, ncol=2)
        U[,,2] <- matrix(c(.25^2, 0.0, 0.0, .25^2), nrow=2, ncol=2)
        U[,,3] <- matrix(c(.5^2, 0.0,  0.0, .5^2), nrow=2, ncol=2)
        U[,,4] <- matrix(c(1, 0.0, 0.0, 1), nrow=2, ncol=2)
        U[,,5] <- matrix(c(2^2, 0.0, 0.0, 2^2), nrow=2, ncol=2)
    
    } else if(scenario == "bignormal"){
    
        # mixture proportions and means
        w <- c(pi0, pi1)
        mu <- array(0, dim=c(2, 1, 2))
    
        # covariance matricies
        U <- array(dim=c(2, 2, 2))
        U[,,1] <- matrix(0, nrow=2, ncol=2)
        U[,,2] <- matrix(c(4^2, 0.0, 0.0, 4^2), nrow=2, ncol=2)
    
    } else if(scenario == "nearnormal"){
    
        # mixture proportions and means
        w <- c(pi0, pi1 * (2/3), pi1 * (1/3))
        mu <- array(0, dim=c(2, 1, 3))
    
        # covariance matricies
        U <- array(dim=c(2, 2, 3))
        U[,,1] <- matrix(0, nrow=2, ncol=2)
        U[,,2] <- matrix(c(1, 0.0, 0.0, 1), nrow=2, ncol=2)
        U[,,3] <- matrix(c(2^2, 0.0, 0.0, 2^2), nrow=2, ncol=2)
    
    } else if(scenario == "bimodal"){ 
    
        # mixture proportions and means
        w <- c(pi0, pi1 * .5, pi1 * .5)
        mu <- array(dim=c(2, 1, 3))     
        mu[,,1] <- c(0, 0)
        mu[,,2] <- c(-2, -2)
        mu[,,3] <- c(2, 2)
    
        # covariance matricies
        U <- array(dim=c(2, 2, 3))
        U_alt <- matrix(c(1, 0.0, 0.0, 1), nrow=2, ncol=2)
        U[,,1] <- matrix(0, nrow=2, ncol=2)
        U[,,2] <- U_alt
        U[,,3] <- U_alt
    
    } else if(scenario == "flattop"){
    
        # mixture proportions and means
        pi1_alt <- pi1 * (1/7)
        w <- c(pi0, pi1_alt, pi1_alt, pi1_alt, pi1_alt, pi1_alt, pi1_alt, pi1_alt)
        mu <- array(dim=c(2, 1, 8))
        mu[,,1] <- c(0, 0)
        mu[,,2] <- c(-1.5, -1.5)
        mu[,,3] <- c(-1, -1)
        mu[,,4] <- c(-.5, -.5)
        mu[,,5] <- c(0, 0)
        mu[,,6] <- c(.5, .5)
        mu[,,7] <- c(1, 1)
        mu[,,8] <- c(1.5, 1.5)
    
        # covariance matricies
        U <- array(dim=c(2, 2, 8))
        U_alt <- matrix(c(.5^2, 0.0, 0.0, .5^2), nrow=2, ncol=2)
        U[,,1] <- matrix(0, nrow=2, ncol=2)
        U[,,2] <- U_alt
        U[,,3] <- U_alt
        U[,,4] <- U_alt
        U[,,5] <- U_alt
        U[,,6] <- U_alt
        U[,,7] <- U_alt
        U[,,8] <- U_alt
    
    } else if(scenario == "skew"){
    
        # mixture proportions and means
        w <- c(pi0, pi1 * (1/4), pi1 * (1/4), pi1 * (1/3), pi1 * (1/6))
        mu <- array(dim=c(2, 1, 5))
        mu[,,1] <- c(0, 0)
        mu[,,2] <- c(-2, -2)
        mu[,,3] <- c(-1, -1)
        mu[,,4] <- c(0, 0)
        mu[,,5] <- c(1, 1)
    
        # covariance matricies
        U <- array(dim=c(2, 2, 5))
        U[,,1] <- matrix(0, nrow=2, ncol=2)
        U[,,2] <- matrix(c(2^2, 0.0, 0.0, 2^2), nrow=2, ncol=2)
        U[,,3] <- matrix(c(1.5^2, 0.0, 0.0, 1.5^2), nrow=2, ncol=2)
        U[,,4] <- matrix(c(1, 0.0, 0.0, 1), nrow=2, ncol=2)
        U[,,5] <- matrix(c(1, 0.0, 0.0, 1), nrow=2, ncol=2)
  
    }

    prior_params <- list(w=w, mu=mu, U=U)

    return(prior_params)

}


#' @title Simulate data from causal model
#'
#' @description simulates effects sizes from the bvash causal factor model
#'
#' @param n integer number of observations to simulate
#' @param lambda_12 double 
#' @param lambda_21 double
#' @param rho double
#' @param prior_params list of mixture proportions, covariance matricies, and 
#'                     means output from get_prior_params
#'
#' @return sim_res list of simulated betahats Y and std. errors S
#' @export
simulate_data_causal <- function(n, lambda_12, lambda_21, rho, prior_params){
  
    # prior params
    w <- prior_params$w
    mu <- prior_params$mu
    U <- prior_params$U
  
    # number of components
    K <- length(w)   

    # matrix of loadings
    Lambda <- matrix(c(1, lambda_12, lambda_21, 1), byrow=TRUE, nrow=2, ncol=2) 
  
    # matrix of factors
    B <- matrix(NA, nrow=n, ncol=2)
    
    # matrix of beta_hats
    Y <- matrix(NA, nrow=n, ncol=2) 

    # matrix of std. errors
    S <- matrix(NA, nrow=n, ncol=2) 

    for(i in 1:n){

        # simulate latent component
        k <- which(rmultinom(1, 1, w) == 1)

        # simulate true beta
        B[i, ] <- mvtnorm::rmvnorm(n=1, mean=mu[,,k], sigma=U[,,k])
    
        # matrix of std. errors
        S[i, ] <- 1 / rgamma(2, 5, 5) 

        # matrix errors
        cov_i <- rho * S[i, 1] * S[i, 2]
        E_i <- matrix(c(S[i, 1]^2, cov_i, cov_i, S[i, 2]^2), nrow=2, ncol=2)
    
        # simulate beta_hat
        Y[i, ] <- mvtnorm::rmvnorm(n=1, mean=Lambda %*% B[i, ], sigma=E_i)
    }

    sim_res <- list(Y=Y, S=S)

    return(sim_res)

}

