#' @useDynLib bvashr
#' @import ashr stats4

#' @title Estimate the mixture proportions
#'
#' @description estimates the mixture proportions via convex optimzation condtional
#'              on the current estimate of the loadings matrix
#' 
#' @param Y matrix of betahats n x 2 
#' @param S matrix of std. errors n x 2
#' @param U array prior diagonal 2 x 2 covariance matricies
#' @param lambda_12 double causal effect of trait 2 --> trait 1
#' @param lambda_21 double causal effect of trait 1 --> trait 2
#' @param rho double correlation of the estimated effects
#'
#' @return w_next vector next estimate of mixture proportions
#' @export
estimate_w_given_lambda_rho <- function(Y, S, U, lambda_12, lambda_21, rho){

    K <- dim(U)[3]
    L <- get_causal_lik_mat(Y, S, U, lambda_12, lambda_12, rho)
    w_res <- ashr::mixIP(L, rep(1.0 / K, K))
    w_next <- w_res$pihat

    return(w_next)

}

#' @title Estimate lambda and rho  
#'
#' @description estimate lambda and rho condtional on 
#'              knowing the mixture proportions using 
#'              a numerical optimizer
#' 
#' @param Y matrix of betahats n x 2 
#' @param S matrix of std. errors n x 2
#' @param U array prior diagonal 2 x 2 covariance matricies
#' @param w vector prior K mixture proportions (weights)
#'
#' @return lambda_rho list of estimated parameters
#' @export
estimate_lambda_rho_given_w <- function(Y, S, U, w){

    minus_log_likelihood <- function(lambda_12, lambda_21, z){
        rho <- tanh(z)
        minus_ll <- -1.0 * get_causal_log_lik(lambda_12=lambda_12, lambda_21=lambda_21, rho=rho, 
                                              Y=Y, U=U, w=w, S=S)
        return(minus_ll)
    }

    # numerically optimize log-likelihood
    fit <- stats4::mle(minuslogl=minus_log_likelihood, 
                       start=list(lambda_12=0.0, lambda_21=0.0, z=0.0) 
                      )

    # TODO: write some assertions about covergence

    # reformat optimized parameters
    mle_res <- list(lambda_12=as.numeric(fit@coef[1]), 
                    lambda_21=as.numeric(fit@coef[2]), 
                    rho=as.numeric(tanh(fit@coef[3])))

    return(mle_res)

}

#' @title Estimate w, lambda, and rho in 2-step iterative algorithim
#'
#' @description estimates w, lambda and rho in a 2-step iterative algorithim where
#'              in each iteration we first estimate w given lambda and rho by convex
#'              optimazation and then estimate lambda and rho given w using black box 
#'              numerical optimization
#'
#' @param Y matrix of betahats n x 2 
#' @param S matrix of std. errors n x 2
#' @param U array prior diagonal 2 x 2 covariance matricies
#'
#' @return w_lambda_rho_res list of optimized parameters 
#' @export
estimate_w_lambda_rho <- function(Y, S, U, max_iter=10){

    # number of components
    K <- dim(U)[3]

    # intialize params
    w_l <- rep(1.0 / K, K)
    lambda_12_l <- 0.0
    lambda_21_l <- 0.0
    rho_l <- 0.0 
    ll <- rep(NA, max_iter)

    # run 2-step iterative algorithim
    for(l in 1:max_iter){

        # estimate the mixture proportions given the loadings and rho
        w_l <- estimate_w_given_lambda_rho(Y, S, U, lambda_12_l, lambda_21_l, rho) 

        # estimate the loadings and rho given the mixture proportions
        lambda_rho_res <- estimate_lambda_rho_given_w(Y, S, U, w_l)
        lambda_12_l <- lambda_rho_res$lambda_12
        lambda_21_l <- lambda_rho_res$lambda_21
        rho_l <- lambda_rho_res$rho

        # compute log-likelihood at current iteration
        ll[l] <- get_causal_log_lik(Y, S, U, w_l, lambda_12_l, lambda_21_l, rho_l)
           
        # print update of current iteration
        print(l)

    }
   
    # list of optimized parameters
    w_lambda_rho_res <- list(lambda_12=lambda_12_l, lambda_21=lambda_21_l, w=w_l, rho=rho_l)

    return(w_lambda_rho_res)

}


