#' @useDynLib bvashr
#' @import ashr stats4

#' @title Estimate the mixture proportions
#'
#' @description estimates the mixture proportions via convex optimzation condtional
#'              on the current estimate of the loadings matrix
#' 
#' @param Y matrix of betahats n x 2 
#' @param lambda_12 double causal effect of trait 2 --> trait 1
#' @param lambda_21 double causal effect of trait 1 --> trait 2
#' @param U array prior diagonal 2 x 2 covariance matricies
#' @param rho double correlation of the estimated effects
#' @param S matrix of std. errors n x 2
#'
#' @return w_next vector next estimate of mixture proportions
#' @export
estimate_w_given_lambda_rho <- function(Y, lambda_12, lambda_21, U, rho, S){

    K <- dim(U)[3]
    L <- get_causal_lik_mat(Y, lambda_12, lambda_12, U, rho, S)
    w_res = ashr::mixIP(L, rep(1.0 / K, K))
    w_next = w_res$pihat

    return(w_next)

}

#' @title Estimate lambda and rho  
#'
#' @description estimate lambda and rho condtional on 
#'              knowing the mixture proportions using 
#'              a numerical optimizer
#' 
#' @param Y matrix of betahats n x 2 
#' @param w vector prior K mixture proportions (weights)
#' @param U array prior diagonal 2 x 2 covariance matricies
#' @param S matrix of std. errors n x 2
#'
#' @return lambda_rho list of estimated parameters
#' @export
estimate_lambda_rho_given_w <- function(Y, w, U, S){

    minus_log_likelihood <- function(lambda_12, lambda_21, rho){
        -1 * get_causal_log_lik(lambda_12=lambda_12, lambda_21=lambda_21, rho=rho, 
                                Y=Y, U=U, w=w, S=S)
    }

    res_mle <- stats4::mle(minuslogl = minus_log_likelihood, 
                           start=list(lambda_12=0.0, lambda_21=0.0, rho=0.0), 
                          )

    return(res_mle)

}
