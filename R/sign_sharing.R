#' @import ashr

#' @title Compute matrix of sign probabilities 
#' 
#' @description computes matrix of of sign probabilities for 
#' a pair of traits
#'
#' @param Y matrix of betahats n x 2 
#' @param S matrix of std. errors n x 2
#'
#' @return M matrix n x 6 of sign probabilities
#' @export
get_sign_prob_mat <- function(Y, S){

    n <- nrow(Y) # number of snps

    # run ash on both traits
    print('running ash on trait 1')
    ash_res_1 <- ashr::ash(Y[, 1], S[, 1], mixcompdist='normal')
    
    print('running ash on trait 2')
    ash_res_2 <- ashr::ash(Y[, 2], S[, 2], mixcompdist='normal')

    # prior prob of 0
    pi0_1 <- ash_res_1$fitted_g$pi[1]
    pi0_2 <-ash_res_2$fitted_g$pi[1]

    # prior prob of being +,-
    pi1_1 <- .5 * (1 - pi0_1)
    pi1_2 <- .5 * (1 - pi0_2)

    # matrix of sign probabilities
    M <- matrix(NA, nrow=n, ncol=6)
    
    # trait 1
    M[,1] <- ash_res_1$result$NegativeProb / pi1_1 # P(D_1 | sign = -)
    M[,2] <- ash_res_1$result$lfdr / pi0_1 # P(D_1 | sign = 0)
    M[,3] <- ash_res_1$result$PositiveProb / pi1_1 # P(D_1 | sign = +)

    # trait 2
    M[,4] <- ash_res_2$result$NegativeProb / pi1_2 # P(D_2 | sign = -)
    M[,5] <- ash_res_2$result$lfdr / pi0_2 # P(D_2 | sign = 0)
    M[,6] <- ash_res_2$result$PositiveProb / pi1_2 # P(D_2 | sign = +)

    return(M)

}

#' @title Compute sign mixture proportions
#' 
#' @description TODO: description
#'
#' @param Y matrix of betahats n x 2 
#' @param S matrix of std. errors n x 2
#'
#' @return pi_df data.frame of sign mixture probabilities
#' @export
estimate_sign_pi <- function(Y, S){

    M <- get_sign_prob_mat(Y, S)
    L <- get_lik_mat_ss(M)

    print('estimating mixture proportions')
    
    # TODO: figure out why mixIP isn't working here
    pi_vec_res <- ashr::mixEM(L, rep(1.0 / 9, 9))
    pihat <- pi_vec_res$pihat
    sign_pairs <- c('--', '-0', '-+', '0-', '00', '0+', '+-', '+0', '++')
    pi_df <- data.frame(pi=pihat, sign_pair=sign_pairs)
    
    return(pi_df)

}

