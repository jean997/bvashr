#' @import dplyr ggplot2 knitr ExtremeDeconvolution MCMCpack mvtnorm ellipse

#' @title Convert std. error matrix to array
#'
#' @description converts std. error matrix (n x 2) into 
#'              array (2 x 2 x n) that in includes that null
#'              correlation on the off diagonal
#'
#' @return S_arr array of estimated effect error covariance matricies
#' @export 
se_mat_to_arr <- function(S, rho){

    N <- nrow(S)
    S_arr <- array(NA, dim=c(2, 2, N))
    for(i in 1:N){
        diag(S_arr[,,i]) <- S[i,]^2
        cov_i <- rho * S[i,1] * S[i,2]
        S_arr[1,2,i] <- cov_i
        S_arr[2,1,i] <- cov_i
    }

    return(S_arr)

} 

#' @title Estimate prior
#'
#' @description estimates parameters of a bivariate ash prior using an 
#'              EM algorithim implemented in the ExtremeDeconvolution package
#'
#' @param X n x 2 matrix of estimated effects
#' @param S n x 2 matrix of std. errors
#' @param K int of mixture components to use in EM
#' @param rho double null correlation in beta_hats
#'
#' @return prior_res list from output of ExtremeDeconvolution
#' @export
estimate_prior <- function(X, S, K, rho=NA){
  
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

        # initialze using inverse-wishart to ensure pos-def
        sigma_init[[k]] <- riwish(3, C) 
        sigma_fix[[k]] <- FALSE

    }
  
    # intialize mixture proportions
    pi_init <- rep(1 / K, K)
  
    # get se array if null corr is provided
    if(!is.na(rho)){
        R <- se_mat_to_arr(S, rho)
    } else {
        R <- S^2
    }
  
    # run the EM
    prior_res <- extreme_deconvolution(X, R, 
                                       pi_init,
                                       mu_fix, sigma_init,
                                       fixmean=TRUE,
                                       fixcovar=sigma_fix)
  
    # add / reformt fitted params to list
    prior_res$U0_init <- array(as.numeric(unlist(sigma_init)), dim=c(2, 2, K))
    prior_res$pi_init <- prior_res$xamp
    prior_res$U0 <- array(as.numeric(unlist(prior_res$xcovar)), dim=c(2, 2, K))
    prior_res$pi <- prior_res$xamp
    prior_res <- prior_res[-c(1, 2, 3)]

    return(prior_res)

}

#' @title Get prior component parameters
#'
#' @description computes formatted table (df, or kable) of fitted prior
#'              paramters including mixture propotions, correlation coef, 
#'              and std devs 
#'
#' @param prior_res list output from estimate_prior
#' @param table boolean if true return kable defaults to False
#'
#' @return parameters of prior components
#' @export
get_prior_comp <- function(prior_res, table=FALSE){

    pi_hat <- prior_res$pi
    K <- length(pi_hat)
    rhos <- c(0.0, rep(NA, K-1))
    sds_1 <- c(0.0, rep(NA, K-1))
    sds_2 <- c(0.0, rep(NA, K-1))
    for(k in 2:K){

        U0_k <- prior_res$U0[,,k]
        sds_1[k] <- sqrt(U0_k[1, 1])
        sds_2[k] <- sqrt(U0_k[2, 2]) 
        rhos[k] <- U0_k[1,2] / (sds_1[k] * sds_2[k])

    }
    comp_df <- data.frame(k=1:K - 1, 
                          pi_k=pi_hat, 
                          rho_k=rhos, 
                          sd_1k=sds_1, 
                          sd_2k=sds_2) 
  
    if(table){
        return(kable(comp_df, digits=4))
    } else {
        return(comp_df)
    }

}

#' @title Plot prior
#'
#' @description plots ellipes of fitted prior covariance
#'              matricies of each mixture component
#'
#' @param trait_1 char first trait name
#' @param trait_2 char second trait name
#' @param prior_res list output from estimate_prior
#'
#' @return p ggplot object of prior plot
#' @export
plot_prior <- function(trait_1, trait_2, prior_res){

    K <- length(prior_res$pi)
    pi_hat <- prior_res$pi
    comp_df <- get_prior_comp(prior_res)
  
    # null comp
    comp_df_0 <- comp_df %>% filter(k == 0)
    df_0 <- data.frame(x=0.0, y=0.0, k=0, pi_k=round(pi_hat[1], digits=3))
  
    # non-null comps
    comp_df_1 <- comp_df %>% filter(k != 0)
    df_1 <- data.frame()
    for(i in 1:nrow(comp_df_1)){
    
        # find points on ellipse for ith comp
        ellipse_i <- ellipse(comp_df_1[i, "rho_k"], 
                             scale=c(comp_df_1[i, "sd_1k"], comp_df_1[i, "sd_2k"]), 
                             centre=c(0.0, 0.0))

        # add points to data.frame 
        df_1 <- rbind(df_1, 
                      as.data.frame(ellipse_i) %>% 
                      mutate(k=comp_df_1[i, "k"], pi_k=round(pi_hat[comp_df_1[i, "k"] + 1], digits=3))
                     )

    }
  
    # create plot
    p <- ggplot() +
         geom_point(data=df_0, aes(x=x, y=y, color=paste0(k, ": ", pi_k)), size=2.5) + 
         geom_path(data=df_1, aes(x=x, y=y, color=paste0(k, ": ", pi_k))) + 
         xlab(paste0("Prior Effect (", toupper(trait_1), ")")) + 
         ylab(paste0("Prior Effect (", toupper(trait_2), ")")) +
         guides(color=guide_legend(title="pi_k")) +
         theme_bw()

    return(p)

}
