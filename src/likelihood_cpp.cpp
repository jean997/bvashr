#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <math.h>

using namespace arma;
using namespace Rcpp;

// Helper functions for likelihood computations --------------------

const double log2pi = std::log(2.0 * M_PI);

//' @title Mahalanobis distance
//' 
//' @description from http://gallery.rcpp.org/articles/dmvnorm_arma/
//' 
//' @param x matrix of data
//' @param mean vector of means for MVN density
//' @param sigma variance-covariance matrix for MVN density
//' 
//' @return mahalanobis distance
//' @export
//'
// [[Rcpp::export]]
arma::vec Mahalanobis(arma::mat x, arma::rowvec center, arma::mat cov) {
    
    int n = x.n_rows;
    arma::mat x_cen;
    x_cen.copy_size(x);
    for (int i=0; i < n; i++) {
        x_cen.row(i) = x.row(i) - center;
    }
    
    return sum((x_cen * cov.i()) % x_cen, 1);
    
}

//' @title Multivariate normal density
//' 
//' @description from http://gallery.rcpp.org/articles/dmvnorm_arma/
//' 
//' @param x matrix of data
//' @param mean vector of means for MVN density
//' @param sigma variance-covariance matrix for MVN density
//' 
//' @return logretval evaluated density
//' @export
//'
// [[Rcpp::export]]
arma::vec dmvnorm_arma(arma::mat x, arma::rowvec mean, arma::mat sigma, bool log = false) { 
    
    arma::vec distval = Mahalanobis(x,  mean, sigma);
    double logdet = sum(arma::log(arma::eig_sym(sigma)));
    arma::vec logretval = -( (x.n_cols * log2pi + logdet + distval) / 2  ) ;
    
    if (log) { 
        return(logretval);
    } else { 
        return(exp(logretval));
    }
    
}

//' @title Bivariate normal density
//' 
//' @description density of bivariate normal distribution
//' 
//' @param y vector of length 2 of observations
//' @param mu vector of length 2 of means 
//' @param sigma 2 x 2 variance-covariance matrix
//' @param log_lik bool 
//' 
//' @return evaluated density
//' @export
//'
//'
// [[Rcpp::export]]
double dbvnorm(arma::rowvec y, arma::rowvec mu, arma::mat sigma, bool log_lik=false){
    
    arma::mat x(1, 2);
    x(0, 0) = y(0);
    x(0, 1) = y(1);
    arma::vec lik_vec;
    lik_vec = dmvnorm_arma(x, mu, sigma);
    
    // if nan return 0
    if(isnan(lik_vec[0]) == 1){
        
        return 0.0;  
    
    } else {
        
        // return log-lik
        if(log_lik){
            return std::log(lik_vec[0]); 
        // return lik
        } else {
            return lik_vec[0];
        }
        
    }
}

//' @title Log sum exponential
//'
//' @description inspired by the implementation in https://github.com/dcgerard/updog/blob/master/src/utitility.cpp
//'
//' @param y vector to be log-sum-exponentiated 
//'
//' @return lse doube of log-sum-exponentiated vector 
//' @export
//'
// [[Rcpp::export]]
double log_sum_exp(arma::vec y){
    
    // logsumexp
    double lse;

    // get max value 
    double max_ele = y.max();
        
    // do logsumexp
    if (max_ele == R_NegInf) {
        lse = R_NegInf;
    } else {
        lse = std::log(arma::sum(arma::exp(y - max_ele))) + max_ele;
    }
    
    return lse;
    
}

// Causal factor model functions -----------------------------

//' @title Gets likelihood matrix for causal factor model (cfm)
//' 
//' @description Computes n x K matrix of component likelihoods
//' 
//' @param y matrix n x 2 of observations
//' @param s matrix n x 2 of std. errors
//' @param u cube 2 x 2 x K of prior covariance matricies
//' @param lambda_12 double effect of varialbe 2 on variable 1
//' @param lambda_21 double effect of variable 1 on variable 2
//' @param rho double correlation of estimated effects
//' 
//' @return lik_mat matrix n x K of component likelihoods
//' @export
//'
// [[Rcpp::export]]
arma::mat get_lik_mat_cfm(const arma::mat y, const arma::mat s, 
                          const arma::cube u, const double lambda_12, 
                          const double lambda_21, const double rho){
    
    // number of observations and components
    const int N = y.n_rows;
    const int K = u.n_slices;
    
    // mean of each bivariate normal component
    arma::rowvec mu; 
    mu.zeros(2);
    
    // loadings matrix
    arma::mat lambda(2, 2);
    lambda(0, 0) = 1.0;
    lambda(0, 1) = lambda_21;
    lambda(1, 0) = lambda_12;
    lambda(1, 1) = 1.0;
    
    // error variance componet
    arma::mat e_i(2, 2);
    
    // covariance matrix
    arma::mat sigma_ik(2, 2);
    
    // likelihood matrix
    arma::mat lik_mat(N, K);
    
    for(int i=0; i<N; i++){
        for(int k=0; k<K; k++){
            
            // fill in error variance componet
            e_i(0, 0) = std::pow(s(i, 0), 2);
            e_i(1, 1) = std::pow(s(i, 1), 2);
            e_i(0, 1) = rho * s(i, 0) * s(i, 1);
            e_i(1, 0) = rho * s(i, 0) * s(i, 1);
            
            // covariance matrix of k component
            sigma_ik = lambda * (u.slice(k) * lambda.t()) + e_i;
            
            // fill in likelihood matrix
            lik_mat(i, k) = dbvnorm(y.row(i), mu, sigma_ik);
            
        }
    }
    
    return lik_mat;
    
}

//' @title Gets log likelihood for causal factor model (cfm)
//'
//' @description Computes likelihood under causal factor model
//'
//' @inheritParams get_lik_mat_cfm
//' @param pi_vec vector vector of mixture proportions
//'
//' @return ll double likelihood under causal factor model
//' @export
//'
// [[Rcpp::export]]
double get_log_likelihood_cfm(const arma::mat y, const arma::mat s, const arma::cube u,
                              const double lambda_12, const double lambda_21, 
                              const double rho, const arma::vec pi_vec){

    // number of observations and components
    const int N = y.n_rows;
    const int K = u.n_rows;

    // mean of each bivariate normal component
    arma::rowvec mu;
    mu.zeros(2);

    // loadings matrix
    arma::mat lambda(2, 2);
    lambda(0, 0) = 1.0;
    lambda(1, 1) = 1.0;
    lambda(0, 1) = lambda_21;
    lambda(1, 0) = lambda_12;

    // error variance componet
    arma::mat e_i(2, 2);

    // covariance matrix
    arma::mat sigma_ik(2, 2);
    
    // vector of log-likelihoods for ith observation
    arma::vec ll_i(K);
    
    // log-likelihood
    double ll = 0.0;

    for(int i=0; i<N; i++){
        for(int k=0; k<K; k++){

            // fill in error variance component
            e_i(0, 0) = std::pow(s(i, 0), 2);
            e_i(1, 1) = std::pow(s(i, 1), 2);
            e_i(0, 1) = rho * s(i, 0) * s(i, 1);
            e_i(1, 0) = rho * s(i, 0) * s(i, 1);

            // covariance matrix of ith observation and kth component
            sigma_ik = lambda * (u.slice(k) * lambda.t()) + e_i;

            // log-likelihood
            ll_i(k) = std::log(pi_vec(k)) + dbvnorm(y.row(i), mu, sigma_ik, true);

        }
        
        // add to log-likelihood over components using log-sum-exp trick
        ll += log_sum_exp(ll_i);

    }

    return ll;

}

//' @title Compute likelihood matrix for signing sharing model
//'
//' @description TODO: description
//' 
//' @param m matrix of sign probabilites
//'
//' @export
// [[Rcpp::export]]
arma::mat get_lik_mat_ss(arma::mat m){
    
    // number of observations and components
    const int n = m.n_rows;
    int j;
    
    // likelihood matrix
    arma::mat lik_mat(m.n_rows, 9);
    
    for(int i=0; i<n; i++){

        j = 0;

        for(int k=0; k<3; k++){
            for(int l=0; l<3; l++){
            
                // fill in likelihood matrix
                lik_mat(i, j) = (m(i, k) / m(i, 1)) * (m(i, 3 + l) / m(i, 4));
                j = j + 1;

            }
        }
    }
    
    return lik_mat;

}

