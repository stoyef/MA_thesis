#include "allprobs_cpp.h"

//' C++ version of negative (penalized) log-likelihood of an AR(p)-HMM
//'
//' Compute the negative log-likelihood, using one or several specified distributions.
//' The distributions have to be specified by their commonly known abbreviation in R, 
//' e.g. one of ['gamma', 'vm', 'pois', 'binom',...].
//' The named list of parameters (one value for each parameter and for each state) have to be
//' in suitable form, i.e. a vector of the working parameters.
//' In the likelihood computation, contemporaneous independence is assumed.
//' Includes an optional penalization term \eqn{\lambda} for parameter selection of \eqn{p}.
//' 
//' @param theta.star Vector of parameters in the following order: 1) Off-diagonal entries of TPM,
//'                   2) Distribution parameters for each state (each distribution at a time, i.e. 
//'                   first all parameters of dist1, then all parameters of dist2 etc.), 
//'                   3) Autoregression parameters (each distribution at a time, i.e. 
//'                   first all autoregression parameters of dist1, then all autocorrelation 
//'                   parameters of dist2 etc.) -> degree parameters for each state of each variable.
//' @param dists Vector containing abbreviated names (in R-jargon) of the distributions 
//'              to be considered in the likelihood computation.
//' @param x Data vector or matrix for which the negative log-likelihood should be computed.
//' @param N Number of states.
//' @param p_auto Vector of autoregression degrees, one value for each state of each variable 
//'               (in case of penalization choose upper bound of number of parameters).
//' @param lambda Complexity penalty (≥0) for autoregression parameters \eqn{\phi}.(default: 0 no penalization).
//' @param scale_kappa Default 1, Scaling factor for kappa to avoid numerical issues in optimization for large kappa.
//' @param zero_inf Default FALSE, indicates if the gamma distributed variables should incorporate zero-inflation.
//' @param alt_data Default NULL, provide data here if variable name x is already taken in wrapper function.
//' 
//' @return Negative (penalized) log-likelihood.
// [[Rcpp::export]]
mllk <- function(theta.star, dists, x, N, p_auto, lambda=0, scale_kappa=1, zero_inf=FALSE, alt_data=NULL){
  
// First: Working to natural parameters, list structure for better handling
// We currently only use distributions with 2 parameters, once we use Poisson distribution etc, we need a re-write
  all_params = unstarize(theta.star=theta.star, N=N, p=p_auto, dists=dists, scale_kappa = scale_kappa, zero_inf = zero_inf)
    Gamma = all_params$Gamma
  delta = all_params$delta
  autocor = all_params$autocor
  params = all_params$params
  
// transform data to matrix, if necessary
  if (is.vector(x)) x <- matrix(x, nrow=length(x))
    
// allprobs calculation in separate function allprobs
    allprobs = allprobs(x=x, dists=dists, autocor = autocor, params=params, N=N, p=p_auto)
      
      foo <- delta%*%diag(allprobs[1,])
      l <- log(sum(foo))
      phi <- foo/sum(foo)
      for (t in 2:dim(x)[1]){
        foo <- phi%*%Gamma%*%diag(allprobs[t,]) // here it can happen that 
// allprobs[t,] = c(0,0) due to numeric issues. 
// Then the function fails
        l <- l+log(sum(foo))
        phi <- foo/sum(foo)
      }
//// Lasso penalty for autoregression parameters
      if (lambda>0){
        sum_auto = 0
        auto = c()
        for (dist in 1:length(dists)){
          for (state in 1:N){
//sum_auto = sum_auto + sum(abs(autocor[[dist]][[state]]),na.rm = T)
            auto = c(auto, autocor[[dist]][[state]])
          }
        }
        sum_auto = sum(sqrt((auto+1e-10)^2), na.rm = T) // Approx wie in Marius Paper
        l <- l - lambda * sum_auto
      }
      return(-l)
}