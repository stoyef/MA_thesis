#include "densities.h"

// allprobs for likelihood computation in C++
//
// Used in mllk.R to speed up loop in allprobs calculation
//
// @param x Data
// @param dists Vector of variable names
// @param autocor Matrix of autoregression coefficients
// @param params List of lists of distribution parameters
// @param p Vector of autoregression degrees
// 
// @return Matrix of allprobs.
// [[Rcpp::export]]
arma::mat allprobs_cpp(arma::mat x, std::vector<std::string> dists, arma::mat autocor,
                    List params, IntegerVector p)
{
  int N = x.n_cols;
  int nObs = x.n_rows;
  arma::mat probs(nObs, N, arma::fill::ones);
  return probs;
}


// Forward algorithm for likelihood computation in C++
//
// Used in mllk.R to speed up loop in forward algorithm
//
// @param allprobs Output from allprobs_cpp
// @param delta Stationary distribution
// @param Gamma Transition probability matrix
// @param autocor Matrix of autoregression coefficients
// @param nObs Number of observations
// @param nVars Number of variables
// @param lambda Complexity penalty
// 
// @return Negative (penalized) log-likelihood.
// [[Rcpp::export]]
double forward_cpp(arma::mat allprobs, arma::rowvec delta, arma::mat Gamma,
                   arma::mat autocor, int nObs, int nVars, double lambda)
  {
  int N = allprobs.n_cols;
  double mllk_scale = 0;
  return -mllk_scale;
}


