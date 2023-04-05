#include "densities.h"

// Compute negative (penalized) log-likelihood of an AR(p)-HMM in C++
//
// Compute the negative log-likelihood, using one or several specified distributions.
// The distributions have to be specified by their commonly known abbreviation in R, 
// e.g. one of ['gamma', 'vm', 'pois', 'binom',...].
// The named list of parameters (one value for each parameter and for each state) have to be
// in suitable form, i.e. a vector of the working parameters.
// In the likelihood computation, contemporaneous independence is assumed.
// Includes an optional penalization term \eqn{\lambda} for parameter selection of \eqn{p}.
// 
// @param theta.star Vector of parameters in the following order: 1) Off-diagonal entries of TPM,
//                   2) Distribution parameters for each state (each distribution at a time, i.e. 
//                   first all parameters of dist1, then all parameters of dist2 etc.), 
//                   3) Autoregression parameters (each distribution at a time, i.e. 
//                   first all autoregression parameters of dist1, then all autocorrelation 
//                   parameters of dist2 etc.) -> degree parameters for each state of each variable.
// @param dists Vector containing abbreviated names (in R-jargon) of the distributions 
//              to be considered in the likelihood computation.
// @param x Data vector or matrix for which the negative log-likelihood should be computed.
// @param N Number of states.
// @param p_auto Vector of autoregression degrees, one value for each state of each variable 
//               (in case of penalization choose upper bound of number of parameters).
// @param lambda Complexity penalty (≥0) for autoregression parameters \eqn{\phi}.(default: 0 no penalization).
// @param scale_kappa Default 1, Scaling factor for kappa to avoid numerical issues in optimization for large kappa.
// @param zero_inf Default FALSE, indicates if the gamma distributed variables should incorporate zero-inflation.
// @param alt_data Default NULL, provide data here if variable name x is already taken in wrapper function.
// 
// @return Negative (penalized) log-likelihood.
// [[Rcpp::export]]
double mllk_cpp(){
  double mllk_scale = 0;
  return -mllk_scale;
}


