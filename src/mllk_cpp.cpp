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
                    arma::rowvec params_1, arma::rowvec params_2, 
                    arma::rowvec p)
{
  
  // map the functions names with the actual functions
  // (the type FunPtr and the density functions are defined in densities.h)
  map<std::string,FunPtr> funMap;
  funMap["gamma"] = dgamma_rcpp;
  funMap["weibull"] = dweibull_rcpp;
  funMap["lnorm"] = dlnorm_rcpp;
  funMap["exp"] = dexp_rcpp;
  funMap["vm"] = dvm_rcpp;
  funMap["wrpcauchy"] = dwrpcauchy_rcpp;
  
  
  int N = x.n_cols;
  int nObs = x.n_rows;
  arma::mat probs(nObs, N, arma::fill::ones);
  
  //for(unsigned int dist=0; dist<dists.size(); dist++){ // loop through variables
  //  string current_dist = dists[dist];
  //  if(arma::any(p.subvec(N*(N-1), N*(N-1)+N)>0)){ // check if variable has autoregression
  //    
  //    for (unsigned int j=0; j<N; j++){ // loop through states
  //      if(p.at(N*(N-1)+j-1)>0){ // check if current state has autoregression
  //        arma::uvec ind = arma::find_finite(x.col(dist)); // build vector of indices
  //        ind = ind.subvec(p.at((dist-1)*N+j), ind.size()-1);
  //        arma::mat autocor_ind(ind.size(), p.at((dist-1)*N+j-1));
  //        
  //        for (int i = 0; i < p.at((dist-1)*N+j-1); i++) { // fill autocorrelation index matrix
  //          autocor_ind.col(i) = ind - p.at((dist-1)*N)+j-1) + i;
  //        }
  //        
  //        // Replace missing values with column means
  //        arma::mat x_wo_na = clone(x);
  //        arma::uvec na_rows = arma::find_nonfinite(x_wo_na.col(dist));
  //        x_wo_na(na_rows, dist).fill(arma::mean(x_wo_na.col(dist)));
  //        
  //        // Compute autocorrelation values
  //        arma::mat autocor_values = arma::zeros(autocor_ind.n_rows, autocor_ind.n_cols);
  //        for (int i = 0; i < autocor_ind.n_cols; i++) {
  //          arma::uvec a = autocor_ind.col(i);
  //          autocor_values.col(i) = x_wo_na(a, dist);
  //        }
  //        
  //        stepProb = funMap[dists[dist]](x.col(dist),params_1.at((dist-1)*N+j-1),
  //                                       params_2.at((dist-1)*N+j-1));

          /// t.b.c
          
          
          
          // multiply entries in allprobs (row is observation, column is state) 
          // with density of current variable
          
  //      } else{ // else, current state has no autoregression
  //
  //        // loop through states, extract parameters, multiply matrix with density
  //        
  //      }
  //    }
  //  } else{ // current variable has no autoregression
  //    // calculate indices
  //    for(unsigned int j=0; j<N; j++){ // loop through states
  //      // extract parameters, multiply matrix with density
  //    }
  //  }
  //}
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
// 
// @return Negative (penalized) log-likelihood.
// [[Rcpp::export]]
double forward_cpp(arma::mat allprobs, arma::rowvec delta, arma::mat Gamma,
                   int nObs, int nVars)
  {
  
  int N = allprobs.n_cols;
  arma::rowvec foo(N);

  foo = delta % allprobs.row(0);
  double mllk_scale = log(sum(foo));
  arma::rowvec phi = foo/sum(foo);
  for (unsigned int i=1; i<allprobs.n_rows; i++){
    foo = (phi*Gamma) % allprobs.row(i);
    mllk_scale = mllk_scale + log(sum(foo));
    phi = foo/sum(foo);
  }
  return -mllk_scale;
}

// [[Rcpp::export]]
int test_cpp(arma::rowvec p){
  Rcout << "Hi" << endl;
  //arma::rowvec p(5);
  Rcout << p.at(1) << endl;
  
  return(0);
}


