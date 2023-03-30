#include "densities_cpp.h"

// From moveHMM: https://github.com/TheoMichelot/moveHMM/blob/master/src/nLogLike.cpp
// map the functions names with the actual functions
// (the type FunPtr and the density functions are defined in densities.h)
map<std::string,FunPtr> funMap;
funMap["gamma"] = dgamma_rcpp;
funMap["weibull"] = dweibull_rcpp;
funMap["lnorm"] = dlnorm_rcpp;
funMap["exp"] = dexp_rcpp;
funMap["vm"] = dvm_rcpp;
funMap["wrpcauchy"] = dwrpcauchy_rcpp;


//' C++ version of matrix of all probabilities
//' 
//' Calculate matrix of all probabilities for an AR(p)-HMM with or without within-state autoregression.
//' 
//' @param x Data matrix the model was fitted to.
//' @param dists Vector of distributions in R-jargon.
//' @param autocor list of autoregression vectors, respecting order of dists, one vector for each state.
//' @param params List of optimized parameters, returned by fitting the HMM.
//' @param N Number of states.
//' @param p Vector of degree of autoregression for each distribution (0 = no autoregression).
//' 
//' @return Matrix of all probabilities.
 // [[Rcpp::export]]

 arma::mat allprobs_cpp(arma::mat x, CharacterVector dists,
                        List autocor, List params, int N, IntegerVector p) {

   arma::mat probs(x.n_rows, N, arma::fill::ones);
   
   for (int dist = 0; dist < dists.size(); dist++) {
     if (any(p.subvec(dist * N, (dist + 1) * N - 1) > 0)) {
       for (int j = 0; j < N; j++) {
         if (p(dist * N + j) > 0) {
           
           arma::uvec ind = find_nonfinite(x.col(dist));
           ind.shed_rows(0, p(dist * N + j) - 1);
           
           arma::mat autocor_ind(ind.n_elem, p(dist * N + j));
           for (int i = 0; i < p(dist * N + j); i++) {
             autocor_ind.col(i) = ind-p(((dist-1)*N) + j)+i-1;
           }
           
           arma::mat x_wo_na = x;
           x_wo_na[which(is.na(x_wo_na[,dist])),dist] = mean(x_wo_na.col(dist).accu()/
             autocor_ind <- apply(autocor_ind, 2, function(a)x_wo_na[a,dist]) 
           
               
   
   
   
   
 
  
  return probs;
}
 
 
 allprobs <- function(x, dists, autocor, params, N, p){
  if (is.vector(x)) x <- matrix(x, nrow=length(x))
    allprobs <- matrix(1,dim(x)[1],N)
    
    for (dist in 1:length(dists)){
      if (any(p[((dist-1)*N) + 1:N]>0)){
        
        for (j in 1:N){ 
          if (p[(dist-1)*N + j] > 0){ # check if current state has autoregression
          
          ind <- which(!is.na(x[,dist]))[-c(1:p[((dist-1)*N) + j])] 
          autocor_ind <- matrix(NA,nrow=length(ind),ncol=p[((dist-1)*N) + j]) 
            
            for (i in 1:p[((dist-1)*N) + j]){
              autocor_ind[,i] <- ind-p[((dist-1)*N) + j]+i-1
            }
            
            x_wo_na = x
            x_wo_na[which(is.na(x_wo_na[,dist])),dist] = mean(x_wo_na[,dist],na.rm = TRUE)
              autocor_ind <- apply(autocor_ind, 2, function(a)x_wo_na[a,dist]) 
              
              theta_j <- params[[dist]]
            for (i in names(theta_j)) theta_j[i][[1]] = theta_j[i][[1]][j] 
            
            allprobs[ind,j] <- allprobs[ind,j] * 
              match.fun(paste('dens_', dists[dist], sep=""))(x[ind,dist], theta_j, 
                        autocor_ind, 
                        autocor[[dist]][[j]],
                                       p[((dist-1)*N) + j])
          } else{ # current state has no autoregression
            ind <- which(!is.na(x[,dist]))
            for (j in 1:N){
              
              theta_j <- params[[dist]]
              for (i in names(theta_j)) theta_j[i][[1]] = theta_j[i][[1]][j] 
              
              allprobs[ind,j] <- allprobs[ind,j] * 
                match.fun(paste('dens_', dists[dist], sep=""))(x[ind,dist], theta_j)
            }
          }
        }
      } else{ 
        
        ind <- which(!is.na(x[,dist]))
        
        for (j in 1:N){
          theta_j <- params[[dist]]
// theta_j consists the parameters of state j for each parameter in theta$params
          for (i in names(theta_j)) theta_j[i][[1]] = theta_j[i][[1]][j] 
          
          allprobs[ind,j] <- allprobs[ind,j] * 
            match.fun(paste('dens_', dists[dist], sep=""))(x[ind,dist], theta_j)
        }
      }
    }
    return(allprobs)
}




