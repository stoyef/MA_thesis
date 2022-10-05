# 2022-06-03
# Functions to calculate the negative log-Likelihood of different HMMs



###### general function for arbitrary distribution


#' Compute negative log-likelihood of an AR(p)-HMM
#'
#' Compute the negative log-likelihood, using one or several specified distributions.
#' The distributions have to be specified by their commonly known abbreviation in R, 
#' e.g. one of ['gamma', 'vm', 'pois', 'binom',...].
#' The named list of parameters (one value for each parameter and for each state) have to be
#' in suitable form, i.e. a vector of the working parameters.
#' In the likelihood computation, contemporaneous independence is assumed.
#' 
#' @param theta.star Vector of parameters in the following order: 1) Off-diagonal entries of TPM,
#'                   2) Distribution parameters for each state (each distribution at a time, i.e. 
#'                   first all parameters of dist1, then all parameters of dist2 etc.), 
#'                   3) Autoregression parameters (each distribution at a time, i.e. 
#'                   first all autoregression parameters of dist1, then all autocorrelation 
#'                   parameters of dist2 etc.).
#' @param dists Vector containing abbreviated names (in R-jargon) of the distributions 
#'              to be considered in the likelihood computation.
#' @param x Data vector or matrix for which the negative log-likelihood should be computed.
#' @param N Number of states.
#' @param p_auto Vector of degree of autoregression for each distribution, 0 = no autoregression
#' @param scale_kappa Default 1, Scaling factor for kappa to avoid numerical issues in optimization for large kappa.
#' @param zero_inf Default FALSE, indicates if the gamma distributed variables should incorporate zero-inflation.
#' 
#' @return Negative log-likelihood.
#' 
#' @export
#' @rdname mllk
mllk <- function(theta.star, dists, x, N, p_auto, scale_kappa=1, zero_inf=FALSE){
  
  # First: Working to natural parameters, list structure for better handling
  # We currently only use distributions with 2 parameters, once we use Poisson distribution etc, we need a re-write
  all_params = unstarize(theta.star=theta.star, N=N, p=p_auto, dists=dists, scale_kappa = scale_kappa, zero_inf = zero_inf)
  Gamma = all_params$Gamma
  delta = all_params$delta
  autocor = all_params$autocor
  params = all_params$params
  
  # transform data to matrix, if necessary
  if (is.vector(x)) x <- matrix(x, nrow=length(x))
  allprobs <- matrix(1,dim(x)[1],N)
  
  for (dist in 1:length(dists)){ # for each distribution (= column of x) to consider

    if (p_auto[dist]>0){
      autocor_m <- matrix(autocor[[dist]], ncol=p_auto[dist], byrow=TRUE) # matrix for easier handling later on
      
      ind <- which(!is.na(x[,dist]))[-c(1:p_auto[dist])] # change: we omit first p steps 
      # in order to always have the step in t-p
      autocor_ind <- matrix(NA,nrow=length(ind),ncol=p_auto[dist]) # matrix for indices of autocor data
      
      for (i in 1:p_auto[dist]){
        autocor_ind[,i] <- ind-p_auto[dist]+i-1
      }
      
      x_wo_na = x
      x_wo_na[which(is.na(x_wo_na[,dist])),dist] = mean(x_wo_na[,dist],na.rm = TRUE)
      autocor_ind <- apply(autocor_ind, 2, function(a)x_wo_na[a,dist]) # substitute indices with values
      # replace NA values in x that are put into autocor_ind with mean value 
      # (so that the data that has NA in previous time steps does not have to be deleted)
      
      for (j in 1:N){
        # here comes the autoregression! -> computed inside the dens_<...> functions, considering the 
        # autocorrelation!
        theta_j <- params[[dist]]
        # theta_j consists the parameters of state j for each parameter in theta$params
        for (i in names(theta_j)) theta_j[i][[1]] = theta_j[i][[1]][j] 
        
        # computation of the different densities is outsourced to the respective 
        # functions with name "dens_<name_of_distribution>, e.g. dens_gamma for gamma distribution
        allprobs[ind,j] <- allprobs[ind,j] * 
                              match.fun(paste('dens_', dists[dist], sep=""))(x[ind,dist], theta_j, autocor_ind, autocor_m[j,], p_auto[dist])
        # here we have to choose mu_auto[ind], because
        # we have an individual mu for each data point
      }
    } else{
      ind <- which(!is.na(x[,dist]))
      
      for (j in 1:N){
        
        theta_j <- params[[dist]]
        # theta_j consists the parameters of state j for each parameter in theta$params
        for (i in names(theta_j)) theta_j[i][[1]] = theta_j[i][[1]][j] 
        
        allprobs[ind,j] <- allprobs[ind,j] * 
                              match.fun(paste('dens_', dists[dist], sep=""))(x[ind,dist], theta_j)
      }
    }
      
  }
  
  foo <- delta%*%diag(allprobs[1,])
  l <- log(sum(foo))
  phi <- foo/sum(foo)
  for (t in 2:dim(x)[1]){
    foo <- phi%*%Gamma%*%diag(allprobs[t,]) # here it can happen that 
                                            # allprobs[t,] = c(0,0) due to numeric issues. 
                                            # Then the function fails
    l <- l+log(sum(foo))
    phi <- foo/sum(foo)
  }
  return(-l)
}





