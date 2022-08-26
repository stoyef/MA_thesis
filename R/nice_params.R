# 2022-06-08
# Functions to manipulate HMM parameters

#' Construct parameter object for HMM likelihood evaluation
#'
#' Construct a nice named list that contains all relevant parameters of an HMM.
#' This gets used in the Log-Likelihood computation to generalize to arbitrary 
#' distributions.
#' 
#' @param Gamma Full TPM.
#' @param autocor Vector or matrix of autocorrelation coefficients in suitable format.
#' @param ... List of other parameters, specified in "param" (e.g. mu, sigma, depending on distribution). 
#'            Here, notational conventions have to be considered when supplying the distributional parameters
#' 
#' @return Named list of all parameters
#' 
#' @export
#' @rdname nice_params
nice_params <- function(Gamma, autocor, ...){
  args = list(Gamma=Gamma,autocor=autocor,...)
  return(args)
}


#' Natural parameters to working parameters (theta->theta.star)
#'
#' Transform the natural parameters of HMMs to their working parameters for optimization (theta -> theta.star).
#' Currently only works for distributions with exactly two parameters.
#' 
#' @param theta Vector of natural parameters.
#' @param N Number of states of the HMM.
#' @param p Vector of degree of autocorrelation within the distributions.
#' @param dists Vector of the distributions in the HMM.
#' @param scale_kappa Default 1, Scaling factor for kappa to avoid numerical issues in optimization for large kappa.
#' 
#' @return Vector of working parameters of the HMM.
#' 
#' @export
#' @rdname starize
starize <- function(theta,N,p,dists, scale_kappa=1){
  # check for right number of parameters
  if (length(theta) != (N*(N-1)+2*N*length(dists)+sum(p)*N)){
    return("ERROR: Wrong number of parameters supplied.")
  } else{
    Gamma <- theta[1:(N*(N-1))]
    param_count <- N*(N-1)
    
    params <- c()
    for (dist in 1:length(dists)){
      if (dists[dist] == 'gamma'){
        params <- c(params, log(theta[param_count+1:(2*N)]))
        param_count <- param_count+2*N
      } else if (dists[dist] == 'vm'){
        params <- c(params, (theta[param_count+N+1:N] * cos(theta[param_count+1:N]))/scale_kappa)
        params <- c(params, (theta[param_count+N+1:N] * sin(theta[param_count+1:N]))/scale_kappa)
        param_count <- param_count+2*N
      } else if (dists[dist] == 'norm'){
        params <- c(params, theta[param_count+1:N])
        params <- c(params, log(theta[param_count+N+1:N]))
        param_count <- param_count+2*N
      }
    }
    
    if (any(p>0)){
      autocor <- qlogis(theta[param_count+1:(sum(p)*N)])
      theta.star <- c(Gamma, params, autocor)
    } else{
      theta.star <- c(Gamma, params)
    }
    return(theta.star)
  }
}



#' Working parameters to natural parameters (theta.star->theta)
#'
#' Transform the working parameters of HMMs to their natural parameters for optimization (theta.star -> theta).
#' Currently only works for distributions with exactly two parameters.
#' 
#' @param theta.star Vector of working parameters.
#' @param N Number of states of the HMM.
#' @param p Vector of degree of autocorrelation within the distributions.
#' @param dists Vector of the distributions in the HMM.
#' @param scale_kappa Default 1, Scaling factor for kappa to avoid numerical issues in optimization for large kappa.
#' 
#' @return List of natural parameters of the HMM.
#' 
#' @export
#' @rdname unstarize
unstarize <- function(theta.star,N,p,dists, scale_kappa=1){
  all_params = list()
  
  ### TPM
  Gamma <- diag(N)
  Gamma[!Gamma] <- exp(theta.star[1:((N-1)*N)])
  all_params$Gamma <- Gamma/rowSums(Gamma) 
  all_params$delta <- solve(t(diag(N)-all_params$Gamma+1),rep(1,N))
  
  ### distribution parameters
  n_dists <- length(dists)
  n_dist_params <- 2*n_dists*N
  counter = N*(N-1)
  params = list()
  for (dist in 1:n_dists){
    if (dists[dist]=='gamma'){
      params[[dist]] = list(mu=exp(theta.star[counter+1:N]),
                            sigma=exp(theta.star[counter+N+1:N]))
    } else if (dists[dist]=='vm'){
      params[[dist]] = list(mu = Arg((theta.star[counter+1:N]*scale_kappa)+1i*(theta.star[counter+N+1:N]*scale_kappa)),
                            kappa = sqrt((theta.star[counter+1:N]*scale_kappa)^2+(theta.star[counter+N+1:N]*scale_kappa)^2)
                            )
      
    } else if (dists[dist]=='norm'){
      params[[dist]] = list(mu=theta.star[counter+1:N],
                            sigma=exp(theta.star[counter+N+1:N]))
    } else{
      return(paste("ERROR: The distribution", dists[dist], "is not implemented."))
    }
    counter = counter+2*N
  }
  all_params$params = params
  
  ### Autocorrelation parameters
  if (any(p>0)){
    auto_params = theta.star[-c(1:counter)]
    autocor = list()
    counter = 0
    for (dist in 1:n_dists){
      autocor[[dist]] = plogis(auto_params[counter+1:(p[dist]*N)])
      counter = counter+p[dist]*N
    }
    all_params$autocor = autocor
  }
  
  return(all_params)
}




#' Generate starting parameters for a distribution
#'
#' Generate meaningful starting parameters for the optimization of a specified distribution using 
#' common sense and the supplied data.
#' 
#' @param data Vector of natural parameters.
#' @param dist Distribution in R-jargon.
#' @param N Number of states of the HMM.
#' 
#' @return Vector of starting parameters.
#' 
#' @export
#' @rdname starting_params_opt
starting_params_opt <- function(data,dist,N){
  params <- c()
  
  # Distribution parameters
  q_1 = seq(0,1,length=N+2) # quantiles to separate data
  q_2 = seq(0,1,length=N+1)
  if (dist=='gamma'){ 
      params <- c(params, as.numeric(quantile(data,probs=q_1[2:(N+1)]))) # mu
      s <- c()
      for (state in 1:N){
        s<-c(s,as.numeric(sd(data[data>quantile(data,q_2[state]) & data<quantile(data,q_2[state+1])]))) # sigma
      }
      params <- c(params,s)
  } else if (dist=='vm'){
    require(CircStats)
    params <- c(params, as.numeric(quantile(data,probs=q_1[2:(N+1)]))) # mu
    k <- c()
    for (state in 1:N){
      k<-c(k,as.numeric(est.kappa(data[data>quantile(data,q_2[state]) & data<quantile(data,q_2[state+1])]))) # kappa
    }
    params <- c(params,k)
  } else if (dist=='norm'){
    params <- c(params, as.numeric(quantile(data,probs=q_1[2:(N+1)]))) # mu
    s <- c()
    for (state in 1:N){
      s<-c(s,as.numeric(sd(data[data>quantile(data,q_2[state]) & data<quantile(data,q_2[state+1])]))) # sigma
    }
    params <- c(params,s)
  }
  return(params)
}


