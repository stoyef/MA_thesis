# 2022-06-08
# Functions to manipulate HMM parameters

#' Construct parameter object for AR(p)-HMM likelihood evaluation
#'
#' Construct a nice named list that contains all relevant parameters of an AR(p)-HMM.
#' This gets used in the log-likelihood computation to generalize to arbitrary 
#' distributions.
#' 
#' @param Gamma Full TPM.
#' @param autocor Vector or matrix of autoregression parameters in suitable format.
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
#' Transform the natural parameters of AR(p)-HMMs to their working parameters for optimization (theta -> theta.star).
#' Currently only works for distributions with exactly two parameters.
#' 
#' @param theta Vector of natural parameters.
#' @param N Number of states of the HMM.
#' @param p Vector of degree of autoregression within the distributions (one value for each state).
#' @param dists Vector of the distributions in the HMM.
#' @param scale_kappa Default 1, Scaling factor for kappa to avoid numerical issues in optimization for large kappa.
#' @param zero_inf Default FALSE, indicates if the gamma distributed variables should incorporate zero-inflation.
#' 
#' @return Vector of working parameters of the AR(p)-HMM.
#' 
#' @export
#' @rdname starize
starize <- function(theta,N,p,dists, scale_kappa=1, zero_inf=FALSE){
  # check for right number of parameters
  if (zero_inf){
    right_number = N*(N-1)+2*N*length(dists)+sum(p)+sum(dists=='gamma')*N
  } else{
    right_number = N*(N-1)+2*N*length(dists)+sum(p)
  }
  if (length(theta) != right_number){
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
      autocor <- qlogis(theta[param_count+1:sum(p)])
      theta.star <- c(Gamma, params, autocor)
    } else{
      theta.star <- c(Gamma, params)
    }
    if (zero_inf){
      zero_inf <- qlogis(theta[-(1:(param_count+sum(p)*N))])
      theta.star <- c(theta.star, zero_inf)
    }
    
    return(theta.star)
  }
}




#' Working parameters to natural parameters (theta.star->theta)
#'
#' Transform the working parameters of AR(p)-HMMs to their natural parameters for optimization (theta.star -> theta).
#' Currently only works for distributions with exactly two parameters.
#' 
#' @param theta.star Vector of working parameters.
#' @param N Number of states of the HMM.
#' @param p Vector of degree of autoregression within the distributions 
#'          (one value of each state for each distribution).
#' @param dists Vector of the distributions in the HMM.
#' @param scale_kappa Default 1, Scaling factor for kappa to avoid numerical issues in optimization for large kappa.
#' @param zero_inf Default FALSE, indicates if the gamma distributed variables should incorporate zero-inflation.
#' 
#' @return List of natural parameters of the AR(p)-HMM.
#' 
#' @export
#' @rdname unstarize
unstarize <- function(theta.star,N,p,dists, scale_kappa=1, zero_inf=FALSE){
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
  counter_zero_inf = 0
  params = list()
  for (dist in 1:n_dists){
    if (dists[dist]=='gamma'){
      if (zero_inf){
        params[[dist]] = list(mu=exp(theta.star[counter+1:N]),
                              sigma=exp(theta.star[counter+N+1:N]),
                              zero_inf=plogis(theta.star[(N*(N-1)+n_dist_params+sum(p)*N+counter_zero_inf+1:N)]))
        counter_zero_inf=counter_zero_inf+N # workaround for number of zero inflation params
      } else{
        params[[dist]] = list(mu=exp(theta.star[counter+1:N]),
                              sigma=exp(theta.star[counter+N+1:N]))
      }
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
    auto_params = theta.star[counter+1:(sum(p))] # no sum(p)*N anymore
    autocor = list()
    counter_param = 0
    counter_p = 1
    for (dist in 1:n_dists){
      autocor[[dist]] = list()
      for (state in 1:N){
        if (p[counter_p]>0){ # check if there is autocorrelation in the state
          autocor[[dist]][[state]] = plogis(auto_params[counter_param+1:p[counter_p]])
          counter_param = counter_param + p[counter_p]
          counter_p = counter_p+1
        } else{
          autocor[[dist]][[state]] = NA
          counter_p = counter_p+1
        }
      }
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


