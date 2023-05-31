## 2022-09-14
# Pseudo residuals for a fitted AR(p)-HMM

#' Forward log-probabilities
#' 
#' Computes foward log-probabilities for an AR(p)-HMM. 
#' Adapted from moveHMM. 
#' Attention: Only works for 2-dim with 1st variable gamma distributed step lengths 
#' and 2nd variable von Mises distributed turning angles!
#' 
#' @param mod Fitted AR(p)-HMM object.
#' @param data Data the model is fitted on.
#' @param N Number of states.
#' @param p Vector of degree of autoregression for each distribution 
#'          (one value for each state) (0 = no autoregression).
#' 
#' @return Matrix of forward log-probabilities.
#' 
#' @export
#' @rdname logAlpha
logAlpha <- function(mod, data, N, p){
  nbObs = nrow(data)
  lalpha = matrix(NA,nbObs,N)
  delta = mod$delta
  trMat = mod$Gamma
  
  if(any(p>0)){
    probs = allprobs(x=data, dists=c('gamma','vm'), autocor=mod$autocorrelation, params=mod$params, N=N, p=p)
  } else{
    probs = allprobs(x=data, dists=c('gamma','vm'), autocor=0, params=mod$params, N=N, p=p)
  }
  foo <- delta*probs[1,]
  lscale = log(sum(foo))
  foo <- foo/sum(foo)
  lalpha[1,] <- log(foo)+lscale
  for (i in 2:nbObs){
    foo <- (foo%*%trMat)*probs[i,]
    lscale = lscale+log(sum(foo))
    foo <- foo/sum(foo)
    lalpha[i,] <- log(foo)+lscale
  }
  
  return(lalpha)
}



#' Compute pseudo residuals
#' 
#' Compute the pseudo residuals of an AR(p)-HMM.
#' Adapted from moveHMM.
#' Attention: Only works for 2-dim with 1st variable gamma distributed step lengths 
#' and 2nd variable von Mises distributed turning angles!
#' 
#' @param mod Fitted AR(p)-HMM object.
#' @param data Data the model is fitted on.
#' @param N Number of states.
#' @param p Vector of degree of autoregression for each distribution 
#'          (one value for every state) (0 = no autoregression).
#' 
#' @return Pseudo residuals for step length and turning angles.
#' 
#' @export
#' @rdname pseudores_arp
#' @import CircStats
pseudores_arp <- function(mod, data, N, p){
  require(CircStats)
  nbObs = nrow(data)
  stepfun = 'pgamma'
  anglefun = 'dvm'
  trMat = mod$Gamma
  autocor = mod$autocorrelation
  
  # forward log-probabilities
  la <- logAlpha(mod=mod,data=data,N=N,p=p)
  
  stepRes <- rep(NA,nbObs)
  pStepMat <- matrix(NA,nbObs,N)
  angleRes <- rep(NA,nbObs)
  pAngleMat <- matrix(NA,nbObs,N)
  
  if (any(p>0)){
    
    #autocor_m_step <- matrix(autocor[[1]], ncol=p[1], byrow=TRUE) # autoregression matrix for easier handling later on
    #autocor_m_angle <- matrix(autocor[[2]], ncol=p[2], byrow=TRUE) # autoregression matrix for easier handling later on
    
    x_wo_na_step = data[,1]
    x_wo_na_angle = data[,2]
    x_wo_na_step[which(is.na(x_wo_na_step))] = mean(x_wo_na_step,na.rm = TRUE)
    x_wo_na_angle[which(is.na(x_wo_na_angle))] = mean(x_wo_na_angle,na.rm = TRUE)
    # replace NA values in x that are put into autocor_ind with mean value 
    # (so that the data that has NA in previous time steps does not have to be deleted)
  }
  for(state in 1:N) {
    # define lists of parameters -> global step
    stepArgs <- list(data[1,1])
    stepArgs[[2]] = mod$params[[1]][[1]][state]
    stepArgs[[3]] = mod$params[[1]][[2]][state]
    zeromass = 0
    cv = stepArgs[[3]] / stepArgs[[2]]
    
    angleArgs <- list(anglefun,-pi,data[1,2]) # to pass to function "integrate" below
    angleArgs[[4]] =  mod$params[[2]][[1]][state]
    angleArgs[[5]] =  mod$params[[2]][[2]][state]
    
    for(i in 1:nbObs) {
      # steps
      if(!is.na(data[i,1])) {
        if (p[state]>0){
          if (i > p[state]){
            mu_auto <- (1-sum(autocor[[1]][[state]]))*stepArgs[[2]] + 
              as.vector(x_wo_na_step[(i-p[state]):(i-1)]%*%autocor[[1]][[state]])
            sigma_auto <- cv*mu_auto
            pStepMat[i,state] <- zeromass+(1-zeromass)*pgamma(data[i,1],
                                                              shape=mu_auto^2/sigma_auto^2,
                                                              scale=sigma_auto^2/mu_auto)
          } else{
            pStepMat[i,state] <- zeromass+(1-zeromass)*pgamma(data[i,1],
                                                              shape=stepArgs[[2]]^2/stepArgs[[3]]^2,
                                                              scale=stepArgs[[3]]^2/stepArgs[[2]])
          }
        } else{
          pStepMat[i,state] <- zeromass+(1-zeromass)*pgamma(data[i,1],
                                                            shape=stepArgs[[2]]^2/stepArgs[[3]]^2,
                                                            scale=stepArgs[[3]]^2/stepArgs[[2]])
        }
      }
      
      # angles
      if(!is.na(data[i,2])) {
        # angle==pi => residual=Inf
        if (p[N+state]>0){
          if (i > p[N+state]){
            if(data[i,2]!=pi) {
              mu_auto <- Arg(
                (1-sum(autocor[[2]][[state]]))*exp(1i*angleArgs[[4]]) + 
                  as.vector(exp(1i*x_wo_na_angle[(i-p[N+state]):(i-1)])%*%autocor[[2]][[state]])
              )
              pAngleMat[i,state] <- integrate(dvm,angleArgs[[2]],data[i,2],
                                              mu_auto,angleArgs[[5]])$value
            }
          } else{
            if(data[i,2]!=pi) {
              pAngleMat[i,state] <- integrate(dvm,angleArgs[[2]],data[i,2],
                                              angleArgs[[4]],angleArgs[[5]])$value
            }
          }
        } else{
          if(data[i,2]!=pi) {
            pAngleMat[i,state] <- integrate(dvm,angleArgs[[2]],data[i,2],
                                            angleArgs[[4]],angleArgs[[5]])$value
          }
        }
      }
    }
  }
  
  # pseudores
  # 1st observation
  stepRes[1] <- qnorm(t(mod$delta)%*%pStepMat[1,])
  if (!is.na(data[1,2])){
    angleRes[1] <- qnorm(t(mod$delta)%*%pAngleMat[1,])
  }
  for (i in 2:nbObs){
    if(!is.na(data[i,1])){
      c = max(la[i-1,])
      a = exp(la[i-1,]-c)
      stepRes[i] <- qnorm(t(a)%*%(trMat/sum(a))%*%pStepMat[i,])
      
      if(!is.na(data[i,2])){
        angleRes[i] <- qnorm(t(a)%*%(trMat/sum(a))%*%pAngleMat[i,])
      }
    }
  }
  
  return(list(stepRes=stepRes, angleRes=angleRes))
}

