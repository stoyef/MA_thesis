# HMM WASP

# compute -logL (from lecture code)

mllk<-function(theta.star,x,N){
  Gamma <- diag(N)
  Gamma[!Gamma] <- exp(theta.star[1:((N-1)*N)])
  Gamma <- Gamma/rowSums(Gamma)
  delta <- solve(t(diag(N)-Gamma+1),rep(1,N))
  mu <- theta.star[(N-1)*N+1:N]
  sigma <- exp(theta.star[(N-1)*N+(N+1):(2*N)])
  allprobs <- matrix(1,length(x),N)
  ind <- which(!is.na(x))
  for (j in 1:N){
    allprobs[ind,j] <- dnorm(x[ind],mu[j],sigma[j]) # in this case Gaussian
  }
  foo <- delta%*%diag(allprobs[1,])
  l <- log(sum(foo))
  phi <- foo/sum(foo)
  for (t in 2:length(x)){
    foo <- phi%*%Gamma%*%diag(allprobs[t,])
    l <- l+log(sum(foo))
    phi <- foo/sum(foo)
  }
  return(-l)
}


###
### 1st aim: Clone model from Lawler et. al (2019)
###


