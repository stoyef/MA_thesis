# 2022-08-16 
# Global minimization for the negative log likelihood

library(MasterThesis)
data = sample_arp(n_samples=1000, 
                  delta=c(0.5,0.5),
                  Gamma=matrix(c(0.9,0.1,0.1,0.9),ncol=2),
                  N=2,
                  params=c(20,40,2,3.5, 0,0,10,20),
                  autocor=list(matrix(c(0.1,0.1,0.1,0.1),ncol=2,byrow=TRUE),
                               matrix(c(0.1,0.1,0.1,0.1),ncol=2,byrow=TRUE)),
                  p=c(2,2),
                  dists=c('gamma','vm'))
data = sample_arp(n_samples=1000, 
                  delta=c(0.5,0.5),
                  Gamma=matrix(c(0.9,0.1,0.1,0.9),ncol=2),
                  N=2,
                  params=c(20,40,2,3.5, 0,0,10,20),
                  autocor=0,#list(matrix(c(0.1,0.1,0.1,0.1),ncol=2,byrow=TRUE),
                            #   matrix(c(0.1,0.1,0.1,0.1),ncol=2,byrow=TRUE)),
                  p=c(0,0),
                  dists=c('gamma','vm'))

theta = c(0.9,0.9,15,40,3,3,0,0,5,30,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1)
theta.star=starize(theta, N=2,p=c(2,2),dists=c('gamma','vm'))

mod_optim = fit_arp_model(mllk, data=data$data, theta.star = theta.star, N=2, p_auto=c(0,0),dists=c('gamma','vm'))
mod_optim

# nlm usually doesn't work
mod_nlm = fit_arp_model(mllk, data=data$data, theta.star = theta.star, N=2, p_auto=c(0,0),dists=c('gamma','vm'),
                        opt_fun = 'nlm')
mod_nlm


theta = c(0.9,0.9,15,40,3,3,0,0,5,25)
theta.star=starize(theta, N=2,p=c(0,0),dists=c('gamma','vm'))

# genoud works often but is slower
library(rgenoud)
mod <- genoud(fn=mllk, nvars=length(theta.star), pop.size=1000, max.generations=100,
              starting.values = theta.star,
              N=2,p_auto=c(0,0),
              dists=c('gamma','vm'),x=data$data)
unstarize(mod$par, N=2,p=c(0,0),dists=c('gamma','vm'))


