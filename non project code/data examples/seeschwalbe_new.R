# 2022-08-15
# Analyzing seeschwalbe data using package


#setwd(getSrcDirectory()[1]) # set wd to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # RStudio way

source("../tsa_functions_schwalbe.R")

schwalben = read_schwalbe("/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/MA_thesis/datasets/seeschwalbe/slick2-h1-h2.csv")

head(schwalben)

count=1
for (i in unique(schwalben$ID)){
  assign(paste("schwalbe_", count, sep=""), select_schwalbe(i, count))
  count=count+1
}


# die 3. nehmen wir

library(moveHMM)

schwalbe_3 = prepData(schwalbe_3)
head(schwalbe_3)
plot(schwalbe_3$x, schwalbe_3$y, type='l')

schwalbe_80 = prepData(schwalbe_80)
head(schwalbe_80)
plot(schwalbe_80$x, schwalbe_80$y, type='l')


library(MasterThesis)

mod_move=fitHMM(schwalbe_80,
       2,
       c(15,30,10,5),
       c(0,0,500,1000))
plot(mod_move)

# no autocorrelation
theta=c(0.90,0.90,15,30,9,3,0,0,500,1000)
theta.star=starize(theta, N=2,p=c(0,0),dists=c('gamma','vm'))
mod=fit_arp_model(mllk=mllk, data=cbind(schwalbe_80$step, schwalbe_80$angle), 
              theta.star=theta.star, N=2, p=c(0,0),dists=c('gamma','vm'))
mod
plot_fitted_dist(data=schwalbe_80$step,
                 dist='gamma', param=list(mu=mod$params[[1]]$mu,sigma=mod$params[[1]]$sigma),
                 N=2,delta=mod$delta, title='schwalbe_80 w/o autocorrelation')
plot_fitted_dist(data=schwalbe_80$angle,
                 dist='vm', param=list(mu=mod$params[[2]]$mu,kappa=mod$params[[2]]$kappa),
                 N=2,delta=mod$delta, title='schwalbe_80 w/o autocorrelation')

# AR(1) in step length (not turning angle)
theta=c(0.9,0.9,15,28,7,3,0,0,500,1000,0.3,0.3)
theta.star=starize(theta, N=2,p=c(1,0),dists=c('gamma','vm'))
mod=fit_arp_model(mllk=mllk, data=cbind(schwalbe_80$step, schwalbe_80$angle), 
              theta.star=theta.star, N=2, p=c(1,0),dists=c('gamma','vm'))
mod
plot_fitted_dist(data=schwalbe_80$step,
                 dist='gamma', param=list(mu=mod$params[[1]]$mu,sigma=mod$params[[1]]$sigma),
                 N=2,delta=mod$delta, title='schwalbe_80 w/o autocorrelation')
plot_fitted_dist(data=schwalbe_80$angle,
                 dist='vm', param=list(mu=mod$params[[2]]$mu,kappa=mod$params[[2]]$kappa),
                 N=2,delta=mod$delta, title='schwalbe_80 w/o autocorrelation')

# AR(1) in turning angle (not step length)
theta=c(0.9,0.9,15,30,7,3,0,0,500,1000,0.3,0.3)
theta.star=starize(theta, N=2,p=c(0,1),dists=c('gamma','vm'))
mod=fit_arp_model(mllk=mllk, data=cbind(schwalbe_80$step, schwalbe_80$angle), 
              theta.star=theta.star, N=2, p=c(0,1),dists=c('gamma','vm'))
mod
plot_fitted_dist(data=schwalbe_80$step,
                 dist='gamma', param=list(mu=mod$params[[1]]$mu,sigma=mod$params[[1]]$sigma),
                 N=2,delta=mod$delta, title='schwalbe_80 w/o autocorrelation')
plot_fitted_dist(data=schwalbe_80$angle,
                 dist='vm', param=list(mu=mod$params[[2]]$mu,kappa=mod$params[[2]]$kappa),
                 N=2,delta=mod$delta, title='schwalbe_80 w/o autocorrelation')

# AR(1)
theta=c(0.9,0.9,15,30,7,3,0,0,500,1000,0.2,0.3,0.2,0.3)
theta.star=starize(theta, N=2,p=c(1,1),dists=c('gamma','vm'))
mod=fit_arp_model(mllk=mllk, data=cbind(schwalbe_80$step, schwalbe_80$angle), 
              theta.star=theta.star, N=2, p=c(1,1),dists=c('gamma','vm'))
mod
plot_fitted_dist(data=schwalbe_80$step,
                 dist='gamma', param=list(mu=mod$params[[1]]$mu,sigma=mod$params[[1]]$sigma),
                 N=2,delta=mod$delta, title='schwalbe_80 w/o autocorrelation')
plot_fitted_dist(data=schwalbe_80$angle,
                 dist='vm', param=list(mu=mod$params[[2]]$mu,kappa=mod$params[[2]]$kappa),
                 N=2,delta=mod$delta, title='schwalbe_80 w/o autocorrelation')

# AR(2)
theta=c(0.9,0.9,15,30,20,10,0,0,500,1500,0.2,0.1,0.2,0.1,0.05,0.05,0.05,0.05)
theta.star=starize(theta, N=2,p=c(2,2),dists=c('gamma','vm'))
mod=fit_arp_model(mllk=mllk, data=cbind(schwalbe_80$step, schwalbe_80$angle), 
              theta.star=theta.star, N=2, p=c(2,2),dists=c('gamma','vm'))
mod
plot_fitted_dist(data=schwalbe_80$step,
                 dist='gamma', param=list(mu=mod$params[[1]]$mu,sigma=mod$params[[1]]$sigma),
                 N=2,delta=mod$delta, title='schwalbe_80 w/o autocorrelation')
plot_fitted_dist(data=schwalbe_80$angle,
                 dist='vm', param=list(mu=mod$params[[2]]$mu,kappa=mod$params[[2]]$kappa),
                 N=2,delta=mod$delta, title='schwalbe_80 w/o autocorrelation')

