## 2022-09-13
# Code to prduce analyses of sea tern data for chapter 4.2

# set wd to current file directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # RStudio way

# read in data
muskox = read.table("../datasets/muskox/muskoxdata.txt",
                    header = TRUE)
head(muskox)
library(moveHMM)
muskox = prepData(muskox, type='UTM')
head(muskox)


## N=2
mod_move=fitHMM(muskox,
                2,
                c(25,150,20,200,rep(0.01,2)), # account for zeromass 
                c(pi,0,0.5,0.5),
                stationary = TRUE) 
plot(mod_move)
mod_move
pseudo_residuals = pseudoRes(mod_move)
qqnorm(pseudo_residuals$stepRes)



# no autocorrelation
theta=c(rep(-2,2),25,150,20,200,pi,0,0.5,0.5) 
theta.star=starize(theta, N=2,p=c(0,0),dists=c('gamma','vm'), scale_kappa = 1)
mod_ar0=fit_arp_model(mllk=mllk, data=cbind(muskox$step, muskox$angle), 
                      theta.star=theta.star, N=2, p_auto=c(0,0),dists=c('gamma','vm'),scale_kappa = 1)
mod_ar0
par(mfrow=c(1,2))
plot_fitted_dist(data=muskox$step,
                 dist='gamma', param=list(mu=mod_ar0$params[[1]]$mu,sigma=mod_ar0$params[[1]]$sigma),
                 N=2,delta=mod_ar0$delta, title='none',breaks=100,
                 xlim=c(0,500),legend=F,xlab='step length')
plot_fitted_dist(data=muskox$angle,
                 dist='vm', param=list(mu=mod_ar0$params[[2]]$mu,kappa=mod_ar0$params[[2]]$kappa),
                 N=2,delta=mod_ar0$delta, title='none',breaks=20,
                 legend=F,xlab='turning angle')
legend('topleft', c(paste("State",1:N),"Total"), bty='n', lwd=2,
       col=pal,inset=c(-0.95,0),xpd='NA',horiz=T)# N=3: inset=c(0.25,-0.1)


# AR(1)
theta=c(rep(-2,2),25,150,20,200,pi,0,0.5,0.5, rep(0.01,4)) 
theta.star=starize(theta, N=2,p=c(1,1),dists=c('gamma','vm'), scale_kappa = 1)
mod_ar1both=fit_arp_model(mllk=mllk, data=cbind(muskox$step, muskox$angle), 
                          theta.star=theta.star, N=2, p_auto=c(1,1),dists=c('gamma','vm'), scale_kappa = 1)
mod_ar1both
plot_fitted_dist(data=muskox$step,
                 dist='gamma', param=list(mu=mod_ar1both$params[[1]]$mu,sigma=mod_ar1both$params[[1]]$sigma),
                 N=2,delta=mod_ar1both$delta, title='none',breaks=100,
                 xlim=c(0,500),legend=F,xlab='step length')
plot_fitted_dist(data=muskox$angle,
                 dist='vm', param=list(mu=mod_ar1both$params[[2]]$mu,kappa=mod_ar1both$params[[2]]$kappa),
                 N=2,delta=mod_ar1both$delta, title='none',breaks=20,
                 legend=F,xlab='turning angle')
legend('topleft', c(paste("State",1:N),"Total"), bty='n', lwd=2,
       col=pal,inset=c(-0.95,0),xpd='NA',horiz=T)# N=3: inset=c(0.25,-0.1)
# Good results!!!! :D

# AR(2)
theta=c(rep(-2,2),10,150,8,165,pi,0,0.5,0.5, rep(0.1,8)) 
theta.star=starize(theta, N=2,p=c(2,2),dists=c('gamma','vm'), scale_kappa = 1)
mod_ar2both=fit_arp_model(mllk=mllk, data=cbind(muskox$step, muskox$angle), 
                          theta.star=theta.star, N=2, p_auto=c(2,2),dists=c('gamma','vm'), scale_kappa = 1)
mod_ar2both
plot_fitted_dist(data=muskox$step,
                 dist='gamma', param=list(mu=mod_ar2both$params[[1]]$mu,sigma=mod_ar2both$params[[1]]$sigma),
                 N=2,delta=mod_ar2both$delta, title='none',breaks=100,
                 xlim=c(0,500),legend=F,xlab='step length')
plot_fitted_dist(data=muskox$angle,
                 dist='vm', param=list(mu=mod_ar2both$params[[2]]$mu,kappa=mod_ar2both$params[[2]]$kappa),
                 N=2,delta=mod_ar2both$delta, title='none',breaks=20,
                 legend=F,xlab='turning angle')
legend('topleft', c(paste("State",1:N),"Total"), bty='n', lwd=2,
       col=pal,inset=c(-0.95,0),xpd='NA',horiz=T)# N=3: inset=c(0.25,-0.1)
# works as well


# All those models take some time to compute 
# For 50 models, parallelize!!
# t.b.c
