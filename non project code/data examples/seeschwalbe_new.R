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


# die 77. nehmen wir

library(moveHMM)

schwalbe_3 = prepData(schwalbe_3)
head(schwalbe_3)
plot(schwalbe_3$x, schwalbe_3$y, type='l')

schwalbe_80 = prepData(schwalbe_80)
head(schwalbe_80)
plot(schwalbe_80$x, schwalbe_80$y, type='l')

schwalbe_77 = prepData(schwalbe_77)
head(schwalbe_77)
plot(schwalbe_77$x, schwalbe_77$y, type='l')


library(MasterThesis)

mod_move=fitHMM(schwalbe_77,
       3,
       c(15,30,35,9,4,3),
       c(0,0,0,20,50,200))
plot(mod_move)

# Modification for seeschwalbe:
# change resolution of kappa parameter in von Mises distribution to avoid numerical 
# issues in optimization function
# -> theta.star for kappa in input of optim function should be smaller

# no autocorrelation
theta=c(rep(-2,6),15,25,30,9,4,3,0,0,0,200,500,1000)
theta.star=starize(theta, N=3,p=c(0,0),dists=c('gamma','vm'), scale_kappa = 20)
mod_ar0=fit_arp_model(mllk=mllk, data=cbind(schwalbe_77$step, schwalbe_77$angle), 
              theta.star=theta.star, N=3, p_auto=c(0,0),dists=c('gamma','vm'),scale_kappa = 20)
mod_ar0
plot_fitted_dist(data=schwalbe_77$step,
                 dist='gamma', param=list(mu=mod_ar0$params[[1]]$mu,sigma=mod_ar0$params[[1]]$sigma),
                 N=3,delta=mod_ar0$delta, title='schwalbe_77 w/o autocorrelation',breaks=100)
plot_fitted_dist(data=schwalbe_77$angle,
                 dist='vm', param=list(mu=mod_ar0$params[[2]]$mu,kappa=mod_ar0$params[[2]]$kappa),
                 N=3,delta=mod_ar0$delta, title='schwalbe_77 w/o autocorrelation',breaks=100)

# nur turning angle
theta=c(rep(-2,6),0,0,0,100,300,800)
theta.star=starize(theta, N=3,p=c(0),dists=c('vm'),scale_kappa = 10)
mod=fit_arp_model(mllk=mllk, data=schwalbe_77$angle, 
                  theta.star=theta.star, N=3, p_auto=c(0),dists=c('vm'),scale_kappa = 10)
mod
plot_fitted_dist(data=schwalbe_77$angle,
                 dist='vm', param=list(mu=mod$params[[1]]$mu,kappa=mod$params[[1]]$kappa),
                 N=3,delta=mod$delta, title='schwalbe_77 w/ AR(1)',breaks=200)


# AR(1) in step length (not turning angle) -> doesn't work
theta=c(rep(-2,6),15,25,30,9,4,3,0,0,0,100,300,800,0.1,0.1,0.1)
theta.star=starize(theta, N=3,p=c(1,0),dists=c('gamma','vm'))
mod=fit_arp_model(mllk=mllk, data=cbind(schwalbe_77$step, schwalbe_77$angle), 
              theta.star=theta.star, N=3, p_auto=c(1,0),dists=c('gamma','vm'))
mod
plot_fitted_dist(data=schwalbe_77$step,
                 dist='gamma', param=list(mu=mod$params[[1]]$mu,sigma=mod$params[[1]]$sigma),
                 N=2,delta=mod$delta, title='schwalbe_80 w/o autocorrelation')
plot_fitted_dist(data=schwalbe_77$angle,
                 dist='vm', param=list(mu=mod$params[[2]]$mu,kappa=mod$params[[2]]$kappa),
                 N=2,delta=mod$delta, title='schwalbe_80 w/o autocorrelation')

# AR(1) in turning angle (not step length)
theta=c(rep(-2,6),15,25,30,9,4,3,0,0,0,200,500,1000,0.1,0.1,0.1)
theta.star=starize(theta, N=3,p=c(0,1),dists=c('gamma','vm'))
mod_ar1turn=fit_arp_model(mllk=mllk, data=cbind(schwalbe_77$step, schwalbe_77$angle), 
              theta.star=theta.star, N=3, p_auto=c(0,1),dists=c('gamma','vm'))
mod_ar1turn
# works
plot_fitted_dist(data=schwalbe_77$step,
                 dist='gamma', param=list(mu=mod_ar1turn$params[[1]]$mu,sigma=mod_ar1turn$params[[1]]$sigma),
                 N=3,delta=mod_ar1turn$delta, title='schwalbe_77 w/ AR(1)')
plot_fitted_dist(data=schwalbe_77$angle,
                 dist='vm', param=list(mu=mod_ar1turn$params[[2]]$mu,kappa=mod_ar1turn$params[[2]]$kappa),
                 N=3,delta=mod_ar1turn$delta, title='schwalbe_77 w/ AR(1)')

# AR(1)
theta=c(rep(-2,6),15,25,30,9,4,3,0,0,0,200,500,1000,rep(0.1,6))
theta.star=starize(theta, N=3,p=c(1,1),dists=c('gamma','vm'), scale_kappa = 1)
mod_ar1both=fit_arp_model(mllk=mllk, data=cbind(schwalbe_77$step, schwalbe_77$angle), 
              theta.star=theta.star, N=3, p_auto=c(1,1),dists=c('gamma','vm'), scale_kappa = 1)
mod_ar1both
plot_fitted_dist(data=schwalbe_77$step,
                 dist='gamma', param=list(mu=mod_ar1both$params[[1]]$mu,sigma=mod_ar1both$params[[1]]$sigma),
                 N=3,delta=mod_ar1both$delta, title='schwalbe_77 w/o autocorrelation')
plot_fitted_dist(data=schwalbe_77$angle,
                 dist='vm', param=list(mu=mod_ar1both$params[[2]]$mu,kappa=mod_ar1both$params[[2]]$kappa),
                 N=3,delta=mod$delta, title='schwalbe_77 w/o autocorrelation')
## AR(1) in step length doesn't work


# AR(2) in turning angle (not step length)
theta=c(rep(-2,6),15,25,30,9,4,3,0,0,0,200,500,1000,rep(0.1,6))
theta.star=starize(theta, N=3,p=c(0,2),dists=c('gamma','vm'))
mod_ar2turn=fit_arp_model(mllk=mllk, data=cbind(schwalbe_77$step, schwalbe_77$angle), 
                          theta.star=theta.star, N=3, p_auto=c(0,2),dists=c('gamma','vm'))
mod_ar2turn
# works
plot_fitted_dist(data=schwalbe_77$step,
                 dist='gamma', param=list(mu=mod_ar2turn$params[[1]]$mu,sigma=mod_ar2turn$params[[1]]$sigma),
                 N=3,delta=mod_ar2turn$delta, title='schwalbe_77 w/ AR(2) in turning angle')
plot_fitted_dist(data=schwalbe_77$angle,
                 dist='vm', param=list(mu=mod_ar2turn$params[[2]]$mu,kappa=mod_ar2turn$params[[2]]$kappa),
                 N=3,delta=mod_ar2turn$delta, title='schwalbe_77 w/ AR(2) in turning angle')


# AR(2)
theta=c(0.9,0.9,15,30,20,10,0,0,500,1500,0.2,0.1,0.2,0.1,0.05,0.05,0.05,0.05)
theta.star=starize(theta, N=2,p=c(2,2),dists=c('gamma','vm'))
mod=fit_arp_model(mllk=mllk, data=cbind(schwalbe_80$step, schwalbe_80$angle), 
              theta.star=theta.star, N=2, p_auto=c(2,2),dists=c('gamma','vm'))
mod
plot_fitted_dist(data=schwalbe_80$step,
                 dist='gamma', param=list(mu=mod$params[[1]]$mu,sigma=mod$params[[1]]$sigma),
                 N=2,delta=mod$delta, title='schwalbe_80 w/o autocorrelation')
plot_fitted_dist(data=schwalbe_80$angle,
                 dist='vm', param=list(mu=mod$params[[2]]$mu,kappa=mod$params[[2]]$kappa),
                 N=2,delta=mod$delta, title='schwalbe_80 w/o autocorrelation')



##
## Tortuosity
##

plot(schwalbe_77$tortuosity,type='l')

## 2 states

# AR(0)
theta=c(rep(-2,2),1.05,1.5,1,3)
theta.star=starize(theta, N=2,p=c(0),dists=c('gamma'))
mod_ar0=fit_arp_model(mllk=mllk, data=schwalbe_77$tortuosity, 
                      theta.star=theta.star, N=2, p_auto=0,dists=c('gamma'))
mod_ar0
plot_fitted_dist(data=schwalbe_77$tortuosity,
                 dist='gamma', param=list(mu=mod_ar0$params[[1]]$mu,sigma=mod_ar0$params[[1]]$sigma),
                 N=2,delta=mod_ar0$delta, title='schwalbe_77 tortuosity AR(0)')

# AR(1)
theta=c(rep(-2,2),1.05,1.5,1,3,rep(0.2,2))
theta.star=starize(theta, N=2,p=c(1),dists=c('gamma'))
mod_ar1=fit_arp_model(mllk=mllk, data=schwalbe_77$tortuosity, 
                      theta.star=theta.star, N=2, p_auto=1,dists=c('gamma'))
mod_ar1
plot_fitted_dist(data=schwalbe_77$tortuosity,
                 dist='gamma', param=list(mu=mod_ar1$params[[1]]$mu,sigma=mod_ar1$params[[1]]$sigma),
                 N=2,delta=mod_ar1$delta, title='schwalbe_77 tortuosity AR(1)')

# AR(2)
theta=c(rep(-2,2),1.05,1.5,1,3,rep(0.1,4))
theta.star=starize(theta, N=2,p=c(2),dists=c('gamma'))
mod_ar2=fit_arp_model(mllk=mllk, data=schwalbe_77$tortuosity, 
                      theta.star=theta.star, N=2, p_auto=2,dists=c('gamma'))
mod_ar2
plot_fitted_dist(data=schwalbe_77$tortuosity,
                 dist='gamma', param=list(mu=mod_ar2$params[[1]]$mu,sigma=mod_ar2$params[[1]]$sigma),
                 N=2,delta=mod_ar2$delta, title='schwalbe_77 tortuosity AR(2)')

# AR(3) -> doesn't work 
theta=c(rep(-2,2),1.05,1.4,1.1,1.6,rep(0.2,6))
theta.star=starize(theta, N=2,p=c(3),dists=c('gamma'))
mod_ar3=fit_arp_model(mllk=mllk, data=schwalbe_77$tortuosity, 
                      theta.star=theta.star, N=2, p_auto=3,dists=c('gamma'))
mod_ar3
plot_fitted_dist(data=schwalbe_77$tortuosity,
                 dist='gamma', param=list(mu=mod_ar3$params[[1]]$mu,sigma=mod_ar3$params[[1]]$sigma),
                 N=2,delta=mod_ar3$delta, title='schwalbe_77 tortuosity AR(3)')


##
## Next try: downsampling to 1Hz
##

# read in again, to get raw data
schwalben = read_schwalbe("/Users/stoye/sciebo/Studium/31-M-Thesis Master's Thesis/MA_thesis/datasets/seeschwalbe/slick2-h1-h2.csv")

head(schwalben)

count=1
for (i in unique(schwalben$ID)){
  assign(paste("schwalbe_", count, sep=""), select_schwalbe(i, count))
  count=count+1
}

head(schwalbe_77)
head(schwalbe_77$t,40) # 30Hz Daten
library(zoo)
schwalbe_77_1hz = rollapply(schwalbe_77[,2:dim(schwalbe_77)[2]], 6, mean)
head(schwalbe_77_1hz)
schwalbe_77_1hz = as.data.frame(schwalbe_77_1hz[seq(1, nrow(schwalbe_77_1hz), 6),]) # make data 5hz
plot(schwalbe_77_1hz$x,schwalbe_77_1hz$y,type='l')
schwalbe_77_1hz = prepData(schwalbe_77_1hz)
plot(schwalbe_77_1hz)

# this sis better to handle

##
## N=2
##
mod_move=fitHMM(schwalbe_77_1hz,
                2,
                c(80,175,30,10), # 5Hz
                c(0,0,7,150)) # 5Hz
plot(mod_move)
mod_move
pseudo_residuals = pseudoRes(mod_move)
qqnorm(pseudo_residuals$stepRes)
qqnorm(pseudo_residuals$angleRes)



# no autocorrelation
#theta=c(rep(-2,6),350,700,850,300,200,100,0,0,0,1,2,6) # 1Hz
theta=c(rep(-2,2),80,175,30,10,0,0,7,150) # 5 Hz
theta.star=starize(theta, N=2,p=c(0,0),dists=c('gamma','vm'), scale_kappa = 5)
mod_ar0=fit_arp_model(mllk=mllk, data=cbind(schwalbe_77_1hz$step, schwalbe_77_1hz$angle), 
                      theta.star=theta.star, N=2, p_auto=c(0,0),dists=c('gamma','vm'),scale_kappa = 5)
mod_ar0
plot_fitted_dist(data=schwalbe_77_1hz$step,
                 dist='gamma', param=list(mu=mod_ar0$params[[1]]$mu,sigma=mod_ar0$params[[1]]$sigma),
                 N=2,delta=mod_ar0$delta, title='schwalbe_77 w/o autocorrelation',breaks=20)
plot_fitted_dist(data=schwalbe_77_1hz$angle,
                 dist='vm', param=list(mu=mod_ar0$params[[2]]$mu,kappa=mod_ar0$params[[2]]$kappa),
                 N=2,delta=mod_ar0$delta, title='schwalbe_77 w/o autocorrelation',breaks=20)

# AR(1)
#theta=c(rep(-2,6),350,700,850,300,200,100,0,0,0,1,2,6, rep(0.1,6)) # 1Hz
theta=c(rep(-2,2),80,160,40,15,0,0,7,100, rep(0.1,4)) # 5 Hz
theta.star=starize(theta, N=2,p=c(1,1),dists=c('gamma','vm'), scale_kappa = 5)
mod_ar1both=fit_arp_model(mllk=mllk, data=cbind(schwalbe_77_1hz$step, schwalbe_77_1hz$angle), 
                          theta.star=theta.star, N=2, p_auto=c(1,1),dists=c('gamma','vm'), scale_kappa = 5)
mod_ar1both
plot_fitted_dist(data=schwalbe_77_1hz$step,
                 dist='gamma', param=list(mu=mod_ar1both$params[[1]]$mu,sigma=mod_ar1both$params[[1]]$sigma),
                 N=2,delta=mod_ar1both$delta, title='schwalbe_77 w/ AR(1)', breaks=20)
plot_fitted_dist(data=schwalbe_77_1hz$angle,
                 dist='vm', param=list(mu=mod_ar1both$params[[2]]$mu,kappa=mod_ar1both$params[[2]]$kappa),
                 N=2,delta=mod$delta, title='schwalbe_77 w/ AR(1)', breaks=20)
# 5Hz -> AR(1) for 2 states doesn't work (phi=1)

# AR(2)
# theta=c(rep(-2,6),350,700,850,300,200,100,0,0,0,1,2,6, rep(0.1,12)) # 1Hz
theta=c(rep(-2,2),80,160,40,15,0,0,7,100, rep(0.1,8)) # 5 Hz
theta.star=starize(theta, N=2,p=c(2,2),dists=c('gamma','vm'), scale_kappa = 10)
mod_ar2both=fit_arp_model(mllk=mllk, data=cbind(schwalbe_77_1hz$step, schwalbe_77_1hz$angle), 
                          theta.star=theta.star, N=2, p_auto=c(2,2),dists=c('gamma','vm'), scale_kappa = 10)
mod_ar2both
plot_fitted_dist(data=schwalbe_77_1hz$step,
                 dist='gamma', param=list(mu=mod_ar2both$params[[1]]$mu,sigma=mod_ar2both$params[[1]]$sigma),
                 N=2,delta=mod_ar2both$delta, title='schwalbe_77 w/ AR(2)', breaks=20)
plot_fitted_dist(data=schwalbe_77_1hz$angle,
                 dist='vm', param=list(mu=mod_ar2both$params[[2]]$mu,kappa=mod_ar2both$params[[2]]$kappa),
                 N=2,delta=mod$delta, title='schwalbe_77 w/ AR(2)', breaks=20)
# doesn't work





##
## N=3
##
mod_move=fitHMM(schwalbe_77_1hz,
                3,
                #c(350,700,850,300,200,100), # 1Hz
                c(80,150,175,30,10,10), # 5Hz
                #c(0,0,0,1,2,6)) # 1Hz
                c(0,0,0,7,60,150)) # 5Hz
plot(mod_move)
mod_move
pseudo_residuals = pseudoRes(mod_move)
qqnorm(pseudo_residuals$stepRes)
qqnorm(pseudo_residuals$angleRes)


# no autocorrelation
#theta=c(rep(-2,6),350,700,850,300,200,100,0,0,0,1,2,6) # 1Hz
theta=c(rep(-2,6),80,150,175,30,10,10,0,0,0,7,60,150) # 5 Hz
theta.star=starize(theta, N=3,p=c(0,0),dists=c('gamma','vm'), scale_kappa = 5)
mod_ar0=fit_arp_model(mllk=mllk, data=cbind(schwalbe_77_1hz$step, schwalbe_77_1hz$angle), 
                      theta.star=theta.star, N=3, p_auto=c(0,0),dists=c('gamma','vm'),scale_kappa = 5)
mod_ar0
plot_fitted_dist(data=schwalbe_77_1hz$step,
                 dist='gamma', param=list(mu=mod_ar0$params[[1]]$mu,sigma=mod_ar0$params[[1]]$sigma),
                 N=3,delta=mod_ar0$delta, title='schwalbe_77 w/o autocorrelation',breaks=20)
plot_fitted_dist(data=schwalbe_77_1hz$angle,
                 dist='vm', param=list(mu=mod_ar0$params[[2]]$mu,kappa=mod_ar0$params[[2]]$kappa),
                 N=3,delta=mod_ar0$delta, title='schwalbe_77 w/o autocorrelation',breaks=20)

# AR(1)
#theta=c(rep(-2,6),350,700,850,300,200,100,0,0,0,1,2,6, rep(0.1,6)) # 1Hz
theta=c(rep(-2,6),80,150,175,30,10,10,0,0,0,7,60,150, rep(0.1,6)) # 5 Hz
theta.star=starize(theta, N=3,p=c(1,1),dists=c('gamma','vm'), scale_kappa = 5)
mod_ar1both=fit_arp_model(mllk=mllk, data=cbind(schwalbe_77_1hz$step, schwalbe_77_1hz$angle), 
                          theta.star=theta.star, N=3, p_auto=c(1,1),dists=c('gamma','vm'), scale_kappa = 5)
mod_ar1both
plot_fitted_dist(data=schwalbe_77_1hz$step,
                 dist='gamma', param=list(mu=mod_ar1both$params[[1]]$mu,sigma=mod_ar1both$params[[1]]$sigma),
                 N=3,delta=mod_ar1both$delta, title='schwalbe_77 w/ AR(1)', breaks=20)
plot_fitted_dist(data=schwalbe_77_1hz$angle,
                 dist='vm', param=list(mu=mod_ar1both$params[[2]]$mu,kappa=mod_ar1both$params[[2]]$kappa),
                 N=3,delta=mod$delta, title='schwalbe_77 w/ AR(1)', breaks=20)
# 1Hz -> AR(1) for 3 states leads to overfitting
# 5Hz -> AR(1) for 3 states doesn't work

# AR(2)
# theta=c(rep(-2,6),350,700,850,300,200,100,0,0,0,1,2,6, rep(0.1,12)) # 1Hz
theta=c(rep(-2,6),80,150,175,30,10,10,0,0,0,7,60,150, rep(0.1,12)) # 5Hz
theta.star=starize(theta, N=3,p=c(2,2),dists=c('gamma','vm'), scale_kappa = 10)
mod_ar2both=fit_arp_model(mllk=mllk, data=cbind(schwalbe_77_1hz$step, schwalbe_77_1hz$angle), 
                          theta.star=theta.star, N=3, p_auto=c(2,2),dists=c('gamma','vm'), scale_kappa = 10)
mod_ar2both
plot_fitted_dist(data=schwalbe_77_1hz$step,
                 dist='gamma', param=list(mu=mod_ar2both$params[[1]]$mu,sigma=mod_ar2both$params[[1]]$sigma),
                 N=3,delta=mod_ar2both$delta, title='schwalbe_77 w/ AR(2)', breaks=20)
plot_fitted_dist(data=schwalbe_77_1hz$angle,
                 dist='vm', param=list(mu=mod_ar2both$params[[2]]$mu,kappa=mod_ar2both$params[[2]]$kappa),
                 N=3,delta=mod$delta, title='schwalbe_77 w/ AR(2)', breaks=20)
# doesn't really work



