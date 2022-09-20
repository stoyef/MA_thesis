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

layout(matrix(c(1,3,3, 2,3,3),byrow=TRUE,ncol=3))
hist(muskox$step, breaks=25,prob=T,main="", xlab='step length', cex.lab=1.25)
hist(muskox$angle, breaks=15,prob=T,main="", xlab='turning angle', cex.lab=1.25)
plot(muskox$x,muskox$y, type='l', bty='n',xlab='x',ylab='y', cex.lab=1.25)



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
curve(1*x,-4,4,add=T)


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
theta=c(rep(-2,2),25,150,20,200,pi,0,0.5,0.5, rep(0.1,4)) 
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



library(parallel)
n_cores = detectCores()

N=3
p=c(0,0)
muskox_wrap <- function(iteration){
  # starting values - N=2
  #step_init = c(runif(1,15,30),runif(1,100,200),
  #              runif(1,10,30),runif(1,100,300))
  #angle_init = c(runif(2,-pi,pi),
  #               runif(1,0.35,0.6),runif(1,0.4,0.65))
  # N=3
  step_init = c(runif(1,5,10),runif(1,40,70),runif(1,150,300),
                runif(1,2,5),runif(1,30,70),runif(1,100,300))
  angle_init = c(runif(3,-pi,pi),
                 runif(3,0.2,0.8))
  theta = c(rep(-2,N*(N-1)),
            step_init,
            angle_init,
            rep(0.1,p[1]*N+p[2]*N))
  # working parameters
  theta.star = starize(theta, N=N, p=p, dists=c('gamma','vm'),scale_kappa = 1)
  begin = Sys.time()
  mod = fit_arp_model(mllk, data=cbind(muskox$step, muskox$angle), theta.star = theta.star, N=N, 
                      p_auto=p, dists=c('gamma','vm'), 
                      scale_kappa = 1)
  res = list()
  time = Sys.time()-begin
  res['time'] = time
  res[['mod']] = mod
  res['aic'] = 2*mod$mllk_optim + 2*length(theta)
  res['bic'] = 2*mod$mllk_optim + log(dim(muskox)[1])*length(theta)
  res['mllk'] = mod$mllk_optim
  return(res)
}

iterations = 1:10
results_ar00 <- mclapply(iterations,
                         muskox_wrap,
                         mc.cores = n_cores
)
#results_ar00
#work=c(1,2,3,4,5,6,7,10)
#work=c(work,work+10,work+20,work+30,work+40)
#work

#(1:50)[-work]


comps = rep(NA, length(iterations))
aics = rep(NA, length(iterations))
bics = rep(NA, length(iterations))
mllks = rep(NA, length(iterations))
mods = list()
c=1
for (i in iterations){
  mods[[i]] = results_ar00[i][[1]]$mod
  comps[i] = results_ar00[i][[1]]$time
  aics[i] = results_ar00[i][[1]]$aic
  bics[i] = results_ar00[i][[1]]$bic
  mllks[i] = results_ar00[i][[1]]$mllk
  c=c+1
}

comps
aics
bics
mllks
mods
best_mod_00 = mods[[which.min(bics)]]
best_mod_11 = mods[[which.min(bics)]]


par(mfrow=c(1,2))
plot_fitted_dist(data=muskox$step,
                 dist='gamma', param=list(mu=mods[[4]]$params[[1]]$mu,sigma=mods[[4]]$params[[1]]$sigma),
                 N=3,delta=mods[[4]]$delta, title='none',breaks=100,
                 xlim=c(0,500),legend=F,xlab='step length')
plot_fitted_dist(data=muskox$angle,
                 dist='vm', param=list(mu=mods[[4]]$params[[2]]$mu,kappa=mods[[4]]$params[[2]]$kappa),
                 N=3,delta=mods[[4]]$delta, title='none',breaks=20,
                 legend=F,xlab='turning angle')
pal=brewer.pal(N+1,'Dark2')
legend('topleft', c(paste("State",1:N),"Total"), bty='n', lwd=2,
       col=pal,inset=c(-1.2,#-0.95,
                       -0.1),xpd='NA',horiz=T)



# state decoding

# N=2
# AR(0,0), best model
states_00 = viterbi_arp(x=cbind(muskox$step, muskox$angle),
                        Gamma=matrix(c(0.857,0.143,0.160,0.840),byrow=T,ncol=2),
                        delta=c(0.528,0.472),
                        dists=c('gamma','vm'),
                        autocor=0,
                        params=list(list(mu=c(21.037, 183.914), sigma=c(21.862, 193.879)),
                                    list(mu=c(3.094, -0.060), kappa=c(0.414, 0.503))),
                        N=2,
                        p=c(0,0))
states_00

# AR(1,1), best model
states_11 = viterbi_arp(x=cbind(muskox$step, muskox$angle),
                        Gamma=matrix(c(0.749,0.251,0.200,0.800),byrow=T,ncol=2),
                        delta=c(0.443,0.557),
                        dists=c('gamma','vm'),
                        autocor=list(c(0.663,0.412),c(0.269,0.011)),
                        params=list(list(mu=c(10.965, 154.827), sigma=c(9.001, 173.828)),
                                    list(mu=c(3.119, -0.068), kappa=c(0.567, 0.458))),
                        N=2,
                        p=c(1,1))
states_11

par(mfrow=c(2,1))
plot(muskox$step[1:300],type='l',bty='n',xlab='time',ylab='step length')
segments(x0 = 1:(300 - 1), y0 = muskox$step[1:(300-1)],
         x1 = 2:300, y1 = muskox$step[2:300],
         col = pal[states_00[1:(300-1)]], lwd = 1.5)
plot(muskox$step[1:300],type='l',bty='n',xlab='time',ylab='step length')
segments(x0 = 1:(300 - 1), y0 = muskox$step[1:(300-1)],
         x1 = 2:300, y1 = muskox$step[2:300],
         col = pal[states_11[1:(300-1)]], lwd = 1.5)
legend('top', c(paste("State",1:N)), bty='n', lwd=2,
       col=pal[1:2],inset=c(0,-0.75),xpd='NA',horiz=T)# N=3: inset=c(0.25,-0.1)

1-sum(states_00==states_11)/dim(muskox)[1]


# N=3
N=3
# AR(0,0), best model
states_00 = viterbi_arp(x=cbind(muskox$step, muskox$angle),
                        Gamma=matrix(c(0.705,0.231,0.064,0.122,0.806,0.072,0.063,0.179,0.767),byrow=T,ncol=3),
                        delta=c(0.260,0.511,0.229),
                        dists=c('gamma','vm'),
                        autocor=0,
                        params=list(list(mu=c(5.731,59.756, 287.373), sigma=c(4.023, 52.325,277.459)),
                                    list(mu=c(3.062, -2.899,0.003), kappa=c(0.472, 0.111, 0.912))),
                        N=3,
                        p=c(0,0))
states_00

# AR(1,1), best model
states_11 = viterbi_arp(x=cbind(muskox$step, muskox$angle),
                        Gamma=matrix(c(0.676,0.010,0.229,0.085,0.841,0.073,0.078,0.089,0.834),byrow=T,ncol=3),
                        delta=c(0.199,0.364,0.437),
                        dists=c('gamma','vm'),
                        autocor=list(c(0.007,0.442,0.381),c(0.295,0.262,0.047)),
                        params=list(list(mu=c(4.433, 37.267,193.574), sigma=c(2.822, 30.480,197.814)),
                                    list(mu=c(3.118, -3.111,-0.044), kappa=c(0.527, 0.365, 0.563))),
                        N=3,
                        p=c(1,1))
states_11

par(mfrow=c(2,1))
plot(muskox$step[1:300],type='l',bty='n',xlab='time',ylab='step length')
segments(x0 = 1:(300 - 1), y0 = muskox$step[1:(300-1)],
         x1 = 2:300, y1 = muskox$step[2:300],
         col = pal[states_00[1:(300-1)]], lwd = 1.5)
plot(muskox$step[1:300],type='l',bty='n',xlab='time',ylab='step length')
segments(x0 = 1:(300 - 1), y0 = muskox$step[1:(300-1)],
         x1 = 2:300, y1 = muskox$step[2:300],
         col = pal[states_11[1:(300-1)]], lwd = 1.5)
legend('top', c(paste("State",1:N)), bty='n', lwd=2,
       col=pal[1:N],inset=c(0,-0.75),xpd='NA',horiz=T)# N=3: inset=c(0.25,-0.1)

1-sum(states_00==states_11)/dim(muskox)[1]



par(mfrow=c(2,2))
## Pseudo residuals
pres_00 = pseudores_arp(mod=best_mod_00,data=cbind(muskox$step, muskox$angle),N=3,p=c(0,0))
qqnorm(pres_00$stepRes[pres_00$stepRes!=-Inf], bty='n',pch=19,
       xlim=c(-4.25,4.25),ylim=c(-3,4.5),main='Basic HMM')
abline(a=0,b=1,lwd=2)


pres_11 = pseudores_arp(mod=best_mod_11,data=cbind(muskox$step, muskox$angle),N=3,p=c(1,1))
qqnorm(pres_11$stepRes[pres_11$stepRes!=-Inf], bty='n',pch=19,
       xlim=c(-4.25,4.25),ylim=c(-3,4.5),main='AR(1,1)-HMM')
abline(a=0,b=1,lwd=2)

plot(density(pres_00$stepRes[pres_00$stepRes!=-Inf],na.rm=T),
     bty='n',main='',xlab='x',lwd=1.5,xlim=c(-5,5))
curve(dnorm(x),-5,5,lty=2,add=T,col=4,lwd=1.5)
plot(density(pres_11$stepRes[pres_11$stepRes!=-Inf],na.rm=T),
     bty='n',main='',lwd=1.5,xlab='x',xlim=c(-5,5))
curve(dnorm(x),-5,5,lty=2,add=T,lwd=1.5,col=4)
legend('topleft',c('KDE of pseudo residuals','N(0,1)'),col=c(1,4),xpd='NA',bty='n',lwd=1.5,
       horiz=T,lty=1,inset=c(-0.7,-0.3))

