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
                3,
                c(25,60,190,5,20,200,rep(0.01,3)), # account for zeromass 
                c(pi,0,0,0.5,0.5,0.5),
                stationary = TRUE) 
plot(mod_move)
mod_move
pseudo_residuals = pseudoRes(mod_move)
qqnorm(pseudo_residuals$stepRes)
curve(1*x,-4,4,add=T)


# no autocorrelation
theta=c(rep(-2,2),25,150,20,200,pi,0,0.5,0.5,0.001,0.001) 
theta.star=starize(theta, N=2,p=c(0,0),dists=c('gamma','vm'), scale_kappa = 1, zero_inf = T)
mod_ar0=fit_arp_model(mllk=mllk, data=cbind(muskox$step, muskox$angle), 
                      theta.star=theta.star, N=2, p_auto=c(0,0),dists=c('gamma','vm'),scale_kappa = 1,
                      zero_inf = T)
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
theta=c(rep(-2,2),25,150,20,200,pi,0,0.5,0.5, rep(0.1,4), 0.001,0.001) 
theta.star=starize(theta, N=2,p=c(1,1),dists=c('gamma','vm'), scale_kappa = 1, zero_inf = T)
mod_ar1both=fit_arp_model(mllk=mllk, data=cbind(muskox$step, muskox$angle), 
                          theta.star=theta.star, N=2, p_auto=c(1,1),dists=c('gamma','vm'), scale_kappa = 1,
                          zero_inf = T)
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
theta=c(rep(-2,2),10,150,8,165,pi,0,0.5,0.5, rep(0.1,8), 0.001,0.001) 
theta.star=starize(theta, N=2,p=c(2,2),dists=c('gamma','vm'), scale_kappa = 1,zero_inf = T)
mod_ar2both=fit_arp_model(mllk=mllk, data=cbind(muskox$step, muskox$angle), 
                          theta.star=theta.star, N=2, p_auto=c(2,2),dists=c('gamma','vm'), scale_kappa = 1,
                          zero_inf = T)
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
# For 10 models, parallelize!!



library(parallel)
n_cores = detectCores()

N=3
p=c(3,3)
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
            rep(0.1,p[1]*N+p[2]*N),
            rep(0.001,N))
  # working parameters
  theta.star = starize(theta, N=N, p=p, dists=c('gamma','vm'),scale_kappa = 1, zero_inf = T)
  begin = Sys.time()
  mod = fit_arp_model(mllk, data=cbind(muskox$step, muskox$angle), theta.star = theta.star, N=N, 
                      p_auto=p, dists=c('gamma','vm'), 
                      scale_kappa = 1, zero_inf = T)
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
which.min(bics)
mean(comps)*60
best_mod_00_n2 = mods[[which.min(bics)]]
best_mod_00_n3 = mods[[which.min(bics)]]
best_mod_11_n2 = mods[[which.min(bics)]]
best_mod_11_n3 = mods[[which.min(bics)]]
best_mod_22_n2 = mods[[which.min(bics)]]
best_mod_22_n3 = mods[[which.min(bics)]]
best_mod_33_n2 = mods[[which.min(bics)]]
best_mod_33_n3 = mods[[which.min(bics)]]





par(mfrow=c(1,2))
plot_fitted_dist(data=muskox$step,
                 dist='gamma', param=list(mu=mods[[1]]$params[[1]]$mu,sigma=mods[[1]]$params[[1]]$sigma),
                 N=3,delta=mods[[1]]$delta, title='none',breaks=100,
                 xlim=c(0,500),legend=F,xlab='step length')
plot_fitted_dist(data=muskox$angle,
                 dist='vm', param=list(mu=mods[[1]]$params[[2]]$mu,kappa=mods[[1]]$params[[2]]$kappa),
                 N=3,delta=mods[[1]]$delta, title='none',breaks=20,
                 legend=F,xlab='turning angle')
pal=brewer.pal(N+1,'Dark2')
legend('topleft', c(paste("State",1:N),"Total"), bty='n', lwd=2,
       col=pal,inset=c(#-1.2,#
                       -1.35,
                       -0.1),xpd='NA',horiz=T)



# state decoding

N=2
# AR(0,0), best model
states_00_n2 = viterbi_arp(x=cbind(muskox$step, muskox$angle),
                        Gamma=matrix(c(0.857,0.143,0.160,0.840),byrow=T,ncol=2),
                        delta=c(0.529,0.471),
                        dists=c('gamma','vm'),
                        autocor=0,
                        params=list(list(mu=c(21.056, 183.989), sigma=c(21.886, 193.946),
                                         zero_inf=c(0.000006,0.000015)),
                                    list(mu=c(3.094, -0.060), kappa=c(0.414, 0.503))),
                        N=2,
                        p=c(0,0))
states_00_n2

# AR(1,1), best model
states_11_n2 = viterbi_arp(x=cbind(muskox$step, muskox$angle),
                        Gamma=matrix(c(0.745,0.255,0.203,0.797),byrow=T,ncol=2),
                        delta=c(0.443,0.557),
                        dists=c('gamma','vm'),
                        autocor=list(c(0.650,0.410),c(0.230,0.040)),
                        params=list(list(mu=c(10.739, 154.714), sigma=c(8.790, 866),
                                         zero_inf=c(0.012,0.001)),
                                    list(mu=c(3.130, -0.059), kappa=c(0.538, 0.451))),
                        N=2,
                        p=c(1,1))
states_11_n2

par(mfrow=c(2,1))
plot(muskox$step[1:300],type='l',bty='n',xlab='time',ylab='step length')
segments(x0 = 1:(300 - 1), y0 = muskox$step[1:(300-1)],
         x1 = 2:300, y1 = muskox$step[2:300],
         col = pal[states_00_n2[1:(300-1)]], lwd = 1.5)
plot(muskox$step[1:300],type='l',bty='n',xlab='time',ylab='step length')
segments(x0 = 1:(300 - 1), y0 = muskox$step[1:(300-1)],
         x1 = 2:300, y1 = muskox$step[2:300],
         col = pal[states_11_n2[1:(300-1)]], lwd = 1.5)
legend('top', c(paste("State",1:N)), bty='n', lwd=2,
       col=pal[1:2],inset=c(0,-0.75),xpd='NA',horiz=T)# N=3: inset=c(0.25,-0.1)

1-sum(states_00_n2==states_11_n2)/dim(muskox)[1]


# N=3
N=3
# AR(0,0), best model
states_00_n3 = viterbi_arp(x=cbind(muskox$step, muskox$angle),
                        Gamma=matrix(c(0.704,0.234,0.062,0.122,0.805,0.073,0.061,0.177,0.762),byrow=T,ncol=3),
                        delta=c(0.259,0.515,0.226),
                        dists=c('gamma','vm'),
                        autocor=0,
                        params=list(list(mu=c(5.728,60.009, 289.824), sigma=c(4.016, 52.687,279.007),
                                         zero_inf=c(0.0001,0.00002,0.0001)),
                                    list(mu=c(3.072, -2.898,0.007), kappa=c(0.474, 0.112, 0.925))),
                        N=3,
                        p=c(0,0))
states_00_n3

# AR(1,1), best model
states_11_n3 = viterbi_arp(x=cbind(muskox$step, muskox$angle),
                        Gamma=matrix(c(0.674,0.094,0.231,0.087,0.834,0.079,0.077,0.094,0.829),byrow=T,ncol=3),
                        delta=c(0.200,0.362,0.438),
                        dists=c('gamma','vm'),
                        autocor=list(c(0.007,0.440,0.387),c(0.309,0.269,0.047)),
                        params=list(list(mu=c(4.408, 36.941,198.893), sigma=c(2.803, 30.127,198.893),
                                         zero_inf=c(0.023,0.002,0.002)),
                                    list(mu=c(3.121, -3.114,-0.042), kappa=c(0.450, 0.372, 0.567))),
                        N=3,
                        p=c(1,1))
states_11_n3

par(mfrow=c(2,1))
plot(muskox$step[1:300],type='p',bty='n',xlab='time',ylab='step length',
     pch=19,cex=0.5,col=pal[states_00_n3])
plot(muskox$step[1:300],type='l',bty='n',xlab='time',ylab='step length')
segments(x0 = 1:(300 - 1), y0 = muskox$step[1:(300-1)],
         x1 = 2:300, y1 = muskox$step[2:300],
         col = pal[states_00_n3[1:(300-1)]], lwd = 1.5)
plot(muskox$step[1:300],type='l',bty='n',xlab='time',ylab='step length')
segments(x0 = 1:(300 - 1), y0 = muskox$step[1:(300-1)],
         x1 = 2:300, y1 = muskox$step[2:300],
         col = pal[states_11_n3[1:(300-1)]], lwd = 1.5)
legend('top', c(paste("State",1:N)), bty='n', lwd=2,
       col=pal[1:N],inset=c(0,-0.75),xpd='NA',horiz=T)# N=3: inset=c(0.25,-0.1)

1-sum(states_00_n3==states_11_n3)/dim(muskox)[1]



par(mfrow=c(2,2))
## Pseudo residuals

# N=2
pres_00_n2 = pseudores_arp(mod=best_mod_00_n2,data=cbind(muskox$step, muskox$angle),N=2,p=c(0,0))
qqnorm(pres_00_n2$stepRes[pres_00_n2$stepRes!=-Inf], bty='n',pch=19,
       xlim=c(-4.25,4.25),ylim=c(-3,6),main='Basic HMM')
abline(a=0,b=1,lwd=2)

pres_11_n2 = pseudores_arp(mod=best_mod_11_n2,data=cbind(muskox$step, muskox$angle),N=2,p=c(1,1))
qqnorm(pres_11_n2$stepRes[pres_11_n2$stepRes!=-Inf], bty='n',pch=19,
       xlim=c(-4.25,4.25),ylim=c(-3,6),main='AR(1,1)-HMM')
abline(a=0,b=1,lwd=2)

plot(density(pres_00_n2$stepRes[pres_00_n2$stepRes!=-Inf],na.rm=T),
     bty='n',main='',xlab='x',lwd=1.5,xlim=c(-6,6))
curve(dnorm(x),-5,5,lty=2,add=T,col=4,lwd=1.5)
plot(density(pres_11_n2$stepRes[pres_11_n2$stepRes!=-Inf],na.rm=T),
     bty='n',main='',lwd=1.5,xlab='x',xlim=c(-6,6))
curve(dnorm(x),-5,5,lty=2,add=T,lwd=1.5,col=4)
legend('topleft',c('KDE of pseudo residuals','N(0,1)'),col=c(1,4),xpd='NA',bty='n',lwd=1.5,
       horiz=T,lty=1,inset=c(-0.7,-0.3))

# N=3
pres_00_n3 = pseudores_arp(mod=best_mod_00_n3,data=cbind(muskox$step, muskox$angle),N=3,p=c(0,0))
qqnorm(pres_00_n3$stepRes[pres_00_n3$stepRes!=-Inf], bty='n',pch=19,
       xlim=c(-4.25,4.25),ylim=c(-3,6),main='Basic HMM')
abline(a=0,b=1,lwd=2)

pres_11_n3 = pseudores_arp(mod=best_mod_11_n3,data=cbind(muskox$step, muskox$angle),N=3,p=c(1,1))
qqnorm(pres_11_n3$stepRes[pres_11_n3$stepRes!=-Inf], bty='n',pch=19,
       xlim=c(-4.25,4.25),ylim=c(-3,6),main='AR(1,1)-HMM')
abline(a=0,b=1,lwd=2)

plot(density(pres_00_n3$stepRes[pres_00_n3$stepRes!=-Inf],na.rm=T),
     bty='n',main='',xlab='x',lwd=1.5,xlim=c(-5,5))
curve(dnorm(x),-5,5,lty=2,add=T,col=4,lwd=1.5)
plot(density(pres_11_n3$stepRes[pres_11_n3$stepRes!=-Inf],na.rm=T),
     bty='n',main='',lwd=1.5,xlab='x',xlim=c(-5,5))
curve(dnorm(x),-5,5,lty=2,add=T,lwd=1.5,col=4)
legend('topleft',c('KDE of pseudo residuals','N(0,1)'),col=c(1,4),xpd='NA',bty='n',lwd=1.5,
       horiz=T,lty=1,inset=c(-0.7,-0.3))



# comparison to AR(3,3) for N=3
qqnorm(pres_00_n3$stepRes[pres_00_n3$stepRes!=-Inf], bty='n',pch=19,
       xlim=c(-4.25,4.25),ylim=c(-3,6),main='Basic HMM')
abline(a=0,b=1,lwd=2)
pres_33_n3 = pseudores_arp(mod=best_mod_33_n3,data=cbind(muskox$step, muskox$angle),N=3,p=c(3,3))
qqnorm(pres_33_n3$stepRes[pres_33_n3$stepRes!=-Inf], bty='n',pch=19,
       xlim=c(-4.25,4.25),ylim=c(-3,6),main='AR(3,3)-HMM')
abline(a=0,b=1,lwd=2)

plot(density(pres_00_n3$stepRes[pres_00_n3$stepRes!=-Inf],na.rm=T),
     bty='n',main='',xlab='x',lwd=1.5,xlim=c(-5,5))
curve(dnorm(x),-5,5,lty=2,add=T,col=4,lwd=1.5)
plot(density(pres_33_n3$stepRes[pres_33_n3$stepRes!=-Inf],na.rm=T),
     bty='n',main='',lwd=1.5,xlab='x',xlim=c(-5,5))
curve(dnorm(x),-5,5,lty=2,add=T,lwd=1.5,col=4)
# lassen wir jetzt die Finger von :)


