## 2022-09-12
# Code to produce analyses of sea tern data for chapter 4.1

# set wd to current file directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # RStudio way

source("tsa_functions_schwalbe.R")
# read in data
schwalben = read_schwalbe("../datasets/seeschwalbe/slick2-h1-h2.csv")
count=1
for (i in unique(schwalben$ID)){
  assign(paste("schwalbe_", count, sep=""), select_schwalbe(i, count))
  count=count+1
}

## for our analyses we take Schwalbe 77

# pre-process data
library(moveHMM)
s_77 = prepData(schwalbe_77)
head(s_77)

layout(matrix(c(1,1,2, 1,1,3),byrow=TRUE,ncol=3))
plot(s_77$x,s_77$y, type='l', bty='n',xlab='x',ylab='y', cex.lab=1.25)
hist(s_77$step, breaks=25,prob=T,main="", xlab='step length', cex.lab=1.25)
hist(s_77$angle, breaks=80,prob=T,main="", xlab='turning angle', cex.lab=1.25,
     xlim=c(-1,1))


### fit different models for original data
data_s_77 = cbind(s_77$step,s_77$angle)[2:(nrow(s_77)-1),]
states = c(2,3,4)
ar_deg = c(0,1,2,3)
mods = list()
comp_times = rep(NA,12)
c=0



# AR(0,0)
p=c(0,0)
N=3
theta = c(rep(-2,N*(N-1)),
          15,20,35,9,4,3,
          0,0,0,20,80,1000)
theta.star = starize(theta, N=N, p=p, dists=c('gamma','vm'),scale_kappa = 20)
begin = Sys.time()
mod = fit_arp_model(mllk, data=data_s_77, theta.star = theta.star, N=N, 
                    p_auto=p, dists=c('gamma','vm'), 
                    scale_kappa = 20)
time = Sys.time()-begin
mod
par(mfrow=c(1,2))
plot_fitted_dist(data=s_77$step,
                 dist='gamma', param=list(mu=mod$params[[1]]$mu,sigma=mod$params[[1]]$sigma),
                 N=N,delta=mod$delta, title='none',breaks=50, xlab='step length', legend=FALSE)
plot_fitted_dist(data=s_77$angle,
                 dist='vm', param=list(mu=mod$params[[2]]$mu,kappa=mod$params[[2]]$kappa),
                 N=N,delta=mod$delta, title='none',breaks=75, xlab='turning angle', legend=FALSE)
library(RColorBrewer)
pal=brewer.pal(N+1,'Dark2')
legend('topright', c(paste("State",1:N),"Total"), bty='n', lwd=2,
       col=pal,inset=c(0.65,-0.1),xpd='NA',horiz=T)



# N=3, AR(0,0)
N=2
p=c(0,3)

schwalbe_wrap <- function(iteration){
  # starting values
  # N=2
  step_init = c(runif(1,5,15),runif(1,22,30),
                runif(1,6,15),runif(1,2,6))
  angle_init = c(runif(2,-0.01,0.01),
                 runif(1,20,100),runif(1,500,1500))
  # N=3
  #step_init = c(runif(1,5,15),runif(1,15,25),runif(1,25,30),
  #              runif(1,6,15),runif(1,3,10),runif(1,2,6))
  #angle_init = c(runif(3,-0.01,0.01),
  #               runif(1,3,15),runif(1,100,250),runif(1,500,1500))
  # N=4
  #step_init = c(runif(1,5,15),runif(1,10,20),runif(1,20,25),runif(1,25,30),
  #              runif(1,6,15),runif(1,3,10),runif(1,3,8),runif(1,2,6))
  #angle_init = c(runif(4,-0.01,0.01),
  #               runif(1,3,15),runif(1,20,50),runif(1,100,250),runif(1,500,1500))
  theta = c(rep(-2,N*(N-1)),
            step_init,
            angle_init,
            rep(0.1,p[1]*N+p[2]*N))
  # working parameters
  theta.star = starize(theta, N=N, p=p, dists=c('gamma','vm'),scale_kappa = 20)
  begin = Sys.time()
  mod = fit_arp_model(mllk, data=data_s_77, theta.star = theta.star, N=N, 
                      p_auto=p, dists=c('gamma','vm'), 
                      scale_kappa = 20)
  res = list()
  time = Sys.time()-begin
  res['time'] = time
  res[['mod']] = mod
  res['aic'] = 2*mod$mllk_optim + 2*length(theta)
  res['bic'] = 2*mod$mllk_optim + log(dim(data_s_77)[1])*length(theta)
  res['mllk'] = mod$mllk_optim
  return(res)
}

library(parallel)
n_cores = detectCores()
iterations = 1:15
results_ar00 <- mclapply(iterations,
                    schwalbe_wrap,
                    mc.cores = n_cores
)
#results_ar00


comps = rep(NA, length(iterations))
aics = rep(NA, length(iterations))
bics = rep(NA, length(iterations))
mllks = rep(NA, length(iterations))
mods = list()
c=1
for (i in iterations){
  mods[[i]] = results_ar00[c][[1]]$mod
  comps[i] = results_ar00[c][[1]]$time
  aics[i] = results_ar00[c][[1]]$aic
  bics[i] = results_ar00[c][[1]]$bic
  mllks[i] = results_ar00[c][[1]]$mllk
  c=c+1
}

comps
aics
bics
mllks
mods

par(mfrow=c(1,2))
plot_fitted_dist(data=s_77$step,
                 dist='gamma', param=list(mu=mods[[29]]$params[[1]]$mu,sigma=mods[[29]]$params[[1]]$sigma),
                 N=N,delta=mods[[29]]$delta, title='none',breaks=50, xlab='step length', legend=FALSE)
plot_fitted_dist(data=s_77$angle,
                 dist='vm', param=list(mu=mods[[29]]$params[[2]]$mu,kappa=mods[[29]]$params[[2]]$kappa),
                 N=N,delta=mods[[29]]$delta, title='none',breaks=75, xlab='turning angle', legend=FALSE)
library(RColorBrewer)
pal=brewer.pal(N+1,'Dark2')
legend('topright', c(paste("State",1:N),"Total"), bty='n', lwd=2,
       col=pal,inset=c(0.25,-0.1),xpd='NA',horiz=T)# N=3: inset=c(0.25,-0.1)

# Model selection
# done, autoregression in step length does not work...
# (one try for univariate gamma HMM)
N=2
p=1
theta = c(rep(-2,N*(N-1)),
          10,25,9,4,
          rep(0.1,p*N))
# working parameters
theta.star = starize(theta, N=N, p=p, dists=c('gamma'),scale_kappa = 20)
mod = fit_arp_model(mllk, data=data_s_77[,1], theta.star = theta.star, N=N, 
                    p_auto=p, dists=c('gamma'))
# doesn't work


# downsample data
library(zoo)
# 5Hz
schwalbe_77_5hz = rollapply(schwalbe_77[,2:dim(schwalbe_77)[2]], 6, mean)
schwalbe_77_5hz = as.data.frame(schwalbe_77_5hz[seq(1, nrow(schwalbe_77_5hz), 6),]) # make data 5hz
data_s_77_5hz = prepData(schwalbe_77_5hz)
data_s_77_5hz = cbind(data_s_77_5hz$step,data_s_77_5hz$angle)[2:(nrow(data_s_77_5hz)-1),]

# 1Hz
schwalbe_77_1hz = rollapply(schwalbe_77[,2:dim(schwalbe_77)[2]], 30, mean)
schwalbe_77_1hz = as.data.frame(schwalbe_77_1hz[seq(1, nrow(schwalbe_77_1hz), 30),]) # make data 5hz
data_s_77_1hz = prepData(schwalbe_77_1hz)
data_s_77_1hz = cbind(data_s_77_1hz$step,data_s_77_1hz$angle)[2:(nrow(data_s_77_1hz)-1),]

par(mfrow=c(1,2))
hist(data_s_77_1hz[,1],probability = T,xlab='step length',breaks=10,
     main='')
hist(data_s_77_1hz[,2],probability = T,xlab='turning angle',breaks=10,
     main='')



# fit different models for downsampled data

# 5Hz

N=2
p=c(2,2)
schwalbe_wrap <- function(iteration){
  # starting values
  step_init = c(runif(1,60,100),runif(1,140,170),
                runif(1,20,50),runif(1,7,20))
  angle_init = c(runif(2,-0.01,0.01),
                 runif(1,4,15),runif(1,50,200))
  theta = c(rep(-2,N*(N-1)),
            step_init,
            angle_init,
            rep(0.1,p[1]*N+p[2]*N))
  # working parameters
  theta.star = starize(theta, N=N, p=p, dists=c('gamma','vm'),scale_kappa = 5)
  begin = Sys.time()
  mod = fit_arp_model(mllk, data=data_s_77_5hz, theta.star = theta.star, N=N, 
                      p_auto=p, dists=c('gamma','vm'), 
                      scale_kappa = 5)
  res = list()
  time = Sys.time()-begin
  res['time'] = time
  res[['mod']] = mod
  res['aic'] = 2*mod$mllk_optim + 2*length(theta)
  res['bic'] = 2*mod$mllk_optim + log(dim(data_s_77_5hz)[1])*length(theta)
  res['mllk'] = mod$mllk_optim
  return(res)
}

iterations = 1:50
results_ar00 <- mclapply(iterations,
                         schwalbe_wrap,
                         mc.cores = n_cores
)
#results_ar00
work=c(1,5,8,9)
work=c(work,work+10,work+20,work+30,work+40)
work

(1:50)[-work]


comps = rep(NA, length(iterations))
aics = rep(NA, length(iterations))
bics = rep(NA, length(iterations))
mllks = rep(NA, length(iterations))
mods = list()
c=1
for (i in work){
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

par(mfrow=c(1,2))
plot_fitted_dist(data=data_s_77_5hz[,1],
                 dist='gamma', param=list(mu=mods[[45]]$params[[1]]$mu,sigma=mods[[45]]$params[[1]]$sigma),
                 N=N,delta=mods[[45]]$delta, title='none',breaks=25, xlab='step length', legend=FALSE)
plot_fitted_dist(data=data_s_77_5hz[,2],
                 dist='vm', param=list(mu=mods[[45]]$params[[2]]$mu,kappa=mods[[45]]$params[[2]]$kappa),
                 N=N,delta=mods[[45]]$delta, title='none',breaks=25, xlab='turning angle', legend=FALSE)
pal=brewer.pal(N+1,'Dark2')
legend('topright', c(paste("State",1:N),"Total"), bty='n', lwd=2,
       col=pal,inset=c(0.45,-0.1),xpd='NA',horiz=T)# N=3: inset=c(0.25,-0.1)


# 1Hz

N=2
p=c(1,1)
schwalbe_wrap <- function(iteration){
  # starting values
  step_init = c(runif(1,300,500),runif(1,700,850),
                runif(1,100,300),runif(1,50,100))
  angle_init = c(runif(2,-0.01,0.01),
                 runif(1,1.5,3),runif(1,2.5,5))
  theta = c(rep(-2,N*(N-1)),
            step_init,
            angle_init,
            rep(0.1,p[1]*N+p[2]*N))
  # working parameters
  theta.star = starize(theta, N=N, p=p, dists=c('gamma','vm'),scale_kappa = 1)
  begin = Sys.time()
  mod = fit_arp_model(mllk, data=data_s_77_1hz, theta.star = theta.star, N=N, 
                      p_auto=p, dists=c('gamma','vm'), 
                      scale_kappa = 1)
  res = list()
  time = Sys.time()-begin
  res['time'] = time
  res[['mod']] = mod
  res['aic'] = 2*mod$mllk_optim + 2*length(theta)
  res['bic'] = 2*mod$mllk_optim + log(dim(data_s_77_1hz)[1])*length(theta)
  res['mllk'] = mod$mllk_optim
  return(res)
}

iterations = 1:50
results_ar00 <- mclapply(iterations,
                         schwalbe_wrap,
                         mc.cores = n_cores
)
#results_ar00
#work=c(1,2,4,5,6,7,8,9,10)
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
best_mod_31 = mods[[which.min(bics)]]


par(mfrow=c(1,2))
plot_fitted_dist(data=data_s_77_1hz[,1],
                 dist='gamma', param=list(mu=mods[[36]]$params[[1]]$mu,sigma=mods[[36]]$params[[1]]$sigma),
                 N=N,delta=mods[[36]]$delta, title='none',breaks=20, xlab='step length', legend=FALSE)
plot_fitted_dist(data=data_s_77_1hz[,2],
                 dist='vm', param=list(mu=mods[[36]]$params[[2]]$mu,kappa=mods[[36]]$params[[2]]$kappa),
                 N=N,delta=mods[[36]]$delta, title='none',breaks=15, xlab='turning angle', legend=FALSE)
pal=brewer.pal(N+1,'Dark2')
legend('topright', c(paste("State",1:N),"Total"), bty='n', lwd=2,
       col=pal,inset=c(0.45,-0.1),xpd='NA',horiz=T)# N=3: inset=c(0.25,-0.1)

par(mfrow=c(1,2))
plot_fitted_dist(data=data_s_77_1hz[,1],
                 dist='gamma', param=list(mu=c(429.303, 765.007),sigma=c(179.193, 63.354)),
                 N=N,delta=c(0.278,0.722), title='none',breaks=20, xlab='step length', legend=FALSE)
plot_fitted_dist(data=data_s_77_1hz[,2],
                 dist='vm', param=list(mu=c(-0.037, 0.206),kappa=c(2.833, 15.631)),
                 N=N,delta=c(0.278,0.722), title='none',breaks=15, xlab='turning angle', legend=FALSE)
pal=brewer.pal(N+1,'Dark2')
legend('topright', c(paste("State",1:N),"Marginal"), bty='n', lwd=2,
       col=pal,inset=c(0.55,-0.1),xpd='NA',horiz=T)# N=3: inset=c(0.25,-0.1)


# state decoding

# AR(0,0), best model
states_00 = viterbi_arp(x=data_s_77_1hz,
            Gamma=matrix(c(0.739,0.261,0.052,0.948),byrow=T,ncol=2),
            delta=c(0.165,0.835),
            dists=c('gamma','vm'),
            autocor=0,
            params=list(list(mu=c(429.303, 765.007), sigma=c(189.151, 72.251)),
                     list(mu=c(-0.113, 0.041), kappa=c(2.756, 4.634))),
            N=2,
            p=c(0,0))
states_00

par(mfrow=c(2,1))
plot(data_s_77_1hz[,1],type='p',bty='n',xlab='time',ylab='step length',pch=19,cex=0.5,
     col=pal[states_00])
plot(data_s_77_1hz[,1],type='l',bty='n',xlab='time',ylab='step length')
segments(x0 = 1:(dim(data_s_77_1hz)[1] - 1), y0 = data_s_77_1hz[-dim(data_s_77_1hz)[1],1],
         x1 = 2:dim(data_s_77_1hz)[1], y1 = data_s_77_1hz[-1,1],
         col = pal[states_00[-(dim(data_s_77_1hz)[1])]], lwd = 1.5)
plot(data_s_77_1hz[,2],type='l',bty='n',xlab='time',ylab='turning angle')
segments(x0 = 1:(dim(data_s_77_1hz)[1] - 1), y0 = data_s_77_1hz[-dim(data_s_77_1hz)[1],2],
         x1 = 2:dim(data_s_77_1hz)[1], y1 = data_s_77_1hz[-1,2],
         col = pal[states_00[-(dim(data_s_77_1hz)[1])]], lwd = 1.5)
legend('top', c(paste("State",1:N)), bty='n', lwd=2,
       col=pal[1:2],inset=c(0,-0.75),xpd='NA',horiz=T)# N=3: inset=c(0.25,-0.1)


# AR(3,1), best model
states_31 = viterbi_arp(x=data_s_77_1hz,
                        Gamma=matrix(c(0.676,0.324,0.125,0.875),byrow=T,ncol=2),
                        delta=c(0.278,0.722),
                        dists=c('gamma','vm'),
                        autocor=list(c(0.018,0.013,0.619,0.048,0.002,0.318),c(0.304,0.795)),
                        params=list(list(mu=c(483.283, 776.814), sigma=c(179.193, 63.354)),
                                    list(mu=c(-0.037, 0.206), kappa=c(2.833, 15.631))),
                        N=2,
                        p=c(3,1))
states_31

par(mfrow=c(2,1))
plot(data_s_77_1hz[,1],type='l',bty='n',xlab='time',ylab='step length')
segments(x0 = 1:(dim(data_s_77_1hz)[1] - 1), y0 = data_s_77_1hz[-dim(data_s_77_1hz)[1],1],
         x1 = 2:dim(data_s_77_1hz)[1], y1 = data_s_77_1hz[-1,1],
         col = pal[states_31[-(dim(data_s_77_1hz)[1])]], lwd = 1.5)
plot(data_s_77_1hz[,2],type='l',bty='n',xlab='time',ylab='turning angle')
segments(x0 = 1:(dim(data_s_77_1hz)[1] - 1), y0 = data_s_77_1hz[-dim(data_s_77_1hz)[1],2],
         x1 = 2:dim(data_s_77_1hz)[1], y1 = data_s_77_1hz[-1,2],
         col = pal[states_31[-(dim(data_s_77_1hz)[1])]], lwd = 1.5)
legend('top', c(paste("State",1:N)), bty='n', lwd=2,
       col=pal[1:2],inset=c(0,-0.75),xpd='NA',horiz=T)# N=3: inset=c(0.25,-0.1)


1-sum(states_00==states_31)/dim(data_s_77_1hz)[1]

# Plot that compares step lengths exclusively:
par(mfrow=c(2,1), mar=c(5.1,4.1,4.1,2.1), xpd=TRUE, mai=c(1,0.8,0.25,0.1))
plot(data_s_77_1hz[,1],type='l',bty='n',xlab='time',ylab='step length',
     col='grey')
#segments(x0 = 1:(dim(data_s_77_1hz)[1] - 1), y0 = data_s_77_1hz[-dim(data_s_77_1hz)[1],1],
#         x1 = 2:dim(data_s_77_1hz)[1], y1 = data_s_77_1hz[-1,1],
#         col = pal[states_00[-(dim(data_s_77_1hz)[1])]], lwd = 1.5)
points(x=1:length(data_s_77_1hz[,1]),y=data_s_77_1hz[,1], col=pal[states_00], cex=0.75, pch=19)
plot(data_s_77_1hz[,1],type='l',bty='n',xlab='time',ylab='step length',
     col='grey')
#segments(x0 = 1:(dim(data_s_77_1hz)[1] - 1), y0 = data_s_77_1hz[-dim(data_s_77_1hz)[1],1],
#         x1 = 2:dim(data_s_77_1hz)[1], y1 = data_s_77_1hz[-1,1],
#         col = pal[states_31[-(dim(data_s_77_1hz)[1])]], lwd = 1.5)
points(x=1:length(data_s_77_1hz[,1]),y=data_s_77_1hz[,1], col=pal[states_31], cex=0.75, pch=19)
legend('top', c(paste("State",1:N)), bty='n', lwd=2,
       col=pal[1:2],inset=c(0,-0.35),xpd='NA',horiz=T)# N=3: inset=c(0.25,-0.1)


par(mfrow=c(2,2))
## Pseudo residuals
pres_00 = pseudores_arp(mod=best_mod_00,data=data_s_77_1hz,N=2,p=c(0,0))
qqnorm(pres_00$stepRes, bty='n',pch=19,xlim=c(-2.75,2.75),ylim=c(-2.5,2.5),
       main='Basic HMM')
abline(a=0,b=1,lwd=2)


pres_31 = pseudores_arp(mod=best_mod_31,data=data_s_77_1hz,N=2,p=c(3,1))
qqnorm(pres_31$stepRes, bty='n',pch=19,xlim=c(-2.75,2.75),ylim=c(-2.5,2.5),
       main='AR(3,1)-HMM')
abline(a=0,b=1,lwd=2)


plot(density(pres_00$stepRes[pres_00$stepRes!=-Inf],na.rm=T),
     bty='n',main='',xlab='x',lwd=1.5,xlim=c(-5,5))
curve(dnorm(x),-5,5,lty=2,add=T,col=4,lwd=1.5)
plot(density(pres_31$stepRes[pres_31$stepRes!=-Inf],na.rm=T),
     bty='n',main='',lwd=1.5,xlab='x',xlim=c(-5,5))
curve(dnorm(x),-5,5,lty=2,add=T,lwd=1.5,col=4)
legend('topleft',c('KDE of pseudo residuals','N(0,1)'),col=c(1,4),xpd='NA',bty='n',lwd=1.5,
       horiz=T,lty=1,inset=c(-0.7,-0.3))

