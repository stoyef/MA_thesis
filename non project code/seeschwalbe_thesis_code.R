## 2022-09-12
# Code to prduce analyses of tea tern data for chapter 4.1

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

layout(matrix(c(1,3,3, 2,3,3),byrow=TRUE,ncol=3))
hist(s_77$step, breaks=25,prob=T,main="", xlab='step length', cex.lab=1.25)
hist(s_77$angle, breaks=80,prob=T,main="", xlab='turning angle', cex.lab=1.25,
     xlim=c(-1,1))
plot(s_77$x,s_77$y, type='l', bty='n',xlab='x',ylab='y', cex.lab=1.25)


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
iterations = 1:50
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
for (i in c(1:2,6:12,16:22,26:32,36:42,46:50)){#iterations){
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

# downsample data


# fit different models for downsampled data



