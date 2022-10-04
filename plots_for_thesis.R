# 2022-08-17
# Generation of plots used in the thesis
library(MasterThesis)

## Chapter 1 - Introduction

## Chapter 2 - HMMs

## Chapter 3 - Simulation Study

library(RColorBrewer)
pal=brewer.pal(3,'Dark2')
# Plot simulated data as time series
# AR(0)
data_ar0 = sample_arp(2000,
                      delta=c(0.5,0.5),
                      Gamma=matrix(c(0.9,0.1,0.1,0.9),ncol=2,byrow=TRUE),
                      N=2,
                      params=c(10,30,5,5,0,0,3,20),
                      autocor=0,
                      p=c(0,0),
                      dists=c('gamma','vm'))
par(mfrow=c(2,1), mar=c(5.1,4.1,4.1,2.1), xpd=TRUE, mai=c(1,0.8,0.25,0))
plot_decoded_data(data_ar0$data[1:300,1], data_ar0$states[1:300], col=pal, name='step length', legend=FALSE)
plot_decoded_data(data_ar0$data[1:300,2], data_ar0$states[1:300], col=pal, name='turning angle', legend=FALSE)
legend('top',c('State 1','State 2'), col=pal[1:2],lwd=1.5,bty='n',
       y.intersp = 0.75,
       horiz=TRUE, inset=c(0,-0.35),xpd='NA')

# AR(1)
data = sample_arp(2000,
                  delta=c(0.5,0.5),
                  Gamma=matrix(c(0.9,0.1,0.1,0.9),ncol=2,byrow=TRUE),
                  N=2,
                  params=c(10,30,5,5,0,0,3,20),
                  autocor=list(matrix(c(0.3,0.6),ncol=1),matrix(c(0.3,0.6),ncol=1)),
                  p=c(1,1),
                  dists=c('gamma','vm'))
par(mfrow=c(2,1), mar=c(5.1,4.1,4.1,2.1), xpd=TRUE, mai=c(1,0.8,0.25,0))
plot_decoded_data(data$data[1:300,1], data$states[1:300], col=pal, name='step length', legend=FALSE)
plot_decoded_data(data$data[1:300,2], data$states[1:300], col=pal, name='turning angle', legend=FALSE)
legend('top',c('State 1','State 2'), col=pal[1:2],lwd=1.5,bty='n',
       y.intersp = 0.75,
       horiz=TRUE, inset=c(0,-0.35),xpd='NA')

# AR(2)
data_ar2 = sample_arp(2000,
                      delta=c(0.5,0.5),
                      Gamma=matrix(c(0.9,0.1,0.1,0.9),ncol=2,byrow=TRUE),
                      N=2,
                      params=c(10,30,5,5,0,0,3,20),
                      autocor=list(matrix(c(0.1,0.2,0.2,0.4),ncol=2, byrow=T),
                                   matrix(c(0.1,0.2,0.2,0.4),ncol=2,byrow=T)),
                      p=c(2,2),
                      dists=c('gamma','vm'))
par(mfrow=c(2,1), mar=c(5.1,4.1,4.1,2.1), xpd=TRUE, mai=c(1,0.8,0.25,0))
plot_decoded_data(data_ar2$data[1:300,1], data_ar2$states[1:300], col=pal, name='step length', legend=FALSE)
plot_decoded_data(data_ar2$data[1:300,2], data_ar2$states[1:300], col=pal, name='turning angle', legend=FALSE)
legend('top',c('State 1','State 2'), col=pal[1:2],lwd=1.5,bty='n',
       y.intersp = 0.75,
       horiz=TRUE, inset=c(0,-0.35),xpd='NA')

# AR(3)
data_ar3 = sample_arp(2000,
                  delta=c(0.5,0.5),
                  Gamma=matrix(c(0.9,0.1,0.1,0.9),ncol=2,byrow=TRUE),
                  N=2,
                  params=c(10,30,5,5,0,0,3,20),
                  autocor=list(matrix(c(0.1,0.1,0.1,0.2,0.2,0.2),ncol=3, byrow=T),
                               matrix(c(0.1,0.1,0.1,0.2,0.2,0.2),ncol=3,byrow=T)),
                  p=c(3,3),
                  dists=c('gamma','vm'))
par(mfrow=c(2,1), mar=c(5.1,4.1,4.1,2.1), xpd=TRUE, mai=c(1,0.8,0.25,0))
plot_decoded_data(data_ar3$data[1:300,1], data_ar3$states[1:300], col=pal, name='step length', legend=FALSE)
plot_decoded_data(data_ar3$data[1:300,2], data_ar3$states[1:300], col=pal, name='turning angle', legend=FALSE)
legend('top',c('State 1','State 2'), col=pal[1:2],lwd=1.5,bty='n',
       y.intersp = 0.75,
       horiz=TRUE, inset=c(0,-0.35),xpd='NA')
dev.off()


# Ground truth of the distributions sampled from (no autocorrelation)
par(mfrow=c(1,2))
curve(0.5*dgamma(x,shape=20^2/5^2, scale=5^2/20),col=pal[1],lwd=2,add=F,
      from=0,to=80,n=1000,bty='n',xlab='step length', ylab='Density')
curve(0.5*dgamma(x,shape=40^2/7^2, scale=7^2/40),col=pal[2],lwd=2,add=T,
      from=0,to=80,n=1000)
curve(0.5*dgamma(x,shape=20^2/5^2, scale=5^2/20)+
        0.5*dgamma(x,shape=40^2/7^2, scale=7^2/40), col=pal[3],lwd=2,add=T,
      from=0,to=80,n=1000)
legend('topright',
       c('State 1', 'State 2', 'Marginal'),
       col=pal,
       lwd=2,bty='n')
library(CircStats)
curve(0.5*dvm(x,mu=0, kappa=2),col=pal[1],lwd=2,add=F,
      from=-pi,to=pi,n=1000,bty='n',xlab='turning angle', ylab='Density',
      ylim=c(0,1))
curve(0.5*dvm(x,mu=0, kappa=12),col=pal[2],lwd=2,add=T,
      from=-pi,to=pi,n=1000)
curve(0.5*dvm(x,mu=0, kappa=2)+ 0.5*dvm(x,mu=0, kappa=12), 
      col=pal[3],lwd=2,add=T,
      from=-pi,to=pi,n=1000)
legend('topright',
       c('State 1', 'State 2', 'Marginal'),
       col=pal,
       lwd=2,bty='n')



# Comparison of distributions and histograms
# AR(0,0)
hist(data_ar0$data[,1],breaks=35,xlab='step length', prob=T,
     xlim=c(0,60),ylim=c(0,0.045),main="Basic HMM")
library(RColorBrewer)
pal=brewer.pal(3,'Dark2')
curve(0.5*dgamma(x,shape=10^2/5^2, scale=5^2/10),col=pal[1],lwd=2,add=T,
      from=0,to=60,n=1000)
curve(0.5*dgamma(x,shape=30^2/5^2, scale=5^2/30),col=pal[2],lwd=2,add=T,
      from=0,to=60,n=1000)
curve(0.5*dgamma(x,shape=10^2/5^2, scale=5^2/10)+
        0.5*dgamma(x,shape=30^2/5^2, scale=5^2/30), col=pal[3],lwd=2,add=T,
      from=0,to=60,n=1000)
legend('topright',
       c('State 1', 'State 2', 'Marginal'),
       col=pal,
       lwd=2,bty='n', inset=c(-0.1,0))
# AR(1,1)
hist(data$data[,1],breaks=35,xlab='step length', prob=T,
     xlim=c(0,60),ylim=c(0,0.045),main="AR(1,1)-HMM")
library(RColorBrewer)
pal=brewer.pal(3,'Dark2')
curve(0.5*dgamma(x,shape=10^2/5^2, scale=5^2/10),col=pal[1],lwd=2,add=T,
      from=0,to=60,n=1000)
curve(0.5*dgamma(x,shape=30^2/5^2, scale=5^2/30),col=pal[2],lwd=2,add=T,
      from=0,to=60,n=1000)
curve(0.5*dgamma(x,shape=10^2/5^2, scale=5^2/10)+
        0.5*dgamma(x,shape=30^2/5^2, scale=5^2/30), col=pal[3],lwd=2,add=T,
      from=0,to=60,n=1000)
#legend('topright',
#       c('State 1', 'State 2', 'Marginal'),
#       col=pal,
#       lwd=2,bty='n')



# circular visualization of von Mises distribution
par(mfrow=c(1,1))
pal=brewer.pal(4,'Dark2')
vm_sample=sample_arp(2000,delta=c(0.5,0.5),
                     Gamma=matrix(c(0.9,0.1,0.1,0.9),nrow=2),
                     N=2,
                     params=c(0,0,5,10),
                     autocor=list(matrix(c(0.5,0.6),nrow=2)),
                     p=c(1),
                     dists=c('vm'))
circ_vm_viz(mu=c(0,0),kappa=c(2,12), delta=c(0.5,0.5), data=vm_sample$data, 
            sum_dist=T, leg=T)

hist(vm_sample$data,breaks=20,main='',xlab='turning angle', prob=T,
     xlim=c(-pi,pi),ylim=c(0,0.95))
curve(0.5*dvm(x,mu=0, kappa=2),col=pal[1],lwd=1.5,add=T,
      from=-pi,to=pi,n=1000)
curve(0.5*dvm(x,mu=0, kappa=12),col=pal[2],lwd=1.5,add=T,
      from=-pi,to=pi,n=1000)
curve(0.5*dvm(x,mu=0, kappa=2)+
        0.5*dvm(x,mu=0, kappa=12), col=pal[3],lwd=1.5,add=T,
      from=-pi,to=pi,n=1000)
legend('topright',c('State 1', 'State 2', 'Marginal'),bty='n',
       col=pal[1:3], lty=1,lwd=1.5, inset=c(0.1,0.1))#, xpd='NA',inset=c(-0.45,0))


## Chapter 4 - Case Study 1

# see file

## Chapter 5 - Case Study 2

# see file

## Chapter 6 - Conclusion