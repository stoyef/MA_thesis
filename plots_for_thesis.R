# 2022-08-17
# Generation of plots used in the thesis
library(MasterThesis)

## Chapter 1 - Introduction

## Chapter 2 - HMMs

## Chapter 3 - Simulation Study

data = sample_arp(2000,
                  delta=c(0.5,0.5),
                  Gamma=matrix(c(0.9,0.1,0.1,0.9),ncol=2,byrow=TRUE),
                  N=2,
                  params=c(20,40,5,7,0,0,2,12),
                  autocor=list(matrix(c(0.45,0.55),ncol=1),matrix(c(0.5,0.6),ncol=1)),
                  p=c(1,1),
                  dists=c('gamma','vm'))

# Plot simulated data as time series
# AR(0)
data_ar0 = sample_arp(2000,
                      delta=c(0.5,0.5),
                      Gamma=matrix(c(0.9,0.1,0.1,0.9),ncol=2,byrow=TRUE),
                      N=2,
                      params=c(20,40,5,7,0,0,2,12),
                      autocor=0,
                      p=c(0,0),
                      dists=c('gamma','vm'))
par(mfrow=c(2,1), mar=c(5.1,4.1,4.1,2.1), xpd=TRUE)
plot_decoded_data(data_ar0$data[1:300,1], data_ar0$states[1:300], col=pal, name='step length', legend=FALSE)
plot_decoded_data(data_ar0$data[1:300,2], data_ar0$states[1:300], col=pal, name='turning angle', legend=FALSE)
legend('top',c('state 1','state 2'), col=pal[1:2],lwd=1.5,bty='n',
       y.intersp = 0.75,
       horiz=TRUE, inset=c(0,-0.75))

# AR(1)
par(mfrow=c(2,1), mar=c(5.1,4.1,4.1,2.1), xpd=TRUE)
plot_decoded_data(data$data[1:300,1], data$states[1:300], col=pal, name='step length', legend=FALSE)
plot_decoded_data(data$data[1:300,2], data$states[1:300], col=pal, name='turning angle', legend=FALSE)
legend('top',c('state 1','state 2'), col=pal[1:2],lwd=1.5,bty='n',
       y.intersp = 0.75,
       horiz=TRUE, inset=c(0,-0.75))

# AR(2)
data_ar2 = sample_arp(2000,
                      delta=c(0.5,0.5),
                      Gamma=matrix(c(0.9,0.1,0.1,0.9),ncol=2,byrow=TRUE),
                      N=2,
                      params=c(20,40,5,7,0,0,2,12),
                      autocor=list(matrix(c(0.2,0.3,0.3,0.4),ncol=2, byrow=T),
                                   matrix(c(0.2,0.3,0.3,0.4),ncol=2,byrow=T)),
                      p=c(2,2),
                      dists=c('gamma','vm'))
par(mfrow=c(2,1), mar=c(5.1,4.1,4.1,2.1), xpd=TRUE)
plot_decoded_data(data_ar2$data[1:300,1], data_ar2$states[1:300], col=pal, name='step length', legend=FALSE)
plot_decoded_data(data_ar2$data[1:300,2], data_ar2$states[1:300], col=pal, name='turning angle', legend=FALSE)
legend('top',c('state 1','state 2'), col=pal[1:2],lwd=1.5,bty='n',
       y.intersp = 0.75,
       horiz=TRUE, inset=c(0,-0.75))

# AR(3)
data_ar3 = sample_arp(2000,
                  delta=c(0.5,0.5),
                  Gamma=matrix(c(0.9,0.1,0.1,0.9),ncol=2,byrow=TRUE),
                  N=2,
                  params=c(20,40,5,7,0,0,2,12),
                  autocor=list(matrix(c(0.1,0.1,0.55,0.1,0.1,0.65),ncol=3, byrow=T),
                               matrix(c(0.1,0.1,0.5,0.1,0.1,0.6),ncol=3,byrow=T)),
                  p=c(3,3),
                  dists=c('gamma','vm'))
par(mfrow=c(2,1), mar=c(5.1,4.1,4.1,2.1), xpd=TRUE)
plot_decoded_data(data_ar3$data[1:300,1], data_ar3$states[1:300], col=pal, name='step length', legend=FALSE)
plot_decoded_data(data_ar3$data[1:300,2], data_ar3$states[1:300], col=pal, name='turning angle', legend=FALSE)
legend('top',c('state 1','state 2'), col=pal[1:2],lwd=1.5,bty='n',
       y.intersp = 0.75,
       horiz=TRUE, inset=c(0,-0.75))
dev.off()



# Comparison of distributions and histograms
data$data
hist(data$data[,1],breaks=35,main='',xlab='x', prob=T,
     xlim=c(0,80),ylim=c(0,0.045))
library(RColorBrewer)
pal=brewer.pal(3,'Dark2')
curve(0.5*dgamma(x,shape=20^2/5^2, scale=5^2/20),col=pal[1],lwd=2,add=T,
      from=0,to=80,n=1000)
curve(0.5*dgamma(x,shape=40^2/7^2, scale=7^2/40),col=pal[2],lwd=2,add=T,
      from=0,to=80,n=1000)

## Chapter 4 - Case Study 1

## Chapter 5 - Case Study 2

## Chapter 6 - Conclusion