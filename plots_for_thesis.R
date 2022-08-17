# 2022-08-17
# Generation of plots used in the thesis
library(MasterThesis)

## Chapter 1 - Introduction

## Chapter 2 - HMMs

## Chapter 3 - Simulation Study

# Comparison of distributions and histograms
data = sample_arp(2000,
                  delta=c(0.5,0.5),
                  Gamma=matrix(c(0.9,0.1,0.1,0.9),ncol=2,byrow=TRUE),
                  N=2,
                  params=c(20,40,5,7),
                  autocor=list(matrix(c(0.45,0.55),ncol=1)),
                  p=c(1),
                  dists=c('gamma'))
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