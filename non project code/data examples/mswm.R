# 2022-08-25
# Try package MSwM

library(MSwM)

data=schwalbe_77$step
plot(data,type='l')

form = glm(data~1,family=Gamma(link='log'))
mod = msmFit(form, k=2, sw=rep(TRUE,1), family='gamma')
summary(mod)

# comparison to my function
theta=c(-2,-2,10,25,5,5)
theta.star=starize(theta,2,0,'gamma')
my_mod = fit_arp_model(mllk, data, theta.star, 2, 0, 'gamma')
my_mod

# comparison to moveHMM
library(moveHMM)
movehmm_mod = fitHMM(schwalbe_77, 2, c(10,25,5,5),angleDist = 'none')
movehmm_mod

####
#### Now, AR(1)

form = glm(data~1,family=Gamma(link='log'))
mod = msmFit(form, k=2, p=2, sw=rep(TRUE,3), family='gamma')
summary(mod)
# doesn't work, the autoregression parameters seem to be larger than 1???

# comparison to my function
theta=c(-2,-2,10,25,5,5,0.1,0.1)
theta.star=starize(theta,2,1,'gamma')
my_mod = fit_arp_model(mllk, data, theta.star, 2, 1, 'gamma')
my_mod
decoded_states = viterbi_arp(data, my_mod$Gamma, my_mod$delta, 'gamma', my_mod$autocorrelation, 
            my_mod$params, 2, 1)

library(RColorBrewer)
pal=brewer.pal(3,'Dark2')

plot_decoded_data(data, col=pal)
     