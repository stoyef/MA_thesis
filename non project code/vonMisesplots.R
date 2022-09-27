# 2022-09-04
# Visualization of von Mises distribution

library(CircStats)

mu=0
kappa=10
data.vm <- rvm(1000, mu+pi, 10)-pi
low=-pi
up=pi
n=1000
x=seq(low,up,length=n)
pdf.vm=dvm(x,mu,kappa)

# plot circle
# prepare "circle data"
radius = 1
center_x = 0
center_y = 0
theta = seq(0, 2 * pi, length = 200) # angles for drawing points around the circle
# initialize a plot
#plot(c(-2, 2), c(-2, 2), type = "n",bty='n',xaxt='n',yaxt='n',xlab='',ylab='')
# draw the circle
#lines(x = radius * cos(theta) + center_x, y = radius * sin(theta) + center_y,
#      lwd=3)
#lines(c(-1.1,1.1),c(0,0))
#lines(c(0,0),c(-1.1,1.1))
#text(1.2,0, '0')
#text(-1.2,0, expression(pi))
#text(0, 1.2, expression(pi/2))
#text(0, -1.2, expression(-pi/2))


polar2cart <- function(r, theta) {
  data.frame(x = r * cos(theta), y = r * sin(theta))
}

cart2polar <- function(x, y) {
  data.frame(r = sqrt(x^2 + y^2), theta = atan2(y, x))
}

data.cart = polar2cart(radius,data.vm)
#points(data.cart, pch=19)


vm.pdf.data = polar2cart(radius+pdf.vm,x)
#points(vm.pdf.data, type='l',col=4)

# To avoid skewed bars in a histogram, we plot the kernel density estimate of the sampled data 
# instead
# KDE of data
data.kde = density(data.vm)
#points(polar2cart(radius+data.kde$y,data.kde$x),type='l',col=3)

library(circular)

circ.dens = density.circular(x=data.vm, from=circular(-pi), to=circular(pi), bw=bw.nrd.circular(data.vm),
                             na.rm=TRUE)
#circ.dens$x
#circ.dens$y
#points(polar2cart(radius+circ.dens$y,circ.dens$x),type='l',col=3)


circ_vm_viz(mu=c(0,2),kappa=c(5,12))
circ_vm_viz(mu=c(0,2),kappa=c(5,12), data=rvm(1000, 0+pi, 5)-pi)
circ_vm_viz(mu=c(0,2),kappa=c(5,12), data=rvm(1000, 0+pi, 5)-pi, leg=TRUE)

# now, for our case
vm_sample=sample_arp(2000,delta=c(0.5,0.5),
                     Gamma=matrix(c(0.9,0.1,0.1,0.9),nrow=2),
                     N=2,
                     params=c(0,0,2,12),
                     autocor=list(matrix(c(0.5,0.6),nrow=2)),
                     p=c(1),
                     dists=c('vm'))
circ_vm_viz(mu=c(0,0),kappa=c(2,12), delta=c(0.5,0.5), data=vm_sample$data, 
            sum_dist=T, leg=T)

