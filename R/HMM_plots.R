# 2022-06-03 
# Functions that plot different HMMs


#' Plot states decoded by Viterbi algorithm
#' 
#' This function plots the states of an AR(p)-HMM decoded by the Viterbi algorithm.
#' 
#' @param states Vector of states (natural numbers, range 1-N).
#' @param names optional, specifies the names of the states. If "none", state 1-N will be used.
#' @param title bool, indicator if the plot should have a title.
#' 
#' @export
#' @rdname plot_states
plot_states <- function(states, names="none", title=TRUE){
  N <- length(unique(states))
  require(RColorBrewer)
  pal <- brewer.pal(N,'Dark2')
  
  if (title){
    plot(states, pch=19, bty='n', main='Decoded states using Viterbi',
         xlab='Index', ylab='State', col=pal[states])
  } else{
    plot(states, pch=19, bty='n', main='',
         xlab='Index', ylab='State', col=pal[states])
  }
  if (names[1]=="none"){
    legend(length(states)/1.5, 1.8, paste("State",1:N), col=pal,lwd=1.5,bty='n')
  } else{
    legend(length(states)/1.5, 1.8, names, col=pal,lwd=1.5,bty='n')
  }
}


#' Density plot of estimated distributions of AR(p)-HMM
#'
#' Plots the estimated distributions (only the state-dependent stationary values are considered) 
#' of a fitted AR(p)-HMM along with a histogram of the data.
#' 
#' @param data Input data for histogram.
#' @param dist string, indicates the kind of distribution, one of ['gamma', 'von Mises'].
#' @param param Vector of distribution parameters, in customary form, e.g. for 2 state gamma
#'              HMM: c(\eqn{\mu1,\mu2,\sigma1,\sigma2}).
#' @param N int, number of states (or weighted distributions) to be plotted.
#' @param delta Input vector for parameter delta.
#' @param title optional, title of the plot, if "none", then no title.
#' @param breaks optional, specify number of breaks of the histogram.
#' @param xlab optional, label of x-axis.
#' @param legend optional, plot legend (default: TRUE).
#' @param xlim optional, set a limit for x-axis for better visibility of the data.
#' 
#' @export
#' @rdname plot_fitted_dist
plot_fitted_dist <- function(data, dist, param, N, delta, title="none", breaks=30, xlab='x', legend=TRUE,
                             xlim='none'){
  require(CircStats)
  require(RColorBrewer)
  pal <- brewer.pal(N+1, 'Dark2')
  dens <- matrix(NA, nrow=10000,ncol=N)
  x <- seq(min(data[data!=0], na.rm=TRUE),max(data, na.rm=TRUE),length=10000)
  
  if (dist=='gamma'){
    mu <- param[1:N]
    sigma <- param[N+1:N]
    for (i in 1:N){
      dens[,i] = delta[i]*dgamma(x, shape=param$mu[i]^2/param$sigma[i]^2, scale=param$sigma[i]^2/param$mu[i])
    }
  } else if (dist=='vm'){
    mu <- param[1:N]
    kappa <- param[N+1:N]
    for (i in 1:N){
      dens[,i] = delta[i]*dvm(x, param$mu[i], param$kappa[i])
    }
  }
  total_dist <- apply(dens,1,sum)
  
  if (title=="none"){
    if (xlim[1] != 'none'){
      hist(data, breaks=breaks, probability = TRUE,
         main="", xlab=xlab, xlim=xlim, ylim=c(0,max(max(hist(data,plot=F,breaks=breaks)$density), 
                                         max(total_dist, na.rm=TRUE))))
    } else{
      hist(data, breaks=breaks, probability = TRUE,
           main="", xlab=xlab, ylim=c(0,max(max(hist(data,plot=F,breaks=breaks)$density), 
                                            max(total_dist, na.rm=TRUE))))
    }
  }else{
    if (xlim[1] != 'none'){
      hist(data, breaks=breaks, probability = TRUE,
           main="", xlab=xlab, xlim=xlim, ylim=c(0,max(max(hist(data,plot=F,breaks=breaks)$density), 
                                                       max(total_dist, na.rm=TRUE))))
    } else{
      hist(data, breaks=breaks, probability = TRUE,
         main=title, xlab=xlab, ylim=c(0,max(max(hist(data,plot=F,breaks=breaks)$density), 
                                            max(total_dist, na.rm=TRUE))))
    }
  }
  for (i in 1:N){
    lines(x,dens[,i],col=pal[i], lwd=2)
  }
  lines(x,total_dist,col=pal[N+1],lwd=2)
  
  if(legend){
  legend('topright', c(paste("State",1:N),"Total"), bty='n', lwd=2,
         col=pal)
  }
}


#' Plot decoded data 
#'
#' Plot a data vector, colored with regard to the Viterbi decoded states.
#' 
#' @param data Data vector to be plotted.
#' @param states Vector of decoded states.
#' @param col Color of the states.
#' @param name optional, y-axis label, if "none" than "data".
#' @param title optional, title of the plot, if "none", then no title.
#' @param legend optional, bool, specifies if there should be a legend.
#' 
#' @export
#' @rdname plot_decoded_data
plot_decoded_data <- function(data, states, col, name="none", title="none", legend=TRUE){
  if (name=="none"){
    plot(NULL,xlim=c(0,length(data)),ylim=c(min(data,na.rm=TRUE),max(data,na.rm=TRUE)),
         ylab='data',xlab='time',bty='n')
  } else{
    plot(data, xlim=c(0,length(data)),ylim=c(min(data,na.rm=TRUE),max(data,na.rm=TRUE)),
         ylab=name,xlab='time',bty='n', col='grey', type='l')
  }
  points(x=1:length(data),y=data, col=col[states], cex=0.75, pch=19)
  if (title=="none"){
    title(main="")
  } else{
    title(main=title)
  }
  if (legend){
    legend('topright',
           c('state 1', 'state 2'),
           col=col,lwd=1.5,
           bty='n')
  }  
}


#' Boxplot of different simulation runs
#'
#' Generate a boxplot of several fitted values, generated by fitting a model
#' to data. Includes a comparison to the true value.
#' 
#' @param data Simulated values of the parameter.
#' @param name Name of the parameter.
#' @param true_value optional, true parameter value, 0 if no horizontal line should 
#' be drawn.
#'                        
#' @export
#' @rdname boxplot_params
#' 
boxplot_params <- function(data, name,true_value,...){
  args=list(x=data,xlab=name,...)
  do.call(boxplot,args)
  if (true_value){
    abline(h=true_value,lwd=2,col=2)
  }
}


#' Plot data 
#'
#' Plot a data vector.
#' 
#' @param data Data vector to be plotted.
#' @param name optional, y-axis label, if "none" than "data".
#' @param title optional, title of the plot, if "none", then no title.
#' 
#' @export
#' @rdname plot_data
plot_data <- function(data, name="none", title="none"){
  if (name=="none"){
    plot(data, typ='l',xlab="Index", ylab="Data",bty='n')
  } else{
    plot(data, typ='l',xlab="Index", ylab=name,bty='n')
  }
  if (title=="none"){
    title(main="")
  } else{
    title(main=title)
  }
  
}


#' Full simulation boxplots of estimated parameters
#'
#' Generate a boxplot of all fitted parameters, generated by fitting an AR(p)-HMM
#' to data. Includes a comparison to the true value. Currently only works
#' for a simulation containing one gamma distributed and one von Mises distributed
#' variable.
#' 
#' @param gamma_mu Fitted values of the parameter mu of the gamma distribution.
#' @param gamma_sigma Fitted values of the parameter sigma of the gamma distribution.
#' @param vm_mu Fitted values of the parameter mu of the von Mises distribution.
#' @param vm_kappa Fitted values of the parameter kappa of the von Mises distribution.
#' @param gamma_auto Fitted values of the autoregression parameters of the gamma distribution. 
#'                   0 if no autoregression was fitted.
#' @param vm_auto Fitted values of the autoregression parameters of the von Mises distribution.
#'                   0 if no autoregression was fitted.
#' @param true_values True parameter values, order: [gamma_mu, gamma_sigma, vm_mu, vm_kappa, gamma_auto, vm_auto]
#' @param nstates Number of states fitted (default: 2). Works only for 2 states!
#' @param vm_head_pos optional, determines position of heading of von Mises distribution. 
#'                    Only gets used for fitted models without autoregression
#'                        
#' @export
#' @rdname boxplot_sim_params
#' 
boxplot_sim_params <- function(gamma_mu, gamma_sigma, vm_mu, vm_kappa, gamma_auto=0,
                               vm_auto=0, true_values, nstates=2, vm_head_pos=0){
  
  if (length(true_values)==(4*nstates)){
    layout(matrix(c(1,2,3,4),ncol=2,byrow=TRUE))
  } else{
    layout(matrix(c(1,2,3,3),ncol=2,byrow=TRUE))
  }
  
  # boxplot gamma_mu
  boxplot(x=gamma_mu, xlab=expression(mu), cex.lab=2,xaxt='n',cex.axis=1.5,
          pch=19, ylim=c(min(min(true_values[1:nstates]),min(gamma_mu)),
                         max(max(true_values[1:nstates]),max(gamma_mu))
                        )
          )
  axis(1, at=1:nstates, labels=rep("",nstates))
  mtext(c(expression(mu[1]),expression(mu[2])),side=1,at=1:nstates,line=1.25,cex=1.5)
  abline(h=true_values[1:nstates],col=2,lwd=2)
  
  # boxplot gamma_sigma
  boxplot(x=gamma_sigma, xlab=expression(sigma), cex.lab=2,xaxt='n',cex.axis=1.5,
          pch=19, ylim=c(min(min(true_values[nstates+1:nstates]),min(gamma_sigma)),
                          max(max(true_values[nstates+1:nstates]),max(gamma_sigma))
          ))
  axis(1, at=1:nstates, labels=rep("",nstates))
  mtext(c(expression(sigma[1]),expression(sigma[2])),side=1,at=1:nstates,line=1.25,cex=1.5)
  abline(h=true_values[nstates+1:nstates],col=2,lwd=2)
  
  # boxplot gamma_auto
  if (length(true_values)>(4*nstates)){
    boxplot(gamma_auto,xlab=expression(phi),cex.lab=2,cex.axis=1.5,ylim=c(0,1),xaxt='n',
            pch=19)
    axis(1, at=1:((length(true_values)-4*nstates)/2), labels=rep("",(length(true_values)-4*nstates)/2))
    if ((length(true_values)-4*nstates)/(nstates*2) == 1){ # AR(1)
      mtext(c(expression(phi["1"]),expression(phi["2"])),
            side=1,at=1:nstates,line=1.25,cex=1.5)
      abline(h=true_values[4*nstates+1:nstates],col=2,lwd=2)
    } else if ((length(true_values)-4*nstates)/(nstates*2) == 2){ # AR(2)
      mtext(c(expression(phi["1,1"]),expression(phi["1,2"]),expression(phi["2,1"]),expression(phi["2,2"])),
            side=1,at=1:(2*nstates),line=1.25,cex=1.5)
      abline(h=true_values[4*nstates+1:(2*nstates)],col=2,lwd=2)
    } else if ((length(true_values)-4*nstates)/(nstates*2) == 3){ # AR(3)
      mtext(c(expression(phi["1,1"]),expression(phi["1,2"]),expression(phi["1,3"]),expression(phi["2,1"]),expression(phi["2,2"]),expression(phi["2,3"])),
            side=1,at=1:(3*nstates),line=1.25,cex=1.5)
      abline(h=true_values[4*nstates+1:(3*nstates)],col=2,lwd=2)
    } else{
      return("ERROR: Wrong length of true_values.")
    }
  }
  title("Estimated parameters of gamma distribution",outer=TRUE,line=-2,cex.main=2)
  
  # boxplot vm_mu
  boxplot(x=vm_mu, xlab=expression(mu), cex.lab=2,xaxt='n',cex.axis=1.5,
          pch=19, ylim=c(min(min(true_values[2*nstates+1:nstates]),min(vm_mu)),
                         max(max(true_values[2*nstates+1:nstates]),max(vm_mu))
          ))
  axis(1, at=1:nstates, labels=rep("",nstates))
  mtext(c(expression(mu[1]),expression(mu[2])),side=1,at=1:nstates,line=1.25,cex=1.5)
  abline(h=true_values[2*nstates+1:nstates],col=2,lwd=2)
  
  # boxplot vm_kappa
  boxplot(x=vm_kappa, xlab=expression(kappa), cex.lab=2,xaxt='n',cex.axis=1.5,
          pch=19, ylim=c(min(min(true_values[3*nstates+1:nstates]),min(vm_kappa)),
                         max(max(true_values[3*nstates+1:nstates]),max(vm_kappa))
          ))
  axis(1, at=1:nstates, labels=rep("",nstates))
  mtext(c(expression(kappa[1]),expression(kappa[2])),side=1,at=1:nstates,line=1.25,cex=1.5)
  abline(h=true_values[3*nstates+1:nstates],col=2,lwd=2)
  
  # boxplot vm_auto
  if (length(true_values)>(4*nstates)){
    boxplot(vm_auto,xlab=expression(phi),cex.lab=2,cex.axis=1.5,ylim=c(0,1),xaxt='n',
            pch=19)
    axis(1, at=1:((length(true_values)-4*nstates)/2), labels=rep("",(length(true_values)-4*nstates)/2))
    if ((length(true_values)-4*nstates)/(nstates*2) == 1){ # AR(1)
      mtext(c(expression(phi["1"]),expression(phi["2"])),
            side=1,at=1:nstates,line=1.25,cex=1.5)
      abline(h=true_values[4*nstates+nstates+1:nstates],col=2,lwd=2)
    } else if ((length(true_values)-4*nstates)/(nstates*2) == 2){ # AR(2)
      mtext(c(expression(phi["1,1"]),expression(phi["1,2"]),expression(phi["2,1"]),expression(phi["2,2"])),
            side=1,at=1:(2*nstates),line=1.25,cex=1.5)
      abline(h=true_values[4*nstates+2*nstates+1:(2*nstates)],col=2,lwd=2)
    } else if ((length(true_values)-4*nstates)/(nstates*2) == 3){ # AR(3)
      mtext(c(expression(phi["1,1"]),expression(phi["1,2"]),expression(phi["1,3"]),expression(phi["2,1"]),expression(phi["2,2"]),expression(phi["2,3"])),
            side=1,at=1:(3*nstates),line=1.25,cex=1.5)
      abline(h=true_values[4*nstates+3*nstates+1:(3*nstates)],col=2,lwd=2)
    } else{
      return("ERROR: Wrong length of true_values.")
    }
  title("Estimated parameters of von Mises distribution",outer=TRUE,line=-2,cex.main=2)
  } else{
    title("Estimated parameters of von Mises distribution",outer=TRUE,line=vm_head_pos,cex.main=2)
  }
  
}



#' Circular KDE and density of von Mises distribution
#'
#' Generate a circular visualization of the density of one or several von Mises
#' distributed variables. Also, if data is supplied, generate the kernel density
#' estimate of the data for comparison.
#' 
#' @param mu Vector of parameter mu for von Mises distribution.
#' @param kappa Vector of parameter mu for von Mises distribution.
#' @param delta Optional, weights of the distributions (equal to delta in HMMs).
#' @param data Optional, provide data vector, KDE of this data will be plotted.
#' @param sum_dist Optional, sum weighted distributions.
#' @param leg Optional, include legend.
#'                        
#' @export
#' @rdname circ_vm_viz
#' 
circ_vm_viz <- function(mu, kappa, delta=1, data=NULL, sum_dist=FALSE, leg=FALSE){
  require(CircStats)
  require(RColorBrewer)
  
  # create unit circle
  radius = 1
  low=-pi
  up=pi
  n=1000
  x=seq(low,up,length=n)
  
  # transform polar coordinates to cartesian coordinates
  polar2cart <- function(r, theta) {
    data.frame(x = r * cos(theta), y = r * sin(theta))
  }
  
  # von mises pdf
  number_dists = length(mu)
  vm_pdfs = matrix(NA,nrow=n,ncol=number_dists)
  vm_pdf.cart.data.x = matrix(NA,nrow=n,ncol=number_dists)
  vm_pdf.cart.data.y = matrix(NA,nrow=n,ncol=number_dists)
  for (dist in 1:number_dists){
    vm_pdf=dvm(x,mu[dist],kappa[dist])
    vm_pdfs[,dist] = delta[dist]*vm_pdf
    polar.data = polar2cart(radius+vm_pdfs[,dist],x)
    vm_pdf.cart.data.x[,dist]=polar.data[,1]
    vm_pdf.cart.data.y[,dist]=polar.data[,2]
  }
  total_dist = rowSums(vm_pdfs)
  total.data = polar2cart(radius+total_dist,x)

  # plot unit circle
  center_x = 0
  center_y = 0
  theta = seq(0, 2 * pi, length = 200)
  plot(c(min(-1.3,min(total.data[,1])), max(1.3,max(total.data[,1]))), c(min(-1.3,min(total.data[,2])), max(1.3,max(total.data[,2]))), 
       type = "n",bty='n',xaxt='n',yaxt='n',xlab='',ylab='')
  lines(x = radius * cos(theta) + center_x, y = radius * sin(theta) + center_y,
        lwd=3)
  lines(c(-1.1,1.1),c(0,0))
  lines(c(0,0),c(-1.1,1.1))
  text(1.2,0, '0')
  text(-1.2,0, expression(pi))
  text(0, 1.2, expression(pi/2))
  text(0, -1.2, expression(-pi/2))
  
  # plot densities
  cols=brewer.pal(number_dists+2, 'Dark2')
  for (dist in 1:number_dists){
    points(vm_pdf.cart.data.x[,dist], vm_pdf.cart.data.y[,dist], type='l',col=cols[dist],
           lwd=1.5)
  }
  
  if (sum_dist){
    points(total.data, type='l',col=cols[number_dists+1],
           lwd=1.5)
  }
  
  if (!is.null(data)){
    require(circular)
    # kde of data
    circ.dens = density.circular(x=data, from=circular(-pi), to=circular(pi), bw=bw.nrd.circular(data),
                                 na.rm=TRUE)
    points(polar2cart(radius+circ.dens$y,circ.dens$x),type='l',col=tail(cols,1),
           lwd=1.5)
  }
  
  # legend
  if (leg){
    if (sum_dist){
      legend('topright',
             c(paste('State',1:number_dists),'Marginal', 'KDE'),
             lwd=1.5,bty='n',
             col=cols)
    } else{
      legend('topright',
             c(paste('State',1:number_dists),'Sampled data'),
             lwd=1.5,bty='n',
             col=cols)
    }
  }
}




