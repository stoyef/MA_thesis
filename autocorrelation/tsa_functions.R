#setwd(getSrcDirectory()[1]) # set wd to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # RStudio way

## read in a data frame of the seeschwalbe data
read_schwalbe <- function(data_path){
  schwalbe = read.csv(data_path, header = T)
  return(schwalbe)
}

## function to select a schwalbe set from the data
select_schwalbe <- function(set_name, count, verbose=T){
  schwalbe = schwalben[schwalben$ID== set_name,]
  if (verbose){
    cat("Schwalbe", count,":", dim(schwalbe)[1], "Observations\n")
    plot(schwalbe$x, schwalbe$y, type='l', bty='n', main=paste("Schwalbe",count))
  }
  return(schwalbe)
}

# function to cluster a schwalbe data set with k-means++
# create new column where the cluster number is denoted
cluster_schwalbe <- function(data, no_clusters){
  clusters = kmeans(data[2:(nrow(data)-1),c(2,3,8,10)], no_clusters)
  data$cluster = c(NA, clusters$cluster, NA)
  return(data)
}

# function to extract the different parts in the data where a cluster occurs
create_cluster_batches <- function(data){
  rows = as.numeric(rownames(data))
  data$batch = c(1,1+cumsum(as.numeric(diff(as.numeric(rownames(data)))>1)))
  return(data)
}

# function to extract a batch from a cluster and to generate ACF and PACF
# of this batch
# (for tortuosity)
batch_autocor <- function(data, batch_no){
  batch = data[which(data$batch == batch_no),]
  batch_acf = acf(batch$tortuosity, lag.max = 25, plot=FALSE)
  batch_pacf = pacf(batch$tortuosity, lag.max = 25, plot=FALSE)
  return(list(batch_acf = batch_acf, batch_pacf = batch_pacf))
}

# function to plot ACF and PACF of a batch of a cluster
plot_autocor <- function(batch_acf, batch_pacf){
  par(mfrow=c(2,1))
  plot(batch_acf, ylim=c(-1,1), type='h')
  plot(batch_pacf, ylim=c(-1,1), type='h')
  par(mfrow=c(1,1))
}


