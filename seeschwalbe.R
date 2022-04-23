#setwd(getSrcDirectory()[1]) # set wd to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # RStudio way

source("tsa_functions.R")

schwalben = read_schwalbe('data/seeschwalbe/slick2-h1-h2.csv')

head(schwalben)

count=1
for (i in unique(schwalben$ID)){
  assign(paste("schwalbe_", count, sep=""), select_schwalbe(i, count))
  count=count+1
}


# die 3. nehmen wir

library(moveHMM)

schwalbe_3 = prepData(schwalbe_3)
head(schwalbe_3)
schwalbe_3_hmm_2_states = fitHMM(schwalbe_3, 2, c(25, 30 , 5, 10), c(0,0,2,5))
schwalbe_3_hmm_2_states
plot(schwalbe_3_hmm_2_states)

plot(schwalbe_3$tortuosity, type='l', bty='n')


####
#### 2 states clustern
####
# mit k-means++
schwalbe_3 = cluster_schwalbe(schwalbe_3, 2)
schwalbe_3[which(schwalbe_3$cluster==1),]
plot(schwalbe_3$x, schwalbe_3$y, col=schwalbe_3$cluster+2)
plot(schwalbe_3$tortuosity, col=schwalbe_3$cluster+2, type='h')

cluster_1 = schwalbe_3[which(schwalbe_3$cluster==1),]
cluster_2 = schwalbe_3[which(schwalbe_3$cluster==2),]
plot(cluster_1$x, cluster_1$y,pch=19) # works
plot(cluster_2$x, cluster_2$y,pch=19) # works

mean(cluster_2$tortuosity)
mean(cluster_1$tortuosity)

cluster_1 = create_cluster_batches(cluster_1)
cluster_2 = create_cluster_batches(cluster_2)


####
#### AR-Prozess
####


# cluster 1 -> das interessante Cluster
## wie modelliert man die unterbrochenen Strecken? 
cluster_1_longest = batch_autocor(data = cluster_1, batch_no = 4)
#cluster_1_longest
# 1. Ansatz: Separate Modellierung -- wir nehmen nur die größte zusammenhängende
# Strecke -> Datenpunkte

plot_autocor(cluster_1_longest$batch_acf, cluster_1_longest$batch_pacf)
# points to AR(2)

# cluster 2 -> das andere Cluster
cluster_2_longest = batch_autocor(data = cluster_2, batch_no = 7)
plot_autocor(cluster_2_longest$batch_acf, cluster_2_longest$batch_pacf)

# keine systematische Autokorrelation zu erkennen


####
#### Schwalbe 77
####

# Longest data set is schwalbe_77
plot(schwalbe_77$x, schwalbe_77$y, type='l', bty='n', main="Schwalbe 77")

## Prep, 2-state-HMM
schwalbe_77 = prepData(schwalbe_77)
head(schwalbe_77)
schwalbe_77_hmm_2_states = fitHMM(schwalbe_77, 2, c(25, 30 , 5, 10), c(0,0,2,5))
schwalbe_77_hmm_2_states
plot(schwalbe_77_hmm_2_states)
plot(schwalbe_77$tortuosity, type='l', bty='n')

####
#### 2 states clustern
####
# mit k-means++
schwalbe_77 = cluster_schwalbe(schwalbe_77, 2)
schwalbe_77[which(schwalbe_77$cluster==1),]
plot(schwalbe_77$x, schwalbe_77$y, col=schwalbe_77$cluster+2)
plot(schwalbe_77$tortuosity, col=schwalbe_77$cluster+2, type='h')

cluster_1 = schwalbe_77[which(schwalbe_77$cluster==1),]
cluster_2 = schwalbe_77[which(schwalbe_77$cluster==2),]
plot(cluster_1$x, cluster_1$y,pch=19) # works
plot(cluster_2$x, cluster_2$y,pch=19) # works

mean(cluster_2$tortuosity)
mean(cluster_1$tortuosity)

mean(cluster_1$speed)
mean(cluster_2$speed)

cluster_1 = create_cluster_batches(cluster_1) # high tortuosity
cluster_2 = create_cluster_batches(cluster_2)
# klappt

####
#### AR-Prozess
####


# cluster 1 -> das interessante Cluster
cluster_1_longest = batch_autocor(data = cluster_1, batch_no = 1)
#cluster_1_longest
# 1. Ansatz: Separate Modellierung -- wir nehmen nur die größte zusammenhängende
# Strecke -> Datenpunkte

plot_autocor(cluster_1_longest$batch_acf, cluster_1_longest$batch_pacf)
# Nicht ganz eindeutig welcher AR(x) aber auf jeden Fall Autokorrelation

# cluster 2 -> das andere Cluster
cluster_2_longest = batch_autocor(data = cluster_2, batch_no = 8)
plot_autocor(cluster_2_longest$batch_acf, cluster_2_longest$batch_pacf)
# schwierig


