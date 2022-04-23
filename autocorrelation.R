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


schwalbe_77_tort = ts(schwalbe_77$tortuosity)

plot(schwalbe_77_tort,type='l')
pacf(schwalbe_77_tort)
lag1_schwalbe_77 <- lag(schwalbe_77_tort, -1)
plot(lag1_schwalbe_77, schwalbe_77_tort)
lagdata <- ts.intersect(schwalbe_77_tort, lag1_schwalbe_77, dframe=T)
mod <- lm(schwalbe_77_tort ~ lag1_schwalbe_77, data=lagdata)
summary(mod)
plot(mod)
plot(schwalbe_77_tort,type='l')
lines(predict(mod, lagdata), col=2, lwd=1)
plot(residuals(mod, response=T), type='h')
pacf(residuals(mod, response=T))
# AR(1) doesn't do the job

lag2_schwalbe_77 <- lag(schwalbe_77_tort, -2)
lagdata <- ts.intersect(schwalbe_77_tort, lag1_schwalbe_77, lag2_schwalbe_77,
                        dframe=T)
mod <- lm(schwalbe_77_tort ~ lag1_schwalbe_77 + lag2_schwalbe_77, data=lagdata)
summary(mod)
plot(mod)
plot(schwalbe_77_tort,type='l')
lines(predict(mod, lagdata), col=2, lwd=1)
plot(residuals(mod, response=T), type='h')
pacf(residuals(mod, response=T))
# AR(2) auch nicht besser

lag3_schwalbe_77 <- lag(schwalbe_77_tort, -3)
lagdata <- ts.intersect(schwalbe_77_tort, lag1_schwalbe_77, lag2_schwalbe_77,
                        lag3_schwalbe_77, dframe=T)
mod <- lm(schwalbe_77_tort ~ lag1_schwalbe_77 + lag2_schwalbe_77 + lag3_schwalbe_77, 
          data=lagdata)
summary(mod)
plot(mod)
plot(schwalbe_77_tort,type='l')
lines(predict(mod, lagdata), col=2, lwd=1)
plot(residuals(mod, response=T), type='h')
pacf(residuals(mod, response=T))
# Manno, auch nicht besser

####
#### Besser: alles auf einmal
####
for (order_ar in 1:6){
    crit = arima(schwalbe_77_tort, c(order_ar,0,0))$aic
    cat("AIC of AR(", order_ar, ") model: ", crit, "\n")
  }
# AIC findet AR(1) hier am besten


## ARIMA
for (order_ar in 1:6){
  for (order_ma in 0:5){
    crit = arima(schwalbe_77_tort, c(order_ar,0,order_ma))$aic
    cat("AIC of AR(", order_ar, ") and MA(", order_ma, ") model: ", crit, "\n")
  }
}
# AIC ist immer noch für AR(1)


#### Was ist, wenn wir uns nicht den gesameten Datensatz angucken, 
# sondern nur einen Batch aus cluster 1?
library(moveHMM)
schwalbe_77 = prepData(schwalbe_77)
schwalbe_77 = cluster_schwalbe(schwalbe_77, 2)
cluster_1 = schwalbe_77[which(schwalbe_77$cluster==1),]
cluster_2 = schwalbe_77[which(schwalbe_77$cluster==2),]
cluster_1 = create_cluster_batches(cluster_1) # high tortuosity
cluster_2 = create_cluster_batches(cluster_2)
cluster_1_longest = batch_autocor(data = cluster_1, batch_no = 1)
schwalbe_77_batch1 = cluster_1[which(cluster_1$batch==1),]
tort = ts(schwalbe_77_batch1$tortuosity)
plot(tort, type='l')

## ARIMA
for (order_ar in 0:6){
  for (order_ma in 0:6){
    crit = arima(tort, c(order_ar,0,order_ma))$aic
    cat("AIC of AR(", order_ar, ") and MA(", order_ma, ") model: ", crit, "\n")
  }
}
# AIC will sehr komplexes Modell
for (order_ar in 0:6){
  for (order_ma in 0:6){
    crit = BIC(arima(tort, c(order_ar,0,order_ma)))
    cat("BIC of AR(", order_ar, ") and MA(", order_ma, ") model: ", crit, "\n")
  }
}
# BIC will AR(4) & MA(2) Modell

# -> Da ist auf jeden Fall was dran, dass man für die Zustände separat
# AR(p)-Koeffizienten schätzen sollte




