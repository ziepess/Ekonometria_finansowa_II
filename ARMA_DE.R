rm(list = ls())

require(eurostat)
require(zoo)
require(urca)
source("Block1Functions.R")
require(ggplot2)
require(knitr)
require(forecast)
require(urca)
require(quantmod)
require(zoo)
require(xts)
require(lubridate)
require(gridExtra)
require(grid)
require(tidyr)

#### przygotowanie szeregu ####
gdp_de <-get_eurostat("namq_10_gdp",
                      filters = list(geo = "DE",
                                     na_item = "B1GQ",
                                     unit = "CLV_PCH_SM",
                                     s_adj = "NSA"))

Y <- na.omit(zoo(gdp_de$values,gdp_de$time))

# roczne tempo wzrostu 
plot(Y, xlab = "Time", ylab = "GDP_yearly_change", main = "Roczna zmiana PKB Niemiec")
abline(h = mean(Y))


diff_Y <- diff(Y)

# Test ADF  - roczne tempo wzrostu 
summary(ur.df(Y, type="drift", lags=1))

summary(ur.df(diff_Y, type="drift", lags=1))


#### Wybranie specyfikacji modelu ####

Y_ts <- ts(coredata(Y),frequency = 4,start =c(1992,1) )

LagSel <- function(x, Pmax=3, Qmax=3, d=0, crit="BIC"){
  IC <- matrix(NA, Pmax+1, Qmax+1)
  for(p in 0:Pmax){
    if(p==Pmax){x0=x} else {x0 = x[-(1:(Pmax-p))]} 
    for(q in 0:Qmax){
      if(crit == "AIC"){ IC[p+1,q+1] <- AIC(Arima(x0,order=c(p,d,q)), k=2) }
      if(crit == "BIC"){ IC[p+1,q+1] <- AIC(Arima(x0,order=c(p,d,q)), k=log(length(x)-p)) }
      if(crit == "HQ"){  IC[p+1,q+1] <- AIC(Arima(x[- (1:(Pmax-p+1))],order=c(p,d,q)), k=2*log(log(length(x)-p))) }
    }}
  rownames(IC) <- paste('ar',0:Pmax, sep=""); colnames(IC) <- paste('ma',0:Qmax, sep="")
  return(IC)
}

# BIC
kable(LagSel(Y_ts, crit="BIC"),digits=2) 

# AIC
kable(LagSel(Y_ts, crit="AIC"),digits=2) 

arma_0 <- arima(Y_ts, order = c(3,0,3))  # AIC
arma_1 <- arima(Y_ts, order = c(3,0,1)) # BIC 


#### Wybor modelu ####

# ARMA(3,0,3) AIC
ar.coef_0 <- arma_0[["coef"]][1:3]
ma.coef_0 <- arma_0[["coef"]][3:6]

### wybrany model ARMA (3,1) BIC
ar.coef_1 <- arma_1[["coef"]][1:3]
ma.coef_1 <- arma_1[["coef"]][4:4]



roots=polyroot(c(1,-ar.coef_1))
1/roots

plot(arma_1)


#### IRF ####
Kmax = 35

IRF <- ARMAtoMA(ar = ar.coef_1,ma = ma.coef_1,Kmax); IRF = c(1,IRF)
IRFtab <- data.frame(H = 0:Kmax,IRF = IRF)

MyPlot(IRFtab, main="IRF from ARMA for Germany")



#### Prognoza ARMA(3,1) ####                  
fcst <- forecast(arma_1, h=12)
head(fcst)
plot(fcst, include=120, shaded=T, shadecols=c('grey','GhostWhite'), main="GDP Germany diff", ylab="", xlab="", bty="l")
abline(h=mean(Y),col="red")


means <- fcst$mean
means


# 2021 
mean(means[1:4])

# 2022 
mean(means[5:8])

# 2023 
mean(means[9:12])
