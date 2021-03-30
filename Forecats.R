rm(list=ls())
require(eurostat)
require(zoo)
require(xts)
require(lubridate)
require(forecast)
require(vars)
require(knitr)
require(ggplot2)
require(grid)
require(gridExtra)
require(tidyr)
source("Block1Functions.R")

temp <- get_eurostat("namq_10_gdp",
                     filters = list(geo = c("DE","FR"),
                                    na_item = "B1GQ",
                                    unit = "CLV_PCH_SM",
                                    s_adj = "SCA"))
tempF <- subset(temp, geo == "FR")
tempH <- subset(temp, geo == "DE")


pF  <- na.omit(zoo(tempF$values,order.by=tempF$time))
pH  <- na.omit(zoo(tempH$values,order.by=tempH$time))
# together
z0 <- merge(pF,pH,all=FALSE)
z  <- z0

plot(z, screens=c(1,1),col=1:2,main="Inflation in DE and FR", type = "l", ylab="",xlab="", bty="l")
legend("bottomleft",legend=c("UK","DE"),col=1:2,lty=1, bty="n")

# 2. Calculating ex-post forecast
#####################

# settings
T          <- dim(z)[1] 	    # number of observations
M          <- 20              # Forecast evaluluation sample size
T1         <- T-M             # sample split (1st fct for T1+1)
horiz      <- 4              # forecast horizon
varN       <- 2               # which variable 

# Calculating forecasts and realizations
#########################################

act      <- Actfct(z[,varN], horiz, T1)
fctRW    <- RWfct( z[,varN], horiz, T1)
fctARMA  <- ARIMAfct(z[,varN], horiz, T1, 1, 0, 1)
fctVAR   <- VARfct(z, horiz, T1, varN, 2)

# 3. Sequential forecast
#####################################

# Choose model "RW", "ARIMA", VAR"
fct   = fctARMA

# make a graph
datesF = c(index(z),last(index(z)) + months(1:horiz))
inf_data    = data.frame(t=index(z[,varN]),y=coredata(z[,varN]))
temp_plot   <- MyPlot(inf_data,main="Sequential forecast",xlab="",hline=mean(z[,varN]))
for (n in 0:(T-T1-1)) {
  fct_data = data.frame(yf = fct[n+1,], tf = datesF[(T1+n+1):(T1+n+horiz)])
  temp_plot <- temp_plot + geom_line(data=fct_data, aes(y=yf,x=tf), colour="red", size=1)
}
temp_plot + 
xlim(as.Date("2010-01-01"),as.Date("2022-01-01"))

# 4. Effectiveness test
###################################

h     = 12                         # Horizon
fct   = fctVAR                      # Choose model "RW", "ARIMA", "VAR"

ylevel  <- act[1:(T-T1-h+1),h] 
xlevel  <- fct[1:(T-T1-h+1),h] 
ychange <- act[1:(T-T1-h+1),h] - z[T1:(T-h),varN]
xchange <- fct[1:(T-T1-h+1),h] - z[T1:(T-h),varN]

par(mfrow=c(1,2), cex = 0.7, bty="l")
a1 <- min(cbind(xlevel,ylevel))-0.01; a2 <- max(cbind(xlevel,ylevel))+0.01; 
plot(xlevel,ylevel, main="level", xlab="forecast", ylab="realization", pch=19, col="red",  xlim=c(a1,a2), ylim = c(a1,a2))
abline(0,1,lwd=2)
a1 <- min(cbind(xchange,ychange))-0.01; a2 <- max(cbind(xchange,ychange))+0.01; 
plot(xchange,ychange, main="change", xlab="foreast", ylab="realization", pch=19, col = "blue", xlim=c(a1,a2), ylim = c(a1,a2))
abline(0,1,lwd=2)

fit <- lm(ylevel ~ xlevel)
summary(fit)
fit <- lm(ychange ~ xchange)
summary(fit)

require(car)
linearHypothesis(fit, hypothesis.matrix = diag(1,2), 
                 rhs=c(0,1), test="F") 


# 5. Table with statistics
##############################

errRW     <- act - fctRW 
errARMA   <- act - fctARMA  
errVAR    <- act - fctVAR 

statRW      <- FctErrStats(errRW)
statARMA    <- FctErrStats(errARMA)
statVAR     <- FctErrStats(errVAR)

dmRW       <- rep(NA,h)
dmARMA     <- DMstat(errARMA, errRW, type="two.sided", power=2)
dmVAR      <- DMstat(errVAR, errRW, type="two.sided", power=2)

ME   <- rbind(statRW$ME,statARMA$ME,statVAR$ME )
RMSE <- rbind(statRW$RMSE,statARMA$RMSE,statVAR$RMSE)
DM   <- rbind(dmRW,dmARMA,dmVAR)

models <- c("RW", "ARIMA", "VAR")
rownames(ME)   <- models; colnames(ME)   <- paste(1:horiz)
rownames(RMSE) <- models; colnames(RMSE) <- paste(1:horiz)
rownames(DM)   <- models; colnames(DM)   <- paste(1:horiz)

Hsel = c(1,2,3,6,9,12)
#Hsel = c(1,2,3,4)



options(digits=3)
ME[,Hsel]
RMSE[,Hsel]
DM[,Hsel]

RMSE1 = RMSE / matrix(rep(RMSE[1,],3),3,h,byrow=TRUE)
RMSE1[1,]=RMSE[1,]
RMSE1[,Hsel]

# Diebold-Mariano test 

h   = 12
e1h = errRW[,h]
e2h = errVAR[,h]
dh  = na.omit(e1h^2 - e2h^2)
plot.ts(dh)
abline(h=0,col=2)
abline(h=mean(dh),col=3)

dhat = mean(dh)
s    = sd(dh)
Th   = length(dh)

require(sandwich)
sig2lr = lrvar(dh, type = "Newey-West")
sig2   = var(dh)/Th

DM = dhat/sqrt(sig2lr)
DM
pnorm(-abs(DM))
