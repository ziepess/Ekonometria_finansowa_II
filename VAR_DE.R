#########################################
# Meeting 3.                            # 
# VAR models                            #
#########################################
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

#############################
# 1. Loading data           #
#############################
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
z <- merge(pF,pH,all=FALSE)

y  <- ts(coredata(z), frequency=4,start=c(year(z[1]),month(z[1])))

plot(z, screens=c(1,1),col=1:2,main="GDP growth in DE and FR", type = "l", ylab="",xlab="", bty="l")
legend("bottomleft",legend=c("EA","DE"),col=1:2,lty=1, bty="n")
#############################
# 2. VAR model              #
#############################

VARselect(z, lag.max=6, type="const")
p    = 1
VARE <- VAR(z, p=p, type="const")

# Testing for autocorrelation: test Ljunga-Boxa
serial.test(VARE, lags.pt = 6, type="PT.adjusted")

par(mfrow=c(2,2), cex = 0.7, bty="l", mar=c(3,3,3,3))
Acf(residuals(VARE)[,1], xlab="",main="e1 vs lags of e1")
Ccf(residuals(VARE)[,1],residuals(VARE)[,2],xlab="",main="e1 vs lags of e2", xlim=c(0,24))
Acf(residuals(VARE)[,2],xlab="",main="e2 vs lags of e2")
Ccf(residuals(VARE)[,2],residuals(VARE)[,1],xlab="",main="e2 vs lags of e1", xlim=c(0,24))

##################################
# 3. Recursive structuralization #
##################################

# causality
causality(VARE, cause="pF")
causality(VARE, cause="pH")

# recursive structuralization 
bmat <- matrix(0,2,2)
bmat[lower.tri(bmat, diag=T)] <- NA
SVARE <- SVAR(VARE, Amat = NULL, Bmat = bmat)
SVARE
B <- SVARE$B

#####################
# 4. structural IRF #
#####################
K    <- 48
SVMA <- Phi(SVARE, nstep=K)

# short run-impact matrix
SVMA[,,1]

# Plot in Vars package
SIRF <- irf(SVARE, n.ahead=K, cum = F, boot = F)
plot(SIRF)

# Single IRF
IRFtab = data.frame(H=0:K,IRF=SVMA[1,1,])
MyPlot(IRFtab, main= "IRF from uEA to infEA")

# Panel of IRFs  
require(grid)
require(gridExtra)
IRFtab$IRF = SVMA[1,1,]; irf11 <- MyPlot(IRFtab, main = "IRF from uFR to infFR")
IRFtab$IRF = SVMA[1,2,]; irf12 <- MyPlot(IRFtab, main = "IRF from uDE to infFR")
IRFtab$IRF = SVMA[2,1,]; irf21 <- MyPlot(IRFtab, main = "IRF from uFR to infDE")
IRFtab$IRF = SVMA[2,2,]; irf22 <- MyPlot(IRFtab, main = "IRF from uDE to infDE")
grid.arrange(irf11,irf12,irf21,irf22, ncol=2)

# where do we see these values?
B

#####################
# 5. FEVD           #
#####################
i      = 2            # which variables
SIRFi  = SVMA[i,,]    # IRF for variable i
temp   = SIRFi^2      # contribution of shock[t] to forecast variance at [t+h]
FVARi  = apply(temp,1,cumsum) # contribution of shocks[t:t+h] to forecast variance at [t+h]
FEVD0  <- prop.table(FVARi,1); colnames(FEVD0) = c("uFR","uDE")
FEVD1  <- data.frame(H = 0:K, eEA = FEVD0[,1], eDE = FEVD0[,2])

kable(FEVD1[c(1:5,13,25,49),],digits=3,row.names = FALSE)

require(tidyr)
FEVD2 <- gather(data = FEVD1, key = shocks, value = value, -c(H))
ggplot(FEVD2, aes(fill=shocks, y=value, x=H)) + 
  geom_bar( stat="identity") +   
  theme_light()+
  labs(title="FEVD for GDP growth in Germany", y="", x="", caption="")+
  scale_fill_manual(values=c("grey70", "grey30"))

# function in the package
temp <- fevd(SVARE,n.ahead=60)
temp$pH
plot(temp)


#################################
# 6. Historical decomposition   #
#################################

i      = 2            # which variables
# a. Structural shocks
e <- residuals(VARE)
B <- SVARE$B
u <- t(solve(B)%*%t(e))
T <- dim(u)[1]
# cor(u) 

# b. IRFs
SVMA <- Phi(SVARE, nstep=T)
SIRF  = t(SVMA[i,,])    # IRF for variable i

# c. Historical decompositoion
HistDec             <- matrix(NA,T,2)
colnames(HistDec)   <- c("uFR","uDE")

for(t in 1:T){
  junk1 <- as.matrix(u[1:t,])
  junk2 <- as.matrix(SIRF[t:1,])
  HistDec[t,] <- colSums(junk1*junk2)
}

# the impact of initial conditions
InitCond <- tail(coredata(z[,i]),T) - colSums(t(HistDec)) 
mu       <- as.numeric(tail(InitCond,1))
InitCond <- InitCond - mu
HistDec  <- cbind(HistDec, InitCond)

HD1   <- as.data.frame(HistDec); HD1$t   <- tail(index(z[,i]),T)
HD2   <- gather(data = HD1, key = shocks, value = value, -c(t))
Y     <- data.frame(y = tail(coredata(z[,i]),T)-mu,t=tail(index(z[,i]),T))

ggplot() + 
  geom_bar(data=HD2, aes(fill=shocks, y=value, x=t), stat="identity",width=50) + # 
  scale_fill_manual(values=c("grey80", "red", "green"))+
  geom_line(data=Y, aes(y=y, x=t), size=1.5)+
  theme_light()+
  labs(title="Historical decomposition for GDP growth in Germany", y="", x="", caption="")

#################
# 7. A forecast #
#################

H   <- 20
T   <- dim(z)[1]

fct <- predict(VARE, n.ahead=H)
fct
fanchart(fct, nc=1, xlim=c(T-60,T+H),bty="l",mar=c(3,3,1,1))

# fan chart in ggplot
datesf <- last(index(z))+months(1:H)
f      <- zoo(fct$fcst$pH[,1],datesf)

i = 2                                                       # select variable
inf_data  = data.frame(t=index(z),y = coredata(z[,i]))
inf_fct   = data.frame(t=datesf  ,y = fct$fcst[[i]][,1]) 
point_fct = rbind(inf_data,inf_fct)
interval_fct = data.frame(t=datesf,lower = fct$fcst[[i]][,2],upper = fct$fcst[[i]][,3]) 
MyPlot(point_fct,main="Forecast for GDP growth from VAR",xlab="",hline=mean(z[,i]))+
  geom_ribbon(data= interval_fct, aes(x=t,ymin=lower, ymax=upper), alpha=0.2) +
  xlim(as.Date("2010-01-01"),as.Date("2026-01-01"))

# Annual forecasts
require(xts)
require(lubridate)
datesf <- seq.Date(from = as.Date("2021-01-01"), to =as.Date("2025-10-01"),by = "quarter"  ))
temp     <- zoo(fct$fcst$pH[,1],datesf)
xf = c(z[,2],temp)

kable(apply.yearly(xf,mean),digits=2)
