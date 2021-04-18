rm(list=ls())
# Stochastic Volatility Model
#  Heston's model
# This model assumes random volatility
# and returns the stochastic volatility and stock price of each simulated path
myHeston <- function(M,N,t,mu,alpha,b,volvol,S0,V0){ 
    dt <- t/N
    sqdt <- sqrt(dt)
    S <- V <- matrix(NA,ncol=N+1,nrow=M)
    S[,1] <- S0
    V[,1] <- V0
    for (i in 1:N){
      Z1 <- matrix(rnorm(M),ncol=1) 
      Z2 <- matrix(rnorm(M),ncol=1)
      S[,i+1] <- S[,i] + mu*S[,i]*dt + sqrt(V[,i])*S[,i]*sqdt*Z1
      V[,i+1] <- V[,i] + alpha*(b - V[,i])*dt + volvol*sqrt(V[,i])*sqdt*Z2
    }
    out <- list("S"=S,"V"=V) 
    return(out)
}
# parameters:
M <- 10000
N <- 252
t <- 1
mu <- 0.1 
alpha <- 2
b <- 0.16 
volvol <-0.3 
V0 <- 0.4
S0 <- 100
H <- myHeston(M,N,t,mu,alpha,b,volvol,S0,V0)

time <- as.matrix(seq(0,t,length=N+1),ncol=1) 
S1 <- t(H$S)
V1 <- t(H$V)
df1 <- as.data.frame(cbind(time,S1))
df2 <- as.data.frame(cbind(time,V1))
# Plot some paths:
mP <- 100 #limit the number of paths being plotted 
Anames <- rep("",mP)
for (i in 1:mP){
  Anames[i] <- paste("A",i,sep="") 
}
names(df1) <- c("time", Anames)
names(df2) <- c("time", Anames)

# This creates a new data frame with columns x, variable and value 
# x is the id, variable holds each of our timeseries designation 
library("reshape")
library("ggplot2")
df1.melted <- melt(df1[,c(1:mP)], id = "time")
df2.melted <- melt(df2[,c(1:mP)], id = "time")
plot.stocks <- ggplot() +
  geom_line(data = df1.melted, aes(x = time, y = value, factor=variable),color="red") + xlab('time') +
  ylab('Stock Price')

print(plot.stocks)


plot.vol <- ggplot() +
  geom_line(data = df2.melted, aes(x = time, y = value, factor=variable),color="black") + xlab('time') +
  ylab('Volatility')




print(plot.vol)

hist(H$S[,N+1])
mean(H$S[,N+1])





