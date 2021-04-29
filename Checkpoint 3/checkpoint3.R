rm(list=ls())

library(quantmod)
library(ggplot2)
library(reshape2)
library(ggpubr)
# install.packages("quantmod")
# install.packages("ggplot2")
# install.packages("reshape2")
# install.packages("ggpubr")

# ===============
# Get market data
# ===============
tickers <- c("SPY") # Get S&P 500 data
P.list <- lapply(tickers, function(x) get(getSymbols(x, from = "2001-04-28")))

sapply(P.list,nrow)

# Get the adjusted prices into a single object
P.adj <- lapply(P.list, function(p) p[,6])

# Merge the elements of the list
P <- Reduce(merge, P.adj)
names(P) <- tickers

# Use mean return to compute expected return 
# Calculate daily return 
returns <- P / lag(P) - 1

# Calculate mean return for all non-empty observations (We remove the the NA)
mu <- apply(returns, 2, function(x) mean(x, na.rm = TRUE)) 

# Calculate standard deviation
sigma <- apply(returns, 2, function(x) sd(x, na.rm = TRUE)) 

# Annualize the data 
mu <- (1 + mu) ** 252 - 1
sigma <- 252 * sigma

# Calculate volvol 
tickers <- c("VIX") # Get S&P 500 data
P.list <- lapply(tickers, function(x) get(getSymbols(x, from = "2001-04-28")))

sapply(P.list,nrow)

# Get the adjusted prices into a single object
P.adj <- lapply(P.list, function(p) p[,6])

# Merge the elements of the list
P <- Reduce(merge, P.adj)
names(P) <- tickers

# Use mean return to compute expected return 
# Calculate daily return 
returns <- P / lag(P) - 1

volvol <- apply(returns, 2, function(x) sd(x, na.rm = TRUE))


# ================
# Delta Hedge
# ================

delta.hedge.CRN <- function(M,N,S0,K,r,sigma,t,mu,alpha,b,volvol,V0,call, barrier = 53){
  # barrier is the argument for the barrier option. It is the "barrier" for when the option is in the money 
  # You don't need it other wise
  print(N)
  
  # Simulate the paths and deltas...also want to save Z's:
  #heston
  X <- V <- Zmat <- matrix(NA,ncol=N+1,nrow=M)
  X[,1] <- S0
  Zmat[,1] <- 0
  dt <- t/N
  sqdt <- sqrt(dt)
  V[,1] <- V0
  for (i in 1:N){
    Z1 <- matrix(rnorm(M),ncol=1) 
    Z2 <- matrix(rnorm(M),ncol=1)
    
    Zmat[, i + 1] <- Z1
    
    
    X[,i+1] <- X[,i] + mu*X[,i]*dt + sqrt(V[,i])*X[,i]*sqdt*Z1
    V[,i+1] <- V[,i] + alpha*(b - V[,i])*dt + volvol*sqrt(V[,i])*sqdt*Z2
    
  }
  
  #end Heston
  
  S <- X
  # european option price:
  getLastElement <- function(vec){
    return(vec[length(vec)])
  }
  
  Sbar.Euro <- apply(X, 1, getLastElement) - K # Get last value do K - apply(X, 1, getLastElement) for put 
  tmp <- which(Sbar.Euro > 0)
  Sbar.Euro[-tmp] <- 0
  P.Euro <- mean(Sbar.Euro)
  
  # Barrier option 
  barrierPrice <- function(vec, call = TRUE) { # CHANGE THIS FOR PUT
    if(call) {
      
      Smax <- max(vec) # If put Smin <- apply(S, 1, min)
      if(Smax >= barrier)
        return(0)
      return(vec[length(vec)])
    }
    else {
      Smin <- min(vec) # If put Smin <- apply(S, 1, min)
      if(Smin <= barrier)
        return(0)
      return(vec[length(vec)])
    }
  }
  
  Sbar.Barrier <- apply(X, 1, barrierPrice) - K # Get last value do K - apply(X, 1, barrier) for put 
  tmp <- which(Sbar.Barrier > 0)
  Sbar.Barrier[-tmp] <- 0
  P.Barrier <- mean(Sbar.Barrier)
  
  
  # ======================
  # Common Random Numbers
  # ======================
  sumS <- t(apply(S,1,cumsum))
  # 3. Payout function:
  if (call == 1){
    f <- function(S,K){
      f <- pmax(S - K,0)
    }
    premium<-P.Euro
  }else if(call == 2){  # Up and out barrier Euro call option { 
    Smax <- apply(S,1,max)
    f <- function(ST, K) {
      f <- pmax(ST - K,0) # First, treat like regular call
      f[Smax >= barrier] = 0 # Remove if Smax exceeds b
      return(f)
    }
  }else if(call == 3) {
    Smin <- apply(S, 1, min)
    f <- function(ST, K) {
      f <- pmax(K - ST, 0) # First, treat like regular put
      f[Smin <= barrier] = 0 # Remove is Smin is less than b
      return(f)
    }
  }else { # Euro Pu
    f <- function(S,K){
      f <- pmax(K - S,0)
    }
  }
  
  #4 now that the pay off function is set we can use a binomial lattice to calculate the options prices
    #calibrate the lattice with CRR
 # u <- exp(sigma*sqrt(dt))
 # d <- 1/u
 # p <- (exp(r*dt) - d) / (u - d)
 # q<-1-p
 # disc <- exp(-r*dt)
 # vec <- rep(NA,length=(2*N+1))
  
 # # Populate the terminal payoffs:
 # NN <- length(vec)
 # nu <- matrix(seq(N,0),ncol=1)
 # nd <- matrix(seq(0,N),ncol=1)
 # vec[seq(1,NN,by=2)] <- f(S0 * u ^ nu * d ^ nd,K)
  
 # # Solve by backward induction:
 # for (i in 1:N){
 #   opts <- vec[seq(i,NN-i+1,by=2)]
 #   vec[seq(1+i,NN-i,by=2)] <- disc*(p*opts[1:(length(opts)-1)] +
 #                                      q*opts[2:(length(opts))])
 # }
 # premium<-vec[N+1]
 # 
 # 4. Compute discount factors (for the deltas):
 #disc <- exp(-r*dt*rev(c(0:N)))
    
  # 5. Finally compute the deltas over time:
  deltas.A <- matrix(NA,ncol=N+1,nrow=M)
  
  delta_CRN <- function(S,V,K,r,sigma){
    fd<-rep(0,M)
    sqdt <- sqrt(dt)
    for (j in 1:M){
      ss <- sample(100000000,1)
      set.seed(ss)
      ST <- (S[j]-h) + mu*S[j]*dt + sqrt(V[j])*S[j]*sqdt*rnorm(M)
      set.seed(ss)
      STh <- (S[j]+h) + mu*S[j]*dt + sqrt(V[j])*S[j]*sqdt*rnorm(M)
      f0 <- f(ST,K)
      f0h <- f(STh,K)
      fd[j] <- exp(-r*t)*mean((f0h - f0) / (2*h))
    }
    return(fd)
  }
  for (i in 1:N){
    deltas.A[,i] <- delta_CRN(X[,i],V[,i],K,r,sigma)
  }
  
  # Fill in terminal deltas (1/0):
  temp<-f(X[,N+1],K)
  for(j in 1:M){
    if(temp[j]>1e-6){
      temp[j]<-1
    }else{
      temp[j]<-0
    }
  }
  deltas.A[,N+1] <- temp
  
  
  
  # Generate a matrix of positions:
  CF <- matrix(NA,ncol=N+1,nrow=M)
  CF[,1] <- -deltas.A[,1]*S0 + premium
  CF[,1] <- CF[,1] - abs(CF[i] * 0.01) # Account for trading costs (1%)
  for (i in 1:(N-1)){
    CF[,i+1] <- -1*(deltas.A[,i+1] - deltas.A[,i])*X[,i+1]
    CF[,i+1] <- CF[,i+1] - abs(CF[,i+1] * 0.01) # Account for trading costs (1%)
  }
  
  Xbar <- apply(X,1,mean)
  IN <- which(Xbar > K)
  
  CF[IN,N+1] <- K - Xbar[IN] + deltas.A[IN,N]*X[IN,N+1]
  CF[-IN,N+1] <- deltas.A[-IN,N]*X[-IN,N+1]
  
  
  # 3. sum the costs:
  disc <- matrix(exp(-r*seq(0,t,length=N+1)),ncol=1)
  PV <- CF%*%disc 
  
  
  # compute performance
  H.perf <- sqrt(var(PV))/P.Euro
  outlist <- list("H.perf"=H.perf,"PV"=PV,"P.Euro"=P.Euro,
                  "deltas.A"=deltas.A,"CF"=CF)
  return(outlist)
}


# ===============
# Test the function:
# M <- 500
# N <- 12
# t <- 20/52
# mu <- 0.13
# alpha<-2
# b<-0.16
# volvol<-2
# V0<-0.4
# S0 <- 49
# K <- 50
# r <- 0.05
# sigma <- 0.20
# h <- 0.1
# set.seed(2500)
S0 <- 49
K <- 50
r <- 0.05
t <- 20/52
M <- 500
N <- 12
h <- 0.1
alpha <- 2
b <- 0.16
V0 <- 0.4
set.seed(2021)
CRN.1 <- delta.hedge.CRN(M,N,S0,K,r,sigma,t,mu,alpha,b,volvol,V0,call=1)
CRN.means <- apply(CRN.1$deltas.A,2,mean)
CRN.sd <- sqrt(apply(CRN.1$deltas.A,2,var))
print(CRN.means)
print(CRN.sd)

# for testing code
N <- c(4,5,10,20,40,80,100)
H <- c(NA)
PV <- c(NA)
for (j in 1:length(N)){
  tmp <-  delta.hedge.CRN(M,N[j],S0,K,r,sigma,t,mu,alpha,b,volvol,V0,call=1)
  H[j] <- tmp$H.perf
  PV[j] <- mean(tmp$PV)
}
print(H)
print(PV)
