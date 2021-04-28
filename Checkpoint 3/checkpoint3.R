rm(list=ls())

library(quantmod)
library(ggplot2)
library(reshape2)
library(ggpubr)
# install.packages("quantmod")
# install.packages("ggplot2")
# install.packages("reshape2")
# install.packages("ggpubr")




# ================
# Delta Hedge
# ================

delta.hedge.CRN <- function(M,N,S0,K,r,sigma,t,mu,alpha,b,volvol,V0,call, barrier = 53){
  # barrier is the argument for the barrier option. It is the "barrier" for when the option is in the money 
  # You don't need it other wise
  print(N)
  
  # Simulate the paths and deltas...also want to save Z's:
  #heston
  X <- V <- deltas <- Zmat <- matrix(NA,ncol=N+1,nrow=M)
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
  
  
  Sbar.Euro <- apply(X, 1, getLastElement) - K # Get last value
  tmp <- which(Sbar.Euro > 0)
  Sbar.Euro[-tmp] <- 0
  P.Euro <- mean(Sbar.Euro)
  
  
  
  
  # Sbar.A <- apply(X,1,mean) - K
  # tmp <- which(Sbar.A > 0)
  # Sbar.A[-tmp] <- 0
  # P.Asian <- mean(Sbar.A)
  
  
  # ======================
  # Common Random Numbers
  # ======================
  sumS <- t(apply(S,1,cumsum))
  # 3. Payout function:
  if (call == 1){
    f <- function(S,K){
      f <- pmax(S - K,0)
    }
  }
  else if(call == 2){  # Up and out barrier Euro call option { 
    Smax <- apply(S,1,max)
    f <- function(ST, Smax, K) {
      f <- pmax(ST - K,0) # First, treat like regular call
      f[Smax >= barrier] = 0 # Remove if Smax exceeds b
      return(f)
    }
  }
  else if(call == 3) {
    Smin <- apply(S, 1, min)
    f <- function(ST, Smin, K) {
      f <- pmax(K - ST, 0) # First, treat like regular put
      f[Smin <= barrier] = 0 # Remove is Smin is less than b
      return(f)
    }
  }
  else { # Euro Pu
    f <- function(S,K){
      f <- pmax(K - S,0)
    }
  }
  
  # 4. Compute discount factors (for the deltas):
  disc <- exp(-r*dt*rev(c(0:N)))
  
  # 5. Finally compute the deltas over time:
  deltas.A <- matrix(NA,ncol=N+1,nrow=M)
  for (i in 1:N){
    
    # 6. Save a matrix of the cumulative sums of Z from t_j to T
    if (i < N){
      message(print(i))
      Zmat.sum <- matrix(NA,nrow=M,ncol=N+1)
      Zmat.sum[,c((i+1):(N+1))] <- t(apply(Zmat[,c((i+1):(N+1))],1,cumsum))
    }
    else{
      Zmat.sum <- matrix(NA,nrow=M,ncol=N+1)
      Zmat.sum[,c((i+1):(N+1))] <- Zmat[,N+1]
    }
    
    # 7. Save the innovatinos (last two terms of Eq. (2))
    RHS <- matrix(NA,nrow=M,ncol=N+1)
    for (k in i:N){
      RHS[,k+1] <- (k-i+1)*(r - 0.5*sigma*sigma)*dt + sigma*sqrt(dt)*Zmat.sum[,k+1]
    }
    
    # 8. For each row, use the RHS matrix to compute the forward prices
    # and combine that with sumS[] to get Sbar.
    for (j in 1:M){
      S2 <- S2h <- matrix(NA,nrow=M,ncol=N+1)
      S2 <- exp(RHS + log(S[j,i]))
      S2[,i] <- S[j,i]
      Sbar <- (apply(S2[,c(i:(N+1))],1,sum) + sumS[j,i] - S[j,i])/(N+1)
      p1 <- mean(disc[i]*f(Sbar,K))
      
      S2h <- exp(RHS + log(S[j,i]+h))
      S2h[,i] <- S[j,i] + h
      Sbarh <- (apply(S2h[,c(i:(N+1))],1,sum) + sumS[j,i] - S[j,i])/(N+1)
      p2 <- mean(disc[i]*f(Sbarh,K))
      
      deltas.A[j,i] <- (p2 - p1)/h
    }
  }
  
  
  
  # Generate a matrix of positions:
  CF <- matrix(NA,ncol=N+1,nrow=M)
  CF[,1] <- -deltas.A[,1]*S0
  for (i in 1:(N-1)){
    CF[,i+1] <- -1*(deltas.A[,i+1] - deltas.A[,i])*X[,i+1]
  }
  
  Xbar <- apply(X,1,mean)
  IN <- which(Xbar > K)
  
  CF[IN,N+1] <- K - Xbar[IN] + deltas.A[IN,N]*X[IN,N+1]
  CF[-IN,N+1] <- deltas.A[-IN,N]*X[-IN,N+1]
  
  # Account for trading costs (1%)
  for(i in (1:length(CF)))
    CF[i] <- CF[i] - CF[i] * 0.01
  
  
  # 3. sum the costs:
  disc <- matrix(exp(-r*seq(0,t,length=N+1)),ncol=1)
  PV <- CF%*%disc 
  
  
  # compute performace
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
sigma <- 0.20
t <- 20/52
mu <- 0.13
M <- 500
N <- 12
h <- 0.1
alpha <- 2
b <- 0.16
V0 <- 0.4
volvol <- 0.3 
set.seed(2021)
CRN.1 <- delta.hedge.CRN(M,N,S0,K,r,sigma,t,mu,alpha,b,volvol,V0,call=1)
CRN.means <- apply(CRN.1$deltas.A,2,mean)
CRN.sd <- sqrt(apply(CRN.1$deltas.A,2,var))
print(CRN.means)
print(CRN.sd)