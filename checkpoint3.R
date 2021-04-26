rm(list=ls())

library(quantmod)
library(ggplot2)
library(reshape2)
library(ggpubr)

# ================
# Delta Hedge
# ================
delta.hedge.PW <- function(M,N,S0,K,r,sigma,t,mu){
  print(N)
  
  # Simulate the paths and deltas...also want to save Z's:
  
  # Heston
  X <- V <- deltas <- Zmat <- matrix(NA,ncol=N+1,nrow=M)
  Zmat[,1] <- 0
  X[, 1] <- S0
  dt <- t / N
  sqdt <- sqrt(dt)
  X[,1] <- S0
  V[,1] <- V0
  for (i in 1:N)
  {
    Z1 <- matrix(rnorm(M),ncol=1) 
    Z2 <- matrix(rnorm(M),ncol=1)
    
    Zmat[, i + 1] <- Z1
    
    X[,i+1] <- X[,i] + mu*X[,i]*dt + sqrt(V[,i])*X[,i]*sqdt*Z1
    V[,i+1] <- V[,i] + alpha*(b - V[,i])*dt + volvol*sqrt(V[,i])*sqdt*Z2
  }
  # End Heston
  
  # Asian option price:
  Sbar.A <- apply(X,1,mean) - K
  tmp <- which(Sbar.A > 0)
  Sbar.A[-tmp] <- 0
  P.Asian <- mean(Sbar.A)
  
  
  # ======================
  # Pathwise Estimator
  # ======================
    sumS <- t(apply(X,1,cumsum))
    S <- X
    # 3. Payout function:
    f <- function(S,K)
    {
      f <- pmax(S - K,0)
      return(f)
    }
    
    # 4. Compute discount factors (for the deltas):
    disc <- exp(-r*dt*rev(c(0:N)))
    
    # 5. Finally compute the deltas over time:
    deltas.A <- matrix(NA,ncol=N+1,nrow=M)
    for (i in 1:N)
    {
      if (i < N)
      {
        Zmat.sum <- matrix(NA,nrow=M,ncol=N+1)
        Zmat.sum[,c((i+1):(N+1))] <- t(apply(Zmat[,c((i+1):(N+1))],1,cumsum))
      }
      else
      {
        Zmat.sum <- matrix(NA,nrow=M,ncol=N+1)
        Zmat.sum[,c((i+1):(N+1))] <- Zmat[,N+1]
      }
      RHS <- matrix(NA,nrow=M,ncol=N+1)
      for (k in i:N)
      {
        RHS[,k+1] <- (k-i+1)*(r - 0.5*sigma*sigma)*dt + sigma*sqrt(dt)*Zmat.sum[,k+1]
      }
      for (j in 1:M)
      {
        S2 <- matrix(NA,nrow=M,ncol=N+1)
        S2 <- exp(RHS + log(S[j,i])) # Matrix of stock prices going forwward
        S2[,i] <- S[j,i] # Initial price
        Sbar <- (apply(S2[,c(i:(N+1))],1,sum) + sumS[j,i] - S[j,i])/(N+1)
        dSdJ <- (apply(S2[,c(i:(N+1))],1,sum))/(N+1)
        I <- rep(0,length=M)
        tmp <- which(Sbar > K)
        I[tmp] <- 1
        deltas.A[j,i] <- disc[i]*mean(I*dSdJ/S2[j,i])
      }
    }

  
  # Generate a matrix of positions:
  CF <- matrix(NA,ncol=N+1,nrow=M)
  CF[,1] <- -deltas.A[,1]*S0
  for (i in 1:(N-1))
  {
    CF[,i+1] <- -1*(deltas.A[,i+1] - deltas.A[,i])*X[,i+1]
  }
  
  Xbar <- apply(X,1,mean)
  IN <- which(Xbar > K)
  
  CF[IN,N+1] <- K - Xbar[IN] + deltas.A[IN,N]*X[IN,N+1]
  CF[-IN,N+1] <- deltas.A[-IN,N]*X[-IN,N+1]
  
  # 3. sum the costs:
  disc <- matrix(exp(-r*seq(0,t,length=N+1)),ncol=1)
  PV <- CF%*%disc
  
  # compute performace
  H.perf <- sqrt(var(PV))/P.Asian
  outlist <- list("H.perf"=H.perf,"PV"=PV,"P.Asian"=P.Asian,
                  "deltas.A"=deltas.A,"CF"=CF)
  return(outlist)
}


# ===============
# Test the function:


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
PW.1 <- delta.hedge.PW(M,N,S0,K,r,sigma,t,mu)
PW.means <- apply(PW.1$deltas.A,2,mean)
PW.sd <- sqrt(apply(PW.1$deltas.A,2,var))

print(PW.means)
print(PW.sd)