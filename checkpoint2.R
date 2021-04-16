# Checkpoint 2

rm(list=ls())

library(quantmod)
library(ggplot2)
library(reshape2)
library(ggpubr)


callPayoff <- function(s, k)
  return(pmax(s- k, 0))

putPayoff <- function(s, k)
  return(pmax(k - s, 0))

set.seed(2021)

# Define 10 companies by ticker (Change these later)
tickers <- c("PLUG", "BP", "DISH", "T", "CRM", "XOM", "GE", "TSM", "MUFG", "ASML", "AAPL")
P.list <- lapply(tickers, function(x) get(getSymbols(x, from = "2010-01-01")))

sapply(P.list,nrow)

# Get the adjusted prices into a single object
P.adj <- lapply(P.list, function(p) p[,6])

# Merge the elements of the list
P <- Reduce(merge, P.adj)
names(P) <- tickers


# Use mean return to compute expected return 
# Calculate daily return 
returns <- P / lag(P) - 1
# Chop off top NA
returns <- returns[-1, ]

# Calculate mean return for all non-empty observations (We remove the the NA)
mu <- apply(returns, 2, function(x) mean(x, na.rm = TRUE)) 

# Calculate standard deviation
sigma <- apply(returns, 2, function(x) sd(x, na.rm = TRUE))

# Calculate the covariance matrix. (Thee built in functions make this so easy...)
Sigma <- var(returns, use = "pairwise")

# Annualize the data 
mu <- (1 + mu) ** 252 - 1
Sigma <- 252 * Sigma

# Calculate correlation 
rho <- cor(returns)


GBMD <- function(M, N, d, t, mu ,X0, sigma, rho)
{
  X <- array(NA, dim = c(M, N+1, d))
  
  X[, 1, ] <- X0
  B <- t(chol(Sigma))
  dt <- t/N
  muMat <- matrix(rep(mu - (1 / 2) * sigma * sigma, M), ncol = d, byrow = TRUE)
  
  for (i in 1:N)
  {
    Zmat <- matrix(rnorm(d*M),ncol=M)
    X[, i+1, ] <- X[, i,] * exp(muMat * dt + sqrt(dt) * t(B %*% Zmat))
  }
  
  return(X)
}

M <- 10000
N <- 52
d <- length(tickers)
t <- 1
X0 <- 100

X <- GBMD(M,N, d, t, mu, X0, sigma, rho)

# X.means <- matrix(NA,ncol=d+4,nrow=(N+1)) 
# X.vars <- matrix(NA,ncol=d+4,nrow=(N+1)) 
# X.means[,1] <- seq(0,t,length=(N+1)) 
# X.vars[,1] <- seq(0,t,length=(N+1))
# tvec <- seq(0,t,length=(N+1))
# for(i in 1:d)
# {
#   X.means[,i+1] <-apply(X[,,i],2,mean) 
#   X.vars[,i+1] <-apply(X[,,i],2,var)
#   X.means[,i+d+1] <- matrix(X0*exp(mu[i]*tvec),ncol=1)
#   X.vars[,i+d+1] <- matrix(X0^2*exp(2*mu[i]*tvec)*(exp(sigma[i]^2*tvec)-1),ncol=1)
# }
# X.means <- as.data.frame(X.means)
# names(X.means) <- c("time","d1","d2","d3","d1true","d2true","d3true") 
# dat.means <- melt(X.means, id.vars="time")
# # Everything on the same plot
# ggplot(dat.means, aes(time,value, col=variable)) + geom_point() +
#   stat_smooth()
