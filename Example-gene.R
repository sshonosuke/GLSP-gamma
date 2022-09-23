###------------------------------------------------------###
###        R code for variance estimation of             ###
###               gene expression data                   ###
###------------------------------------------------------###
rm(list=ls())
set.seed(1)

## load dataset and R functions 
load("prostatedata.RData")
source("Gamma-GLSP-MCMC.R")

## data 
Y <- apply(prostatedata[,1:50], 1, var)
plot(Y)
Delta <- 49/2

## tuning parameters
a <- 2
b <- 1/2
mu <- 1
a_ta <- 1
b_ta <- 1
de <- 10^(-8) 
mc <- 3000
burn <- 1000

## SB prior
a_be <- 0.1
b_be <- 0.1
mcmc.SB <- SB_prior(y=Y, A_y=Delta, B_y=Delta, a_ta, b_ta, a_be, b_be, a, b, mc=mc, burn=burn)
est1 <- colMeans(1/mcmc.SB$La)
CI1 <- apply(1/mcmc.SB$La, 2, quantile, prob=c(0.025, 0.975))


## IRB prior
a_be <- 0.1
b_be <- 0.1
mcmc.SB <- IRB_prior(y=Y, A_y=Delta, B_y=Delta, a_ta, b_ta, a_be, b_be, a, b, mc=mc, burn=burn)
est2 <- colMeans(1/mcmc.SB$La)
CI2 <- apply(1/mcmc.SB$La, 2, quantile, prob=c(0.025, 0.975))


## Global (gamma) prior
mcmc.G <- global_prior(y=Y, A_y=Delta, B_y=Delta, a_ta, b_ta, a_be, b_be, mc=mc, burn=burn)[["La"]]   
est3 <- colMeans(mcmc.G)
CI3 <- apply(mcmc.G, 2, quantile, prob=c(0.025, 0.975))


## VASH method 
library(vashr)
fit <- vash(Y, df=2*Delta) 
est5 <- fit$sd.post



## Figure 
pdf("gene-hist-control.pdf", height=4, width=16, pointsize=18)
par(mfcol=c(1,4))
hist(Y, col=grey(0.9), main="SB", xlab="variance", ylim=c(0, 3000))
hist(est1, add=T, col="blue")
hist(Y, col=grey(0.9), main="IRB", xlab="variance", ylim=c(0, 3000))
hist(est2, add=T, col="blue")
hist(Y, col=grey(0.9), main="GL", xlab="variance", ylim=c(0, 3000))
hist(est3, add=T, col="blue")
hist(Y, col=grey(0.9), main="VS", xlab="variance", ylim=c(0, 3000))
hist(est5, add=T, col="blue")
dev.off()




## interval lengths 
mean(apply(CI1, 2, diff))
mean(apply(CI2, 2, diff))
mean(apply(CI3, 2, diff))


