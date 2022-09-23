library(GIGrvg)
###--------------------------------------------------------------###
###        R code for implementing Scaled beta (SB) prior        ###
###            and inverse rescaled beta (IRB) prior             ###
###--------------------------------------------------------------###

##     Function of Miller's sampling algorithm
## (used in MCMC algorithm in the proposed priors)
miller <- function(x, mu, a0, b0, ep, M){
  n <- length(x)
  R <- sum(log(x))
  S <- sum(x)
  T. <- S / mu - R + n * log(mu) - n
  A <- a0 + n / 2; B <- b0 + T.
  for(j in 1:M){
    a <- A / B
    A <- a0 - n * a + n * a^2 * trigamma(a)
    B <- b0 + (A - a0) / a - n * log(a) + n * digamma(a) + T.
    if(abs(a / (A / B) - 1) < ep){
      return(c(A = A, B = B, j = j))
    }
  }
  return(c(A = A, B = B, j = M + 1))
}



###    Inverse Rescaled Beta (IRB) prior    ###
## model settings 
# y \sim Ga(A_y, B_y/lambda)
# lambda \sim Ga(1+tau*u, beta*tau*u)
# u \sim IRSB(b, a)
# tau \sim Ga(a_ta, b_ta)
# beta \sim Ga(a_be, b_be)
## IMPUT
# y: sequence of observed values
# A_y, B_y: fixed values in the sampling model 
# (a_ta, b_ta): hyperparameter in prior of tau
# (a_be, b_be): hyperparameter in prior of beta
# a, b: shape parameters of IRB prior 
# ep, M: tuning parameters in the Miller's algorithm
# mc: MCMC length 
# burn: burn-in period
## OUTPUT
# posterior samples of parameters 

IRB_prior <- function(y, A_y, B_y, a_ta, b_ta, a_be, b_be, a, b, ep=10^(-8), M=10, mc, burn){
  # preparation
  m <- length(y)
  # initial value
  la <- A_y / (B_y * y)
  v <- rep(NA, m)
  w <- rep(NA, m)
  z <- rep(NA, m)
  eta <- rep(1, m)
  ta <- 1
  be <- 1
  
  # objects to store posterior samples
  La <- matrix(NA, mc, m)
  V <- matrix(NA, mc, m)
  W <- matrix(NA, mc, m)
  Z <- matrix(NA, mc, m)
  Eta <- matrix(NA, mc, m)
  Ta <- rep(NA, mc)
  Be <- rep(NA, mc)
  
  # MCMC 
  for(iota in 1:mc){
    # la
    la <- rgamma(m, shape = A_y + 1 + eta, rate = B_y * y + be * eta)
    La[iota, ] <- la
    
    # (v, w, z)
    v <- rgamma(m, shape = 1 - b, rate = log(1 + ta / eta))
    w <- rgamma(m, shape = b + a, rate = 1 + log(1 + ta / eta))
    z <- rgamma(m, v + w + 1, rate = 1 + eta / ta)
    V[iota, ] <- v
    W[iota, ] <- w
    Z[iota, ] <- z
    
    # be
    be <- rgamma(1, shape = a_be + sum(eta + 1), rate = b_be + sum(eta * la))
    Be[iota] <- be
    
    # ta
    ta <- rgig(1, lambda = a_ta - sum(v + w), chi = 2 * sum(z * eta), psi = 2 * b_ta)
    Ta[iota] <- ta
    
    # eta
    ABJ <- apply(rbind(la, v + w, z / ta), 2, function(s) miller(x = s[1], mu = 1 / be, a0 = s[2], b0 = s[3], ep = ep, M = M))
    A <- ABJ[1, ]
    B <- ABJ[2, ]
    eta_proposal <- rgamma(m, shape = A, rate = B)
    uniform <- runif(m, min = 0, max = 1)
    logratio <- ((v + w - 1) * log(eta_proposal) - z * eta_proposal / ta + eta_proposal * log(be * eta_proposal) - lgamma(eta_proposal) + eta_proposal * log(la) - la * be * eta_proposal - (A - 1) * log(eta_proposal) + B * eta_proposal) - ((v + w - 1) * log(eta) - z * eta / ta + eta * log(be * eta) - lgamma(eta) + eta * log(la) - la * be * eta - (A - 1) * log(eta) + B * eta)
    eta <- ifelse(test = (log(uniform) <= logratio), yes = eta_proposal, no = eta)
    Eta[iota, ] <- eta
  }
  
  # summary
  Res <- list(La=1/La[-(1:burn), , drop=FALSE], Be=Be[-(1:burn)], Ta=Ta[-(1:burn)], 
              V=V[-(1:burn), , drop=FALSE], W=W[-(1:burn), , drop=FALSE], Z=Z[-(1:burn), , drop=FALSE], 
              Eta=Eta[-(1:burn), , drop=FALSE])
  return(Res)
}




###    Inverse Rescaled Beta (IRB) prior    ###
## model settings 
# y \sim Ga(A_y, B_y/lambda)
# lambda \sim Ga(1+tau*u, beta*tau*u)
# u \sim SB(a, b)
# tau \sim Ga(a_ta, b_ta)
# beta \sim Ga(a_be, b_be)
## IMPUT
# y: sequence of observed values
# A_y, B_y: fixed values in the sampling model 
# (a_ta, b_ta): hyperparameter in prior of tau
# (a_be, b_be): hyperparameter in prior of beta
# a, b: shape parameters of SB prior 
# ep, M: tuning parameters in the Miller's algorithm
# mc: MCMC length 
# burn: burn-in period
## OUTPUT
# posterior samples of parameters 

SB_prior <- function(y, A_y, B_y, a_ta, b_ta, a_be, b_be, a, b, de=10^(-8), ep=10^(-8), M=10, mc, burn){
  # preparation
  m <- length(y)
  # initial values
  la <- A_y / (B_y * y)
  t. <- rep(1, m)
  eta <- rep(1, m)
  ta <- 1
  be <- 1
  
  # objects to store posterior samples
  La <- matrix(NA, mc, m)
  T. <- matrix(NA, mc, m)
  Eta <- matrix(NA, mc, m)
  Ta <- rep(NA, mc)
  Be <- rep(NA, mc)
  
  # MCMC
  for(iota in 1:mc){
    # la
    la <- rgamma(m, shape = A_y + 1 + eta, rate = B_y * y + be * eta)
    La[iota, ] <- la
    
    # t
    t._proposal <- sapply(2 * (1 + eta / ta) + 2 * de, function(twicerate) rgig(1, lambda = a + b, chi = 2 * de, psi = twicerate))
    uniform <- runif(m, min = 0, max = 1)
    logratio <- dgamma(t._proposal, shape = a + b, rate = 1 + eta / ta, log = TRUE) - dgamma(t., shape = a + b, rate = 1 + eta / ta, log = TRUE) - ((a + b - 1) * log(t._proposal) - ((2 * (1 + eta / ta) + 2 * de) * t._proposal + (2 * de) / t._proposal) / 2 - (a + b - 1) * log(t.) + ((2 * (1 + eta / ta) + 2 * de) * t. + (2 * de) / t.) / 2)
    t. <- ifelse(test = (log(uniform) <= logratio), yes = t._proposal, no = t.)
    T.[iota, ] <- t.
    
    # be
    be <- rgamma(1, shape = a_be + sum(eta + 1), rate = b_be + sum(eta * la))
    Be[iota] <- be
    
    # ta
    ta <- rgig(1, lambda = a_ta - m * a, chi = 2 * sum(t. * eta), psi = 2 * b_ta)
    Ta[iota] <- ta
    
    # eta
    ABJ <- apply(rbind(la, a, t. / ta), 2, function(s) miller(x = s[1], mu = 1 / be, a0 = s[2], b0 = s[3], ep = ep, M = M))
    A <- ABJ[1, ]
    B <- ABJ[2, ]
    eta_proposal <- rgamma(m, shape = A, rate = B)
    uniform <- runif(m, min = 0, max = 1)
    logratio <- ((a - 1) * log(eta_proposal) - t. * eta_proposal / ta + eta_proposal * log(be * eta_proposal) - lgamma(eta_proposal) + eta_proposal * log(la) - la * be * eta_proposal - (A - 1) * log(eta_proposal) + B * eta_proposal) - ((a - 1) * log(eta) - t. * eta / ta + eta * log(be * eta) - lgamma(eta) + eta * log(la) - la * be * eta - (A - 1) * log(eta) + B * eta)
    eta <- ifelse(test = (log(uniform) <= logratio), yes = eta_proposal, no = eta)
    Eta[iota, ] <- eta
  }
  
  # summary
  Res <- list(La=1/La[-(1:burn), , drop=FALSE], Be=Be[-(1:burn)], Ta=Ta[-(1:burn)], 
              T.=T.[-(1:burn), , drop=FALSE], Eta=Eta[-(1:burn), , drop=FALSE])
  return(Res)
}





###    Global (gamma) prior     ###
## model settings 
# y \sim Ga(A_y, B_y/lambda)
# lambda \sim Ga(1+tau, beta*tau)
# tau \sim Ga(a_ta, b_ta)
# beta \sim Ga(a_be, b_be)
## IMPUT
# y: sequence of observed values
# A_y, B_y: fixed values in the sampling model 
# (a_ta, b_ta): hyperparameter in prior of tau
# (a_be, b_be): hyperparameter in prior of beta
## OUTPUT
# posterior samples of parameters 

global_prior <- function(y, A_y, B_y, a_ta, b_ta, a_be, b_be, ep=10^(-8), M=10, mc, burn){
  # preparation
  m <- length(y)
  # initial values
  la <- A_y / (B_y * y)
  ta <- 1
  be <- 1
  
  # objercts to store posterior samples
  La <- matrix(NA, mc, m)
  Ta <- rep(NA, mc)
  Be <- rep(NA, mc)
  
  # MCMC
  for(iota in 1:mc){
    # la
    la <- rgamma(m, shape = A_y + 1 + ta, rate = B_y * y + be * ta)
    La[iota, ] <- la
    
    # be
    be <- rgamma(1, shape = a_be + m * (ta + 1), rate = b_be + ta * sum(la))
    Be[iota] <- be
    
    # ta
    ABJ <- miller(x = la, mu = 1 / be, a0 = a_ta, b0 = b_ta, ep = ep, M = M)
    A <- ABJ[1]
    B <- ABJ[2]
    ta_proposal <- rgamma(1, shape = A, rate = B)
    uniform <- runif(1, min = 0, max = 1)
    logratio <- ((a_ta - 1) * log(ta_proposal) - b_ta * ta_proposal + m * ta_proposal * log(be * ta_proposal) - m * lgamma(ta_proposal) + ta_proposal * sum(log(la)) - be * ta_proposal * sum(la) - (A - 1) * log(ta_proposal) + B * ta_proposal) - ((a_ta - 1) * log(ta) - b_ta * ta + m * ta * log(be * ta) - m * lgamma(ta) + ta * sum(log(la)) - be * ta * sum(la) - (A - 1) * log(ta) + B * ta)
    if(log(uniform) <= logratio){
      ta <- ta_proposal
    }else{
      ta <- ta
    }
    Ta[iota] <- ta
  }
  
  # summary
  Res <- list(La=1/La[-(1:burn), , drop=FALSE], Be=Be[-(1:burn)], Ta=Ta[-(1:burn)])
  return(Res)
}



##     DasGupta's method
## (only provides point estimation)
DasGupta <- function(y, A_y, B_y, const=NA){
  m <- length(y)
  y_tilde <- B_y * y
  if(is.na(const)){
    const <- ((m - 1) / m) * (gamma(mean(A_y) + 1 / m) / gamma(mean(A_y) + 2 / m))^m
  }
  est_DasGupta <- (y_tilde / (A_y + 1)) * (1 + (const / y_tilde) * (prod(y_tilde))^(1 / m))
  return( est_DasGupta )
}

