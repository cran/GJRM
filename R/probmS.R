probmS <- function(eta, margin){ 

#library(Deriv); library(numDeriv)
#expr <- expression(   -((1 - exp(eta)) * exp(-exp(eta)) * exp(eta))   )
#Simplify(D(D(expr, "eta"), "eta")) 
#Simplify(D(expr, "eta")) 
#func0 <- function(eta){ -((1 - exp(eta)) * exp(-exp(eta)) * exp(eta))  }
#grad(func0 , eta)

epsilon <- 0.0000001 
  

if( margin == "probit" ){
 
  pr  <- pnorm(-eta)
  dS  <- -dnorm(-eta)
  d2S <- eta * dnorm(-eta)
  d3S <- (1 - eta^2) * dnorm(-eta)
  
}



if( margin == "logit" ){ # to change
 
  pr  <- 1/(1+exp(eta))    # plogis(-eta)
  dS  <- -(exp(eta)/(1 + exp(eta))^2)
  d2S <- -((1 - 2 * (exp(eta)/(1 + exp(eta)))) * exp(eta)/(1 + exp(eta))^2)
  d3S <- -((1 - (2 + 2 * (1 - 2 * (exp(eta)/(1 + exp(eta)))) + 2 * (1 - 
    exp(eta)/(1 + exp(eta)))) * exp(eta)/(1 + exp(eta))) * exp(eta)/(1 + 
    exp(eta))^2)
}



if( margin == "cloglog" ){ # to change
 
  pr  <- exp(-exp(eta))
  dS  <- -(exp(-exp(eta)) * exp(eta))
  d2S <- -((1 - exp(eta)) * exp(-exp(eta)) * exp(eta))
  d3S <- -((1 - (3 - exp(eta)) * exp(eta)) * exp(-exp(eta)) * exp(eta))
}


pr <- mm(pr) 
dS <- abs(dS)
dS <- -ifelse(dS < epsilon, epsilon, dS )

list(pr = pr, dS = dS, d2S = d2S, d3S = d3S)  
 
}    

