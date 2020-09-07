probmS <- function(eta, margin, min.dn, min.pr, max.pr){ 

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


pr <- mm(pr, min.pr = min.pr, max.pr = max.pr) 
dS <- abs(dS)
dS <- -ifelse(dS < min.dn, min.dn, dS )



list(pr = pr, dS = dS, d2S = d2S, d3S = d3S)  
 
}    

