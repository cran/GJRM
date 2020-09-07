probm <- function(eta, margin, only.pr = TRUE, bc = FALSE, tau = NULL, CLM = FALSE,
                  min.dn, min.pr, max.pr){ # bc stands for binary continuous case
   
  derp1.dereta1 <- der2p1.dereta1eta1 <- d.n <- der2p.dereta <- tauetaIND <- NULL
 
 
if( margin == "GEVlink" ){
 
  taueta    <- tau*eta 
  tauetaIND <- ifelse(taueta < -1, TRUE, FALSE) 
  taueta[tauetaIND] <- -1
  
  pr  <- exp(-(1+taueta)^(-1/tau))
  
  if(only.pr == FALSE){
  
  d.n <- exp(-(1 + taueta)^-(1/tau))/(1 + taueta)^(1 + 1/tau)
  der2p.dereta <- -(exp(-(1 + taueta)^-(1/tau)) * (tau * (1 + 1/tau)/(1 + taueta)^(1/tau + 2) - 1/(1 + taueta)^(2 * (1 + 1/tau))))

                      }
  
  if(bc == TRUE){ 
  
  
  derp1.dereta1      <- -(exp(-(1 + taueta)^-(1/tau))/(1 + taueta)^(1 + 1/tau))   
  der2p1.dereta1eta1 <- exp(-(1 + taueta)^-(1/tau)) * (tau * (1 + 1/tau)/(1 + taueta)^(1/tau + 2) - 1/(1 + taueta)^(2 * (1 + 1/tau))) 
  
                }  
      
} 
 
 
if( margin == "probit" ){
 
  pr  <- pnorm(eta)
  
  if(only.pr == FALSE){
  
  d.n <- dnorm(eta)
  der2p.dereta <- -(eta * dnorm(eta))          # second deriv
  
                      }
  
  if(bc == TRUE){ 
  
  # RECALL that this is deriv of Y = 0 wrt eta1 and not the usual case of Y = 1!!
  derp1.dereta1      <- -dnorm(-eta)         # First derivative of 1-prob(eta1) respect to eta1
  
  der2p1.dereta1eta1 <- eta * dnorm(-eta)  ## This is the second derivative of 1 - p1 respect to eta1 
  
                }  
      
}


if( margin == "logit" ){
 
 

if(CLM == FALSE){
 
  pr  <- plogis(eta)
  
  if(only.pr == FALSE){
  
  d.n <- dlogis(eta)
  der2p.dereta <- -((1 - 2 * (exp(-eta)/(1 + exp(-eta)))) * exp(-eta)/(1 + exp(-eta))^2)
  
  }
  
  if(bc == TRUE){
  derp1.dereta1    <- -((1 - exp(-eta)/(1 + exp(-eta))) * exp(-eta)/(1 + exp(-eta))) # First derivative of 1-prob(eta1) respect to eta1.
      
  der2p1.dereta1eta1 <- (1 - (3 - 2 * (exp(-eta)/(1 + exp(-eta)))) * exp(-eta)/(1 + exp(-eta))) * 
                           exp(-eta)/(1 + exp(-eta))

  }


}

if(CLM == TRUE){



  ex_eta <- exp(-eta) ; ex_eta[ex_eta == Inf] <- 1e+25 # FD: CLM addition


  pr  <- plogis(eta)
  
  if(only.pr == FALSE){
  
  d.n <- dlogis(eta)
  der2p.dereta <- -((1 - 2 * (ex_eta/(1 + ex_eta))) * ex_eta/(1 + ex_eta)^2)
  
  }
  
  if(bc == TRUE){
  derp1.dereta1    <- -((1 - ex_eta/(1 + ex_eta)) * ex_eta/(1 + ex_eta)) # First derivative of 1-prob(eta1) respect to eta1.
      
  der2p1.dereta1eta1 <- (1 - (3 - 2 * (ex_eta/(1 + ex_eta))) * ex_eta/(1 + ex_eta)) * 
                           ex_eta/(1 + ex_eta)

  }
  





}


  
}






if( margin == "cloglog" ){
 
  pr  <- 1-exp(-exp(eta))
  
  
  if(only.pr == FALSE){
  
  d.n <- exp(-exp(eta)) * exp(eta)
  der2p.dereta <- (1 - exp(eta)) * exp(-exp(eta)) * exp(eta) 
  
  
  }
  
  
  if(bc == TRUE){
  
  derp1.dereta1    <-  -(exp(-exp(eta)) * exp(eta))
                   
  der2p1.dereta1eta1 <- -((1 - exp(eta)) * exp(-exp(eta)) * exp(eta))
 
                }
 
}




if( margin == "cauchit" ){
 
  pr  <- 1 / pi * atan(eta) + 0.5
  
  
  if(only.pr == FALSE){
  
  d.n <- 1 / (pi * (1 + eta^2))
  der2p.dereta <- -(2 * (eta/(pi * (1 + eta^2)^2)))
  
  }
  
  
  if(bc == TRUE){
  
  derp1.dereta1    <- -(1/(pi * (1 + eta^2)))
         
  der2p1.dereta1eta1 <- 2 * (eta/(pi * (1 + eta^2)^2))


  }
  
  
}








if( margin == "log" ){
 
  pr  <- exp(eta)
  
  
  if(only.pr == FALSE){
  
  d.n <- exp(eta)
  der2p.dereta <- exp(eta)
  
  
  }
  
  
  if(bc == TRUE){
  
  derp1.dereta1    <- -exp(eta)  
                   
  der2p1.dereta1eta1 <- -exp(eta)
 
                }
 
}


 if(CLM == FALSE){
 
 pr  <- mm(pr, min.pr = min.pr, max.pr = max.pr) 
 d.n <- ifelse(d.n < min.dn, min.dn, d.n)

  
 }
 
 
 
    
 list(pr = pr, d.n = d.n, der2p.dereta = der2p.dereta, 
      derp1.dereta1 = derp1.dereta1, der2p1.dereta1eta1 = der2p1.dereta1eta1, 
      tauetaIND = tauetaIND )  
 
}    