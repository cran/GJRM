distrHsAT <- function(y2, eta2, sigma2, nu, margin2,
                      min.dn, min.pr, max.pr, left.trunc = 0){




# sigma2 is sigma


sigma <- sigma2


if(margin2 == "probit"){

  pdf2          <- dnorm(y2, 0, 1)
    p2          <- pnorm(y2, 0, 1)
     
}


if( margin2 == "logit" ){
 
  p2  <- plogis(y2)
  pdf2 <- dlogis(y2)
  
  
}



if( margin2 == "cloglog" ){
 
  p2  <- 1-exp(-exp(y2))
  pdf2 <- exp(-exp(y2)) * exp(y2)
  
}


if( margin2 == "cauchit" ){
 
  p2  <- 1 / pi * atan(y2) + 0.5
  pdf2 <- 1 / (pi * (1 + y2^2))
  
  
}



if(margin2 == "N"){


  pdf2          <- dnorm(y2, mean = eta2, sd = sigma)
    p2          <- pnorm(y2, mean = eta2, sd = sigma)
     
}




if(margin2 %in% c("tN")){


  ltr <- rep(left.trunc, length(eta2))

  pdf2          <- dnorm(y2, mean = eta2, sd = sigma)/(1 - pnorm(ltr, mean = eta2, sd = sigma))   
    p2          <- (pnorm(y2, mean = eta2, sd = sigma) - pnorm(ltr, mean = eta2, sd = sigma))/(1 - pnorm(ltr, mean = eta2, sd = sigma))
   
              

}





if(margin2 == "LN"){

  pdf2          <- dlnorm(y2, meanlog = eta2, sdlog = sigma)
    p2          <- plnorm(y2, meanlog = eta2, sdlog = sigma)

}


if(margin2 == "WEI"){

  pdf2          <- sigma/exp(eta2)*(y2/exp(eta2))^(sigma-1) * exp(-(y2/exp(eta2))^sigma)
                   
    p2          <-  1-exp(-(y2/exp(eta2))^sigma) 

}

if(margin2 == "GO"){

  pdf2          <- exp(eta2)*exp(sigma*y2)*exp( ( 1 - exp(sigma*y2) )*exp(eta2)/sigma )
                   
    p2          <- 1 - exp( ( 1 - exp(sigma*y2) )*exp(eta2)/sigma )  

}


if(margin2 == "FISK"){

  pdf2          <- sigma*y2^(sigma-1) / (exp(eta2)^sigma*(1+(y2/exp(eta2))^sigma)^2)
                   
    p2          <-  1/(1+(y2/exp(eta2))^-sigma)

}




if(margin2 == "IG"){


#sigma2 <- ifelse(sigma2 < 0.001, 0.001, sigma2)


  pdf2          <- exp(-0.5 * log(2 * pi) - log(sigma) - (3/2) * log(y2) - 
                   ((y2 - exp(eta2))^2)/(2 * sigma^2 * (exp(eta2)^2) * y2))
                   
    p2          <-  pnorm(((y2/exp(eta2)) - 1)/(sigma * sqrt(y2))) + 
                    exp(2/(exp(eta2)*sigma^2))* pnorm(-((y2/exp(eta2)) + 1)/(sigma * sqrt(y2)))
                
}



if(margin2 == "LO"){

  pdf2          <- dlogis(y2, eta2, sigma) # exp(-(y2-eta2)/sqrt(sigma2))/(sqrt(sigma2)*(1+exp(-(y2-eta2)/sqrt(sigma2)))^2)
    p2          <- plogis(y2, eta2, sigma) #1/(1+exp(-(y2-eta2)/sqrt(sigma2)))
                
}


if(margin2 == "rGU"){

  pdf2          <- 1/sigma*exp(-((y2-eta2)/sigma+exp(-((y2-eta2)/sigma))))
    p2          <- exp(-(exp(-(y2-eta2)/sigma)))
                
}



if(margin2 == "GU"){


  
  bit <- (exp((y2 - eta2)/sigma) * (1/sigma))

  bit <- ifelse(bit == "Inf", 1e+290, bit)
  
  pdf2          <- exp(-exp((y2 - eta2)/sigma)) * bit
  
    p2          <- 1 - exp(-exp((y2 - eta2)/sigma))
    

}


if(margin2 %in% c("GA")){

sigma    <- ifelse(sigma < 0.006, 0.006, sigma) # related to gamma function

 pdf2          <-  dgamma(y2, shape = 1/sigma^2, scale = exp(eta2) * sigma^2)
                   
    p2         <-  pgamma(y2, shape = 1/sigma^2, scale = exp(eta2) * sigma^2)
    
    }
    
    
    
    
#if(margin2 %in% c("GA2")){
#
#sigma2    <- ifelse(sigma2 < 0.006, 0.006, sigma2) # related to gamma function
#
#pdf2 <- exp(eta2)^(sqrt(sigma2))  * y2^(sqrt(sigma2) -1) * exp(-(y2*exp(eta2))) / gamma(sqrt(sigma2))          
#                   
#p2  <- pgamma(y2, shape =  sqrt(sigma2), rate= exp(eta2))
#    
#    }
    
    
    
#if(margin2 %in% c("GGA")){
#
#
#pdf2 <- (sqrt(sigma2) / gamma(nu))*(y2^(sqrt(sigma2)*nu -1) / exp(eta2)^(sqrt(sigma2)*nu))* exp(-(y2/exp(eta2))^sqrt(sigma2))           
#                   
#p2  <- pgamma(exp((log(y2) - log(exp(eta2))) * sqrt(sigma2)), shape =  nu)
#    
#    }    
    
    
    
    
    
    
    
    
if(margin2 %in% c("GAi")){

sigma <- ifelse(sigma < 0.006, 0.006, sigma) # related to gamma function
eta2   <- ifelse(eta2 < 1e-07, 1e-07, eta2)

 pdf2  <-  dgamma(y2, shape = 1/sigma^2, scale = eta2 * sigma^2)       
    p2 <-  pgamma(y2, shape = 1/sigma^2, scale = eta2 * sigma^2)
    
    }    
    


if(margin2 == "DAGUM"){

pdf2 <- sigma*nu/y2*( ((y2/exp(eta2))^(sigma*nu)) /  ( (y2/exp(eta2))^sigma + 1 )^(nu+1) )            

    p2  <- ( 1 + (y2/exp(eta2))^-sigma )^-nu 


}




if(margin2 == "TW"){




    
   
aTW <- 1.001  
bTW <- 1.999

mu         <- exp(eta2)
#nu.stTW    <- log(nu)
#sigma.stTW <- log( (sigma - aTW) / (bTW - sigma) ) 

nu.stTW    <- log( (nu - aTW) / (bTW - nu) )
sigma.stTW <- log(sigma) 









#TWob <- ldTweedie(y2, mu = mu, p = NA, phi = NA, rho = nu.stTW, theta = sigma.stTW, all.derivs = TRUE) 


TWob <- ldTweedie(y2, mu = mu, p = NA, phi = NA, rho = sigma.stTW, theta = nu.stTW, all.derivs = TRUE) 



pdf2 <- exp(TWob[, 1])   


p2 <- NA


if(length(sigma) == 1) sigma <- rep(sigma, length(y2))
if(length(nu) == 1)    nu    <- rep(nu, length(y2))
if(length(mu) == 1)    mu    <- rep(mu, length(y2))


 
for(i in 1:length(y2)) p2[i] <- pTweed(y = y2[i], mu = mu[i], phi = sigma[i], p = nu[i])$d0 


    

}





if(margin2 == "SM"){

pdf2 <- sigma*nu*y2^(sigma-1)*(exp(eta2)^sigma*(1+(y2/exp(eta2))^sigma)^(nu+1) )^-1            

    p2  <- 1 - (1+(y2/exp(eta2))^sigma)^-nu 


}



if(margin2 == "BE"){

pdf2 <- dbeta(y2, shape1 = plogis(eta2) * (1 - sigma^2)/(sigma^2), shape2 = (1-plogis(eta2))*(1 - sigma^2)/(sigma^2))          

    p2  <- pbeta(y2, shape1 = plogis(eta2) * (1 - sigma^2)/(sigma^2), shape2 = (1-plogis(eta2))*(1 - sigma^2)/(sigma^2))


}






if(margin2 == "GP"){

indx <- (1 + eta2*y2/sigma) > 0 
pdf2 <- suppressWarnings(  1/sigma*(1 +  eta2*y2/sigma)^(-1/eta2-1)   )
 p2  <- suppressWarnings(  1 - (1 +  eta2*y2/sigma)^(-1/eta2)   )


pdf2 <- ifelse( indx == TRUE, pdf2, 0)
p2   <- ifelse( indx == TRUE, p2, 0)


}



if(margin2 == "GPII"){


mu2  <- exp(eta2) - 0.5 # exi

indx <- (1 + mu2*y2/sigma) > 0 
pdf2 <- suppressWarnings(  1/sigma*(1 +  mu2*y2/sigma)^(-1/mu2-1)   )
 p2  <- suppressWarnings(  1 - (1 +  mu2*y2/sigma)^(-1/mu2)   )


pdf2 <- ifelse( indx == TRUE, pdf2, 0)
p2   <- ifelse( indx == TRUE, p2, 0)


}



if(margin2 == "GPo"){

mu2  <- exp(eta2) - 0.5 # exi

indx <- (1 + mu2*y2/(sigma/(1+mu2)  )  ) > 0 
pdf2 <- suppressWarnings(  1/(sigma/(1+mu2))*(1 +  mu2*y2/(sigma/(1+mu2)))^(-1/mu2-1)   )
 p2  <- suppressWarnings(  1 - (1 +  mu2*y2/(sigma/(1+mu2)))^(-1/mu2)   )

pdf2 <- ifelse( indx == TRUE, pdf2, 0)
p2   <- ifelse( indx == TRUE, p2, 0)

}




pdf2 <- ifelse(pdf2 < min.dn, min.dn, pdf2)
p2   <- mm(p2, min.pr = min.pr, max.pr = max.pr) 



list(pdf2 = ifef(pdf2), p2 = ifef(p2))     


}




     























