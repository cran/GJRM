distrHsAT1 <- function(y2.st, eta2, sigma2, nu, margin2){
 

 
# TW is missing here 



sigma <- sigma2
 
if(margin2 %in% c("N")){

y2 <- y2.st

mu2 <- eta2

der2p2.dery22  <- -(1/sqrt(2 * pi * sigma^2) * (exp(-0.5 * (y2 - mu2)^2/sigma^2) * (0.5 * (2 * (y2 - mu2))/sigma^2)))
dery.dery.st   <- 1
der2y.dery.st2 <- 0


}


if(margin2 %in% c("LN")){

y2 <- exp(y2.st)
dery.dery.st <- exp(y2.st)
der2y.dery.st2 <- exp(y2.st)

mu2 <- exp(eta2)

der2p2.dery22  <- -(((log(y2) - mu2)/(sigma^3 * y2^2 * sqrt(2 * pi)) + sigma * 
    sqrt(2 * pi)/(sigma * y2 * sqrt(2 * pi))^2) * exp(-(0.5 * 
    ((log(y2) - mu2)^2/sigma^2))))
    
    
    #1/(y2*sigma*sqrt(2*pi))*exp(-0.5*(log(y2)-mu2)^2/sigma^2)



}

###################

if(margin2 == "DAGUM"){

y2 <- exp(y2.st)
dery.dery.st <- exp(y2.st)
der2y.dery.st2 <- exp(y2.st)

mu2 <- exp(eta2)

der2p2.dery22  <- nu * sigma * (sigma * (nu * (y2/mu2)^(nu * sigma - 1)/((y2/mu2)^sigma + 
    1)^(1 + nu) - ((y2/mu2)^sigma + 1)^(nu - 2 * (1 + nu)) * 
    (1 + nu) * (y2/mu2)^(sigma * (1 + nu) - 1))/mu2 - (y2/mu2)^(nu * 
    sigma)/(y2 * ((y2/mu2)^sigma + 1)^(1 + nu)))/y2 

         
}  

        
        



####



if(margin2 == "SM"){


y2 <- exp(y2.st)
dery.dery.st <- exp(y2.st)
der2y.dery.st2 <- exp(y2.st)

mu2 <- exp(eta2)

der2p2.dery22  <- nu * sigma * (y2^(sigma - 2) * (sigma - 1)/(mu2^sigma * ((y2/mu2)^sigma + 
    1)^(1 + nu)) - mu2^(sigma - 1) * sigma * y2^(sigma - 1) * 
    ((y2/mu2)^sigma + 1)^nu * (1 + nu) * (y2/mu2)^(sigma - 1)/(mu2^sigma * 
    ((y2/mu2)^sigma + 1)^(1 + nu))^2) 



}



 


###


if(margin2 == "FISK"){


y2 <- exp(y2.st)
dery.dery.st <- exp(y2.st)
der2y.dery.st2 <- exp(y2.st)

mu2 <- exp(eta2)

der2p2.dery22  <- sigma * (y2^(sigma - 2) * (sigma - 1)/(mu2^sigma * ((y2/mu2)^sigma + 
    1)^2) - 2 * (mu2^(sigma - 1) * sigma * y2^(sigma - 1) * ((y2/mu2)^sigma + 
    1) * (y2/mu2)^(sigma - 1)/(mu2^sigma * ((y2/mu2)^sigma + 
    1)^2)^2)) 



}




if(margin2 %in% c("GP","GPII")){


y2 <- exp(y2.st)
dery.dery.st <- exp(y2.st)
der2y.dery.st2 <- exp(y2.st)


if(margin2 == "GP")   mu2 <- eta2
if(margin2 == "GPII") mu2 <- exp(eta2) - 0.5


der2p2.dery22  <- suppressWarnings(  -(mu2 * (1 + 1/mu2)/(sigma^2 * (1 + mu2 * y2/sigma)^(1/mu2 + 2)))     )

# done but exclusions not implemented



}





if(margin2 %in% c("GPo")){


y2 <- exp(y2.st)
dery.dery.st <- exp(y2.st)
der2y.dery.st2 <- exp(y2.st)

mu2 <- exp(eta2) - 0.5

der2p2.dery22  <- suppressWarnings(  -(mu2 * (1 + 1/mu2)/((sigma/(1+mu2)  )^2 * (1 + mu2 * y2/(sigma/(1+mu2)  ))^(1/mu2 + 2)))     )

# done but exclusions not implemented



}






#library(Deriv); library(numDeriv)
#expr <- expression( 1 - (1 + mu2*y2/ (sigma/(1+mu2)) )^(-1/mu2) )
#Simplify( D(D(expr, "y2"),"y2") )
#func0 <- function(y2){ 1 - (1 + mu2*y2/ (sigma/(1+mu2)) )^(-1/mu2)  }
#hessian(func0 , y2)





####



if(margin2 == "WEI"){

y2 <- exp(y2.st)
dery.dery.st <- exp(y2.st)
der2y.dery.st2 <- exp(y2.st)

mu2 <- exp(eta2)

der2p2.dery22  <- sigma * ((sigma - 1) * (y2/mu2)^(sigma - 2) - sigma * (y2/mu2)^(2 * (sigma - 1))) * exp(-(y2/mu2)^sigma)/mu2^2


}



#if(margin2 == "GO"){
#
#y2 <- exp(y2.st)
#dery.dery.st <- exp(y2.st)
#der2y.dery.st2 <- exp(y2.st)
#
#mu2 <- exp(eta2)
#
#der2p2.dery22 <- mu2 * exp(mu2 * (1 - exp(y2 * sqrt(sigma2)))/sqrt(sigma2)) * 
#    exp(y2 * sqrt(sigma2)) * (sigma2/sqrt(sigma2) - mu2 * exp(y2 * 
#    sqrt(sigma2))) 
#
#}





if(margin2 == "BE"){


y2  <- plogis(y2.st)
mu2 <- plogis(eta2)

a <- mu2 * (1 - sigma^2)/(sigma^2)
b <- a * (1 - mu2)/mu2

dery.dery.st   <- (1 - exp(y2.st)/(1 + exp(y2.st))) * exp(y2.st)/(1 + exp(y2.st))
der2y.dery.st2 <- (1 - (3 - 2 * (exp(y2.st)/(1 + exp(y2.st)))) * exp(y2.st)/(1 + exp(y2.st))) * exp(y2.st)/(1 + exp(y2.st))

derp2.dery2 <- dbeta(y2, shape1 = mu2 * (1 - sigma^2)/(sigma^2), shape2 = (1-mu2)*(1 - sigma^2)/(sigma^2))

der2p2.dery22  <- derp2.dery2  * ( (a + b - 2) * y2 - (a - 1) ) / ( ( y2 - 1 ) * y2 ) # from wikipedia




  }
  





####

if(margin2 == "iG"){

y2 <- exp(y2.st)
dery.dery.st <- exp(y2.st)
der2y.dery.st2 <- exp(y2.st)

mu2 <- exp(eta2)


der2p2.dery22  <- -(exp(-0.5 * log(2 * pi) - log(sigma) - (3/2) * log(y2) - ((y2 - 
    mu2)^2)/(2 * sigma^2 * (mu2^2) * y2)) * ((3/2) * (1/y2) + 
    (2 * (y2 - mu2)/(2 * sigma^2 * (mu2^2) * y2) - ((y2 - mu2)^2) * 
        (2 * sigma^2 * (mu2^2))/(2 * sigma^2 * (mu2^2) * y2)^2))) 




}


####

if(margin2 == "LO"){

y2 <- y2.st

mu2 <- eta2


der2p2.dery22  <- (2 * (exp(-((y2 - mu2)/sigma))/(1 + exp(-((y2 - mu2)/sigma)))) - 
    1) * exp(-((y2 - mu2)/sigma))/(sigma^2 * (1 + exp(-((y2 - 
    mu2)/sigma)))^2) 


dery.dery.st   <- 1
der2y.dery.st2 <- 0



}


##########################


if(margin2 == "rGU"){

y2 <- y2.st

mu2 <- eta2

der2p2.dery22  <- -((1 - exp(-((y2 - eta2)/sigma))) * exp(-((y2 - eta2)/sigma + 
    exp(-((y2 - eta2)/sigma))))/sigma^2) 

dery.dery.st   <- 1
der2y.dery.st2 <- 0
  

 
}


######################################


if(margin2 == "GU"){

y2 <- y2.st

mu2 <- eta2

  der2p2.dery22 <- (1 - exp((y2 - eta2)/sigma)) * exp(-exp((y2 - eta2)/sigma)) * 
    exp((y2 - eta2)/sigma)/sigma^2 
         
dery.dery.st   <- 1
der2y.dery.st2 <- 0





}




if(margin2 == "GA"){

sigma2    <- ifelse(sigma2 < 0.006, 0.006, sigma2) 

y2 <- exp(y2.st)
dery.dery.st <- exp(y2.st)
der2y.dery.st2 <- exp(y2.st)

mu2 <- exp(eta2)

  der2p2.dery22 <- ((1/y2 - 1/mu2)/sigma^2 - 1/y2) * exp((log(y2) - (2 * log(sigma) + 
    log(mu2) + y2/mu2))/sigma^2 - (lgamma(1/sigma^2) + log(y2)))   



#expr <- expression( exp( (1/sigma^2) * log(y2/(mu2 * sigma^2)) - y2/(mu2 * sigma^2) - log(y2) - lgamma(1/sigma^2) )  )
#D(expr, "y2")
#Simplify( D(expr, "y2") )
#
#
#func0 <- function(y2){ exp( (1/sigma^2) * log(y2/(mu2 * sigma^2)) - y2/(mu2 * sigma^2) - log(y2) - lgamma(1/sigma^2) )}
#grad(func0 , y2)
#dgamma(y2, shape = 1/sigma^2, scale = mu2 * sigma^2)  

    }
    
    
    


list(der2p2.dery22 = ifef(der2p2.dery22), dery.dery.st = ifef(dery.dery.st), der2y.dery.st2 = ifef(der2y.dery.st2))     


}
 




    