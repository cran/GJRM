distrHsAT1 <- function(y2.st, eta2, sigma2, nu, margin2){

 
if(margin2 %in% c("N")){

y2 <- y2.st

mu2 <- eta2

der2p2.dery22  <- -(1/sqrt(2 * pi * sigma2) * (exp(-0.5 * (y2 - mu2)^2/sigma2) * (0.5 * (2 * (y2 - mu2))/sigma2)))
dery.dery.st   <- 1
der2y.dery.st2 <- 0


}


if(margin2 %in% c("LN")){

y2 <- exp(y2.st)
dery.dery.st <- exp(y2.st)
der2y.dery.st2 <- exp(y2.st)

mu2 <- exp(eta2)

der2p2.dery22  <- -(1/(y2 * sqrt(sigma2) * sqrt(2 * pi)) * (exp(-0.5 * (log(y2) - mu2)^2/sigma2) * (0.5 * (2 * (1/y2 * (log(y2) - mu2)))/sigma2)) + 
    sqrt(sigma2) * sqrt(2 * pi)/(y2 * sqrt(sigma2) * sqrt(2 * pi))^2 * exp(-0.5 * (log(y2) - mu2)^2/sigma2))#D(expression(1/(y2*sqrt(sigma2)*sqrt(2*pi))*exp(-0.5*(log(y2)-mu2)^2/sigma2)),"y2")


}

###################

if(margin2 == "DAGUM"){

y2 <- exp(y2.st)
dery.dery.st <- exp(y2.st)
der2y.dery.st2 <- exp(y2.st)

mu2 <- exp(eta2)

der2p2.dery22  <- sqrt(sigma2) * nu/y2 * ((y2/mu2)^((sqrt(sigma2) * nu) - 1) * 
    ((sqrt(sigma2) * nu) * (1/mu2))/((y2/mu2)^sqrt(sigma2) + 
    1)^(nu + 1) - ((y2/mu2)^(sqrt(sigma2) * nu)) * (((y2/mu2)^sqrt(sigma2) + 
    1)^((nu + 1) - 1) * ((nu + 1) * ((y2/mu2)^(sqrt(sigma2) - 
    1) * (sqrt(sigma2) * (1/mu2)))))/(((y2/mu2)^sqrt(sigma2) + 
    1)^(nu + 1))^2) - sqrt(sigma2) * nu/y2^2 * (((y2/mu2)^(sqrt(sigma2) * 
    nu))/((y2/mu2)^sqrt(sigma2) + 1)^(nu + 1))

         
}  

        
        



####



if(margin2 == "SM"){


y2 <- exp(y2.st)
dery.dery.st <- exp(y2.st)
der2y.dery.st2 <- exp(y2.st)

mu2 <- exp(eta2)

der2p2.dery22  <- sqrt(sigma2) * nu * (y2^((sqrt(sigma2) - 1) - 1) * (sqrt(sigma2) - 
    1)) * (mu2^sqrt(sigma2) * (1 + (y2/mu2)^sqrt(sigma2))^(nu + 
    1))^-1 - sqrt(sigma2) * nu * y2^(sqrt(sigma2) - 1) * ((mu2^sqrt(sigma2) * 
    (1 + (y2/mu2)^sqrt(sigma2))^(nu + 1))^-(1 + 1) * (mu2^sqrt(sigma2) * 
    ((1 + (y2/mu2)^sqrt(sigma2))^((nu + 1) - 1) * ((nu + 1) * 
        ((y2/mu2)^(sqrt(sigma2) - 1) * (sqrt(sigma2) * (1/mu2)))))))


}



 
 



###


if(margin2 == "FISK"){


y2 <- exp(y2.st)
dery.dery.st <- exp(y2.st)
der2y.dery.st2 <- exp(y2.st)

mu2 <- exp(eta2)

der2p2.dery22  <- sqrt(sigma2) * (y2^((sqrt(sigma2) - 1) - 1) * (sqrt(sigma2) - 
    1))/(mu2^sqrt(sigma2) * (1 + (y2/mu2)^sqrt(sigma2))^2) - 
    sqrt(sigma2) * y2^(sqrt(sigma2) - 1) * (mu2^sqrt(sigma2) * 
        (2 * ((y2/mu2)^(sqrt(sigma2) - 1) * (sqrt(sigma2) * (1/mu2)) * 
            (1 + (y2/mu2)^sqrt(sigma2)))))/(mu2^sqrt(sigma2) * 
        (1 + (y2/mu2)^sqrt(sigma2))^2)^2


}




if(margin2 == "GP"){


y2 <- exp(y2.st)
dery.dery.st <- exp(y2.st)
der2y.dery.st2 <- exp(y2.st)

mu2 <- eta2

der2p2.dery22  <- suppressWarnings(   -(mu2 * (1 + 1/mu2)/((1 + mu2 * y2/sqrt(sigma2))^(1/mu2 + 2) * sigma2^1))    )

# done but exclusions not implemented

}







####



if(margin2 == "WEI"){

y2 <- exp(y2.st)
dery.dery.st <- exp(y2.st)
der2y.dery.st2 <- exp(y2.st)

mu2 <- exp(eta2)

der2p2.dery22  <- sqrt(sigma2)/mu2 * ((y2/mu2)^((sqrt(sigma2) - 1) - 1) * ((sqrt(sigma2) - 
    1) * (1/mu2))) * exp(-(y2/mu2)^sqrt(sigma2)) - sqrt(sigma2)/mu2 * 
    (y2/mu2)^(sqrt(sigma2) - 1) * (exp(-(y2/mu2)^sqrt(sigma2)) * 
    ((y2/mu2)^(sqrt(sigma2) - 1) * (sqrt(sigma2) * (1/mu2))))

}



if(margin2 == "GO"){

y2 <- exp(y2.st)
dery.dery.st <- exp(y2.st)
der2y.dery.st2 <- exp(y2.st)

mu2 <- exp(eta2)

der2p2.dery22 <- mu2 * exp(mu2 * (1 - exp(y2 * sqrt(sigma2)))/sqrt(sigma2)) * 
    exp(y2 * sqrt(sigma2)) * (sigma2/sqrt(sigma2) - mu2 * exp(y2 * 
    sqrt(sigma2))) 

}



if(margin2 == "BE"){

#########################

y2  <- plogis(y2.st)
mu2 <- plogis(eta2)

a <- mu2 * (1 - sigma2)/(sigma2)
b <- a * (1 - mu2)/mu2

derp2.dery2 <- dbeta(y2, shape1 = mu2 * (1 - sigma2)/(sigma2), shape2 = (1-mu2)*(1 - sigma2)/(sigma2))

der2p2.dery22  <- derp2.dery2  * ( (a + b - 2) * y2 - (a - 1) ) / ( ( y2 - 1 ) * y2 ) # from wikipedia

dery.dery.st   <- (1 - exp(y2.st)/(1 + exp(y2.st))) * exp(y2.st)/(1 + exp(y2.st))
der2y.dery.st2 <- (1 - (3 - 2 * (exp(y2.st)/(1 + exp(y2.st)))) * exp(y2.st)/(1 + exp(y2.st))) * exp(y2.st)/(1 + exp(y2.st))


  }
  



####

if(margin2 == "iG"){

y2 <- exp(y2.st)
dery.dery.st <- exp(y2.st)
der2y.dery.st2 <- exp(y2.st)

mu2 <- exp(eta2)


der2p2.dery22  <- -(exp(-0.5 * log(2 * pi) - log(sqrt(sigma2)) - (3/2) * log(y2) - 
    ((y2 - mu2)^2)/(2 * sigma2 * (mu2^2) * y2)) * ((3/2) * (1/y2) + 
    (2 * (y2 - mu2)/(2 * sigma2 * (mu2^2) * y2) - ((y2 - mu2)^2) * 
        (2 * sigma2 * (mu2^2))/(2 * sigma2 * (mu2^2) * y2)^2)))

}


####

if(margin2 == "LO"){

y2 <- y2.st

mu2 <- eta2



der2p2.dery22  <- exp((y2 - mu2)/sqrt(sigma2)) * (1/sqrt(sigma2))/(sqrt(sigma2) * 
    (1 + exp((y2 - mu2)/sqrt(sigma2)))^2) - exp((y2 - mu2)/sqrt(sigma2)) * 
    (sqrt(sigma2) * (2 * (exp((y2 - mu2)/sqrt(sigma2)) * (1/sqrt(sigma2)) * 
        (1 + exp((y2 - mu2)/sqrt(sigma2))))))/(sqrt(sigma2) * 
    (1 + exp((y2 - mu2)/sqrt(sigma2)))^2)^2


dery.dery.st   <- 1
der2y.dery.st2 <- 0


}


##########################


if(margin2 == "rGU"){

y2 <- y2.st

mu2 <- eta2

der2p2.dery22  <- -(1/sqrt(sigma2) * (exp(-((y2 - mu2)/sqrt(sigma2) + exp(-((y2 - 
    mu2)/sqrt(sigma2))))) * (1/sqrt(sigma2) - exp(-((y2 - mu2)/sqrt(sigma2))) * 
    (1/sqrt(sigma2)))))

dery.dery.st   <- 1
der2y.dery.st2 <- 0
  
 
}


######################################


if(margin2 == "GU"){

y2 <- y2.st

mu2 <- eta2

  der2p2.dery22 <- exp(-exp((y2 - mu2)/sqrt(sigma2))) * (exp((y2 - mu2)/sqrt(sigma2)) * 
    (1/sqrt(sigma2)) * (1/sqrt(sigma2))) - exp(-exp((y2 - mu2)/sqrt(sigma2))) * 
    (exp((y2 - mu2)/sqrt(sigma2)) * (1/sqrt(sigma2))) * (exp((y2 - 
    mu2)/sqrt(sigma2)) * (1/sqrt(sigma2)))
         
dery.dery.st   <- 1
der2y.dery.st2 <- 0

}




if(margin2 == "GA"){

sigma2    <- ifelse(sigma2 < 0.006, 0.006, sigma2) 

y2 <- exp(y2.st)
dery.dery.st <- exp(y2.st)
der2y.dery.st2 <- exp(y2.st)

mu2 <- exp(eta2)

  der2p2.dery22 <- (y2^((1/sigma2 - 1) - 1) * (1/sigma2 - 1) * exp(-y2/(mu2 * 
    sigma2)) - y2^(1/sigma2 - 1) * (exp(-y2/(mu2 * sigma2)) * 
    (1/(mu2 * sigma2))))/(gamma(1/sigma2)*(mu2 * sigma2))^(1/sigma2)


    }
    
    
if(margin2 == "GA2"){

sigma2    <- ifelse(sigma2 < 0.006, 0.006, sigma2) 

y2 <- exp(y2.st)
dery.dery.st <- exp(y2.st)
der2y.dery.st2 <- exp(y2.st)

mu2 <- exp(eta2)

  der2p2.dery22 <- exp(-(mu2 * y2)) * (mu2^sqrt(sigma2) * y2^(sqrt(sigma2) - 2) * 
    (sqrt(sigma2) - 1) - mu2^(1 + sqrt(sigma2)) * y2^(sqrt(sigma2) - 
    1))/gamma(sqrt(sigma2))

    }    
    
    
    
if(margin2 == "GGA"){

y2 <- exp(y2.st)
dery.dery.st <- exp(y2.st)
der2y.dery.st2 <- exp(y2.st)

mu2 <- exp(eta2)

der2p2.dery22 <- exp(-(y2/mu2)^sqrt(sigma2)) * (y2^(nu * sqrt(sigma2) - 2) * (nu * 
    sqrt(sigma2) - 1) * sqrt(sigma2)/mu2^(nu * sqrt(sigma2)) - 
    sigma2 * y2^(nu * sqrt(sigma2) - 1) * (y2/mu2)^(sqrt(sigma2) - 
        1)/mu2^(1 + nu * sqrt(sigma2)))/gamma(nu)


    }       
    
 
 
    
    
    
if(margin2 %in% c("GAi")){

sigma2 <- ifelse(sigma2 < 0.006, 0.006, sigma2) # related to gamma function
eta2   <- ifelse(eta2 < 0.0000001, 0.0000001, eta2)

y2 <- exp(y2.st)
dery.dery.st <- exp(y2.st)
der2y.dery.st2 <- exp(y2.st)

mu2 <- eta2
      

 der2p2.dery22 <- (y2^((1/sigma2 - 1) - 1) * (1/sigma2 - 1) * exp(-y2/(mu2 * sigma2)) - 
    y2^(1/sigma2 - 1) * (exp(-y2/(mu2 * sigma2)) * (1/(mu2 * 
        sigma2))))/ (gamma(1/sigma2)*(mu2 * sigma2))^(1/sigma2)

    
    } 
    
    

ifef <- function(dv){

epsilon <- 0.0000001 
dv <- ifelse(is.na(dv), epsilon, dv ) 
dv <- ifelse(dv == Inf ,  8.218407e+20, dv )
dv <- ifelse(dv == -Inf ,  -8.218407e+20, dv )
dv

}

der2p2.dery22  <- ifef(der2p2.dery22 )
dery.dery.st   <- ifef(dery.dery.st  )
der2y.dery.st2 <- ifef(der2y.dery.st2)


list(der2p2.dery22 = der2p2.dery22, dery.dery.st = dery.dery.st, der2y.dery.st2 = der2y.dery.st2)     


}
 




    