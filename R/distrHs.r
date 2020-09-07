distrHs <- function(y2, eta2, sigma2, sigma2.st, nu, nu.st, margin2, naive = FALSE,
                    min.dn, min.pr, max.pr){

# sigma2 here means sigma, it was just a typo

sigma <- sigma2

p2 <- derp2.dersigma.st <- derp2.dereta2 <- der2p2.dereta2eta2 <- der2p2.dersigma2.st2 <- der2p2.dereta2dersigma2.st <- indx <- 1

der2pdf2.dereta2dernu.st    = 1
der2pdf2.sigma2.st2dernu.st = 1
derpdf2.dernu.st            = 1
der2pdf2.dernu.st2          = 1
derp2.nu.st                 = 1
der2p2.dernu.st2            = 1
der2p2.dereta2dernu.st      = 1
der2p2.dersigma2.stdernu.st = 1

cont2par <- c("WEI","iG","LO","rGU","GU","GA","GAi","BE","FISK", "N","LN","GP","GPII","GPo") 
cont3par <- c("DAGUM", "SM", "TW")

############################################################################
# remember that eta2 will have to disappear if we change default link on mu2
# this only applies to cases in which mu2 must be positive
# otherwise things are fine # library(Deriv); library(numDeriv)
############################################################################


if(margin2 %in% c("N","LN")){

                   mu2 <- eta2
                 sigma <- sigma2
                                      
dermu2.dereta2         <- 1
der2mu2.dereta2eta2    <- 0 

dersigma2.dersigma2.st  <- exp(sigma2.st)   
dersigma2.dersigma2.st2 <- exp(sigma2.st)   

  pdf2                <- dnorm(y2, mean = mu2, sd = sigma)  # 1/(sqrt(2*sigma^2*pi))*exp(-0.5*(y2-mu2)^2/sigma^2)
derpdf2.dermu2        <- exp(-(0.5 * ((y2 - mu2)^2/sigma^2))) * (y2 - mu2)/(sigma^2 * sqrt(2 * (pi * sigma^2)))   
der2pdf2.dermu2       <- ((y2 - mu2)^2/sigma^2 - 1) * exp(-(0.5 * ((y2 - mu2)^2/sigma^2)))/(sigma^2 * sqrt(2 * (pi * sigma^2))) 

derpdf2.sigma2        <- -((1 - (y2 - mu2)^2/sigma^2) * exp(-(0.5 * ((y2 - mu2)^2/sigma^2)))/(sigma * sqrt(2 * (pi * sigma^2)))) 
der2pdf2.dersigma22   <- -(((3 - (y2 - mu2)^2/sigma^2) * (y2 - mu2)^2/(sigma^4 * sqrt(2 * 
    (pi * sigma^2))) - (1 - (y2 - mu2)^2/sigma^2) * (2 * (pi * 
    sigma^2/sqrt(2 * (pi * sigma^2))) + sqrt(2 * (pi * sigma^2)))/(sigma * 
    sqrt(2 * (pi * sigma^2)))^2) * exp(-(0.5 * ((y2 - mu2)^2/sigma^2))))
der2pdf2.mu2dersigma2 <- ((y2 - mu2)^2/(sigma^5 * sqrt(2 * (pi * sigma^2))) - sigma * 
    (2 * (pi * sigma^2/sqrt(2 * (pi * sigma^2))) + 2 * sqrt(2 * 
        (pi * sigma^2)))/(sigma^2 * sqrt(2 * (pi * sigma^2)))^2) * 
    exp(-(0.5 * ((y2 - mu2)^2/sigma^2))) * (y2 - mu2)  
  

if(naive == FALSE){  
  
    p2                 <- pnorm(y2, mean = mu2, sd = sigma)
derp2.dermu2           <- -pdf2                
 der2p2.dermu22        <- -(exp(-(0.5 * ((y2 - mu2)^2/sigma^2))) * (y2 - mu2)/(sigma^2 * sqrt(2 * (pi * sigma^2)))) 
 
derp2.dersigma2        <- -dnorm((y2-mu2)/sigma)*(y2-mu2)/sigma^2  
der2p2.dersigma22      <- -(((y2 - mu2)^2/sigma^2 - 2) * dnorm((y2 - mu2)/sigma) * (y2 - mu2)/sigma^3) 
der2p2.derdermu2sigma2 <- (1 - (y2 - mu2)^2/sigma^2) * exp(-(0.5 * ((y2 - mu2)^2/sigma^2)))/(sigma * sqrt(2 * (pi * sigma^2))) 
                
                
                }




}











###################

if(margin2 == "DAGUM"){

mu2 <- exp(eta2)
dermu2.dereta2 <- exp(eta2)
der2mu2.dereta2eta2 <- exp(eta2) 

dersigma2.dersigma2.st  <- exp(sigma2.st)   
dersigma2.dersigma2.st2 <- exp(sigma2.st)   

dernu.dernu.st  <- exp(nu.st)
dernu.dernu.st2 <- exp(nu.st)


pdf2 <- sigma*nu/y2*( ((y2/mu2)^(sigma*nu))/( (y2/mu2)^sigma + 1 )^(nu+1) )            
  
derpdf2.dermu2 <- -(nu * sigma^2 * (nu * (y2/mu2)^(nu * sigma - 1)/((y2/mu2)^sigma + 
    1)^(1 + nu) - ((y2/mu2)^sigma + 1)^(nu - 2 * (1 + nu)) * 
    (1 + nu) * (y2/mu2)^(sigma * (1 + nu) - 1))/mu2^2) 
    
derpdf2.sigma2 <- nu * ((y2/mu2)^(nu * sigma)/((y2/mu2)^sigma + 1)^(1 + nu) + sigma * 
    (log(y2) - log(mu2)) * (nu * (y2/mu2)^(nu * sigma)/((y2/mu2)^sigma + 
    1)^(1 + nu) - ((y2/mu2)^sigma + 1)^(nu - 2 * (1 + nu)) * 
    (1 + nu) * (y2/mu2)^(sigma * (1 + nu))))/y2  
    
derpdf2.nu <- sigma * ((y2/mu2)^(nu * sigma)/((y2/mu2)^sigma + 1)^(1 + nu) + 
    nu * (sigma * (log(y2) - log(mu2)) * (y2/mu2)^(nu * sigma)/((y2/mu2)^sigma + 
        1)^(1 + nu) - ((y2/mu2)^sigma + 1)^(1 + nu - 2 * (1 + 
        nu)) * log1p((y2/mu2)^sigma) * (y2/mu2)^(nu * sigma)))/y2

der2pdf2.dermu2 <- nu * sigma^2 * (2 * (nu * (y2/mu2)^(nu * sigma - 1)/((y2/mu2)^sigma + 
    1)^(1 + nu) - ((y2/mu2)^sigma + 1)^(nu - 2 * (1 + nu)) * 
    (1 + nu) * (y2/mu2)^(sigma * (1 + nu) - 1)) + y2 * (nu * 
    ((nu * sigma - 1) * (y2/mu2)^(nu * sigma - 2)/((y2/mu2)^sigma + 
        1)^(1 + nu) - sigma * ((y2/mu2)^sigma + 1)^(nu - 2 * 
        (1 + nu)) * (1 + nu) * (y2/mu2)^(sigma * (1 + nu) - 2)) - 
    (((y2/mu2)^sigma + 1)^(nu - 2 * (1 + nu)) * (sigma * (1 + 
        nu) - 1) * (y2/mu2)^(sigma * (1 + nu) - 2) + sigma * 
        ((y2/mu2)^sigma + 1)^(nu - (1 + 2 * (1 + nu))) * (nu - 
        2 * (1 + nu)) * (y2/mu2)^(sigma * (2 + nu) - 2)) * (1 + 
        nu))/mu2)/mu2^3 


der2pdf2.dersigma22 <- nu * (2 * (nu * (y2/mu2)^(nu * sigma)/((y2/mu2)^sigma + 1)^(1 + 
    nu)) + sigma * (log(y2) - log(mu2)) * (nu * (nu * (y2/mu2)^(nu * 
    sigma)/((y2/mu2)^sigma + 1)^(1 + nu) - ((y2/mu2)^sigma + 
    1)^(nu - 2 * (1 + nu)) * (1 + nu) * (y2/mu2)^(sigma * (1 + 
    nu))) - (((y2/mu2)^sigma + 1)^(nu - (1 + 2 * (1 + nu))) * 
    (nu - 2 * (1 + nu)) * (y2/mu2)^(sigma * (2 + nu)) + ((y2/mu2)^sigma + 
    1)^(nu - 2 * (1 + nu)) * (1 + nu) * (y2/mu2)^(sigma * (1 + 
    nu))) * (1 + nu)) - 2 * (((y2/mu2)^sigma + 1)^(nu - 2 * (1 + 
    nu)) * (1 + nu) * (y2/mu2)^(sigma * (1 + nu)))) * (log(y2) - 
    log(mu2))/y2 


der2pdf2.dernu2 <- sigma * (2 * (sigma * (log(y2) - log(mu2)) * (y2/mu2)^(nu * sigma)/((y2/mu2)^sigma + 
    1)^(1 + nu)) + nu * (sigma * (log(y2) - log(mu2)) * (sigma * 
    (log(y2) - log(mu2)) * (y2/mu2)^(nu * sigma)/((y2/mu2)^sigma + 
    1)^(1 + nu) - ((y2/mu2)^sigma + 1)^(1 + nu - 2 * (1 + nu)) * 
    log1p((y2/mu2)^sigma) * (y2/mu2)^(nu * sigma)) - log1p((y2/mu2)^sigma) * 
    (sigma * ((y2/mu2)^sigma + 1)^(1 + nu - 2 * (1 + nu)) * (log(y2) - 
        log(mu2)) * (y2/mu2)^(nu * sigma) - ((y2/mu2)^sigma + 
        1)^(1 + nu - 2 * (1 + nu)) * log1p((y2/mu2)^sigma) * 
        (y2/mu2)^(nu * sigma))) - 2 * (((y2/mu2)^sigma + 1)^(1 + 
    nu - 2 * (1 + nu)) * log1p((y2/mu2)^sigma) * (y2/mu2)^(nu * 
    sigma)))/y2
    
    

der2pdf2.dersigma2dernu <- ((y2/mu2)^(nu * sigma)/((y2/mu2)^sigma + 1)^(1 + nu) + nu * (sigma * 
    (((y2/mu2)^(nu * sigma) + nu * sigma * (log(y2) - log(mu2)) * 
        (y2/mu2)^(nu * sigma))/((y2/mu2)^sigma + 1)^(1 + nu) + 
        (y2/mu2)^(nu * sigma)/((y2/mu2)^sigma + 1)^(1 + nu) - 
        ((((y2/mu2)^sigma + 1)^(nu - 2 * (1 + nu)) - ((y2/mu2)^sigma + 
            1)^(nu - 2 * (1 + nu)) * (1 + nu) * log1p((y2/mu2)^sigma)) * 
            (y2/mu2)^(sigma * (1 + nu)) + nu * ((y2/mu2)^sigma + 
            1)^(1 + nu - 2 * (1 + nu)) * log1p((y2/mu2)^sigma) * 
            (y2/mu2)^(nu * sigma) + sigma * ((y2/mu2)^sigma + 
            1)^(nu - 2 * (1 + nu)) * (1 + nu) * (log(y2) - log(mu2)) * 
            (y2/mu2)^(sigma * (1 + nu)))) * (log(y2) - log(mu2)) - 
    ((y2/mu2)^sigma + 1)^(1 + nu - 2 * (1 + nu)) * log1p((y2/mu2)^sigma) * 
        (y2/mu2)^(nu * sigma)) + sigma * (log(y2) - log(mu2)) * 
    (nu * (y2/mu2)^(nu * sigma)/((y2/mu2)^sigma + 1)^(1 + nu) - 
        ((y2/mu2)^sigma + 1)^(nu - 2 * (1 + nu)) * (1 + nu) * 
            (y2/mu2)^(sigma * (1 + nu))))/y2  
 
 
 
der2pdf2.mu2dernu <- -(sigma^2 * (nu * (((y2/mu2)^(nu * sigma - 1) + nu * sigma * 
    (log(y2) - log(mu2)) * (y2/mu2)^(nu * sigma - 1))/((y2/mu2)^sigma + 
    1)^(1 + nu) + (y2/mu2)^(nu * sigma - 1)/((y2/mu2)^sigma + 
    1)^(1 + nu) - ((((y2/mu2)^sigma + 1)^(nu - 2 * (1 + nu)) - 
    ((y2/mu2)^sigma + 1)^(nu - 2 * (1 + nu)) * (1 + nu) * log1p((y2/mu2)^sigma)) * 
    (y2/mu2)^(sigma * (1 + nu) - 1) + nu * ((y2/mu2)^sigma + 
    1)^(1 + nu - 2 * (1 + nu)) * log1p((y2/mu2)^sigma) * (y2/mu2)^(nu * 
    sigma - 1) + sigma * ((y2/mu2)^sigma + 1)^(nu - 2 * (1 + 
    nu)) * (1 + nu) * (log(y2) - log(mu2)) * (y2/mu2)^(sigma * 
    (1 + nu) - 1))) - ((y2/mu2)^sigma + 1)^(nu - 2 * (1 + nu)) * 
    (1 + nu) * (y2/mu2)^(sigma * (1 + nu) - 1))/mu2^2) 
    
    
 
der2pdf2.mu2dersigma2 <- -(nu * sigma * (2 * (nu * (y2/mu2)^(nu * sigma - 1)/((y2/mu2)^sigma + 
    1)^(1 + nu) - ((y2/mu2)^sigma + 1)^(nu - 2 * (1 + nu)) * 
    (1 + nu) * (y2/mu2)^(sigma * (1 + nu) - 1)) + sigma * (log(y2) - 
    log(mu2)) * (nu * (nu * (y2/mu2)^(nu * sigma - 1)/((y2/mu2)^sigma + 
    1)^(1 + nu) - ((y2/mu2)^sigma + 1)^(nu - 2 * (1 + nu)) * 
    (1 + nu) * (y2/mu2)^(sigma * (1 + nu) - 1)) - (((y2/mu2)^sigma + 
    1)^(nu - (1 + 2 * (1 + nu))) * (nu - 2 * (1 + nu)) * (y2/mu2)^(sigma * 
    (2 + nu) - 1) + ((y2/mu2)^sigma + 1)^(nu - 2 * (1 + nu)) * 
    (1 + nu) * (y2/mu2)^(sigma * (1 + nu) - 1)) * (1 + nu)))/mu2^2)
        
        
        
         
        
        
if(naive == FALSE){   
 
 
    p2  <- ( 1 + (y2/mu2)^-sigma )^-nu 


derp2.dermu2 <- -(nu * sigma * y2/(mu2^2 * (1 + 1/(y2/mu2)^sigma)^(1 + nu) * 
    (y2/mu2)^(1 + sigma)))  
               
derp2.dersigma2 <- nu * (log(y2) - log(mu2))/((1 + 1/(y2/mu2)^sigma)^(1 + nu) * 
    (y2/mu2)^sigma)        
    

derp2.dernu <- -(log1p(1/(y2/mu2)^sigma)/(1 + 1/(y2/mu2)^sigma)^nu)



der2p2.dermu22 <- nu * sigma * y2 * ((2 * (mu2 * (1 + 1/(y2/mu2)^sigma)^(1 + nu)) + 
    sigma * y2 * (1 + 1/(y2/mu2)^sigma)^nu * (1 + nu)/(y2/mu2)^(1 + 
        sigma)) * (y2/mu2)^(1 + sigma) - y2 * (1 + 1/(y2/mu2)^sigma)^(1 + 
    nu) * (1 + sigma) * (y2/mu2)^sigma)/(mu2^2 * (1 + 1/(y2/mu2)^sigma)^(1 + 
    nu) * (y2/mu2)^(1 + sigma))^2  



der2p2.dersigma22 <- -(nu * ((1 + 1/(y2/mu2)^sigma)^(1 + nu) * (y2/mu2)^sigma - (1 + 
    1/(y2/mu2)^sigma)^nu * (1 + nu)) * (log(y2) - log(mu2))^2/((1 + 
    1/(y2/mu2)^sigma)^(1 + nu) * (y2/mu2)^sigma)^2) 
   
   
   
der2p2.dernu2 <- log1p(1/(y2/mu2)^sigma)^2/(1 + 1/(y2/mu2)^sigma)^nu 
    
    
        
der2p2.dersigma2dernu <- (1/((1 + 1/(y2/mu2)^sigma)^(1 + nu) * (y2/mu2)^sigma) - nu * 
    (1 + 1/(y2/mu2)^sigma)^(1 + nu) * log1p(1/(y2/mu2)^sigma) * 
    (y2/mu2)^sigma/((1 + 1/(y2/mu2)^sigma)^(1 + nu) * (y2/mu2)^sigma)^2) * 
    (log(y2) - log(mu2))
    
    

der2p2.dermu2dernu <- -(sigma * y2 * (1/(mu2^2 * (1 + 1/(y2/mu2)^sigma)^(1 + nu) * 
    (y2/mu2)^(1 + sigma)) - mu2^2 * nu * (1 + 1/(y2/mu2)^sigma)^(1 + 
    nu) * log1p(1/(y2/mu2)^sigma) * (y2/mu2)^(1 + sigma)/(mu2^2 * 
    (1 + 1/(y2/mu2)^sigma)^(1 + nu) * (y2/mu2)^(1 + sigma))^2)) 


der2p2.derdermu2sigma2 <- -(nu * y2 * (1/(mu2^2 * (1 + 1/(y2/mu2)^sigma)^(1 + nu) * (y2/mu2)^(1 + 
    sigma)) - mu2^2 * sigma * ((1 + 1/(y2/mu2)^sigma)^(1 + nu) * 
    (y2/mu2)^(1 + sigma) - y2 * (1 + 1/(y2/mu2)^sigma)^nu * (1 + 
    nu)/mu2) * (log(y2) - log(mu2))/(mu2^2 * (1 + 1/(y2/mu2)^sigma)^(1 + 
    nu) * (y2/mu2)^(1 + sigma))^2)) 
 

    
                                        }


}


##############


if(margin2 == "TW"){

# sigma must not touch the boundaries
# sigma.stTW beween -14.5 and 20 as in eta.tr

aTW <- 1.001  
bTW <- 1.999

mu2        <- exp(eta2)
nu.stTW    <- log(nu)
sigma.stTW <- log( (sigma - aTW) / (bTW - sigma) ) # sigma = (aTW + bTW*exp(sigma.stTW))/(1 + exp(sigma.stTW))


TWob <- ldTweedie(y2, mu = mu2, p = NA, phi = NA, rho = nu.stTW, theta = sigma.stTW, all.derivs = TRUE) 


pdf2           <- exp(TWob[, 1])            
derpdf2.dermu2 <- TWob[, 7]*pdf2
derpdf2.sigma2 <- TWob[, 4]*pdf2 
derpdf2.nu     <- TWob[, 2]*pdf2 

der2pdf2.dermu2     <- pdf2*TWob[, 8] + derpdf2.dermu2^2/pdf2 
der2pdf2.dersigma22 <- pdf2*TWob[, 5] + derpdf2.sigma2^2/pdf2 
der2pdf2.dernu2     <- pdf2*TWob[, 3] +     derpdf2.nu^2/pdf2 
    
der2pdf2.mu2dersigma2   <- pdf2*TWob[, 9]  + derpdf2.dermu2*derpdf2.sigma2/pdf2  
der2pdf2.mu2dernu       <- pdf2*TWob[, 10] + derpdf2.dermu2*derpdf2.nu/pdf2  
der2pdf2.dersigma2dernu <- pdf2*TWob[, 6]  + derpdf2.nu*derpdf2.sigma2/pdf2 
        
        
        
             
if(naive == FALSE){ 
 
if(length(sigma) == 1) sigma <- rep(sigma, length(mu2))
if(length(nu) == 1)    nu    <- rep(nu, length(mu2))



###

dermu2.dereta2      <- 1
der2mu2.dereta2eta2 <- 0   
  
dersigma2.dersigma2.st  <- (bTW - (aTW + bTW * exp(sigma2.st))/(1 + exp(sigma2.st))) * exp(sigma2.st)/(1 + exp(sigma2.st))  
dersigma2.dersigma2.st2 <- (bTW * (1 - exp(sigma2.st)/(1 + exp(sigma2.st))) - ((2 * bTW - 
    2 * ((aTW + bTW * exp(sigma2.st))/(1 + exp(sigma2.st)))) * 
    exp(sigma2.st) + aTW)/(1 + exp(sigma2.st))) * exp(sigma2.st)/(1 + 
    exp(sigma2.st))   

dernu.dernu.st <- dernu.dernu.st2 <- exp(nu.st)  

###



p2 <- derp2.dermu2 <- derp2.dersigma2 <- derp2.dernu <- der2p2.dermu22 <- der2p2.dersigma22 <- der2p2.dernu2 <- der2p2.dersigma2dernu <- der2p2.dermu2dernu <- der2p2.derdermu2sigma2 <- NA 
  
 
for(i in 1:length(mu2)){ 
 
p2[i] <- pTweed(y = y2[i], mu = mu2[i], phi = nu[i], p = sigma[i])$d0 

scTW <- pTweed(y = y2[i], mu = mu2[i], phi = nu[i], p = sigma[i], deriv = 1)$d1
heTW <- pTweed(y = y2[i], mu = mu2[i], phi = nu[i], p = sigma[i], deriv = 2)$d2

derp2.dermu2[i]    <- scTW[1]

derp2.dersigma2[i] <- scTW[3]
derp2.dernu[i]     <- scTW[2]

der2p2.dermu22[i]    <- heTW[1, 1] 
der2p2.dersigma22[i] <- heTW[3, 3]
der2p2.dernu2[i]     <- heTW[2, 2]
    
der2p2.dersigma2dernu[i]  <- heTW[3, 2] 
der2p2.dermu2dernu[i]     <- heTW[2, 1]
der2p2.derdermu2sigma2[i] <- heTW[3, 1]

}



# I attach link function parametrisation here (but need to do that for mu as done for the pdf)


derp2.dereta2                <- derp2.dermu2*dermu2.dereta2
der2p2.dereta2eta2           <- der2p2.dermu22*dermu2.dereta2^2 + derp2.dermu2*der2mu2.dereta2eta2   

derp2.dersigma.st            <- derp2.dersigma2 *  dersigma2.dersigma2.st 
der2p2.dersigma2.st2         <- der2p2.dersigma22 * dersigma2.dersigma2.st^2 + derp2.dersigma2 * dersigma2.dersigma2.st2

der2p2.dereta2dersigma2      <- der2p2.derdermu2sigma2* dermu2.dereta2    
der2p2.dereta2dersigma2.st   <- der2p2.dereta2dersigma2 *  dersigma2.dersigma2.st  

der2p2.dereta2dernu          <- der2p2.dermu2dernu* dermu2.dereta2 
der2p2.dereta2dernu.st       <- der2p2.dereta2dernu * dernu.dernu.st 
der2p2.dersigma2.stdernu.st  <- der2p2.dersigma2dernu * dersigma2.dersigma2.st * dernu.dernu.st  

derp2.nu.st                  <- derp2.dernu *  dernu.dernu.st 
der2p2.dernu.st2             <- der2p2.dernu2 * dernu.dernu.st^2 + derp2.dernu * dernu.dernu.st2


# now change names


derp2.dermu2      <- derp2.dereta2              
der2p2.dermu22    <- der2p2.dereta2eta2         

derp2.dersigma2   <- derp2.dersigma.st          
der2p2.dersigma22 <- der2p2.dersigma2.st2     

derp2.dernu       <- derp2.nu.st                
der2p2.dernu2     <- der2p2.dernu.st2 

der2p2.dermu2dernu     <- der2p2.dereta2dernu.st     
der2p2.derdermu2sigma2 <- der2p2.dereta2dersigma2.st 
der2p2.dersigma2dernu  <- der2p2.dersigma2.stdernu.st

                    }
                   
dermu2.dereta2      <- exp(eta2)
der2mu2.dereta2eta2 <- exp(eta2) 
                   
                    
dersigma2.dersigma2.st  <- 1  
dersigma2.dersigma2.st2 <- 0   

dernu.dernu.st  <- 1
dernu.dernu.st2 <- 0                    
                    


}


##############







if(margin2 == "SM"){


mu2                 <- exp(eta2) 
dermu2.dereta2      <- exp(eta2) 
der2mu2.dereta2eta2 <- exp(eta2) 

dersigma2.dersigma2.st  <- exp(sigma2.st) 
dersigma2.dersigma2.st2 <- exp(sigma2.st)   


dernu.dernu.st  <- exp(nu.st)
dernu.dernu.st2 <- exp(nu.st)

 

pdf2 <- sigma*nu*y2^(sigma-1)*(mu2^sigma*(1+(y2/mu2)^sigma)^(nu+1) )^-1


derpdf2.dermu2 <- -(nu * sigma^2 * y2^(sigma - 1) * (mu2^(sigma - 1) * ((y2/mu2)^sigma + 
    1)^(1 + nu) - mu2^(sigma - 2) * y2 * ((y2/mu2)^sigma + 1)^nu * 
    (1 + nu) * (y2/mu2)^(sigma - 1))/(mu2^sigma * ((y2/mu2)^sigma + 
    1)^(1 + nu))^2)

   
derpdf2.sigma2 <- nu * ((sigma * y2^(sigma - 1) * log(y2) + y2^(sigma - 1))/(mu2^sigma * 
    ((y2/mu2)^sigma + 1)^(1 + nu)) - sigma * y2^(sigma - 1) * 
    (mu2^sigma * ((y2/mu2)^sigma + 1)^(1 + nu) * log(mu2) + mu2^sigma * 
        ((y2/mu2)^sigma + 1)^nu * (1 + nu) * (log(y2) - log(mu2)) * 
        (y2/mu2)^sigma)/(mu2^sigma * ((y2/mu2)^sigma + 1)^(1 + 
    nu))^2)


derpdf2.nu <- sigma * (y2^(sigma - 1)/(mu2^sigma * ((y2/mu2)^sigma + 1)^(1 + 
    nu)) - mu2^sigma * nu * y2^(sigma - 1) * ((y2/mu2)^sigma + 
    1)^(1 + nu) * log1p((y2/mu2)^sigma)/(mu2^sigma * ((y2/mu2)^sigma + 
    1)^(1 + nu))^2)


der2pdf2.dermu2 <- -(nu * sigma^2 * (y2^(sigma - 1) * (mu2^(sigma - 2) * ((y2/mu2)^sigma + 
    1)^(1 + nu) * (sigma - 1) - y2 * ((mu2^(sigma - 3) * ((y2/mu2)^sigma + 
    1)^nu * (sigma - 2) - mu2^(sigma - 4) * nu * sigma * y2 * 
    ((y2/mu2)^sigma + 1)^(nu - 1) * (y2/mu2)^(sigma - 1)) * (y2/mu2)^(sigma - 
    1) + mu2^(sigma - 3) * sigma * ((y2/mu2)^sigma + 1)^nu * 
    (y2/mu2)^(sigma - 1) - mu2^(sigma - 4) * y2 * ((y2/mu2)^sigma + 
    1)^nu * (sigma - 1) * (y2/mu2)^(sigma - 2)) * (1 + nu)) - 
    2 * (mu2^sigma * sigma * y2^(sigma - 1) * ((y2/mu2)^sigma + 
        1)^(1 + nu) * (mu2^(sigma - 1) * ((y2/mu2)^sigma + 1)^(1 + 
        nu) - mu2^(sigma - 2) * y2 * ((y2/mu2)^sigma + 1)^nu * 
        (1 + nu) * (y2/mu2)^(sigma - 1))^2/(mu2^sigma * ((y2/mu2)^sigma + 
        1)^(1 + nu))^2))/(mu2^sigma * ((y2/mu2)^sigma + 1)^(1 + 
    nu))^2)


der2pdf2.dersigma22 <- nu * (log(y2) * (sigma * y2^(sigma - 1) * log(y2) + y2^(sigma - 
    1) + y2^(sigma - 1))/(mu2^sigma * ((y2/mu2)^sigma + 1)^(1 + 
    nu)) - ((mu2^sigma * ((y2/mu2)^sigma + 1)^(1 + nu) * log(mu2) + 
    mu2^sigma * ((y2/mu2)^sigma + 1)^nu * (1 + nu) * (log(y2) - 
        log(mu2)) * (y2/mu2)^sigma) * (sigma * (2 * (y2^(sigma - 
    1) * log(y2)) - 2 * (mu2^sigma * y2^(sigma - 1) * ((y2/mu2)^sigma + 
    1)^(1 + nu) * (mu2^sigma * ((y2/mu2)^sigma + 1)^(1 + nu) * 
    log(mu2) + mu2^sigma * ((y2/mu2)^sigma + 1)^nu * (1 + nu) * 
    (log(y2) - log(mu2)) * (y2/mu2)^sigma)/(mu2^sigma * ((y2/mu2)^sigma + 
    1)^(1 + nu))^2)) + y2^(sigma - 1) + y2^(sigma - 1)) + sigma * 
    y2^(sigma - 1) * (((mu2^sigma * ((y2/mu2)^sigma + 1)^nu * 
    log(mu2) + mu2^sigma * nu * ((y2/mu2)^sigma + 1)^(nu - 1) * 
    (log(y2) - log(mu2)) * (y2/mu2)^sigma) * (y2/mu2)^sigma + 
    mu2^sigma * ((y2/mu2)^sigma + 1)^nu * (log(y2) - log(mu2)) * 
        (y2/mu2)^sigma) * (1 + nu) * (log(y2) - log(mu2)) + log(mu2) * 
    (mu2^sigma * ((y2/mu2)^sigma + 1)^(1 + nu) * log(mu2) + mu2^sigma * 
        ((y2/mu2)^sigma + 1)^nu * (1 + nu) * (log(y2) - log(mu2)) * 
        (y2/mu2)^sigma)))/(mu2^sigma * ((y2/mu2)^sigma + 1)^(1 + 
    nu))^2)


der2pdf2.dernu2 <- -(sigma * log1p((y2/mu2)^sigma) * (mu2^sigma * y2^(sigma - 1) * 
    ((y2/mu2)^sigma + 1)^(1 + nu) + mu2^sigma * y2^(sigma - 1) * 
    ((y2/mu2)^sigma + 1)^(1 + nu) + nu * log1p((y2/mu2)^sigma) * 
    (mu2^sigma * y2^(sigma - 1) * ((y2/mu2)^sigma + 1)^(1 + nu) - 
        2 * (mu2^(3 * sigma) * y2^(sigma - 1) * ((y2/mu2)^sigma + 
            1)^(1 + 2 * (1 + nu) + nu)/(mu2^sigma * ((y2/mu2)^sigma + 
            1)^(1 + nu))^2)))/(mu2^sigma * ((y2/mu2)^sigma + 
    1)^(1 + nu))^2)
    
       
der2pdf2.dersigma2dernu <- sigma * (y2^(sigma - 1) * log(y2)/(mu2^sigma * ((y2/mu2)^sigma + 
    1)^(1 + nu)) - (nu * ((((y2/mu2)^sigma + 1)^(1 + nu) * (mu2^sigma * 
    y2^(sigma - 1) * log(mu2) + mu2^sigma * y2^(sigma - 1) * 
    log(y2)) + mu2^sigma * y2^(sigma - 1) * ((y2/mu2)^sigma + 
    1)^nu * (1 + nu) * (log(y2) - log(mu2)) * (y2/mu2)^sigma - 
    2 * (mu2^(2 * sigma) * y2^(sigma - 1) * ((y2/mu2)^sigma + 
        1)^(2 * (1 + nu)) * (mu2^sigma * ((y2/mu2)^sigma + 1)^(1 + 
        nu) * log(mu2) + mu2^sigma * ((y2/mu2)^sigma + 1)^nu * 
        (1 + nu) * (log(y2) - log(mu2)) * (y2/mu2)^sigma)/(mu2^sigma * 
        ((y2/mu2)^sigma + 1)^(1 + nu))^2)) * log1p((y2/mu2)^sigma) + 
    mu2^sigma * y2^(sigma - 1) * ((y2/mu2)^sigma + 1)^nu * (log(y2) - 
        log(mu2)) * (y2/mu2)^sigma) + y2^(sigma - 1) * (mu2^sigma * 
    ((y2/mu2)^sigma + 1)^(1 + nu) * log(mu2) + mu2^sigma * ((y2/mu2)^sigma + 
    1)^nu * (1 + nu) * (log(y2) - log(mu2)) * (y2/mu2)^sigma))/(mu2^sigma * 
    ((y2/mu2)^sigma + 1)^(1 + nu))^2) + y2^(sigma - 1)/(mu2^sigma * 
    ((y2/mu2)^sigma + 1)^(1 + nu)) - mu2^sigma * nu * y2^(sigma - 
    1) * ((y2/mu2)^sigma + 1)^(1 + nu) * log1p((y2/mu2)^sigma)/(mu2^sigma * 
    ((y2/mu2)^sigma + 1)^(1 + nu))^2 
 
 
der2pdf2.mu2dernu <- -(sigma^2 * (nu * (log1p((y2/mu2)^sigma) * (mu2^(sigma - 1) * 
    y2^(sigma - 1) * ((y2/mu2)^sigma + 1)^(1 + nu) - (2 * (mu2^(2 * 
    sigma) * y2^(sigma - 1) * ((y2/mu2)^sigma + 1)^(2 * (1 + 
    nu)) * (mu2^(sigma - 1) * ((y2/mu2)^sigma + 1)^(1 + nu) - 
    mu2^(sigma - 2) * y2 * ((y2/mu2)^sigma + 1)^nu * (1 + nu) * 
        (y2/mu2)^(sigma - 1))/(mu2^sigma * ((y2/mu2)^sigma + 
    1)^(1 + nu))^2) + mu2^(sigma - 2) * y2^sigma * ((y2/mu2)^sigma + 
    1)^nu * (1 + nu) * (y2/mu2)^(sigma - 1))) - mu2^(sigma - 
    2) * y2^sigma * ((y2/mu2)^sigma + 1)^nu * (y2/mu2)^(sigma - 
    1)) + y2^(sigma - 1) * (mu2^(sigma - 1) * ((y2/mu2)^sigma + 
    1)^(1 + nu) - mu2^(sigma - 2) * y2 * ((y2/mu2)^sigma + 1)^nu * 
    (1 + nu) * (y2/mu2)^(sigma - 1)))/(mu2^sigma * ((y2/mu2)^sigma + 
    1)^(1 + nu))^2)

    
der2pdf2.mu2dersigma2 <-  -(nu * sigma * ((mu2^(sigma - 1) * ((y2/mu2)^sigma + 1)^(1 + 
    nu) - mu2^(sigma - 2) * y2 * ((y2/mu2)^sigma + 1)^nu * (1 + 
    nu) * (y2/mu2)^(sigma - 1)) * (sigma * (y2^(sigma - 1) * 
    log(y2) - 2 * (mu2^sigma * y2^(sigma - 1) * ((y2/mu2)^sigma + 
    1)^(1 + nu) * (mu2^sigma * ((y2/mu2)^sigma + 1)^(1 + nu) * 
    log(mu2) + mu2^sigma * ((y2/mu2)^sigma + 1)^nu * (1 + nu) * 
    (log(y2) - log(mu2)) * (y2/mu2)^sigma)/(mu2^sigma * ((y2/mu2)^sigma + 
    1)^(1 + nu))^2)) + y2^(sigma - 1)) + y2^(sigma - 1) * (((sigma * 
    (log(y2) - log(mu2)) * (mu2^(sigma - 1) * ((y2/mu2)^sigma + 
    1)^nu - mu2^(sigma - 2) * nu * y2 * ((y2/mu2)^sigma + 1)^(nu - 
    1) * (y2/mu2)^(sigma - 1)) - mu2^(sigma - 1) * ((y2/mu2)^sigma + 
    1)^nu) * (y2/mu2)^sigma - mu2^(sigma - 2) * sigma * y2 * 
    ((y2/mu2)^sigma + 1)^nu * (log(y2) - log(mu2)) * (y2/mu2)^(sigma - 
    1)) * (1 + nu) + mu2^(sigma - 1) * ((y2/mu2)^sigma + 1)^(1 + 
    nu) + sigma * log(mu2) * (mu2^(sigma - 1) * ((y2/mu2)^sigma + 
    1)^(1 + nu) - mu2^(sigma - 2) * y2 * ((y2/mu2)^sigma + 1)^nu * 
    (1 + nu) * (y2/mu2)^(sigma - 1))))/(mu2^sigma * ((y2/mu2)^sigma + 
    1)^(1 + nu))^2) 
        
        
                
                
                
                
                
        
if(naive == FALSE){   
 
 
p2  <-  1 - (1+(y2/mu2)^sigma)^-nu

derp2.dermu2 <- -(nu * sigma * y2 * (y2/mu2)^(sigma - 1)/(mu2^2 * ((y2/mu2)^sigma + 
    1)^(1 + nu)))

                                    
derp2.dersigma2 <- nu * (log(y2) - log(mu2)) * (y2/mu2)^sigma/((y2/mu2)^sigma + 
    1)^(1 + nu)         
    

derp2.dernu <- log1p((y2/mu2)^sigma)/((y2/mu2)^sigma + 1)^nu


der2p2.dermu22 <- nu * sigma * y2 * ((2 * (mu2 * ((y2/mu2)^sigma + 1)^(1 + nu)) - 
    sigma * y2 * ((y2/mu2)^sigma + 1)^nu * (1 + nu) * (y2/mu2)^(sigma - 
        1)) * (y2/mu2)^(sigma - 1)/(mu2^2 * ((y2/mu2)^sigma + 
    1)^(1 + nu))^2 + y2 * (sigma - 1) * (y2/mu2)^(sigma - 2)/(mu2^4 * 
    ((y2/mu2)^sigma + 1)^(1 + nu)))
 


der2p2.dersigma22 <- nu * ((y2/mu2)^sigma/((y2/mu2)^sigma + 1)^(1 + nu) - ((y2/mu2)^sigma + 
    1)^(nu - 2 * (1 + nu)) * (1 + nu) * (y2/mu2)^(2 * sigma)) * 
    (log(y2) - log(mu2))^2

   
   
der2p2.dernu2 <- -(log1p((y2/mu2)^sigma)^2/((y2/mu2)^sigma + 1)^nu)


  
der2p2.dersigma2dernu <- ((y2/mu2)^sigma/((y2/mu2)^sigma + 1)^(1 + nu) - nu * log1p((y2/mu2)^sigma) * 
    (y2/mu2)^sigma/((y2/mu2)^sigma + 1)^(1 + nu)) * (log(y2) - 
    log(mu2))

    
der2p2.dermu2dernu <- -(sigma * y2 * ((y2/mu2)^(sigma - 1)/((y2/mu2)^sigma + 1)^(1 + 
    nu) - nu * log1p((y2/mu2)^sigma) * (y2/mu2)^(sigma - 1)/((y2/mu2)^sigma + 
    1)^(1 + nu))/mu2^2)

 


der2p2.derdermu2sigma2 <- -(nu * (((y2/mu2)^sigma + sigma * y2 * (log(y2) - log(mu2)) * 
    (y2/mu2)^(sigma - 1)/mu2)/((y2/mu2)^sigma + 1)^(1 + nu) - 
    sigma * y2 * ((y2/mu2)^sigma + 1)^(nu - 2 * (1 + nu)) * (1 + 
        nu) * (log(y2) - log(mu2)) * (y2/mu2)^(2 * sigma - 1)/mu2)/mu2)


 
                   }


}





if(margin2 %in% c("GP","GPII")){



if(margin2 == "GP"){

                   mu2  <- eta2 # exi
dermu2.dereta2          <- 1
der2mu2.dereta2eta2     <- 0 

}




if(margin2 == "GPII"){

                   mu2  <- exp(eta2) - 0.5 # exi
dermu2.dereta2          <- exp(eta2)
der2mu2.dereta2eta2     <- exp(eta2)

}


dersigma2.dersigma2.st  <- exp(sigma2.st)  # mu
dersigma2.dersigma2.st2 <- exp(sigma2.st)  



indx <- (1 + mu2*y2/sigma) > 0 



pdf2 <- suppressWarnings(   1/sigma*(1 + mu2*y2/sigma)^(-1/mu2-1)    )


derpdf2.dermu2 <- suppressWarnings(  (log1p(mu2 * y2/sigma)/(mu2^2 * (1 + mu2 * y2/sigma)^(1 + 1/mu2)) - 
    y2 * (1 + 1/mu2)/(sigma * (1 + mu2 * y2/sigma)^(1/mu2 + 2)))/sigma   ) 
    
               
derpdf2.sigma2 <- suppressWarnings( -((1/(1 + mu2 * y2/sigma)^(1 + 1/mu2) - mu2 * y2 * (1 + 1/mu2)/(sigma * 
    (1 + mu2 * y2/sigma)^(1/mu2 + 2)))/sigma^2) )
    

der2pdf2.dermu2 <-   suppressWarnings(  (y2 * (2/(mu2^2 * sigma * (1 + mu2 * y2/sigma)^(1/mu2 + 2)) + 
    sigma * (1 + 1/mu2) * (y2 * (1 + mu2 * y2/sigma)^(1 + 1/mu2) * 
        (1/mu2 + 2)/sigma - (1 + mu2 * y2/sigma)^(1/mu2 + 2) * 
        log1p(mu2 * y2/sigma)/mu2^2)/(sigma * (1 + mu2 * y2/sigma)^(1/mu2 + 
        2))^2) - mu2 * (2 * (1 + mu2 * y2/sigma)^(1 + 1/mu2) + 
    mu2 * (y2 * (1 + 1/mu2) * (1 + mu2 * y2/sigma)^(1/mu2)/sigma - 
        (1 + mu2 * y2/sigma)^(1 + 1/mu2) * log1p(mu2 * y2/sigma)/mu2^2)) * 
    log1p(mu2 * y2/sigma)/(mu2^2 * (1 + mu2 * y2/sigma)^(1 + 
    1/mu2))^2)/sigma   )
    
    
      
der2pdf2.dersigma22 <- suppressWarnings(  -((mu2 * y2 * (((1 + mu2 * y2/sigma)^(1/mu2 + 2) - mu2 * y2 * 
    (1 + mu2 * y2/sigma)^(1 + 1/mu2) * (1/mu2 + 2)/sigma)/(sigma * 
    (1 + mu2 * y2/sigma)^(1/mu2 + 2))^2 + (1 + mu2 * y2/sigma)^(1/mu2 - 
    2 * (1 + 1/mu2))/sigma^2) * (1 + 1/mu2) - 2 * ((1/(1 + mu2 * 
    y2/sigma)^(1 + 1/mu2) - mu2 * y2 * (1 + 1/mu2)/(sigma * (1 + 
    mu2 * y2/sigma)^(1/mu2 + 2)))/sigma))/sigma^2)   )
        
     
der2pdf2.mu2dersigma2 <- suppressWarnings(  -(((log1p(mu2 * y2/sigma)/(mu2^2 * (1 + mu2 * y2/sigma)^(1 + 
    1/mu2)) - y2 * (1 + 1/mu2)/(sigma * (1 + mu2 * y2/sigma)^(1/mu2 + 
    2)))/sigma + y2 * ((1/(mu2 * (1 + mu2 * y2/sigma)^(1/mu2 + 
    2)) - mu2^3 * (1 + 1/mu2) * (1 + mu2 * y2/sigma)^(1/mu2) * 
    log1p(mu2 * y2/sigma)/(mu2^2 * (1 + mu2 * y2/sigma)^(1 + 
    1/mu2))^2)/sigma^2 - ((1 + mu2 * y2/sigma)^(1/mu2 + 2) - 
    mu2 * y2 * (1 + mu2 * y2/sigma)^(1 + 1/mu2) * (1/mu2 + 2)/sigma) * 
    (1 + 1/mu2)/(sigma * (1 + mu2 * y2/sigma)^(1/mu2 + 2))^2))/sigma)    )
        


pdf2                  <- ifelse( indx == TRUE, pdf2, 1)                 # 'cause in output log(1) = 0 hence it will not add anything to the lik   
derpdf2.dermu2        <- ifelse( indx == TRUE, derpdf2.dermu2, 0)  
derpdf2.sigma2        <- ifelse( indx == TRUE, derpdf2.sigma2, 0)  
der2pdf2.dermu2       <- ifelse( indx == TRUE, der2pdf2.dermu2, 0)  
der2pdf2.dersigma22   <- ifelse( indx == TRUE, der2pdf2.dersigma22, 0)          
der2pdf2.mu2dersigma2 <- ifelse( indx == TRUE, der2pdf2.mu2dersigma2, 0)  






# not done for cdf for the time being

 
if(naive == FALSE){   
 
 
 
    p2  <- suppressWarnings(    1 - (1 + mu2*y2/sigma)^(-1/mu2)    )
  
  
  
  
    derp2.dermu2 <- suppressWarnings(   -((log1p(mu2 * y2/sigma)/(mu2 * (1 + mu2 * y2/sigma)^(1/mu2)) - 
    y2/(sigma * (1 + mu2 * y2/sigma)^(1 + 1/mu2)))/mu2)    )


                 
derp2.dersigma2 <- suppressWarnings(   -(y2/(sigma^2 * (1 + mu2 * y2/sigma)^(1 + 1/mu2)))    )   
    

der2p2.dermu22 <- suppressWarnings(  -((y2 * (1/(mu2 * sigma * (1 + mu2 * y2/sigma)^(1 + 1/mu2)) + 
    sigma * (y2 * (1 + 1/mu2) * (1 + mu2 * y2/sigma)^(1/mu2)/sigma - 
        (1 + mu2 * y2/sigma)^(1 + 1/mu2) * log1p(mu2 * y2/sigma)/mu2^2)/(sigma * 
        (1 + mu2 * y2/sigma)^(1 + 1/mu2))^2) - (((1 + mu2 * y2/sigma)^(1/mu2) + 
    y2 * (1 + mu2 * y2/sigma)^(1/mu2 - 1)/sigma - (1 + mu2 * 
    y2/sigma)^(1/mu2) * log1p(mu2 * y2/sigma)/mu2) * log1p(mu2 * 
    y2/sigma)/(mu2 * (1 + mu2 * y2/sigma)^(1/mu2))^2 + (log1p(mu2 * 
    y2/sigma)/(mu2 * (1 + mu2 * y2/sigma)^(1/mu2)) - y2/(sigma * 
    (1 + mu2 * y2/sigma)^(1 + 1/mu2)))/mu2))/mu2)  )


der2p2.dersigma22 <- suppressWarnings(  y2 * (2 * (sigma * (1 + mu2 * y2/sigma)^(1 + 1/mu2)) - mu2 * 
    y2 * (1 + 1/mu2) * (1 + mu2 * y2/sigma)^(1/mu2))/(sigma^2 * 
    (1 + mu2 * y2/sigma)^(1 + 1/mu2))^2   )


der2p2.derdermu2sigma2 <- suppressWarnings(   sigma^2 * y2 * (y2 * (1 + 1/mu2) * (1 + mu2 * y2/sigma)^(1/mu2)/sigma - 
    (1 + mu2 * y2/sigma)^(1 + 1/mu2) * log1p(mu2 * y2/sigma)/mu2^2)/(sigma^2 * 
    (1 + mu2 * y2/sigma)^(1 + 1/mu2))^2   )



                                             
                                             
                   }



}







if(margin2 %in% c("GPo")){

                   mu2  <- exp(eta2) - 0.5 # exi
dermu2.dereta2          <- exp(eta2)
der2mu2.dereta2eta2     <- exp(eta2)


dersigma2.dersigma2.st  <- exp(sigma2.st)  # mu
dersigma2.dersigma2.st2 <- exp(sigma2.st)  



indx <- (1 + mu2*y2/(sigma/(1+mu2))) > 0 



pdf2 <- suppressWarnings(   1/(sigma/(1+mu2))*(1 + mu2*y2/(sigma/(1+mu2)))^(-1/mu2-1)    )


derpdf2.dermu2 <- suppressWarnings(  (1 + mu2) * (log1p(mu2 * y2 * (1 + mu2)/sigma)/(mu2^2 * (1 + 
    mu2 * y2 * (1 + mu2)/sigma)^(1 + 1/mu2)) - y2 * ((1 + mu2)/sigma + 
    mu2 * sigma/((1 + mu2)^2 * (sigma/(1 + mu2))^2)) * (1 + 1/mu2)/(1 + 
    mu2 * y2 * (1 + mu2)/sigma)^(1/mu2 + 2))/sigma + sigma/((1 + 
    mu2)^2 * (1 + mu2 * y2 * (1 + mu2)/sigma)^(1 + 1/mu2) * (sigma/(1 + 
    mu2))^2)   ) 
    
               
derpdf2.sigma2 <- suppressWarnings(  -((1/((1 + mu2) * (1 + mu2 * y2 * (1 + mu2)/sigma)^(1 + 1/mu2)) - 
    mu2 * y2 * (1 + 1/mu2)/(sigma * (1 + mu2 * y2 * (1 + mu2)/sigma)^(1/mu2 + 
        2)))/(sigma/(1 + mu2))^2)  )
    

der2pdf2.dermu2 <-   suppressWarnings( (((y2 * ((1 + mu2)/sigma + mu2 * sigma/((1 + mu2)^2 * (sigma/(1 + 
    mu2))^2))/(1 + mu2 * y2 * (1 + mu2)/sigma) - 2 * (log1p(mu2 * 
    y2 * (1 + mu2)/sigma)/mu2))/(1 + mu2 * y2 * (1 + mu2)/sigma)^(1 + 
    1/mu2) + log1p(mu2 * y2 * (1 + mu2)/sigma) * (log1p(mu2 * 
    y2 * (1 + mu2)/sigma)/(mu2^2 * (1 + mu2 * y2 * (1 + mu2)/sigma)^(1 + 
    1/mu2)) - y2 * ((1 + mu2)/sigma + mu2 * sigma/((1 + mu2)^2 * 
    (sigma/(1 + mu2))^2)) * (1 + 1/mu2)/(1 + mu2 * y2 * (1 + 
    mu2)/sigma)^(1/mu2 + 2)))/mu2^2 + y2 * ((((1 + mu2)/sigma + 
    mu2 * sigma/((1 + mu2)^2 * (sigma/(1 + mu2))^2))/mu2^2 - 
    sigma * (1 + 1/mu2) * (2 + mu2 * (2 * (sigma^2/((1 + mu2)^2 * 
        (sigma/(1 + mu2))^2)) - 2)/(1 + mu2))/((1 + mu2)^2 * 
        (sigma/(1 + mu2))^2))/(1 + mu2 * y2 * (1 + mu2)/sigma)^(1/mu2 + 
    2) - ((1 + mu2)/sigma + mu2 * sigma/((1 + mu2)^2 * (sigma/(1 + 
    mu2))^2)) * (1 + 1/mu2) * (log1p(mu2 * y2 * (1 + mu2)/sigma)/(mu2^2 * 
    (1 + mu2 * y2 * (1 + mu2)/sigma)^(1/mu2 + 2)) - y2 * ((1 + 
    mu2)/sigma + mu2 * sigma/((1 + mu2)^2 * (sigma/(1 + mu2))^2)) * 
    (1/mu2 + 2)/(1 + mu2 * y2 * (1 + mu2)/sigma)^(1/mu2 + 3)))) * 
    (1 + mu2)/sigma + sigma * (2 * (log1p(mu2 * y2 * (1 + mu2)/sigma)/(mu2^2 * 
    (1 + mu2 * y2 * (1 + mu2)/sigma)^(1 + 1/mu2))) - ((2 - 2 * 
    (sigma^2/((1 + mu2)^2 * (sigma/(1 + mu2))^2)))/((1 + mu2) * 
    (1 + mu2 * y2 * (1 + mu2)/sigma)^(1 + 1/mu2)) + 2 * (y2 * 
    ((1 + mu2)/sigma + mu2 * sigma/((1 + mu2)^2 * (sigma/(1 + 
        mu2))^2)) * (1 + 1/mu2)/(1 + mu2 * y2 * (1 + mu2)/sigma)^(1/mu2 + 
    2))))/((1 + mu2)^2 * (sigma/(1 + mu2))^2)  )
    
    
      
der2pdf2.dersigma22 <- suppressWarnings( ((2 * (sigma/((1 + mu2) * (1 + mu2 * y2 * (1 + mu2)/sigma)^(1 + 
    1/mu2))) - mu2 * y2 * (1 + 1/mu2)/(1 + mu2 * y2 * (1 + mu2)/sigma)^(1/mu2 + 
    2))/(1 + mu2) + mu2 * y2 * ((mu2 * y2 * (1/mu2 + 2)/(1 + 
    mu2 * y2 * (1 + mu2)/sigma)^(1/mu2 + 3) - 2 * (sigma/((1 + 
    mu2) * (1 + mu2 * y2 * (1 + mu2)/sigma)^(1/mu2 + 2))))/sigma - 
    1/((1 + mu2) * (1 + mu2 * y2 * (1 + mu2)/sigma)^(1/mu2 + 
        2))) * (1 + 1/mu2))/((1 + mu2) * (sigma/(1 + mu2))^4)   )
        
     
der2pdf2.mu2dersigma2 <- suppressWarnings(   -(((log1p(mu2 * y2 * (1 + mu2)/sigma)/(mu2^2 * (1 + mu2 * y2 * 
    (1 + mu2)/sigma)^(1 + 1/mu2)) - ((1 - 2 * (sigma^2/((1 + 
    mu2)^2 * (sigma/(1 + mu2))^2)))/((1 + mu2) * (1 + mu2 * y2 * 
    (1 + mu2)/sigma)^(1 + 1/mu2)) + y2 * ((1 + mu2)/sigma + mu2 * 
    sigma/((1 + mu2)^2 * (sigma/(1 + mu2))^2)) * (1 + 1/mu2)/(1 + 
    mu2 * y2 * (1 + mu2)/sigma)^(1/mu2 + 2)))/(1 + mu2) + y2 * 
    (((1/mu2 - (1 + 1/mu2) * (1 + mu2 * (2 * (sigma^2/((1 + mu2)^2 * 
        (sigma/(1 + mu2))^2)) - 1)/(1 + mu2)))/(1 + mu2 * y2 * 
        (1 + mu2)/sigma)^(1/mu2 + 2) - mu2 * (1 + 1/mu2) * (log1p(mu2 * 
        y2 * (1 + mu2)/sigma)/(mu2^2 * (1 + mu2 * y2 * (1 + mu2)/sigma)^(1/mu2 + 
        2)) - y2 * ((1 + mu2)/sigma + mu2 * sigma/((1 + mu2)^2 * 
        (sigma/(1 + mu2))^2)) * (1/mu2 + 2)/(1 + mu2 * y2 * (1 + 
        mu2)/sigma)^(1/mu2 + 3)))/sigma - mu2 * sigma * (1 + 
        1/mu2)/((1 + mu2)^3 * (1 + mu2 * y2 * (1 + mu2)/sigma)^(1/mu2 + 
        2) * (sigma/(1 + mu2))^2)))/(sigma/(1 + mu2))^2)   )
        


pdf2                  <- ifelse( indx == TRUE, pdf2, 1)                 # 'cause in output log(1) = 0 hence it will not add anything to the lik   
derpdf2.dermu2        <- ifelse( indx == TRUE, derpdf2.dermu2, 0)  
derpdf2.sigma2        <- ifelse( indx == TRUE, derpdf2.sigma2, 0)  
der2pdf2.dermu2       <- ifelse( indx == TRUE, der2pdf2.dermu2, 0)  
der2pdf2.dersigma22   <- ifelse( indx == TRUE, der2pdf2.dersigma22, 0)          
der2pdf2.mu2dersigma2 <- ifelse( indx == TRUE, der2pdf2.mu2dersigma2, 0)  




# not done for cdf for the time being

 
if(naive == FALSE){   
 
    p2  <- suppressWarnings(    1 - (1 + mu2*y2/ (sigma/(1+mu2)) )^(-1/mu2)    )
  

    derp2.dermu2 <- suppressWarnings(  -((log1p(mu2 * y2 * (1 + mu2)/sigma)/(mu2 * (1 + mu2 * y2 * (1 + 
    mu2)/sigma)^(1/mu2)) - y2 * ((1 + mu2)/sigma + mu2 * sigma/((1 + 
    mu2)^2 * (sigma/(1 + mu2))^2))/(1 + mu2 * y2 * (1 + mu2)/sigma)^(1 + 
    1/mu2))/mu2)    )


                 
derp2.dersigma2 <- suppressWarnings(  -(y2/((1 + mu2) * (1 + mu2 * y2 * (1 + mu2)/sigma)^(1 + 1/mu2) * 
    (sigma/(1 + mu2))^2))   )   
    

der2p2.dermu22 <- suppressWarnings( -((((y2 * ((1 + mu2)/sigma + mu2 * sigma/((1 + mu2)^2 * (sigma/(1 + 
    mu2))^2))/(1 + mu2 * y2 * (1 + mu2)/sigma) - 2 * (log1p(mu2 * 
    y2 * (1 + mu2)/sigma)/mu2))/(1 + mu2 * y2 * (1 + mu2)/sigma)^(1/mu2) + 
    log1p(mu2 * y2 * (1 + mu2)/sigma) * (log1p(mu2 * y2 * (1 + 
        mu2)/sigma)/(mu2 * (1 + mu2 * y2 * (1 + mu2)/sigma)^(1/mu2)) - 
        y2 * ((1 + mu2)/sigma + mu2 * sigma/((1 + mu2)^2 * (sigma/(1 + 
            mu2))^2))/(1 + mu2 * y2 * (1 + mu2)/sigma)^(1 + 1/mu2))/mu2)/mu2 + 
    y2 * ((((1 + mu2)/sigma + mu2 * sigma/((1 + mu2)^2 * (sigma/(1 + 
        mu2))^2))/mu2 - sigma * (2 + mu2 * (2 * (sigma^2/((1 + 
        mu2)^2 * (sigma/(1 + mu2))^2)) - 2)/(1 + mu2))/((1 + 
        mu2)^2 * (sigma/(1 + mu2))^2))/(1 + mu2 * y2 * (1 + mu2)/sigma)^(1 + 
        1/mu2) - ((1 + mu2)/sigma + mu2 * sigma/((1 + mu2)^2 * 
        (sigma/(1 + mu2))^2)) * (log1p(mu2 * y2 * (1 + mu2)/sigma)/(mu2^2 * 
        (1 + mu2 * y2 * (1 + mu2)/sigma)^(1 + 1/mu2)) - y2 * 
        ((1 + mu2)/sigma + mu2 * sigma/((1 + mu2)^2 * (sigma/(1 + 
            mu2))^2)) * (1 + 1/mu2)/(1 + mu2 * y2 * (1 + mu2)/sigma)^(1/mu2 + 
        2))))/mu2)  )


der2p2.dersigma22 <- suppressWarnings( -(y2 * (mu2 * y2 * (1 + 1/mu2)/(1 + mu2 * y2 * (1 + mu2)/sigma)^(1/mu2 + 
    2) - 2 * (sigma/((1 + mu2) * (1 + mu2 * y2 * (1 + mu2)/sigma)^(1 + 
    1/mu2))))/((1 + mu2)^2 * (sigma/(1 + mu2))^4))   )


der2p2.derdermu2sigma2 <- suppressWarnings( -(y2 * ((2 * (sigma^2/((1 + mu2)^2 * (sigma/(1 + mu2))^2)) - 
    1)/((1 + mu2) * (1 + mu2 * y2 * (1 + mu2)/sigma)^(1 + 1/mu2)) + 
    log1p(mu2 * y2 * (1 + mu2)/sigma)/(mu2^2 * (1 + mu2 * y2 * 
        (1 + mu2)/sigma)^(1 + 1/mu2)) - y2 * ((1 + mu2)/sigma + 
    mu2 * sigma/((1 + mu2)^2 * (sigma/(1 + mu2))^2)) * (1 + 1/mu2)/(1 + 
    mu2 * y2 * (1 + mu2)/sigma)^(1/mu2 + 2))/((1 + mu2) * (sigma/(1 + 
    mu2))^2))   )



                                             
                                             
                   }



}



































if(margin2 == "FISK"){


mu2                 <- exp(eta2)
dermu2.dereta2      <- exp(eta2)
der2mu2.dereta2eta2 <- exp(eta2) 

dersigma2.dersigma2.st  <- exp(sigma2.st) 
dersigma2.dersigma2.st2 <- exp(sigma2.st)  





pdf2 <- sigma*y2^(sigma-1) / (mu2^sigma*(1+(y2/mu2)^sigma)^2)
  
derpdf2.dermu2 <- -(sigma^2 * y2^(sigma - 1) * ((y2/mu2)^sigma + 1) * (mu2^(sigma - 
    1) * ((y2/mu2)^sigma + 1) - 2 * (mu2^(sigma - 2) * y2 * (y2/mu2)^(sigma - 
    1)))/(mu2^sigma * ((y2/mu2)^sigma + 1)^2)^2)

                
derpdf2.sigma2 <- (sigma * y2^(sigma - 1) * log(y2) + y2^(sigma - 1))/(mu2^sigma * 
    ((y2/mu2)^sigma + 1)^2) - sigma * y2^(sigma - 1) * ((y2/mu2)^sigma + 
    1) * (2 * (mu2^sigma * (log(y2) - log(mu2)) * (y2/mu2)^sigma) + 
    mu2^sigma * ((y2/mu2)^sigma + 1) * log(mu2))/(mu2^sigma * 
    ((y2/mu2)^sigma + 1)^2)^2


der2pdf2.dermu2 <- -(sigma^2 * (((y2/mu2)^sigma + 1) * (y2^(sigma - 1) * (mu2^(sigma - 
    2) * ((y2/mu2)^sigma + 1) * (sigma - 1) - y2 * (2 * (mu2^(sigma - 
    3) * (sigma - 2) * (y2/mu2)^(sigma - 1) - mu2^(sigma - 4) * 
    y2 * (sigma - 1) * (y2/mu2)^(sigma - 2)) + mu2^(sigma - 3) * 
    sigma * (y2/mu2)^(sigma - 1))) - 2 * (mu2^sigma * sigma * 
    y2^(sigma - 1) * ((y2/mu2)^sigma + 1)^3 * (mu2^(sigma - 1) * 
    ((y2/mu2)^sigma + 1) - 2 * (mu2^(sigma - 2) * y2 * (y2/mu2)^(sigma - 
    1)))^2/(mu2^sigma * ((y2/mu2)^sigma + 1)^2)^2)) - sigma * 
    y2^sigma * (mu2^(sigma - 1) * ((y2/mu2)^sigma + 1) - 2 * 
    (mu2^(sigma - 2) * y2 * (y2/mu2)^(sigma - 1))) * (y2/mu2)^(sigma - 
    1)/mu2^2)/(mu2^sigma * ((y2/mu2)^sigma + 1)^2)^2)

            
der2pdf2.dersigma22 <- log(y2) * (sigma * y2^(sigma - 1) * log(y2) + y2^(sigma - 1) + 
    y2^(sigma - 1))/(mu2^sigma * ((y2/mu2)^sigma + 1)^2) - ((((y2/mu2)^sigma + 
    1) * (sigma * y2^(sigma - 1) * log(y2) + y2^(sigma - 1)) + 
    sigma * y2^(sigma - 1) * (log(y2) - log(mu2)) * (y2/mu2)^sigma) * 
    (2 * (mu2^sigma * (log(y2) - log(mu2)) * (y2/mu2)^sigma) + 
        mu2^sigma * ((y2/mu2)^sigma + 1) * log(mu2)) + ((2 * 
    (mu2^sigma * (log(y2) - log(mu2)) * (y2/mu2)^sigma) + mu2^sigma * 
    ((y2/mu2)^sigma + 1) * log(mu2)) * (sigma * y2^(sigma - 1) * 
    log(y2) + y2^(sigma - 1)) + sigma * (y2^(sigma - 1) * (2 * 
    ((log(y2) - log(mu2)) * (mu2^sigma * (log(y2) - log(mu2)) * 
        (y2/mu2)^sigma + mu2^sigma * log(mu2) * (y2/mu2)^sigma)) + 
    log(mu2) * (mu2^sigma * ((y2/mu2)^sigma + 1) * log(mu2) + 
        mu2^sigma * (log(y2) - log(mu2)) * (y2/mu2)^sigma)) - 
    2 * (mu2^sigma * y2^(sigma - 1) * ((y2/mu2)^sigma + 1)^3 * 
        (2 * (mu2^sigma * (log(y2) - log(mu2)) * (y2/mu2)^sigma) + 
            mu2^sigma * ((y2/mu2)^sigma + 1) * log(mu2))^2/(mu2^sigma * 
        ((y2/mu2)^sigma + 1)^2)^2))) * ((y2/mu2)^sigma + 1))/(mu2^sigma * 
    ((y2/mu2)^sigma + 1)^2)^2

         
der2pdf2.mu2dersigma2 <- -(sigma * ((((y2/mu2)^sigma + 1) * (2 * y2^(sigma - 1) + sigma * 
    y2^(sigma - 1) * log(y2)) + sigma * y2^(sigma - 1) * (log(y2) - 
    log(mu2)) * (y2/mu2)^sigma) * (mu2^(sigma - 1) * ((y2/mu2)^sigma + 
    1) - 2 * (mu2^(sigma - 2) * y2 * (y2/mu2)^(sigma - 1))) + 
    sigma * ((y2/mu2)^sigma + 1) * (y2^(sigma - 1) * (mu2^(sigma - 
        1) * ((y2/mu2)^sigma + 1) * log(mu2) + mu2^(sigma - 1) * 
        (log(y2) - log(mu2)) * (y2/mu2)^sigma - 2 * (y2 * (mu2^(sigma - 
        2) * (log(y2) - log(mu2)) * (y2/mu2)^(sigma - 1) + mu2^(sigma - 
        2) * log(mu2) * (y2/mu2)^(sigma - 1)))) - 2 * (mu2^sigma * 
        y2^(sigma - 1) * ((y2/mu2)^sigma + 1)^3 * (2 * (mu2^sigma * 
        (log(y2) - log(mu2)) * (y2/mu2)^sigma) + mu2^sigma * 
        ((y2/mu2)^sigma + 1) * log(mu2)) * (mu2^(sigma - 1) * 
        ((y2/mu2)^sigma + 1) - 2 * (mu2^(sigma - 2) * y2 * (y2/mu2)^(sigma - 
        1)))/(mu2^sigma * ((y2/mu2)^sigma + 1)^2)^2)))/(mu2^sigma * 
    ((y2/mu2)^sigma + 1)^2)^2)





  
 
if(naive == FALSE){   
 
    p2  <- 1/(1+(y2/mu2)^-sigma)
    
    derp2.dermu2 <- -((y2/mu2)^-(sigma + 1) * (sigma * (y2/mu2^2))/(1 + (y2/mu2)^-sigma)^2)


                 
derp2.dersigma2 <- (y2/mu2)^-sigma * log((y2/mu2))/(1 + (y2/mu2)^-sigma)^2 
             
    

der2p2.dermu22 <- -(sigma * y2 * (y2 * ((1 + sigma)/(y2/mu2)^(2 + sigma) - 2 * 
    (sigma/((1 + 1/(y2/mu2)^sigma) * (y2/mu2)^(2 * (1 + sigma)))))/mu2 - 
    2/(y2/mu2)^(1 + sigma))/(mu2^3 * (1 + 1/(y2/mu2)^sigma)^2))


der2p2.dersigma22 <- -((1/(y2/mu2)^sigma - 2/((1 + 1/(y2/mu2)^sigma) * (y2/mu2)^(2 * 
    sigma))) * (log(y2) - log(mu2))^2/(1 + 1/(y2/mu2)^sigma)^2)



der2p2.derdermu2sigma2 <- (sigma * y2 * (1/(y2/mu2)^(1 + sigma) - 2/((1 + 1/(y2/mu2)^sigma) * 
    (y2/mu2)^(1 + 2 * sigma))) * (log(y2) - log(mu2))/mu2 - 1/(y2/mu2)^sigma)/(mu2 * 
    (1 + 1/(y2/mu2)^sigma)^2)

                                               
                                             
                   }



}

####

if(margin2 == "WEI"){


mu2                 <- exp(eta2) 
dermu2.dereta2      <- exp(eta2) 
der2mu2.dereta2eta2 <- exp(eta2) 

dersigma2.dersigma2.st  <- exp(sigma2.st)  
dersigma2.dersigma2.st2 <- exp(sigma2.st)  


pdf2 <- sigma/mu2*(y2/mu2)^(sigma-1) * exp(-(y2/mu2)^sigma)  
  
  
derpdf2.dermu2 <- sigma * exp(-(y2/mu2)^sigma) * (y2 * (sigma * (y2/mu2)^(2 * (sigma - 
    1)) - (sigma - 1) * (y2/mu2)^(sigma - 2))/mu2 - (y2/mu2)^(sigma - 
    1))/mu2^2 
  
                  
derpdf2.sigma2 <- ((y2/mu2)^(sigma - 1) + sigma * ((y2/mu2)^(sigma - 1) - (y2/mu2)^(2 * 
    sigma - 1)) * (log(y2) - log(mu2))) * exp(-(y2/mu2)^sigma)/mu2


der2pdf2.dermu2 <- sigma * exp(-(y2/mu2)^sigma) * (y2 * (sigma * ((y2 * (sigma * 
    (y2/mu2)^(2 * (sigma - 1)) - (sigma - 1) * (y2/mu2)^(sigma - 
    2))/mu2 - (y2/mu2)^(sigma - 1)) * (y2/mu2)^(sigma - 1) - 
    (y2/mu2)^(2 * (sigma - 1))) - (sigma - 1) * (y2 * (2 * (sigma * 
    (y2/mu2)^(2 * (sigma - 1) - 1)) - (sigma - 2) * (y2/mu2)^(sigma - 
    3))/mu2 - 2 * (y2/mu2)^(sigma - 2)))/mu2 - 2 * (y2 * (sigma * 
    (y2/mu2)^(2 * (sigma - 1)) - (sigma - 1) * (y2/mu2)^(sigma - 
    2))/mu2 - (y2/mu2)^(sigma - 1)))/mu2^3
            
            
 
der2pdf2.dersigma22 <- (2 * (y2/mu2)^(sigma - 1) + sigma * ((y2/mu2)^(sigma - 1) - 2 * 
    (y2/mu2)^(2 * sigma - 1)) * (log(y2) - log(mu2)) - (((y2/mu2)^(sigma - 
    1) + sigma * ((y2/mu2)^(sigma - 1) - (y2/mu2)^(2 * sigma - 
    1)) * (log(y2) - log(mu2))) * (y2/mu2)^sigma + (y2/mu2)^(2 * 
    sigma - 1))) * exp(-(y2/mu2)^sigma) * (log(y2) - log(mu2))/mu2
       
       
       
der2pdf2.mu2dersigma2 <- ((1 - sigma * (log(y2) - log(mu2)) * (y2/mu2)^sigma) * (y2 * 
    (sigma * (y2/mu2)^(2 * (sigma - 1)) - (sigma - 1) * (y2/mu2)^(sigma - 
        2))/mu2 - (y2/mu2)^(sigma - 1)) + sigma * (y2 * ((2 * 
    (sigma * (y2/mu2)^(2 * (sigma - 1))) - (sigma - 1) * (y2/mu2)^(sigma - 
    2)) * (log(y2) - log(mu2)) + (y2/mu2)^(2 * (sigma - 1)) - 
    (y2/mu2)^(sigma - 2))/mu2 - (log(y2) - log(mu2)) * (y2/mu2)^(sigma - 
    1))) * exp(-(y2/mu2)^sigma)/mu2^2 

  
 
if(naive == FALSE){   
 
    p2  <-  1-exp(-(y2/mu2)^sigma) 
    
    derp2.dermu2 <- -(sigma * y2 * exp(-(y2/mu2)^sigma) * (y2/mu2)^(sigma - 1)/mu2^2)
                 
derp2.dersigma2 <-  exp(-(y2/mu2)^sigma) * ((y2/mu2)^sigma * log((y2/mu2)))      
    

der2p2.dermu22 <- -(sigma * y2 * exp(-(y2/mu2)^sigma) * (y2 * (sigma * (y2/mu2)^(2 * 
    (sigma - 1)) - (sigma - 1) * (y2/mu2)^(sigma - 2))/mu2 - 
    2 * (y2/mu2)^(sigma - 1))/mu2^3)


der2p2.dersigma22 <-  ((y2/mu2)^sigma - (y2/mu2)^(2 * sigma)) * exp(-(y2/mu2)^sigma) * 
    (log(y2) - log(mu2))^2


der2p2.derdermu2sigma2 <- exp(-(y2/mu2)^sigma) * (sigma * y2 * ((y2/mu2)^(2 * sigma - 1) - 
    (y2/mu2)^(sigma - 1)) * (log(y2) - log(mu2))/mu2 - (y2/mu2)^sigma)/mu2 
   
   
   
                                             }





}


#####################




if(margin2 == "BE"){

#########################

mu2   <- plogis(eta2)
sigma <- sigma2 <- plogis(sigma2.st)  

dermu2.dereta2      <- (1 - exp(eta2)/(1 + exp(eta2))) * exp(eta2)/(1 + exp(eta2))
der2mu2.dereta2eta2 <-  (1 - (3 - 2 * (exp(eta2)/(1 + exp(eta2)))) * exp(eta2)/(1 + exp(eta2))) * 
                        exp(eta2)/(1 + exp(eta2))
         
         

dersigma2.dersigma2.st  <- (1 - exp(sigma2.st)/(1 + exp(sigma2.st))) * exp(sigma2.st)/(1 + exp(sigma2.st))
dersigma2.dersigma2.st2 <- (1 - 2 * (exp(sigma2.st)/(1 + exp(sigma2.st)))) * (1 - exp(sigma2.st)/(1 + exp(sigma2.st))) * exp(sigma2.st)/(1 + exp(sigma2.st))


#########################

#    a <- mu2 * (1 - sigma^2)/(sigma^2)
#    b <- a * (1 - mu2)/mu2
#    a <- mu * (1 - sigma^2)/(sigma^2)
#    b <- a * (1 - mu)/mu
#    dbeta(x, shape1 = a, shape2 = b, ncp = 0, log = log)
#    1/beta(mu2 * (1 - sigma^2)/(sigma^2), (1-mu2)*(1 - sigma^2)/(sigma^2))*y2^( mu2 * (1 - sigma^2)/(sigma^2) -1)*(1-y2)^((1-mu2)*(1 - sigma^2)/(sigma^2)-1) 


pdf2 <- dbeta(y2, shape1 = mu2 * (1 - sigma^2)/(sigma^2), shape2 = (1-mu2)*(1 - sigma^2)/(sigma^2))
  
  
  
derpdf2.dermu2 <- ((-1 + sigma^2)* (1 - y2)^(-1 + ((-1 + mu2)* (-1 + sigma^2))/sigma^2)*
    y2^(-1 + mu2 * (-1 + 1/sigma^2))* (log(1 - y2) - log(y2) + 
     psigamma(mu2* (-1 + 1/sigma^2)) - 
     psigamma(((-1 + mu2)* (-1 + sigma^2))/sigma^2)))/(sigma^2* beta(
    mu2 *(-1 + 1/sigma^2), ((-1 + mu2)* (-1 + sigma^2))/sigma^2))
  
                     
derpdf2.sigma2 <- (-((1 - y2)^(-1 + ((-1 + mu2)* (-1 + sigma^2))/sigma^2)*
     y2^(-1 + mu2 * (-1 + 1/sigma^2)) *(log(1 - y2) - mu2* log(1 - y2) + 
      mu2 * log(y2) + psigamma(-1 + 1/sigma^2) - 
      mu2* psigamma(mu2 *(-1 + 1/sigma^2)) - 
      psigamma(((-1 + mu2) *(-1 + sigma^2))/sigma^2) + 
      mu2* psigamma(((-1 + mu2)* (-1 + sigma^2))/
        sigma^2)))/((sigma^2)^2 *beta(
     mu2* (-1 + 1/sigma^2), ((-1 + mu2)* (-1 + sigma^2))/sigma^2)))*2*sigma 


der2pdf2.dermu2 <- ((-1 + sigma^2)^2* (1 - y2)^(-1 + ((-1 + mu2)* (-1 + sigma^2))/sigma^2)*
    y2^(-1 + 
    mu2* (-1 + 1/sigma^2)) *(log(1 - y2)^2 - 2* log(1 - y2)* log(y2) + 
     log(y2)^2 + psigamma(mu2* (-1 + 1/sigma^2))^2 + 
     2 *psigamma( 
       mu2* (-1 + 1/sigma^2))* (log(1 - y2) - log(y2) - 
        psigamma(((-1 + mu2) *(-1 + sigma^2))/sigma^2)) - 
     2 *(log(1 - y2) - log(y2))* psigamma(
       ((-1 + mu2)* (-1 + sigma^2))/sigma^2) + 
     psigamma(((-1 + mu2)* (-1 + sigma^2))/sigma^2)^2 - 
     psigamma(mu2 *(-1 + 1/sigma^2),1) - 
     psigamma( ((-1 + mu2)* (-1 + sigma^2))/sigma^2,1)))/((sigma^2)^2* beta(
    mu2* (-1 + 1/sigma^2), ((-1 + mu2) *(-1 + sigma^2))/sigma^2)) 
  
  
  
  
der2pdf2.dersigma22 <- (-((1 - y2)^(-1 + ((-1 + mu2) *(-1 + sigma^2))/sigma^2)*
     y2^(-1 + 
     mu2* (-1 + 1/sigma^2))* (-2 *sigma^2* log(1 - y2) + 
      2* mu2* sigma^2* log(1 - y2) - log(1 - y2)^2 + 2* mu2* log(1 - y2)^2 -
       mu2^2* log(1 - y2)^2 - 2* mu2* sigma^2* log(y2) - 
      2 *mu2* log(1 - y2)* log(y2) + 2* mu2^2 *log(1 - y2)* log(y2) - 
      mu2^2* log(y2)^2 - psigamma(-1 + 1/sigma^2)^2 - 
      mu2^2 *psigamma(mu2 *(-1 + 1/sigma^2))^2 + 
      2 *sigma^2* psigamma(((-1 + mu2) *(-1 + sigma^2))/sigma^2) - 
      2 *mu2* sigma^2* psigamma( ((-1 + mu2) *(-1 + sigma^2))/sigma^2) + 
      2* log(1 - y2)* psigamma(((-1 + mu2) *(-1 + sigma^2))/sigma^2) - 
      4 *mu2* log(1 - y2)* psigamma(((-1 + mu2)* (-1 + sigma^2))/
        sigma^2) + 
      2* mu2^2* log(1 - y2)* psigamma(((-1 + mu2)* (-1 + sigma^2))/
        sigma^2) + 
      2 *mu2 *log(y2)* psigamma(((-1 + mu2)* (-1 + sigma^2))/sigma^2) - 
      2 *mu2^2* log(y2)* psigamma( ((-1 + mu2)* (-1 + sigma^2))/
        sigma^2) - psigamma( ((-1 + mu2)* (-1 + sigma^2))/sigma^2)^2 + 
      2* mu2* psigamma(((-1 + mu2)* (-1 + sigma^2))/sigma^2)^2 - 
      mu2^2 *psigamma( ((-1 + mu2)* (-1 + sigma^2))/sigma^2)^2 + 
      2* mu2* psigamma( 
        mu2* (-1 + 1/sigma^2))* (sigma^2 - (-1 + mu2)* log(1 - y2) + 
         mu2 *log(y2) + (-1 + mu2)* psigamma(
           ((-1 + mu2)* (-1 + sigma^2))/sigma^2)) - 
      2* psigamma(-1 + 1/sigma^2)* (sigma^2 + log(1 - y2) - mu2* log(1 - y2) + 
         mu2 *log(y2) - 
         mu2* psigamma( mu2* (-1 + 1/sigma^2)) + (-1 + mu2)* psigamma(((-1 + mu2)* (-1 + sigma^2))/sigma^2)) - 
      psigamma(-1 + 1/sigma^2,1) + 
      mu2^2* psigamma( mu2* (-1 + 1/sigma^2),1) + 
      psigamma( ((-1 + mu2)* (-1 + sigma^2))/sigma^2,1) - 
      2* mu2* psigamma(((-1 + mu2)* (-1 + sigma^2))/sigma^2,1) + 
      mu2^2 *psigamma(((-1 + mu2)* (-1 + sigma^2))/
        sigma^2,1)))/((sigma^2)^4* beta(
     mu2 *(-1 + 1/sigma^2), ((-1 + mu2)* (-1 + sigma^2))/sigma^2)))*(2*sigma)^2 + derpdf2.sigma2/sigma 

      
der2pdf2.mu2dersigma2 <- (((1 - y2)^(-1 + ((-1 + mu2) *(-1 + sigma^2))/sigma^2)*
    y2^(-1 + 
    mu2* (-1 + 1/
       sigma^2))* (-(-1 + sigma^2)* sigma^2* (log(1 - y2) - log(y2) + 
        psigamma(mu2* (-1 + 1/sigma^2)) - 
        psigamma(((-1 + mu2)* (-1 + sigma^2))/sigma^2)) + 
     (sigma^2)^2* (log(1 - y2) - log(y2) + 
        psigamma(mu2 *(-1 + 1/sigma^2)) - 
        psigamma(((-1 + mu2)* (-1 + sigma^2))/sigma^2)) + (-1 + 
        mu2)* (-1 + sigma^2)* log(
       1 - y2)* (log(1 - y2) - log(y2) + 
        psigamma(mu2* (-1 + 1/sigma^2)) - 
        psigamma(((-1 + mu2)* (-1 + sigma^2))/sigma^2)) + 
     mu2* (-1 + sigma^2)* log(
       y2)* (-log(1 - y2) + log(y2) - 
        psigamma(mu2 *(-1 + 1/sigma^2)) + 
        psigamma(((-1 + mu2)* (-1 + sigma^2))/sigma^2)) + (-1 + 
        sigma^2)* (log(1 - y2) - log(y2) + 
        psigamma(mu2* (-1 + 1/sigma^2)) - 
        psigamma(((-1 + mu2)* (-1 + sigma^2))/
         sigma^2))* (-psigamma(-1 + 1/sigma^2) + 
        mu2* psigamma(mu2 *(-1 + 1/sigma^2)) - (-1 + mu2)* psigamma(
          ((-1 + mu2)* (-1 + sigma^2))/sigma^2)) - (-1 + 
        sigma^2)* (mu2* psigamma( 
          mu2 *(-1 + 1/sigma^2),1) + (-1 + mu2)* psigamma(
          ((-1 + mu2) *(-1 + sigma^2))/sigma^2,1))))/((sigma^2)^3* beta(
    mu2 *(-1 + 1/sigma^2), ((-1 + mu2)* (-1 + sigma^2))/sigma^2))  )*2*sigma  
   

   
if(naive == FALSE){   


p2  <-  pbeta(y2, shape1 = mu2 * (1 - sigma^2)/(sigma^2), shape2 = (1-mu2)*(1 - sigma^2)/(sigma^2))

funcD <- function(para) pbeta(y2, shape1 = para * (1 - sigma^2)/sigma^2, shape2 = (1-para)*(1 - sigma^2)/sigma^2)
 
nde <- numgh(funcD, mu2) 
 
derp2.dermu2   <- nde$fi
der2p2.dermu22 <- nde$se

funcD <- function(para) pbeta(y2, shape1 = mu2 * (1 - para^2)/para^2, shape2 = (1-mu2)*(1 - para^2)/para^2 )
 

nde <- numgh(funcD, sigma) 
 
derp2.dersigma2   <- nde$fi
der2p2.dersigma22 <- nde$se

funcD1 <- function(pms1, pms2) pbeta(y2, shape1 = pms1 * (1 - pms2^2)/pms2^2, shape2 = (1-pms1)*(1 - pms2^2)/pms2^2)

der2p2.derdermu2sigma2 <- numch(funcD1, mu2, sigma)


                                             }




}












####

if(margin2 == "iG"){


#
#sigma2    <- ifelse(sigma2 < 0.0001234098, 0.0001234098, sigma2)
#sigma2.st <- log(sigma2) 
#



mu2                 <- exp(eta2)
dermu2.dereta2      <- exp(eta2)
der2mu2.dereta2eta2 <- exp(eta2)

dersigma2.dersigma2.st  <- exp(sigma2.st)    
dersigma2.dersigma2.st2 <- exp(sigma2.st)                


                
pdf2          <- exp(-0.5 * log(2 * pi) - log(sigma) - (3/2) * log(y2) - ((y2 - mu2)^2)/(2 * sigma^2 * (mu2^2) * y2))
  
  
  
derpdf2.dermu2 <-  exp(-0.5 * log(2 * pi) - log(sigma) - (3/2) * log(y2) - ((y2 - 
    mu2)^2)/(2 * sigma^2 * (mu2^2) * y2)) * (2 * (y2 - mu2)/(2 * 
    sigma^2 * (mu2^2) * y2) + ((y2 - mu2)^2) * (2 * sigma^2 * 
    (2 * mu2) * y2)/(2 * sigma^2 * (mu2^2) * y2)^2)
         
           
   derpdf2.sigma2 <- -(exp(-0.5 * log(2 * pi) - log(sigma) - (3/2) * log(y2) - ((y2 - 
    mu2)^2)/(2 * sigma^2 * (mu2^2) * y2)) * (1/sigma - ((y2 - 
    mu2)^2) * (2 * (2 * sigma) * (mu2^2) * y2)/(2 * sigma^2 * 
    (mu2^2) * y2)^2)) 


der2pdf2.dermu2 <-  exp(-0.5 * log(2 * pi) - log(sigma) - (3/2) * log(y2) - ((y2 - 
    mu2)^2)/(2 * sigma^2 * (mu2^2) * y2)) * (2 * (y2 - mu2)/(2 * 
    sigma^2 * (mu2^2) * y2) + ((y2 - mu2)^2) * (2 * sigma^2 * 
    (2 * mu2) * y2)/(2 * sigma^2 * (mu2^2) * y2)^2) * (2 * (y2 - 
    mu2)/(2 * sigma^2 * (mu2^2) * y2) + ((y2 - mu2)^2) * (2 * 
    sigma^2 * (2 * mu2) * y2)/(2 * sigma^2 * (mu2^2) * y2)^2) + 
    exp(-0.5 * log(2 * pi) - log(sigma) - (3/2) * log(y2) - ((y2 - 
        mu2)^2)/(2 * sigma^2 * (mu2^2) * y2)) * ((((y2 - mu2)^2) * 
        (2 * sigma^2 * 2 * y2) - 2 * (y2 - mu2) * (2 * sigma^2 * 
        (2 * mu2) * y2))/(2 * sigma^2 * (mu2^2) * y2)^2 - ((y2 - 
        mu2)^2) * (2 * sigma^2 * (2 * mu2) * y2) * (2 * (2 * 
        sigma^2 * (2 * mu2) * y2 * (2 * sigma^2 * (mu2^2) * y2)))/((2 * 
        sigma^2 * (mu2^2) * y2)^2)^2 - (2/(2 * sigma^2 * (mu2^2) * 
        y2) + 2 * (y2 - mu2) * (2 * sigma^2 * (2 * mu2) * y2)/(2 * 
        sigma^2 * (mu2^2) * y2)^2))

   
der2pdf2.dersigma22 <- exp(-0.5 * log(2 * pi) - log(sigma) - (3/2) * log(y2) - ((y2 - 
    mu2)^2)/(2 * sigma^2 * (mu2^2) * y2)) * (1/sigma^2 + (((y2 - 
    mu2)^2) * (2 * 2 * (mu2^2) * y2)/(2 * sigma^2 * (mu2^2) * 
    y2)^2 - ((y2 - mu2)^2) * (2 * (2 * sigma) * (mu2^2) * y2) * 
    (2 * (2 * (2 * sigma) * (mu2^2) * y2 * (2 * sigma^2 * (mu2^2) * 
        y2)))/((2 * sigma^2 * (mu2^2) * y2)^2)^2)) + exp(-0.5 * 
    log(2 * pi) - log(sigma) - (3/2) * log(y2) - ((y2 - mu2)^2)/(2 * 
    sigma^2 * (mu2^2) * y2)) * (1/sigma - ((y2 - mu2)^2) * (2 * 
    (2 * sigma) * (mu2^2) * y2)/(2 * sigma^2 * (mu2^2) * y2)^2) * 
    (1/sigma - ((y2 - mu2)^2) * (2 * (2 * sigma) * (mu2^2) * 
        y2)/(2 * sigma^2 * (mu2^2) * y2)^2)
 

der2pdf2.mu2dersigma2 <- exp(-0.5 * log(2 * pi) - log(sigma) - (3/2) * log(y2) - ((y2 - 
    mu2)^2)/(2 * sigma^2 * (mu2^2) * y2)) * (((y2 - mu2)^2) * 
    (2 * (2 * sigma) * (2 * mu2) * y2)/(2 * sigma^2 * (mu2^2) * 
    y2)^2 - ((y2 - mu2)^2) * (2 * sigma^2 * (2 * mu2) * y2) * 
    (2 * (2 * (2 * sigma) * (mu2^2) * y2 * (2 * sigma^2 * (mu2^2) * 
        y2)))/((2 * sigma^2 * (mu2^2) * y2)^2)^2 - 2 * (y2 - 
    mu2) * (2 * (2 * sigma) * (mu2^2) * y2)/(2 * sigma^2 * (mu2^2) * 
    y2)^2) - exp(-0.5 * log(2 * pi) - log(sigma) - (3/2) * log(y2) - 
    ((y2 - mu2)^2)/(2 * sigma^2 * (mu2^2) * y2)) * (1/sigma - 
    ((y2 - mu2)^2) * (2 * (2 * sigma) * (mu2^2) * y2)/(2 * sigma^2 * 
        (mu2^2) * y2)^2) * (2 * (y2 - mu2)/(2 * sigma^2 * (mu2^2) * 
    y2) + ((y2 - mu2)^2) * (2 * sigma^2 * (2 * mu2) * y2)/(2 * 
    sigma^2 * (mu2^2) * y2)^2)





                  
if(naive == FALSE){                   
          

p2          <-  pnorm(((y2/mu2) - 1)/(sigma * sqrt(y2))) + exp(2/(mu2*sigma^2))* pnorm(-((y2/mu2) + 1)/(sigma * sqrt(y2)))
                                   
derp2.dermu2 <- exp(2/(mu2 * sigma^2)) * (y2 * dnorm(-((1 + y2/mu2)/(sigma * 
    sqrt(y2))))/(mu2^2 * sigma * sqrt(y2)) - 2 * (sigma^2 * pnorm(-((1 + 
    y2/mu2)/(sigma * sqrt(y2))))/(mu2 * sigma^2)^2)) - y2 * dnorm((y2/mu2 - 
    1)/(sigma * sqrt(y2)))/(mu2^2 * sigma * sqrt(y2))
      
derp2.dersigma2 <-  ((1 + y2/mu2) * dnorm(-((1 + y2/mu2)/(sigma * sqrt(y2)))) * sqrt(y2)/(sigma * 
    sqrt(y2))^2 - 4 * (mu2 * sigma * pnorm(-((1 + y2/mu2)/(sigma * 
    sqrt(y2))))/(mu2 * sigma^2)^2)) * exp(2/(mu2 * sigma^2)) - 
    dnorm((y2/mu2 - 1)/(sigma * sqrt(y2))) * sqrt(y2) * (y2/mu2 - 
        1)/(sigma * sqrt(y2))^2        
    

der2p2.dermu22 <-  exp(2/(mu2 * sigma^2)) * (y2 * ((1 + y2/mu2)/(mu2^4 * sigma^3 * 
    sqrt(y2)) - 2 * (mu2 * sigma * sqrt(y2)/(mu2^2 * sigma * 
    sqrt(y2))^2)) * dnorm(-((1 + y2/mu2)/(sigma * sqrt(y2)))) - 
    sigma * (2 * (sigma * (y2 * dnorm(-((1 + y2/mu2)/(sigma * 
        sqrt(y2))))/(mu2^2 * sigma * sqrt(y2)) - 2 * (sigma^2 * 
        pnorm(-((1 + y2/mu2)/(sigma * sqrt(y2))))/(mu2 * sigma^2)^2))) + 
        2 * (y2 * dnorm(-((1 + y2/mu2)/(sigma * sqrt(y2))))/(mu2^2 * 
            sqrt(y2)) - 2 * (mu2 * sigma^5 * pnorm(-((1 + y2/mu2)/(sigma * 
            sqrt(y2))))/(mu2 * sigma^2)^2)))/(mu2 * sigma^2)^2) - 
    y2 * ((y2/mu2 - 1)/(mu2^4 * sigma^3 * sqrt(y2)) - 2 * (mu2 * 
        sigma * sqrt(y2)/(mu2^2 * sigma * sqrt(y2))^2)) * dnorm((y2/mu2 - 
        1)/(sigma * sqrt(y2)))

  
der2p2.dersigma22 <-  (((1 + y2/mu2)^2/sigma - 2 * (sigma * y2)) * (1 + y2/mu2) * dnorm(-((1 + 
    y2/mu2)/(sigma * sqrt(y2)))) * sqrt(y2)/(sigma * sqrt(y2))^4 - 
    mu2 * (4 * ((1 - 4 * (mu2^2 * sigma^4/(mu2 * sigma^2)^2)) * 
        pnorm(-((1 + y2/mu2)/(sigma * sqrt(y2)))) + sigma * (1 + 
        y2/mu2) * dnorm(-((1 + y2/mu2)/(sigma * sqrt(y2)))) * 
        sqrt(y2)/(sigma * sqrt(y2))^2) + 4 * (sigma * ((1 + y2/mu2) * 
        dnorm(-((1 + y2/mu2)/(sigma * sqrt(y2)))) * sqrt(y2)/(sigma * 
        sqrt(y2))^2 - 4 * (mu2 * sigma * pnorm(-((1 + y2/mu2)/(sigma * 
        sqrt(y2))))/(mu2 * sigma^2)^2))))/(mu2 * sigma^2)^2) * 
    exp(2/(mu2 * sigma^2)) - ((y2/mu2 - 1)^2/sigma - 2 * (sigma * 
    y2)) * dnorm((y2/mu2 - 1)/(sigma * sqrt(y2))) * sqrt(y2) * 
    (y2/mu2 - 1)/(sigma * sqrt(y2))^4



der2p2.derdermu2sigma2 <- (((1 + y2/mu2)^2/sigma^2 - y2) * dnorm(-((1 + y2/mu2)/(sigma * 
    sqrt(y2)))) * sqrt(y2)/(mu2^2 * (sigma * sqrt(y2))^2) - (2 * 
    (sigma^2 * ((1 + y2/mu2) * dnorm(-((1 + y2/mu2)/(sigma * 
        sqrt(y2)))) * sqrt(y2)/(sigma * sqrt(y2))^2 - 4 * (mu2 * 
        sigma * pnorm(-((1 + y2/mu2)/(sigma * sqrt(y2))))/(mu2 * 
        sigma^2)^2))) + 4 * (sigma * (1 - 2 * (mu2^2 * sigma^4/(mu2 * 
    sigma^2)^2)) * pnorm(-((1 + y2/mu2)/(sigma * sqrt(y2)))) + 
    y2 * dnorm(-((1 + y2/mu2)/(sigma * sqrt(y2))))/(mu2 * sqrt(y2))))/(mu2 * 
    sigma^2)^2) * exp(2/(mu2 * sigma^2)) - ((y2/mu2 - 1)^2/sigma^2 - 
    y2) * dnorm((y2/mu2 - 1)/(sigma * sqrt(y2))) * sqrt(y2)/(mu2^2 * 
    (sigma * sqrt(y2))^2)


                                                }
                                                
                                                
}


####

if(margin2 == "LO"){

mu2 <- eta2
dermu2.dereta2 <- 1
der2mu2.dereta2eta2 <- 0 

dersigma2.dersigma2.st  <- exp(sigma2.st) 
dersigma2.dersigma2.st2 <- exp(sigma2.st)  



  pdf2          <- 1/sigma*( exp( -(y2 - mu2)/sigma ) )*( 1 + exp( -(y2 - mu2)/sigma )  )^-2 # dlogis(y2, eta2, sigma) 
  
derpdf2.dermu2 <-  (1 - 2 * (exp(-((y2 - mu2)/sigma))/(1 + exp(-((y2 - mu2)/sigma))))) * 
    exp(-((y2 - mu2)/sigma))/(sigma^2 * (1 + exp(-((y2 - mu2)/sigma)))^2)
  
derpdf2.sigma2 <- ((1 - 2 * (exp(-((y2 - mu2)/sigma))/(1 + exp(-((y2 - mu2)/sigma))))) * 
    (y2 - mu2)/sigma - 1) * exp(-((y2 - mu2)/sigma))/(sigma^2 * 
    (1 + exp(-((y2 - mu2)/sigma)))^2)
  
 
 der2pdf2.dermu2 <- ((1 - (2 + 2 * (1 - exp(-((y2 - mu2)/sigma))/(1 + exp(-((y2 - 
    mu2)/sigma))))) * exp(-((y2 - mu2)/sigma))/(1 + exp(-((y2 - 
    mu2)/sigma))))/(sigma^3 * (1 + exp(-((y2 - mu2)/sigma)))^2) - 
    2 * (sigma * (1 - 2 * (exp(-((y2 - mu2)/sigma))/(1 + exp(-((y2 - 
        mu2)/sigma))))) * (1 + exp(-((y2 - mu2)/sigma))) * exp(-((y2 - 
        mu2)/sigma))/(sigma^2 * (1 + exp(-((y2 - mu2)/sigma)))^2)^2)) * 
    exp(-((y2 - mu2)/sigma))
  
  
der2pdf2.dersigma22 <- (((1 - 2 * (exp(-((y2 - mu2)/sigma))/(1 + exp(-((y2 - mu2)/sigma))))) * 
    (y2 - mu2)/sigma - ((2 * ((1 - exp(-((y2 - mu2)/sigma))/(1 + 
    exp(-((y2 - mu2)/sigma)))) * (y2 - mu2)/sigma) - 2) * exp(-((y2 - 
    mu2)/sigma))/(1 + exp(-((y2 - mu2)/sigma))) + 2)) * (y2 - 
    mu2)/(sigma^4 * (1 + exp(-((y2 - mu2)/sigma)))^2) - ((1 - 
    2 * (exp(-((y2 - mu2)/sigma))/(1 + exp(-((y2 - mu2)/sigma))))) * 
    (y2 - mu2)/sigma - 1) * (1 + exp(-((y2 - mu2)/sigma))) * 
    (2 * (exp(-((y2 - mu2)/sigma)) * (y2 - mu2)) + 2 * (sigma * 
        (1 + exp(-((y2 - mu2)/sigma)))))/(sigma^2 * (1 + exp(-((y2 - 
    mu2)/sigma)))^2)^2) * exp(-((y2 - mu2)/sigma))
  
  
der2pdf2.mu2dersigma2 <-  ((1 - (2 + 2 * (1 - exp(-((y2 - mu2)/sigma))/(1 + exp(-((y2 - 
    mu2)/sigma))))) * exp(-((y2 - mu2)/sigma))/(1 + exp(-((y2 - 
    mu2)/sigma)))) * (y2 - mu2)/(sigma^4 * (1 + exp(-((y2 - mu2)/sigma)))^2) - 
    (1 - 2 * (exp(-((y2 - mu2)/sigma))/(1 + exp(-((y2 - mu2)/sigma))))) * 
        (1 + exp(-((y2 - mu2)/sigma))) * (2 * (exp(-((y2 - mu2)/sigma)) * 
        (y2 - mu2)) + 2 * (sigma * (1 + exp(-((y2 - mu2)/sigma)))))/(sigma^2 * 
        (1 + exp(-((y2 - mu2)/sigma)))^2)^2) * exp(-((y2 - mu2)/sigma))
  
  
  
  
  
if(naive == FALSE){   
  
    p2          <- ( 1 + exp( -(y2 - mu2)/sigma )  )^-1 #  plogis(y2, eta2, sigma) 
             

derp2.dermu2    <- -(exp(-((y2 - mu2)/sigma))/(sigma * (1 + exp(-((y2 - mu2)/sigma)))^2))
                    
         
derp2.dersigma2 <-   -(exp(-((y2 - mu2)/sigma)) * (y2 - mu2)/(sigma^2 * (1 + exp(-((y2 - 
    mu2)/sigma)))^2))                  
    


der2p2.dermu22 <- -((1/(sigma^2 * (1 + exp(-((y2 - mu2)/sigma)))^2) - 2 * ((1 + 
    exp(-((y2 - mu2)/sigma))) * exp(-((y2 - mu2)/sigma))/(sigma * 
    (1 + exp(-((y2 - mu2)/sigma)))^2)^2)) * exp(-((y2 - mu2)/sigma)))


der2p2.dersigma22 <-  -(((y2 - mu2)/(sigma^4 * (1 + exp(-((y2 - mu2)/sigma)))^2) - 
    (1 + exp(-((y2 - mu2)/sigma))) * (2 * (exp(-((y2 - mu2)/sigma)) * 
        (y2 - mu2)) + 2 * (sigma * (1 + exp(-((y2 - mu2)/sigma)))))/(sigma^2 * 
        (1 + exp(-((y2 - mu2)/sigma)))^2)^2) * exp(-((y2 - mu2)/sigma)) * 
    (y2 - mu2))

der2p2.derdermu2sigma2 <-  -((((y2 - mu2)/sigma - 1)/(sigma^2 * (1 + exp(-((y2 - mu2)/sigma)))^2) - 
    2 * (sigma * (1 + exp(-((y2 - mu2)/sigma))) * exp(-((y2 - 
        mu2)/sigma)) * (y2 - mu2)/(sigma^2 * (1 + exp(-((y2 - 
        mu2)/sigma)))^2)^2)) * exp(-((y2 - mu2)/sigma)))
            

                                         }

















}


##########################


if(margin2 == "rGU"){

mu2 <- eta2
dermu2.dereta2 <- 1
der2mu2.dereta2eta2 <- 0 

dersigma2.dersigma2.st  <- exp(sigma2.st) 
dersigma2.dersigma2.st2 <- exp(sigma2.st)    



  pdf2          <- 1/sigma*exp(-((y2-mu2)/sigma+exp(-((y2-mu2)/sigma))))
  
derpdf2.dermu2 <- -(exp(-((y2 - mu2)/sigma + exp(-((y2 - mu2)/sigma)))) * (exp(-((y2 - 
    mu2)/sigma)) - 1)/sigma^2) 
  
derpdf2.sigma2 <- -(((exp(-((y2 - mu2)/sigma)) - 1) * (y2 - mu2)/sigma + 1) * exp(-((y2 - 
    mu2)/sigma + exp(-((y2 - mu2)/sigma))))/sigma^2)


  
 der2pdf2.dermu2 <-  -(exp(-((y2 - mu2)/sigma + exp(-((y2 - mu2)/sigma)))) * (exp(-((y2 - 
    mu2)/sigma)) - (exp(-((y2 - mu2)/sigma)) - 1)^2)/sigma^3)
  
der2pdf2.dersigma22 <- -(((((y2 - mu2)/sigma - 1) * exp(-((y2 - mu2)/sigma)) + 1 - ((exp(-((y2 - 
    mu2)/sigma)) - 1) * (y2 - mu2)/sigma + 1) * (exp(-((y2 - 
    mu2)/sigma)) - 1)) * (y2 - mu2)/sigma - 2 * ((exp(-((y2 - 
    mu2)/sigma)) - 1) * (y2 - mu2)/sigma + 1)) * exp(-((y2 - 
    mu2)/sigma + exp(-((y2 - mu2)/sigma))))/sigma^3) 
  
  

der2pdf2.mu2dersigma2 <-  -(((exp(-((y2 - mu2)/sigma)) - (exp(-((y2 - mu2)/sigma)) - 1)^2) * 
    (y2 - mu2)/sigma - 2 * (exp(-((y2 - mu2)/sigma)) - 1)) * 
    exp(-((y2 - mu2)/sigma + exp(-((y2 - mu2)/sigma))))/sigma^3)
  
  
  
 
  
 
  
if(naive == FALSE){   
  
  
  
    p2          <- exp(-(exp(-(y2-mu2)/sigma)))
                

derp2.dermu2    <- -(exp(-(exp(-(y2 - mu2)/sigma))) * (exp(-(y2 - mu2)/sigma) * 
    (1/sigma)))
                      
                          
derp2.dersigma2 <-     -(exp(-((y2 - mu2)/sigma)) * exp(-exp(-((y2 - mu2)/sigma))) * 
    (y2 - mu2)/sigma^2)                 
    

der2p2.dermu22 <- -((1 - exp(-((y2 - mu2)/sigma))) * exp(-((y2 - mu2)/sigma)) * 
    exp(-exp(-((y2 - mu2)/sigma)))/sigma^2)


der2p2.dersigma22 <- -(((1 - exp(-((y2 - mu2)/sigma))) * (y2 - mu2)/sigma - 2) * exp(-((y2 - 
    mu2)/sigma)) * exp(-exp(-((y2 - mu2)/sigma))) * (y2 - mu2)/sigma^3)

der2p2.derdermu2sigma2 <- -(((1 - exp(-((y2 - mu2)/sigma))) * (y2 - mu2)/sigma - 1) * exp(-((y2 - 
    mu2)/sigma)) * exp(-exp(-((y2 - mu2)/sigma)))/sigma^2)
            

                                   }




}


######################################


if(margin2 == "GU"){

mu2                 <- eta2
dermu2.dereta2      <- 1
der2mu2.dereta2eta2 <- 0 

dersigma2.dersigma2.st  <- exp(sigma2.st) 
dersigma2.dersigma2.st2 <- exp(sigma2.st)  



  pdf2          <- exp(-exp((y2 - mu2)/sigma)) * (exp((y2 - mu2)/sigma) * (1/sigma))
                   
 derpdf2.dermu2 <-   exp(-exp((y2 - mu2)/sigma)) * exp((y2 - mu2)/sigma) * (exp((y2 - 
    mu2)/sigma) - 1)/sigma^2        
                   
derpdf2.sigma2 <-  ((exp((y2 - mu2)/sigma) - 1) * (y2 - mu2)/sigma - 1) * exp(-exp((y2 - 
    mu2)/sigma)) * exp((y2 - mu2)/sigma)/sigma^2

          
der2pdf2.dermu2 <-  ((exp((y2 - mu2)/sigma) - 1)^2 - exp((y2 - mu2)/sigma)) * exp(-exp((y2 - 
    mu2)/sigma)) * exp((y2 - mu2)/sigma)/sigma^3        
                   
der2pdf2.dersigma22 <- ((((exp((y2 - mu2)/sigma) - 2) * (y2 - mu2)/sigma - 2) * exp((y2 - 
    mu2)/sigma) + 2 - (exp((y2 - mu2)/sigma) - 1) * (y2 - mu2)/sigma) * 
    (y2 - mu2)/sigma - 2 * ((exp((y2 - mu2)/sigma) - 1) * (y2 - 
    mu2)/sigma - 1)) * exp(-exp((y2 - mu2)/sigma)) * exp((y2 - 
    mu2)/sigma)/sigma^3             
                   
                                     
der2pdf2.mu2dersigma2 <- (((exp((y2 - mu2)/sigma) - 1)^2 - exp((y2 - mu2)/sigma)) * (y2 - 
    mu2)/sigma - 2 * (exp((y2 - mu2)/sigma) - 1)) * exp(-exp((y2 - 
    mu2)/sigma)) * exp((y2 - mu2)/sigma)/sigma^3    
                   
            
            
            
                 
            
            
if(naive == FALSE){                    
                   
    p2          <- 1 - exp(-exp((y2 - mu2)/sigma))
    
 

derp2.dermu2    <- -(exp(-exp((y2 - mu2)/sigma)) * exp((y2 - mu2)/sigma)/sigma)
                  


                          
derp2.dersigma2 <-  -(exp(-exp((y2 - mu2)/sigma)) * exp((y2 - mu2)/sigma) * (y2 - 
    mu2)/sigma^2)                    
    


der2p2.dermu22 <-  -(exp(-exp((y2 - mu2)/sigma)) * exp((y2 - mu2)/sigma) * (exp((y2 - 
    mu2)/sigma) - 1)/sigma^2)


 der2p2.dersigma22 <- -(((exp((y2 - mu2)/sigma) - 1) * (y2 - mu2)/sigma - 2) * exp(-exp((y2 - 
    mu2)/sigma)) * exp((y2 - mu2)/sigma) * (y2 - mu2)/sigma^3)
 

der2p2.derdermu2sigma2 <-  -(((exp((y2 - mu2)/sigma) - 1) * (y2 - mu2)/sigma - 1) * exp(-exp((y2 - 
    mu2)/sigma)) * exp((y2 - mu2)/sigma)/sigma^2)

            

                                           }
                                           
                                           

}



if(margin2 %in% c("GA","GAi")){

#
sigma <- sigma2 <- ifelse(sigma < 0.006, 0.006, sigma) # related to gamma function
sigma2.st <- log(sigma) 
#



if(margin2 == "GAi"){eta2 <- ifelse(eta2 < 1e-05, 1e-05, eta2); mu2 <- eta2; dermu2.dereta2 <- 1; der2mu2.dereta2eta2 <- 0}

if(margin2 == "GA") { 
   mu2                 <- exp(eta2)
   dermu2.dereta2      <- exp(eta2)
   der2mu2.dereta2eta2 <- exp(eta2) 
                    }

dersigma2.dersigma2.st  <- exp(sigma2.st)   
dersigma2.dersigma2.st2 <- exp(sigma2.st)   




pdf2 <- dgamma(y2, shape = 1/sigma^2, scale = mu2 * sigma^2) # exp( (1/sigma^2) * log(y2/(mu2 * sigma^2)) - y2/(mu2 * sigma^2) - log(y2) - lgamma(1/sigma^2) ) 


derpdf2.dermu2 <- -(sigma^2 * exp((log(y2) - (2 * log(sigma) + log(mu2) + y2/mu2))/sigma^2 - 
    (lgamma(1/sigma^2) + log(y2))) * (mu2 - y2)/(mu2 * sigma^2)^2)  
    
 
derpdf2.sigma2 <- -(((2 * (log(y2) - (2 * log(sigma) + log(mu2))) - 2 * digamma(1/sigma^2))/sigma^3 + 
    mu2 * sigma * (2 * mu2 - 2 * y2)/(mu2 * sigma^2)^2) * exp((log(y2) - 
    (2 * log(sigma) + log(mu2) + y2/mu2))/sigma^2 - (lgamma(1/sigma^2) + 
    log(y2)))) 
 
 der2pdf2.dermu2 <- -(exp((log(y2) - (2 * log(sigma) + log(mu2) + y2/mu2))/sigma^2 - 
    (lgamma(1/sigma^2) + log(y2))) * (sigma^2 - ((1 - y2/mu2)/mu2 + 
    2 * (mu2 * sigma^6/(mu2 * sigma^2)^2)) * (mu2 - y2))/(mu2 * 
    sigma^2)^2) 
            
      
der2pdf2.dersigma22 <- -(exp((log(y2) - (2 * log(sigma) + log(mu2) + y2/mu2))/sigma^2 - 
    (lgamma(1/sigma^2) + log(y2))) * (mu2 * (1 - 4 * (mu2^2 * 
    sigma^4/(mu2 * sigma^2)^2)) * (2 * mu2 - 2 * y2)/(mu2 * sigma^2)^2 - 
    (((2 * (log(y2) - (2 * log(sigma) + log(mu2))) - 2 * digamma(1/sigma^2))/sigma^3 + 
        mu2 * sigma * (2 * mu2 - 2 * y2)/(mu2 * sigma^2)^2) * 
        (2 + 2 * (log(y2) - (2 * log(sigma) + log(mu2) + y2/mu2)) - 
            2 * digamma(1/sigma^2)) + (3 * (2 * (log(y2) - (2 * 
        log(sigma) + log(mu2))) - 2 * digamma(1/sigma^2)) + 4 - 
        4 * (trigamma(1/sigma^2)/sigma^2))/sigma)/sigma^3))  
                    
 der2pdf2.mu2dersigma2 <- -(exp((log(y2) - (2 * log(sigma) + log(mu2) + y2/mu2))/sigma^2 - 
    (lgamma(1/sigma^2) + log(y2))) * (mu2 - y2) * (sigma * (2 - 
    4 * (mu2^2 * sigma^4/(mu2 * sigma^2)^2)) - (2 + 2 * (log(y2) - 
    (2 * log(sigma) + log(mu2) + y2/mu2)) - 2 * digamma(1/sigma^2))/sigma)/(mu2 * 
    sigma^2)^2)







    
if(naive == FALSE){      
              
    p2  <-  pgamma(y2, shape = 1/sigma^2, scale = mu2 * sigma^2)
                                  

   derp2.dermu2 <- -((exp(-(y2/(mu2 *sigma^2)))* y2 *(y2/(mu2* sigma^2))^(-1 + 1/sigma^2))/(mu2^2 *sigma^2*gamma(1/sigma^2)))    # done in Mathematica
   
  
der2p2.dermu22 <- y2 * (2 * (y2/(mu2 * sigma^2))^(1/sigma^2 - 1) + y2 * 
    ((1/sigma^2 - 1) * (y2/(mu2 * sigma^2))^(1/sigma^2 - 2) - (y2/(mu2 * 
        sigma^2))^(1/sigma^2 - 1))/(mu2 * sigma^2)) * exp(-(y2/(mu2 * 
    sigma^2)))/(mu2^3 * sigma^2 * gamma(1/sigma^2))
   
      
der2p2.derdermu2sigma2 <- (y2 * (((1 + log(y2) - (log(mu2) + log(sigma^2)))/sigma^2 - 
    1) * (y2/(mu2 * sigma^2))^(1/sigma^2 - 1) + (y2/(mu2 * sigma^2))^(1/sigma^2 - 
    1) - (psigamma(1/sigma^2, 0) * (y2/(mu2 * sigma^2))^(1/sigma^2 - 
    1) + y2 * (y2/(mu2 * sigma^2))^(1/sigma^2 - 1)/mu2)/sigma^2) * 
    exp(-(y2/(mu2 * sigma^2)))/(mu2^2 * (sigma^2)^2 * gamma(1/sigma^2)))*2*sigma
   

funcD <- function(para) pgamma(y2, shape = 1/para^2, scale = mu2 * para^2)
 
nde <- numgh(funcD, sigma) 
 
derp2.dersigma2   <- nde$fi
der2p2.dersigma22 <- nde$se

                              }



}







##########################################################


if(margin2 %in% c(cont2par,cont3par)){ 

derpdf2.dereta2              <- derpdf2.dermu2*dermu2.dereta2       
der2pdf2.dereta2             <- der2pdf2.dermu2* dermu2.dereta2^2 + derpdf2.dermu2*der2mu2.dereta2eta2     

derpdf2.dersigma2.st         <- derpdf2.sigma2*dersigma2.dersigma2.st   
der2pdf2.dersigma2.st2       <- der2pdf2.dersigma22 * dersigma2.dersigma2.st^2 + derpdf2.sigma2  * dersigma2.dersigma2.st2 


der2pdf2.dereta2dersigma2    <- der2pdf2.mu2dersigma2* dermu2.dereta2
der2pdf2.dereta2dersigma2.st <- der2pdf2.dereta2dersigma2 *  dersigma2.dersigma2.st


if(margin2 %in% cont3par){


der2pdf2.dereta2dernu        <- der2pdf2.mu2dernu* dermu2.dereta2
der2pdf2.dereta2dernu.st     <- der2pdf2.dereta2dernu * dernu.dernu.st

der2pdf2.sigma2.st2dernu.st  <- der2pdf2.dersigma2dernu * dersigma2.dersigma2.st * dernu.dernu.st  

derpdf2.dernu.st             <- derpdf2.nu * dernu.dernu.st 
der2pdf2.dernu.st2           <- der2pdf2.dernu2 * dernu.dernu.st^2 +  derpdf2.nu  * dernu.dernu.st2 



                         }
                         
                         
                  
if(naive == FALSE){  



derp2.dereta2                <- derp2.dermu2*dermu2.dereta2
der2p2.dereta2eta2           <- der2p2.dermu22*dermu2.dereta2^2 + derp2.dermu2*der2mu2.dereta2eta2   

derp2.dersigma.st            <- derp2.dersigma2 *  dersigma2.dersigma2.st 
der2p2.dersigma2.st2         <- der2p2.dersigma22 * dersigma2.dersigma2.st^2 + derp2.dersigma2 * dersigma2.dersigma2.st2

der2p2.dereta2dersigma2      <- der2p2.derdermu2sigma2* dermu2.dereta2    
der2p2.dereta2dersigma2.st   <- der2p2.dereta2dersigma2 *  dersigma2.dersigma2.st  



if(margin2 %in% cont3par){

der2p2.dereta2dernu          <- der2p2.dermu2dernu* dermu2.dereta2 
der2p2.dereta2dernu.st       <- der2p2.dereta2dernu * dernu.dernu.st 
der2p2.dersigma2.stdernu.st  <- der2p2.dersigma2dernu * dersigma2.dersigma2.st * dernu.dernu.st  

derp2.nu.st                  <- derp2.dernu *  dernu.dernu.st 
der2p2.dernu.st2             <- der2p2.dernu2 * dernu.dernu.st^2 + derp2.dernu * dernu.dernu.st2

                         }


                 }
                 
                 
                 

}###############


pdf2 <- ifelse(pdf2 < min.dn, min.dn, pdf2)
p2   <- mm(p2, min.pr = min.pr, max.pr = max.pr) 


list(pdf2                         = ifef(pdf2),
     p2                           = ifef(p2), 
     derpdf2.dereta2              = ifef(derpdf2.dereta2), 
     derpdf2.dersigma2.st         = ifef(derpdf2.dersigma2.st), 
     derp2.dersigma.st            = ifef(derp2.dersigma.st),
     derp2.dereta2                = ifef(derp2.dereta2),
     der2p2.dereta2eta2           = ifef(der2p2.dereta2eta2), 
     der2pdf2.dereta2             = ifef(der2pdf2.dereta2),
     der2p2.dersigma2.st2         = ifef(der2p2.dersigma2.st2), 
     der2pdf2.dersigma2.st2       = ifef(der2pdf2.dersigma2.st2),
     der2p2.dereta2dersigma2.st   = ifef(der2p2.dereta2dersigma2.st),            
     der2pdf2.dereta2dersigma2.st = ifef(der2pdf2.dereta2dersigma2.st),
     der2pdf2.dereta2dernu.st     = ifef(der2pdf2.dereta2dernu.st),   
     der2pdf2.sigma2.st2dernu.st  = ifef(der2pdf2.sigma2.st2dernu.st),
     derpdf2.dernu.st             = ifef(derpdf2.dernu.st),           
     der2pdf2.dernu.st2           = ifef(der2pdf2.dernu.st2),         
     derp2.nu.st                  = ifef(derp2.nu.st),                
     der2p2.dernu.st2             = ifef(der2p2.dernu.st2),           
     der2p2.dereta2dernu.st       = ifef(der2p2.dereta2dernu.st),     
     der2p2.dersigma2.stdernu.st  = ifef(der2p2.dersigma2.stdernu.st), 
     indx = indx)     


}




    