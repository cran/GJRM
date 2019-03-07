distrHsDiscr <- function(y2, eta2, sigma2, sigma2.st, nu, nu.st, margin2, naive = FALSE, y2m){

p2 <- derp2.dersigma.st <- derp2.dereta2 <- der2p2.dereta2eta2 <- der2p2.dersigma2.st2 <- der2p2.dereta2dersigma2.st <- indx <- 1

der2pdf2.dereta2dernu.st    = 1
der2pdf2.sigma2.st2dernu.st = 1
derpdf2.dernu.st            = 1
der2pdf2.dernu.st2          = 1
derp2.nu.st                 = 1
der2p2.dernu.st2            = 1
der2p2.dereta2dernu.st      = 1
der2p2.dersigma2.stdernu.st = 1

cont1par <- c("PO","ZTP")
cont2par <- c("NBI","NBII","NBIa","NBIIa","PIG","PO","ZTP","DGP")
cont3par <- c("DEL","SICHEL")

# library(Deriv); library(numDeriv)
# expr <- expression(  )
# Simplify( D(D(expr, "mu2"),"sigma2") )
# func0 <- function(mu2){   }
# grad(func0 , mu2)
#func <- function(x) pNBI(y2, mu = x[1], sigma = sqrt(x[2])) 
#hessian(func, c(mu2,sigma2))
############################################################################
# remember that eta2 will have to disappear if we change default link on mu2
# this only applies to cases in which mu2 must be positive
# otherwise things are fine
############################################################################

#######################################################################
# Define generic numerical derivative functions
#######################################################################



if(margin2 %in% c("PIG","NBI","NBII","NBIa","NBIIa")){

derpdf2.dermu2FUNC2p <- function(func, y2, mu2, sigma2) numgh(func, mu2)
derpdf2.sigma2FUNC2p <- function(func, y2, mu2, sigma2) numgh(func, sigma2)  

der2pdf2.mu2dersigma2FUNC2p <- function(func, y2, mu2, sigma2) numch(func, mu2, sigma2)

}



#######################################################################

if(margin2 == "PO"){

mu2 <- exp(eta2)
dermu2.dereta2 <- exp(eta2)
der2mu2.dereta2eta2 <- exp(eta2)

dersigma2.dersigma2.st  <- 0  
dersigma2.dersigma2.st2 <- 0


if(max(y2) > 170){

prec <- pmax(53, getPrec(mu2), getPrec(y2))
        
mu2 <- mpfr(mu2, prec)
y2  <- mpfr( y2, prec)        
        
}        
        
        
# exp(y2*log(mu2) - mu2 - log(gamma(y2 + 1)))
# should be more stable but looks the same
# from a few experiments
        
        
pdf2 <-  as.numeric( (exp(-mu2)*mu2^y2)/factorial(y2) )   

derpdf2.dermu2FUNCpo <- function(y2, mu2) exp(-mu2) * (mu2^(y2 - 1) * y2 - mu2^y2)/factorial(y2) 
derpdf2.dermu2       <- as.numeric( derpdf2.dermu2FUNCpo(y2, mu2) )
    
der2pdf2.dermu2FUNCpo <- function(y2, mu2) exp(-mu2) * (mu2^y2 + y2 * (mu2^(y2 - 2) * (y2 - 1) - 2 * mu2^(y2 - 1)))/factorial(y2)  
der2pdf2.dermu2       <- as.numeric( der2pdf2.dermu2FUNCpo(y2, mu2) )
    
derpdf2.sigma2        <- 0
der2pdf2.dersigma22   <- 0
der2pdf2.mu2dersigma2 <- 0



if(naive == FALSE){   # needs y2m 
 
p2  <- pPO(as.numeric(y2), mu = as.numeric(mu2)) 
 
mu2 <- c(mu2)

derp2.dermu2           <- rowSums( matrix(as.numeric(derpdf2.dermu2FUNCpo(y2m, mu2)), dim(y2m)[1],dim(y2m)[2]), na.rm = TRUE ) 
derp2.dersigma2        <- 0
der2p2.dermu22         <- rowSums( matrix(as.numeric(der2pdf2.dermu2FUNCpo(y2m, mu2)),dim(y2m)[1],dim(y2m)[2]), na.rm = TRUE ) 
der2p2.dersigma22      <- 0
der2p2.derdermu2sigma2 <- 0


                      
    
                   }

}



if(margin2 == "ZTP"){

mu2 <- exp(eta2)
dermu2.dereta2 <- exp(eta2)
der2mu2.dereta2eta2 <- exp(eta2) 

dersigma2.dersigma2.st  <- 0
dersigma2.dersigma2.st2 <- 0

if(max(y2) > 170){

prec <- pmax(53, getPrec(mu2), getPrec(y2))
        
mu2 <- mpfr(mu2, prec)
y2  <- mpfr( y2, prec) 

}


pdf2FUNCztp <- function(y2, mu2) mu2^y2/(exp(mu2)-1)*1/factorial(y2)  
pdf2        <- as.numeric( pdf2FUNCztp(y2, mu2) )  

derpdf2.dermu2FUNCztp <- function(y2, mu2) (mu2^(y2 - 1) * y2 - mu2^y2 * exp(mu2)/(exp(mu2) - 1))/(factorial(y2) * (exp(mu2) - 1))  
derpdf2.dermu2        <- as.numeric( derpdf2.dermu2FUNCztp(y2, mu2) )
    
der2pdf2.dermu2FUNCztp <- function(y2, mu2) (y2 * (mu2^(y2 - 2) * (y2 - 1) - mu2^(y2 - 1) * exp(mu2)/(exp(mu2) - 
    1)) - exp(mu2) * (mu2^(y2 - 1) * y2 + mu2^y2 - 2 * (mu2^y2 * 
    exp(mu2)/(exp(mu2) - 1)))/(exp(mu2) - 1))/(factorial(y2) * (exp(mu2) - 
    1))   
der2pdf2.dermu2        <- as.numeric( der2pdf2.dermu2FUNCztp(y2, mu2) ) 
    
derpdf2.sigma2        <- 0
der2pdf2.dersigma22   <- 0
der2pdf2.mu2dersigma2 <- 0


if(naive == FALSE){   #needs y2m

mu2 <- c(mu2)

p2  <- rowSums( matrix(as.numeric(pdf2FUNCztp(y2m, mu2)),dim(y2m)[1],dim(y2m)[2]), na.rm = TRUE )

derp2.dermu2           <- rowSums( matrix(as.numeric(derpdf2.dermu2FUNCztp(y2m, mu2)), dim(y2m)[1],dim(y2m)[2]), na.rm = TRUE )
derp2.dersigma2        <- 0
der2p2.dermu22         <- rowSums( matrix(as.numeric(der2pdf2.dermu2FUNCztp(y2m, mu2)),dim(y2m)[1],dim(y2m)[2]), na.rm = TRUE )
der2p2.dersigma22      <- 0
der2p2.derdermu2sigma2 <- 0
                      
    
                   }

}




####














if(margin2 == "NBIa"){ # all analytical

sigma2    <- ifelse(sigma2 < 4.151334e-05, 4.151334e-05, sigma2) # related to gamma function
sigma2.st <- log(sigma2) 

mu2 <- exp(eta2)
dermu2.dereta2 <- exp(eta2)
der2mu2.dereta2eta2 <- exp(eta2) 

dersigma2.dersigma2.st  <- exp(sigma2.st)  
dersigma2.dersigma2.st2 <- exp(sigma2.st)


pdf2 <- dNBI(y2, mu = mu2, sigma = sqrt(sigma2))

derpdf2.dermu2FUNC <- function(y2, mu2, sigma2) (gamma(y2+1/sqrt(sigma2))/(gamma(1/sqrt(sigma2))*gamma(y2+1))*(sqrt(sigma2)*mu2/(1+sqrt(sigma2)*mu2))^y2*(1/(1+sqrt(sigma2)*mu2))^(1/sqrt(sigma2)))/mu2*(y2*(mu2*sqrt(sigma2))^(-1)*sqrt(sigma2)*mu2-(y2+1/sqrt(sigma2))*(mu2*sqrt(sigma2)+1)^(-1)*sqrt(sigma2)*mu2)    

derpdf2.dermu2     <- derpdf2.dermu2FUNC(y2, mu2, sigma2) 
    
derpdf2.sigma2FUNC <- function(y2, mu2, sigma2) 0.5*(1/sqrt(sigma2))*(gamma(y2+1/sqrt(sigma2))/(gamma(1/sqrt(sigma2))*gamma(y2+1))*(sqrt(sigma2)*mu2/(1+sqrt(sigma2)*mu2))^y2*(1/(1+sqrt(sigma2)*mu2))^(1/sqrt(sigma2)))*(digamma(y2+1/sqrt(sigma2))*(-1/sigma2)+y2*(mu2*sqrt(sigma2))^(-1)*mu2-(digamma(1/sqrt(sigma2))*(-1/sigma2)+(-1/sigma2)*log(mu2*sqrt(sigma2)+1)+(y2+1/sqrt(sigma2))*(1/(mu2*sqrt(sigma2)+1))*mu2))

derpdf2.sigma2     <- derpdf2.sigma2FUNC(y2, mu2, sigma2) 
     
der2pdf2.dermu2FUNC <- function(y2, mu2, sigma2) gamma(1/sqrt(sigma2) + y2) * (y2 * (((mu2 * sqrt(sigma2)/(1 + 
    mu2 * sqrt(sigma2)))^(y2 - 2) * (sqrt(sigma2) - mu2 * sigma2/(1 + 
    mu2 * sqrt(sigma2)))^2 * (y2 - 1) - sigma2 * (2 - 2 * (mu2 * 
    sqrt(sigma2)/(1 + mu2 * sqrt(sigma2)))) * (mu2 * sqrt(sigma2)/(1 + 
    mu2 * sqrt(sigma2)))^(y2 - 1)) * (1/(1 + mu2 * sqrt(sigma2)))^(1/sqrt(sigma2)) - 
    (1/(1 + mu2 * sqrt(sigma2)))^(1/sqrt(sigma2) - 1) * (mu2 * 
        sqrt(sigma2)/(1 + mu2 * sqrt(sigma2)))^(y2 - 1) * (sqrt(sigma2) - 
        mu2 * sigma2/(1 + mu2 * sqrt(sigma2)))/(1 + mu2 * sqrt(sigma2))) - 
    (y2 * (1/(1 + mu2 * sqrt(sigma2)))^(1/sqrt(sigma2) - 1) * 
        (mu2 * sqrt(sigma2)/(1 + mu2 * sqrt(sigma2)))^(y2 - 1) * 
        (sqrt(sigma2) - mu2 * sigma2/(1 + mu2 * sqrt(sigma2))) - 
        ((1/(1 + mu2 * sqrt(sigma2)))^(1/sqrt(sigma2) - 2) * 
            (1/sqrt(sigma2) - 1) * sqrt(sigma2)/(1 + mu2 * sqrt(sigma2)) + 
            2 * (sigma2 * (1/(1 + mu2 * sqrt(sigma2)))^(1/sqrt(sigma2) - 
                1)/sqrt(sigma2))) * (mu2 * sqrt(sigma2)/(1 + 
            mu2 * sqrt(sigma2)))^y2)/(1 + mu2 * sqrt(sigma2)))/((1 + 
    mu2 * sqrt(sigma2))^2 * gamma(1 + y2) * gamma(1/sqrt(sigma2)))
    
der2pdf2.dermu2     <- der2pdf2.dermu2FUNC(y2, mu2, sigma2) 
    
    
der2pdf2.dersigma22FUNC <- function(y2, mu2, sigma2) (gamma(y2 + 1/sqrt(sigma2))/(gamma(1/sqrt(sigma2)) * gamma(y2 + 
    1)) * ((sqrt(sigma2) * mu2/(1 + sqrt(sigma2) * mu2))^((y2 - 
    1) - 1) * ((y2 - 1) * (0.5 * sigma2^-0.5 * mu2/(1 + sqrt(sigma2) * 
    mu2) - sqrt(sigma2) * mu2 * (0.5 * sigma2^-0.5 * mu2)/(1 + 
    sqrt(sigma2) * mu2)^2)) * (y2 * (0.5 * sigma2^-0.5 * mu2/(1 + 
    sqrt(sigma2) * mu2) - sqrt(sigma2) * mu2 * (0.5 * sigma2^-0.5 * 
    mu2)/(1 + sqrt(sigma2) * mu2)^2)) + (sqrt(sigma2) * mu2/(1 + 
    sqrt(sigma2) * mu2))^(y2 - 1) * (y2 * (0.5 * (-0.5 * sigma2^-1.5) * 
    mu2/(1 + sqrt(sigma2) * mu2) - 0.5 * sigma2^-0.5 * mu2 * 
    (0.5 * sigma2^-0.5 * mu2)/(1 + sqrt(sigma2) * mu2)^2 - ((0.5 * 
    sigma2^-0.5 * mu2 * (0.5 * sigma2^-0.5 * mu2) + sqrt(sigma2) * 
    mu2 * (0.5 * (-0.5 * sigma2^-1.5) * mu2))/(1 + sqrt(sigma2) * 
    mu2)^2 - sqrt(sigma2) * mu2 * (0.5 * sigma2^-0.5 * mu2) * 
    (2 * (0.5 * sigma2^-0.5 * mu2 * (1 + sqrt(sigma2) * mu2)))/((1 + 
    sqrt(sigma2) * mu2)^2)^2)))) - (0.5 * sigma2^-0.5/sqrt(sigma2)^2 * 
    (gamma(y2 + 1/sqrt(sigma2)) * digamma(y2 + 1/sqrt(sigma2)))/(gamma(1/sqrt(sigma2)) * 
    gamma(y2 + 1)) - gamma(y2 + 1/sqrt(sigma2)) * (0.5 * sigma2^-0.5/sqrt(sigma2)^2 * 
    (gamma(1/sqrt(sigma2)) * digamma(1/sqrt(sigma2))) * gamma(y2 + 
    1))/(gamma(1/sqrt(sigma2)) * gamma(y2 + 1))^2) * ((sqrt(sigma2) * 
    mu2/(1 + sqrt(sigma2) * mu2))^(y2 - 1) * (y2 * (0.5 * sigma2^-0.5 * 
    mu2/(1 + sqrt(sigma2) * mu2) - sqrt(sigma2) * mu2 * (0.5 * 
    sigma2^-0.5 * mu2)/(1 + sqrt(sigma2) * mu2)^2))) - ((((0.5 * 
    (-0.5 * sigma2^-1.5)/sqrt(sigma2)^2 - 0.5 * sigma2^-0.5 * 
    (2 * (0.5 * sigma2^-0.5 * sqrt(sigma2)))/(sqrt(sigma2)^2)^2) * 
    (gamma(y2 + 1/sqrt(sigma2)) * digamma(y2 + 1/sqrt(sigma2))) - 
    0.5 * sigma2^-0.5/sqrt(sigma2)^2 * (gamma(y2 + 1/sqrt(sigma2)) * 
        (0.5 * sigma2^-0.5/sqrt(sigma2)^2 * trigamma(y2 + 1/sqrt(sigma2))) + 
        0.5 * sigma2^-0.5/sqrt(sigma2)^2 * (gamma(y2 + 1/sqrt(sigma2)) * 
            digamma(y2 + 1/sqrt(sigma2))) * digamma(y2 + 1/sqrt(sigma2))))/(gamma(1/sqrt(sigma2)) * 
    gamma(y2 + 1)) + 0.5 * sigma2^-0.5/sqrt(sigma2)^2 * (gamma(y2 + 
    1/sqrt(sigma2)) * digamma(y2 + 1/sqrt(sigma2))) * (0.5 * 
    sigma2^-0.5/sqrt(sigma2)^2 * (gamma(1/sqrt(sigma2)) * digamma(1/sqrt(sigma2))) * 
    gamma(y2 + 1))/(gamma(1/sqrt(sigma2)) * gamma(y2 + 1))^2 - 
    ((gamma(y2 + 1/sqrt(sigma2)) * (((0.5 * (-0.5 * sigma2^-1.5)/sqrt(sigma2)^2 - 
        0.5 * sigma2^-0.5 * (2 * (0.5 * sigma2^-0.5 * sqrt(sigma2)))/(sqrt(sigma2)^2)^2) * 
        (gamma(1/sqrt(sigma2)) * digamma(1/sqrt(sigma2))) - 0.5 * 
        sigma2^-0.5/sqrt(sigma2)^2 * (gamma(1/sqrt(sigma2)) * 
        (0.5 * sigma2^-0.5/sqrt(sigma2)^2 * trigamma(1/sqrt(sigma2))) + 
        0.5 * sigma2^-0.5/sqrt(sigma2)^2 * (gamma(1/sqrt(sigma2)) * 
            digamma(1/sqrt(sigma2))) * digamma(1/sqrt(sigma2)))) * 
        gamma(y2 + 1)) - 0.5 * sigma2^-0.5/sqrt(sigma2)^2 * (gamma(y2 + 
        1/sqrt(sigma2)) * digamma(y2 + 1/sqrt(sigma2))) * (0.5 * 
        sigma2^-0.5/sqrt(sigma2)^2 * (gamma(1/sqrt(sigma2)) * 
        digamma(1/sqrt(sigma2))) * gamma(y2 + 1)))/(gamma(1/sqrt(sigma2)) * 
        gamma(y2 + 1))^2 + gamma(y2 + 1/sqrt(sigma2)) * (0.5 * 
        sigma2^-0.5/sqrt(sigma2)^2 * (gamma(1/sqrt(sigma2)) * 
        digamma(1/sqrt(sigma2))) * gamma(y2 + 1)) * (2 * (0.5 * 
        sigma2^-0.5/sqrt(sigma2)^2 * (gamma(1/sqrt(sigma2)) * 
        digamma(1/sqrt(sigma2))) * gamma(y2 + 1) * (gamma(1/sqrt(sigma2)) * 
        gamma(y2 + 1))))/((gamma(1/sqrt(sigma2)) * gamma(y2 + 
        1))^2)^2)) * (sqrt(sigma2) * mu2/(1 + sqrt(sigma2) * 
    mu2))^y2 + (0.5 * sigma2^-0.5/sqrt(sigma2)^2 * (gamma(y2 + 
    1/sqrt(sigma2)) * digamma(y2 + 1/sqrt(sigma2)))/(gamma(1/sqrt(sigma2)) * 
    gamma(y2 + 1)) - gamma(y2 + 1/sqrt(sigma2)) * (0.5 * sigma2^-0.5/sqrt(sigma2)^2 * 
    (gamma(1/sqrt(sigma2)) * digamma(1/sqrt(sigma2))) * gamma(y2 + 
    1))/(gamma(1/sqrt(sigma2)) * gamma(y2 + 1))^2) * ((sqrt(sigma2) * 
    mu2/(1 + sqrt(sigma2) * mu2))^(y2 - 1) * (y2 * (0.5 * sigma2^-0.5 * 
    mu2/(1 + sqrt(sigma2) * mu2) - sqrt(sigma2) * mu2 * (0.5 * 
    sigma2^-0.5 * mu2)/(1 + sqrt(sigma2) * mu2)^2))))) * (1/(1 + 
    sqrt(sigma2) * mu2))^(1/sqrt(sigma2)) - (gamma(y2 + 1/sqrt(sigma2))/(gamma(1/sqrt(sigma2)) * 
    gamma(y2 + 1)) * ((sqrt(sigma2) * mu2/(1 + sqrt(sigma2) * 
    mu2))^(y2 - 1) * (y2 * (0.5 * sigma2^-0.5 * mu2/(1 + sqrt(sigma2) * 
    mu2) - sqrt(sigma2) * mu2 * (0.5 * sigma2^-0.5 * mu2)/(1 + 
    sqrt(sigma2) * mu2)^2))) - (0.5 * sigma2^-0.5/sqrt(sigma2)^2 * 
    (gamma(y2 + 1/sqrt(sigma2)) * digamma(y2 + 1/sqrt(sigma2)))/(gamma(1/sqrt(sigma2)) * 
    gamma(y2 + 1)) - gamma(y2 + 1/sqrt(sigma2)) * (0.5 * sigma2^-0.5/sqrt(sigma2)^2 * 
    (gamma(1/sqrt(sigma2)) * digamma(1/sqrt(sigma2))) * gamma(y2 + 
    1))/(gamma(1/sqrt(sigma2)) * gamma(y2 + 1))^2) * (sqrt(sigma2) * 
    mu2/(1 + sqrt(sigma2) * mu2))^y2) * ((1/(1 + sqrt(sigma2) * 
    mu2))^(1/sqrt(sigma2)) * (log((1/(1 + sqrt(sigma2) * mu2))) * 
    (0.5 * sigma2^-0.5/sqrt(sigma2)^2)) + (1/(1 + sqrt(sigma2) * 
    mu2))^((1/sqrt(sigma2)) - 1) * ((1/sqrt(sigma2)) * (0.5 * 
    sigma2^-0.5 * mu2/(1 + sqrt(sigma2) * mu2)^2))) - ((gamma(y2 + 
    1/sqrt(sigma2))/(gamma(1/sqrt(sigma2)) * gamma(y2 + 1)) * 
    ((sqrt(sigma2) * mu2/(1 + sqrt(sigma2) * mu2))^(y2 - 1) * 
        (y2 * (0.5 * sigma2^-0.5 * mu2/(1 + sqrt(sigma2) * mu2) - 
            sqrt(sigma2) * mu2 * (0.5 * sigma2^-0.5 * mu2)/(1 + 
                sqrt(sigma2) * mu2)^2))) - (0.5 * sigma2^-0.5/sqrt(sigma2)^2 * 
    (gamma(y2 + 1/sqrt(sigma2)) * digamma(y2 + 1/sqrt(sigma2)))/(gamma(1/sqrt(sigma2)) * 
    gamma(y2 + 1)) - gamma(y2 + 1/sqrt(sigma2)) * (0.5 * sigma2^-0.5/sqrt(sigma2)^2 * 
    (gamma(1/sqrt(sigma2)) * digamma(1/sqrt(sigma2))) * gamma(y2 + 
    1))/(gamma(1/sqrt(sigma2)) * gamma(y2 + 1))^2) * (sqrt(sigma2) * 
    mu2/(1 + sqrt(sigma2) * mu2))^y2) * ((1/(1 + sqrt(sigma2) * 
    mu2))^(1/sqrt(sigma2)) * (log((1/(1 + sqrt(sigma2) * mu2))) * 
    (0.5 * sigma2^-0.5/sqrt(sigma2)^2)) + (1/(1 + sqrt(sigma2) * 
    mu2))^((1/sqrt(sigma2)) - 1) * ((1/sqrt(sigma2)) * (0.5 * 
    sigma2^-0.5 * mu2/(1 + sqrt(sigma2) * mu2)^2))) + gamma(y2 + 
    1/sqrt(sigma2))/(gamma(1/sqrt(sigma2)) * gamma(y2 + 1)) * 
    (sqrt(sigma2) * mu2/(1 + sqrt(sigma2) * mu2))^y2 * ((1/(1 + 
    sqrt(sigma2) * mu2))^(1/sqrt(sigma2)) * (log((1/(1 + sqrt(sigma2) * 
    mu2))) * (0.5 * (-0.5 * sigma2^-1.5)/sqrt(sigma2)^2 - 0.5 * 
    sigma2^-0.5 * (2 * (0.5 * sigma2^-0.5 * sqrt(sigma2)))/(sqrt(sigma2)^2)^2) - 
    0.5 * sigma2^-0.5 * mu2/(1 + sqrt(sigma2) * mu2)^2/(1/(1 + 
        sqrt(sigma2) * mu2)) * (0.5 * sigma2^-0.5/sqrt(sigma2)^2)) - 
    ((1/(1 + sqrt(sigma2) * mu2))^(1/sqrt(sigma2)) * (log((1/(1 + 
        sqrt(sigma2) * mu2))) * (0.5 * sigma2^-0.5/sqrt(sigma2)^2)) + 
        (1/(1 + sqrt(sigma2) * mu2))^((1/sqrt(sigma2)) - 1) * 
            ((1/sqrt(sigma2)) * (0.5 * sigma2^-0.5 * mu2/(1 + 
                sqrt(sigma2) * mu2)^2))) * (log((1/(1 + sqrt(sigma2) * 
        mu2))) * (0.5 * sigma2^-0.5/sqrt(sigma2)^2)) + ((1/(1 + 
    sqrt(sigma2) * mu2))^((1/sqrt(sigma2)) - 1) * ((1/sqrt(sigma2)) * 
    (0.5 * (-0.5 * sigma2^-1.5) * mu2/(1 + sqrt(sigma2) * mu2)^2 - 
        0.5 * sigma2^-0.5 * mu2 * (2 * (0.5 * sigma2^-0.5 * mu2 * 
            (1 + sqrt(sigma2) * mu2)))/((1 + sqrt(sigma2) * mu2)^2)^2) - 
    0.5 * sigma2^-0.5/sqrt(sigma2)^2 * (0.5 * sigma2^-0.5 * mu2/(1 + 
        sqrt(sigma2) * mu2)^2)) - ((1/(1 + sqrt(sigma2) * mu2))^((1/sqrt(sigma2)) - 
    1) * (log((1/(1 + sqrt(sigma2) * mu2))) * (0.5 * sigma2^-0.5/sqrt(sigma2)^2)) + 
    (1/(1 + sqrt(sigma2) * mu2))^(((1/sqrt(sigma2)) - 1) - 1) * 
        (((1/sqrt(sigma2)) - 1) * (0.5 * sigma2^-0.5 * mu2/(1 + 
            sqrt(sigma2) * mu2)^2))) * ((1/sqrt(sigma2)) * (0.5 * 
    sigma2^-0.5 * mu2/(1 + sqrt(sigma2) * mu2)^2)))))
der2pdf2.dersigma22     <- der2pdf2.dersigma22FUNC(y2, mu2, sigma2) 
   
   
der2pdf2.mu2dersigma2FUNC <- function(y2, mu2, sigma2) gamma(1/sqrt(sigma2) + y2) * (y2 * ((((0.5/sqrt(sigma2) - mu2 * 
    (1.5 - mu2 * sigma2/((1 + mu2 * sqrt(sigma2)) * sqrt(sigma2)))/(1 + 
    mu2 * sqrt(sigma2))) * (mu2 * sqrt(sigma2)/(1 + mu2 * sqrt(sigma2)))^(y2 - 
    1) + mu2 * (0.5/sqrt(sigma2) - 0.5 * (mu2/(1 + mu2 * sqrt(sigma2)))) * 
    (mu2 * sqrt(sigma2)/(1 + mu2 * sqrt(sigma2)))^(y2 - 2) * 
    (sqrt(sigma2) - mu2 * sigma2/(1 + mu2 * sqrt(sigma2))) * 
    (y2 - 1)/(1 + mu2 * sqrt(sigma2)))/(gamma(1 + y2) * gamma(1/sqrt(sigma2))) - 
    (0.5 * (digamma(1/sqrt(sigma2) + y2)/(gamma(1 + y2) * gamma(1/sqrt(sigma2)))) - 
        0.5 * (digamma(1/sqrt(sigma2)) * gamma(1 + y2) * gamma(1/sqrt(sigma2))/(gamma(1 + 
            y2) * gamma(1/sqrt(sigma2)))^2)) * (mu2 * sqrt(sigma2)/(1 + 
        mu2 * sqrt(sigma2)))^(y2 - 1) * (sqrt(sigma2) - mu2 * 
        sigma2/(1 + mu2 * sqrt(sigma2)))/(sigma2 * sqrt(sigma2))) * 
    (1/(1 + mu2 * sqrt(sigma2)))^(1/sqrt(sigma2)) - (0.5 * (mu2 * 
    (1/(1 + mu2 * sqrt(sigma2)))^(1/sqrt(sigma2) - 1)/(1 + mu2 * 
    sqrt(sigma2))^2) - 0.5 * ((1/(1 + mu2 * sqrt(sigma2)))^(1/sqrt(sigma2)) * 
    log1p(mu2 * sqrt(sigma2))/sqrt(sigma2))) * (mu2 * sqrt(sigma2)/(1 + 
    mu2 * sqrt(sigma2)))^(y2 - 1) * (sqrt(sigma2) - mu2 * sigma2/(1 + 
    mu2 * sqrt(sigma2)))/(sigma2 * gamma(1 + y2) * gamma(1/sqrt(sigma2)))) - 
    ((((0.5/sqrt(sigma2) - mu2/(1 + mu2 * sqrt(sigma2)))/sqrt(sigma2) - 
        0.5/sigma2) * (1/(1 + mu2 * sqrt(sigma2)))^(1/sqrt(sigma2) - 
        1) - (0.5 * (mu2 * (1/(1 + mu2 * sqrt(sigma2)))^(1/sqrt(sigma2) - 
        2) * (1/sqrt(sigma2) - 1)/(1 + mu2 * sqrt(sigma2))^2) - 
        0.5 * ((1/(1 + mu2 * sqrt(sigma2)))^(1/sqrt(sigma2) - 
            1) * log1p(mu2 * sqrt(sigma2))/sigma2))/sqrt(sigma2)) * 
        (mu2 * sqrt(sigma2)/(1 + mu2 * sqrt(sigma2)))^y2/(gamma(1 + 
        y2) * gamma(1/sqrt(sigma2))) + (1/(1 + mu2 * sqrt(sigma2)))^(1/sqrt(sigma2) - 
        1) * (mu2 * y2 * (0.5/sqrt(sigma2) - 0.5 * (mu2/(1 + 
        mu2 * sqrt(sigma2)))) * (mu2 * sqrt(sigma2)/(1 + mu2 * 
        sqrt(sigma2)))^(y2 - 1)/((1 + mu2 * sqrt(sigma2)) * gamma(1 + 
        y2) * gamma(1/sqrt(sigma2))) - (0.5 * (digamma(1/sqrt(sigma2) + 
        y2)/(gamma(1 + y2) * gamma(1/sqrt(sigma2)))) - 0.5 * 
        (digamma(1/sqrt(sigma2)) * gamma(1 + y2) * gamma(1/sqrt(sigma2))/(gamma(1 + 
            y2) * gamma(1/sqrt(sigma2)))^2)) * (mu2 * sqrt(sigma2)/(1 + 
        mu2 * sqrt(sigma2)))^y2/(sigma2 * sqrt(sigma2))))/(1 + 
        mu2 * sqrt(sigma2)))/(1 + mu2 * sqrt(sigma2))  
der2pdf2.mu2dersigma2     <- der2pdf2.mu2dersigma2FUNC(y2, mu2, sigma2)

if(naive == FALSE){   
 
p2   <- pNBI(y2, mu = mu2, sigma = sqrt(sigma2))  

ly2 <- length(y2)
if(length(sigma2) == 1) sigma2 <- c(rep(sigma2, ly2))


mu2 <- c(mu2)
sigma2 <- c(sigma2)

derp2.dermu2           <- rowSums( derpdf2.dermu2FUNC(        y2m, mu2, sigma2) , na.rm = TRUE )
derp2.dersigma2        <- rowSums( derpdf2.sigma2FUNC(        y2m, mu2, sigma2) , na.rm = TRUE )
der2p2.dermu22         <- rowSums( der2pdf2.dermu2FUNC(       y2m, mu2, sigma2) , na.rm = TRUE )
der2p2.dersigma22      <- rowSums( der2pdf2.dersigma22FUNC(   y2m, mu2, sigma2) , na.rm = TRUE )
der2p2.derdermu2sigma2 <- rowSums( der2pdf2.mu2dersigma2FUNC( y2m, mu2, sigma2) , na.rm = TRUE )
                      
                        
                        
                   }

}




# seems that numerical is best (in terms of convergence and speed) when 
# sigma approaches 0
# followed by half-half and fully analytical
# full numerical does not give the most accurate results 
# when sigma is not zero, although it tends to be the fastest
# perhaps on occasion it may not convergence
# due to poor gradient and H approximation



# 11 oct 2018, decided to use all numerical for now...

if(margin2 == "NBI"){ # all Numerical - does not need y2m

sigma2    <- ifelse(sigma2 < 4.151334e-06, 4.151334e-06, sigma2) # related to gamma function
sigma2.st <- log(sigma2) 

mu2 <- exp(eta2)
dermu2.dereta2 <- exp(eta2)
der2mu2.dereta2eta2 <- exp(eta2) 

dersigma2.dersigma2.st  <- exp(sigma2.st)  
dersigma2.dersigma2.st2 <- exp(sigma2.st)  


pdf2 <- dNBI(y2, mu = mu2, sigma = sqrt(sigma2))  

derpdf2.dermu2F <- derpdf2.dermu2FUNC2p(function(mu2) dNBI(y2, mu = mu2, sigma = sqrt(sigma2)), y2, mu2, sigma2) 
derpdf2.dermu2  <- derpdf2.dermu2F$fi 
der2pdf2.dermu2 <- derpdf2.dermu2F$se  


derpdf2.sigma2F     <- derpdf2.sigma2FUNC2p(function(sigma2) dNBI(y2, mu = mu2, sigma = sqrt(sigma2)), y2, mu2, sigma2) 
derpdf2.sigma2      <- derpdf2.sigma2F$fi      
der2pdf2.dersigma22 <- derpdf2.sigma2F$se    


der2pdf2.mu2dersigma2 <- der2pdf2.mu2dersigma2FUNC2p(function(mu2, sigma2) dNBI(y2, mu = mu2, sigma = sqrt(sigma2)), y2, mu2, sigma2)


if(naive == FALSE){   
 
p2  <- pNBI(y2, mu = mu2, sigma = sqrt(sigma2))  
 
derp2.dermu2F  <- derpdf2.dermu2FUNC2p(function(mu2) pNBI(y2, mu = mu2, sigma = sqrt(sigma2)), y2, mu2, sigma2)
derp2.dermu2   <- derp2.dermu2F$fi
der2p2.dermu22 <- derp2.dermu2F$se

derp2.dersigma2F <- derpdf2.sigma2FUNC2p(function(sigma2) pNBI(y2, mu = mu2, sigma = sqrt(sigma2)), y2, mu2, sigma2) 
derp2.dersigma2  <- derp2.dersigma2F$fi 
der2p2.dersigma22<- derp2.dersigma2F$se 

der2p2.derdermu2sigma2 <- der2pdf2.mu2dersigma2FUNC2p(function(mu2, sigma2) pNBI(y2, mu = mu2, sigma = sqrt(sigma2)), y2, mu2, sigma2) 

  
                   }

}





if(margin2 == "NBIhh"){ # half analy half numerical - we need y2m, was deemed best choice up until 03/05/2017
                      # in pahse of robust testing the presence of the gamma function created problems
                      # so using all numerical maybe is better but there is no general rule, we will leave the default for now 

sigma2    <- ifelse(sigma2 < 4.151334e-05, 4.151334e-05, sigma2) # related to gamma function
sigma2.st <- log(sigma2) 

mu2 <- exp(eta2)
dermu2.dereta2 <- exp(eta2)
der2mu2.dereta2eta2 <- exp(eta2) 

dersigma2.dersigma2.st  <- exp(sigma2.st)
dersigma2.dersigma2.st2 <- exp(sigma2.st)  



pdf2 <- dNBI(y2, mu = mu2, sigma = sqrt(sigma2))

derpdf2.dermu2FUNC <- function(y2, mu2, sigma2) (gamma(y2+1/sqrt(sigma2))/(gamma(1/sqrt(sigma2))*gamma(y2+1))*(sqrt(sigma2)*mu2/(1+sqrt(sigma2)*mu2))^y2*(1/(1+sqrt(sigma2)*mu2))^(1/sqrt(sigma2)))/mu2*(y2*(mu2*sqrt(sigma2))^(-1)*sqrt(sigma2)*mu2-(y2+1/sqrt(sigma2))*(mu2*sqrt(sigma2)+1)^(-1)*sqrt(sigma2)*mu2)

derpdf2.dermu2     <- derpdf2.dermu2FUNC(y2, mu2, sigma2) 
    
derpdf2.sigma2FUNC <- function(y2, mu2, sigma2) 0.5*(1/sqrt(sigma2))*(gamma(y2+1/sqrt(sigma2))/(gamma(1/sqrt(sigma2))*gamma(y2+1))*(sqrt(sigma2)*mu2/(1+sqrt(sigma2)*mu2))^y2*(1/(1+sqrt(sigma2)*mu2))^(1/sqrt(sigma2)))*(digamma(y2+1/sqrt(sigma2))*(-1/sigma2)+y2*(mu2*sqrt(sigma2))^(-1)*mu2-(digamma(1/sqrt(sigma2))*(-1/sigma2)+(-1/sigma2)*log(mu2*sqrt(sigma2)+1)+(y2+1/sqrt(sigma2))*(1/(mu2*sqrt(sigma2)+1))*mu2))

derpdf2.sigma2     <- derpdf2.sigma2FUNC(y2, mu2, sigma2)  
                                      
der2pdf2.dermu2       <- derpdf2.dermu2FUNC2p(function(mu2)    derpdf2.dermu2FUNC(y2, mu2, sigma2), y2, mu2, sigma2)$fi  
der2pdf2.dersigma22   <- derpdf2.sigma2FUNC2p(function(sigma2) derpdf2.sigma2FUNC(y2, mu2, sigma2), y2, mu2, sigma2)$fi  
der2pdf2.mu2dersigma2 <- derpdf2.sigma2FUNC2p(function(sigma2) derpdf2.dermu2FUNC(y2, mu2, sigma2), y2, mu2, sigma2)$fi       
      
  
       
if(naive == FALSE){   
 
p2  <- pNBI(y2, mu = mu2, sigma = sqrt(sigma2))  
 
ly2 <- length(y2)
if(length(sigma2) == 1) sigma2 <- sigma2 <- c(rep(sigma2, ly2))

mu2 <- c(mu2)
sigma2 <- c(sigma2)

derp2.dermu2           <- rowSums( derpdf2.dermu2FUNC(y2m, mu2, sigma2 ) , na.rm = TRUE )
derp2.dersigma2        <- rowSums( derpdf2.sigma2FUNC(y2m, mu2, sigma2 ) , na.rm = TRUE )


der2p2.dermu22         <- derpdf2.dermu2FUNC2p(function(mu2)    rowSums( derpdf2.dermu2FUNC( y2m, mu2, sigma2 ) , na.rm = TRUE ), y2, mu2, sigma2)$fi
der2p2.dersigma22      <- derpdf2.sigma2FUNC2p(function(sigma2) rowSums( derpdf2.sigma2FUNC( y2m, mu2, sigma2 ) , na.rm = TRUE ), y2, mu2, sigma2)$fi  
der2p2.derdermu2sigma2 <- derpdf2.sigma2FUNC2p(function(sigma2) rowSums( derpdf2.dermu2FUNC( y2m, mu2, sigma2 ) , na.rm = TRUE ), y2, mu2, sigma2)$fi




                   }
}




#######





if(margin2 == "NBII"){ # all Numerical - does not need y2m

sigma2    <- ifelse(sigma2 < 4.151334e-06, 4.151334e-06, sigma2) # related to gamma function
sigma2.st <- log(sigma2) 

mu2 <- exp(eta2)
dermu2.dereta2 <- exp(eta2)
der2mu2.dereta2eta2 <- exp(eta2) 

dersigma2.dersigma2.st  <- exp(sigma2.st)  
dersigma2.dersigma2.st2 <- exp(sigma2.st)  


pdf2 <- dNBII(y2, mu = mu2, sigma = sqrt(sigma2))  

derpdf2.dermu2F <- derpdf2.dermu2FUNC2p(function(mu2) dNBII(y2, mu = mu2, sigma = sqrt(sigma2)), y2, mu2, sigma2) 
derpdf2.dermu2  <- derpdf2.dermu2F$fi 
der2pdf2.dermu2 <- derpdf2.dermu2F$se  


derpdf2.sigma2F     <- derpdf2.sigma2FUNC2p(function(sigma2) dNBII(y2, mu = mu2, sigma = sqrt(sigma2)), y2, mu2, sigma2) 
derpdf2.sigma2      <- derpdf2.sigma2F$fi      
der2pdf2.dersigma22 <- derpdf2.sigma2F$se    


der2pdf2.mu2dersigma2 <- der2pdf2.mu2dersigma2FUNC2p(function(mu2, sigma2) dNBII(y2, mu = mu2, sigma = sqrt(sigma2)), y2, mu2, sigma2)


if(naive == FALSE){   
 
p2  <- pNBII(y2, mu = mu2, sigma = sqrt(sigma2))  
 
derp2.dermu2F  <- derpdf2.dermu2FUNC2p(function(mu2) pNBII(y2, mu = mu2, sigma = sqrt(sigma2)), y2, mu2, sigma2)
derp2.dermu2   <- derp2.dermu2F$fi
der2p2.dermu22 <- derp2.dermu2F$se

derp2.dersigma2F <- derpdf2.sigma2FUNC2p(function(sigma2) pNBII(y2, mu = mu2, sigma = sqrt(sigma2)), y2, mu2, sigma2) 
derp2.dersigma2  <- derp2.dersigma2F$fi 
der2p2.dersigma22<- derp2.dersigma2F$se 

der2p2.derdermu2sigma2 <- der2pdf2.mu2dersigma2FUNC2p(function(mu2, sigma2) pNBII(y2, mu = mu2, sigma = sqrt(sigma2)), y2, mu2, sigma2) 

  
                   }

}



















if(margin2 == "NBIIhahn"){ # half analy half numer - needs y2m

sigma2    <- ifelse(sigma2 < 4.151334e-05, 4.151334e-05, sigma2) # related to gamma function
sigma2.st <- log(sigma2) 

mu2 <- exp(eta2)
dermu2.dereta2 <- exp(eta2)
der2mu2.dereta2eta2 <- exp(eta2) 

dersigma2.dersigma2.st  <- exp(sigma2.st) 
dersigma2.dersigma2.st2 <- exp(sigma2.st)  


pdf2 <- dNBII(y2, mu = mu2, sigma = sqrt(sigma2))

derpdf2.dermu2FUNC <- function(y2, mu2, sigma2) gamma(mu2/sqrt(sigma2) + y2) * (sigma2^(y2/2) * digamma(mu2/sqrt(sigma2) + 
    y2)/((1 + sqrt(sigma2))^(mu2/sqrt(sigma2) + y2) * gamma(1 + 
    y2) * gamma(mu2/sqrt(sigma2))) - sigma2^(y2/2) * ((1 + sqrt(sigma2))^(mu2/sqrt(sigma2) + 
    y2) * digamma(mu2/sqrt(sigma2)) + (1 + sqrt(sigma2))^(mu2/sqrt(sigma2) + 
    y2) * log1p(sqrt(sigma2))) * gamma(1 + y2) * gamma(mu2/sqrt(sigma2))/((1 + 
    sqrt(sigma2))^(mu2/sqrt(sigma2) + y2) * gamma(1 + y2) * gamma(mu2/sqrt(sigma2)))^2)/sqrt(sigma2) 

derpdf2.dermu2     <- derpdf2.dermu2FUNC(y2, mu2, sigma2) 
    
derpdf2.sigma2FUNC <- function(y2, mu2, sigma2) ((0.5 * (sigma2^((y2 - 1)/2) * y2) - 0.5 * (mu2 * sigma2^(y2/2 - 
    1) * digamma(mu2/sqrt(sigma2) + y2)))/((1 + sqrt(sigma2))^(mu2/sqrt(sigma2) + 
    y2) * gamma(1 + y2) * gamma(mu2/sqrt(sigma2))) - sigma2^(y2/2) * 
    (0.5 * ((1 + sqrt(sigma2))^(mu2/sqrt(sigma2) + y2 - 1) * 
        (mu2/sqrt(sigma2) + y2)) - mu2 * (0.5 * ((1 + sqrt(sigma2))^(mu2/sqrt(sigma2) + 
        y2) * digamma(mu2/sqrt(sigma2))) + 0.5 * ((1 + sqrt(sigma2))^(mu2/sqrt(sigma2) + 
        y2) * log1p(sqrt(sigma2))))/sigma2) * gamma(1 + y2) * 
    gamma(mu2/sqrt(sigma2))/((1 + sqrt(sigma2))^(mu2/sqrt(sigma2) + 
    y2) * gamma(1 + y2) * gamma(mu2/sqrt(sigma2)))^2) * gamma(mu2/sqrt(sigma2) + 
    y2)/sqrt(sigma2) 
    
derpdf2.sigma2     <- derpdf2.sigma2FUNC(y2, mu2, sigma2)  
                
                          
der2pdf2.dermu2       <- derpdf2.dermu2FUNC2p(function(mu2)    derpdf2.dermu2FUNC(y2, mu2, sigma2), y2, mu2, sigma2)$fi  
der2pdf2.dersigma22   <- derpdf2.sigma2FUNC2p(function(sigma2) derpdf2.sigma2FUNC(y2, mu2, sigma2), y2, mu2, sigma2)$fi  
der2pdf2.mu2dersigma2 <- derpdf2.sigma2FUNC2p(function(sigma2) derpdf2.dermu2FUNC(y2, mu2, sigma2), y2, mu2, sigma2)$fi       
      
  
       
if(naive == FALSE){   
 
p2  <- pNBII(y2, mu = mu2, sigma = sqrt(sigma2))  
 
ly2 <- length(y2)
if(length(sigma2) == 1) sigma2 <- sigma2 <- c(rep(sigma2, ly2))

mu2 <- c(mu2)
sigma2 <- c(sigma2)

derp2.dermu2           <- rowSums( derpdf2.dermu2FUNC(y2m, mu2, sigma2 ) , na.rm = TRUE )
derp2.dersigma2        <- rowSums( derpdf2.sigma2FUNC(y2m, mu2, sigma2 ) , na.rm = TRUE )


der2p2.dermu22         <- derpdf2.dermu2FUNC2p(function(mu2)    rowSums( derpdf2.dermu2FUNC( y2m, mu2, sigma2 ) , na.rm = TRUE ), y2, mu2, sigma2)$fi
der2p2.dersigma22      <- derpdf2.sigma2FUNC2p(function(sigma2) rowSums( derpdf2.sigma2FUNC( y2m, mu2, sigma2 ) , na.rm = TRUE ), y2, mu2, sigma2)$fi  
der2p2.derdermu2sigma2 <- derpdf2.sigma2FUNC2p(function(sigma2) rowSums( derpdf2.dermu2FUNC( y2m, mu2, sigma2 ) , na.rm = TRUE ), y2, mu2, sigma2)$fi




                   }
}














if(margin2 == "NBIIa"){ # all analytical - needs y2m

sigma2    <- ifelse(sigma2 < 4.151334e-05, 4.151334e-05, sigma2) # related to gamma function
sigma2.st <- log(sigma2) 

mu2 <- exp(eta2) 
dermu2.dereta2 <- exp(eta2)  
der2mu2.dereta2eta2 <- exp(eta2) 

dersigma2.dersigma2.st  <- exp(sigma2.st)
dersigma2.dersigma2.st2 <- exp(sigma2.st)


pdf2 <- (gamma(y2 + mu2/sqrt(sigma2))*sqrt(sigma2)^y2)/(gamma(mu2/sqrt(sigma2))*gamma(y2+1)*(1+sqrt(sigma2))^(y2+mu2/sqrt(sigma2)))    

derpdf2.dermu2FUNC <- function(y2, mu2, sigma2) gamma(mu2/sqrt(sigma2) + y2) * (sigma2^(y2/2) * digamma(mu2/sqrt(sigma2) + 
    y2)/((1 + sqrt(sigma2))^(mu2/sqrt(sigma2) + y2) * gamma(1 + 
    y2) * gamma(mu2/sqrt(sigma2))) - sigma2^(y2/2) * ((1 + sqrt(sigma2))^(mu2/sqrt(sigma2) + 
    y2) * digamma(mu2/sqrt(sigma2)) + (1 + sqrt(sigma2))^(mu2/sqrt(sigma2) + 
    y2) * log1p(sqrt(sigma2))) * gamma(1 + y2) * gamma(mu2/sqrt(sigma2))/((1 + 
    sqrt(sigma2))^(mu2/sqrt(sigma2) + y2) * gamma(1 + y2) * gamma(mu2/sqrt(sigma2)))^2)/sqrt(sigma2) 

derpdf2.dermu2     <- derpdf2.dermu2FUNC(y2, mu2, sigma2) 
    
    
derpdf2.sigma2FUNC <- function(y2, mu2, sigma2) ((0.5 * (sigma2^((y2 - 1)/2) * y2) - 0.5 * (mu2 * sigma2^(y2/2 - 
    1) * digamma(mu2/sqrt(sigma2) + y2)))/((1 + sqrt(sigma2))^(mu2/sqrt(sigma2) + 
    y2) * gamma(1 + y2) * gamma(mu2/sqrt(sigma2))) - sigma2^(y2/2) * 
    (0.5 * ((1 + sqrt(sigma2))^(mu2/sqrt(sigma2) + y2 - 1) * 
        (mu2/sqrt(sigma2) + y2)) - mu2 * (0.5 * ((1 + sqrt(sigma2))^(mu2/sqrt(sigma2) + 
        y2) * digamma(mu2/sqrt(sigma2))) + 0.5 * ((1 + sqrt(sigma2))^(mu2/sqrt(sigma2) + 
        y2) * log1p(sqrt(sigma2))))/sigma2) * gamma(1 + y2) * 
    gamma(mu2/sqrt(sigma2))/((1 + sqrt(sigma2))^(mu2/sqrt(sigma2) + 
    y2) * gamma(1 + y2) * gamma(mu2/sqrt(sigma2)))^2) * gamma(mu2/sqrt(sigma2) + 
    y2)/sqrt(sigma2) 

derpdf2.sigma2     <- derpdf2.sigma2FUNC(y2, mu2, sigma2) 
     
     
der2pdf2.dermu2FUNC <- function(y2, mu2, sigma2) gamma(mu2/sqrt(sigma2) + y2) * (sigma2^(y2/2 - 1) * (digamma(mu2/sqrt(sigma2) + 
    y2)^2 + trigamma(mu2/sqrt(sigma2) + y2))/((1 + sqrt(sigma2))^(mu2/sqrt(sigma2) + 
    y2) * gamma(1 + y2) * gamma(mu2/sqrt(sigma2))) - (((1 + sqrt(sigma2))^(mu2/sqrt(sigma2) + 
    y2) * digamma(mu2/sqrt(sigma2)) + (1 + sqrt(sigma2))^(mu2/sqrt(sigma2) + 
    y2) * log1p(sqrt(sigma2))) * (2 * (sigma2^(y2/2 - 1) * digamma(mu2/sqrt(sigma2) + 
    y2)) - 2 * (sigma2^(y2/2 - 1) * ((1 + sqrt(sigma2))^(mu2/sqrt(sigma2) + 
    y2) * digamma(mu2/sqrt(sigma2)) + (1 + sqrt(sigma2))^(mu2/sqrt(sigma2) + 
    y2) * log1p(sqrt(sigma2))) * (1 + sqrt(sigma2))^(mu2/sqrt(sigma2) + 
    y2) * gamma(1 + y2)^2 * gamma(mu2/sqrt(sigma2))^2/((1 + sqrt(sigma2))^(mu2/sqrt(sigma2) + 
    y2) * gamma(1 + y2) * gamma(mu2/sqrt(sigma2)))^2)) + sigma2^(y2/2 - 
    1) * (((1 + sqrt(sigma2))^(mu2/sqrt(sigma2) + y2) * log1p(sqrt(sigma2)) + 
    2 * ((1 + sqrt(sigma2))^(mu2/sqrt(sigma2) + y2) * digamma(mu2/sqrt(sigma2)))) * 
    log1p(sqrt(sigma2)) + (1 + sqrt(sigma2))^(mu2/sqrt(sigma2) + 
    y2) * (digamma(mu2/sqrt(sigma2))^2 + trigamma(mu2/sqrt(sigma2))))) * 
    gamma(1 + y2) * gamma(mu2/sqrt(sigma2))/((1 + sqrt(sigma2))^(mu2/sqrt(sigma2) + 
    y2) * gamma(1 + y2) * gamma(mu2/sqrt(sigma2)))^2)

der2pdf2.dermu2     <- der2pdf2.dermu2FUNC(y2, mu2, sigma2) 
    
    
der2pdf2.dersigma22FUNC <- function(y2, mu2, sigma2) (gamma(y2 + mu2/sqrt(sigma2)) * (sqrt(sigma2)^((y2 - 1) - 1) * 
    ((y2 - 1) * (0.5 * sigma2^-0.5)) * (y2 * (0.5 * sigma2^-0.5)) + 
    sqrt(sigma2)^(y2 - 1) * (y2 * (0.5 * (-0.5 * sigma2^-1.5)))) - 
    mu2 * (0.5 * sigma2^-0.5)/sqrt(sigma2)^2 * (gamma(y2 + mu2/sqrt(sigma2)) * 
        digamma(y2 + mu2/sqrt(sigma2))) * (sqrt(sigma2)^(y2 - 
        1) * (y2 * (0.5 * sigma2^-0.5))) - (((mu2 * (0.5 * (-0.5 * 
    sigma2^-1.5))/sqrt(sigma2)^2 - mu2 * (0.5 * sigma2^-0.5) * 
    (2 * (0.5 * sigma2^-0.5 * sqrt(sigma2)))/(sqrt(sigma2)^2)^2) * 
    (gamma(y2 + mu2/sqrt(sigma2)) * digamma(y2 + mu2/sqrt(sigma2))) - 
    mu2 * (0.5 * sigma2^-0.5)/sqrt(sigma2)^2 * (gamma(y2 + mu2/sqrt(sigma2)) * 
        (mu2 * (0.5 * sigma2^-0.5)/sqrt(sigma2)^2 * trigamma(y2 + 
            mu2/sqrt(sigma2))) + mu2 * (0.5 * sigma2^-0.5)/sqrt(sigma2)^2 * 
        (gamma(y2 + mu2/sqrt(sigma2)) * digamma(y2 + mu2/sqrt(sigma2))) * 
        digamma(y2 + mu2/sqrt(sigma2)))) * sqrt(sigma2)^y2 + 
    mu2 * (0.5 * sigma2^-0.5)/sqrt(sigma2)^2 * (gamma(y2 + mu2/sqrt(sigma2)) * 
        digamma(y2 + mu2/sqrt(sigma2))) * (sqrt(sigma2)^(y2 - 
        1) * (y2 * (0.5 * sigma2^-0.5)))))/(gamma(mu2/sqrt(sigma2)) * 
    gamma(y2 + 1) * (1 + sqrt(sigma2))^(y2 + mu2/sqrt(sigma2))) - 
    (gamma(y2 + mu2/sqrt(sigma2)) * (sqrt(sigma2)^(y2 - 1) * 
        (y2 * (0.5 * sigma2^-0.5))) - mu2 * (0.5 * sigma2^-0.5)/sqrt(sigma2)^2 * 
        (gamma(y2 + mu2/sqrt(sigma2)) * digamma(y2 + mu2/sqrt(sigma2))) * 
        sqrt(sigma2)^y2) * (gamma(mu2/sqrt(sigma2)) * gamma(y2 + 
        1) * ((1 + sqrt(sigma2))^((y2 + mu2/sqrt(sigma2)) - 1) * 
        ((y2 + mu2/sqrt(sigma2)) * (0.5 * sigma2^-0.5)) - (1 + 
        sqrt(sigma2))^(y2 + mu2/sqrt(sigma2)) * (log((1 + sqrt(sigma2))) * 
        (mu2 * (0.5 * sigma2^-0.5)/sqrt(sigma2)^2))) - mu2 * 
        (0.5 * sigma2^-0.5)/sqrt(sigma2)^2 * (gamma(mu2/sqrt(sigma2)) * 
        digamma(mu2/sqrt(sigma2))) * gamma(y2 + 1) * (1 + sqrt(sigma2))^(y2 + 
        mu2/sqrt(sigma2)))/(gamma(mu2/sqrt(sigma2)) * gamma(y2 + 
        1) * (1 + sqrt(sigma2))^(y2 + mu2/sqrt(sigma2)))^2 - 
    (((gamma(y2 + mu2/sqrt(sigma2)) * (sqrt(sigma2)^(y2 - 1) * 
        (y2 * (0.5 * sigma2^-0.5))) - mu2 * (0.5 * sigma2^-0.5)/sqrt(sigma2)^2 * 
        (gamma(y2 + mu2/sqrt(sigma2)) * digamma(y2 + mu2/sqrt(sigma2))) * 
        sqrt(sigma2)^y2) * (gamma(mu2/sqrt(sigma2)) * gamma(y2 + 
        1) * ((1 + sqrt(sigma2))^((y2 + mu2/sqrt(sigma2)) - 1) * 
        ((y2 + mu2/sqrt(sigma2)) * (0.5 * sigma2^-0.5)) - (1 + 
        sqrt(sigma2))^(y2 + mu2/sqrt(sigma2)) * (log((1 + sqrt(sigma2))) * 
        (mu2 * (0.5 * sigma2^-0.5)/sqrt(sigma2)^2))) - mu2 * 
        (0.5 * sigma2^-0.5)/sqrt(sigma2)^2 * (gamma(mu2/sqrt(sigma2)) * 
        digamma(mu2/sqrt(sigma2))) * gamma(y2 + 1) * (1 + sqrt(sigma2))^(y2 + 
        mu2/sqrt(sigma2))) + (gamma(y2 + mu2/sqrt(sigma2)) * 
        sqrt(sigma2)^y2) * (gamma(mu2/sqrt(sigma2)) * gamma(y2 + 
        1) * (((1 + sqrt(sigma2))^(((y2 + mu2/sqrt(sigma2)) - 
        1) - 1) * (((y2 + mu2/sqrt(sigma2)) - 1) * (0.5 * sigma2^-0.5)) - 
        (1 + sqrt(sigma2))^((y2 + mu2/sqrt(sigma2)) - 1) * (log((1 + 
            sqrt(sigma2))) * (mu2 * (0.5 * sigma2^-0.5)/sqrt(sigma2)^2))) * 
        ((y2 + mu2/sqrt(sigma2)) * (0.5 * sigma2^-0.5)) + (1 + 
        sqrt(sigma2))^((y2 + mu2/sqrt(sigma2)) - 1) * ((y2 + 
        mu2/sqrt(sigma2)) * (0.5 * (-0.5 * sigma2^-1.5)) - mu2 * 
        (0.5 * sigma2^-0.5)/sqrt(sigma2)^2 * (0.5 * sigma2^-0.5)) - 
        (((1 + sqrt(sigma2))^((y2 + mu2/sqrt(sigma2)) - 1) * 
            ((y2 + mu2/sqrt(sigma2)) * (0.5 * sigma2^-0.5)) - 
            (1 + sqrt(sigma2))^(y2 + mu2/sqrt(sigma2)) * (log((1 + 
                sqrt(sigma2))) * (mu2 * (0.5 * sigma2^-0.5)/sqrt(sigma2)^2))) * 
            (log((1 + sqrt(sigma2))) * (mu2 * (0.5 * sigma2^-0.5)/sqrt(sigma2)^2)) + 
            (1 + sqrt(sigma2))^(y2 + mu2/sqrt(sigma2)) * (0.5 * 
                sigma2^-0.5/(1 + sqrt(sigma2)) * (mu2 * (0.5 * 
                sigma2^-0.5)/sqrt(sigma2)^2) + log((1 + sqrt(sigma2))) * 
                (mu2 * (0.5 * (-0.5 * sigma2^-1.5))/sqrt(sigma2)^2 - 
                  mu2 * (0.5 * sigma2^-0.5) * (2 * (0.5 * sigma2^-0.5 * 
                    sqrt(sigma2)))/(sqrt(sigma2)^2)^2)))) - mu2 * 
        (0.5 * sigma2^-0.5)/sqrt(sigma2)^2 * (gamma(mu2/sqrt(sigma2)) * 
        digamma(mu2/sqrt(sigma2))) * gamma(y2 + 1) * ((1 + sqrt(sigma2))^((y2 + 
        mu2/sqrt(sigma2)) - 1) * ((y2 + mu2/sqrt(sigma2)) * (0.5 * 
        sigma2^-0.5)) - (1 + sqrt(sigma2))^(y2 + mu2/sqrt(sigma2)) * 
        (log((1 + sqrt(sigma2))) * (mu2 * (0.5 * sigma2^-0.5)/sqrt(sigma2)^2))) - 
        (((mu2 * (0.5 * (-0.5 * sigma2^-1.5))/sqrt(sigma2)^2 - 
            mu2 * (0.5 * sigma2^-0.5) * (2 * (0.5 * sigma2^-0.5 * 
                sqrt(sigma2)))/(sqrt(sigma2)^2)^2) * (gamma(mu2/sqrt(sigma2)) * 
            digamma(mu2/sqrt(sigma2))) - mu2 * (0.5 * sigma2^-0.5)/sqrt(sigma2)^2 * 
            (gamma(mu2/sqrt(sigma2)) * (mu2 * (0.5 * sigma2^-0.5)/sqrt(sigma2)^2 * 
                trigamma(mu2/sqrt(sigma2))) + mu2 * (0.5 * sigma2^-0.5)/sqrt(sigma2)^2 * 
                (gamma(mu2/sqrt(sigma2)) * digamma(mu2/sqrt(sigma2))) * 
                digamma(mu2/sqrt(sigma2)))) * gamma(y2 + 1) * 
            (1 + sqrt(sigma2))^(y2 + mu2/sqrt(sigma2)) + mu2 * 
            (0.5 * sigma2^-0.5)/sqrt(sigma2)^2 * (gamma(mu2/sqrt(sigma2)) * 
            digamma(mu2/sqrt(sigma2))) * gamma(y2 + 1) * ((1 + 
            sqrt(sigma2))^((y2 + mu2/sqrt(sigma2)) - 1) * ((y2 + 
            mu2/sqrt(sigma2)) * (0.5 * sigma2^-0.5)) - (1 + sqrt(sigma2))^(y2 + 
            mu2/sqrt(sigma2)) * (log((1 + sqrt(sigma2))) * (mu2 * 
            (0.5 * sigma2^-0.5)/sqrt(sigma2)^2))))))/(gamma(mu2/sqrt(sigma2)) * 
        gamma(y2 + 1) * (1 + sqrt(sigma2))^(y2 + mu2/sqrt(sigma2)))^2 - 
        (gamma(y2 + mu2/sqrt(sigma2)) * sqrt(sigma2)^y2) * (gamma(mu2/sqrt(sigma2)) * 
            gamma(y2 + 1) * ((1 + sqrt(sigma2))^((y2 + mu2/sqrt(sigma2)) - 
            1) * ((y2 + mu2/sqrt(sigma2)) * (0.5 * sigma2^-0.5)) - 
            (1 + sqrt(sigma2))^(y2 + mu2/sqrt(sigma2)) * (log((1 + 
                sqrt(sigma2))) * (mu2 * (0.5 * sigma2^-0.5)/sqrt(sigma2)^2))) - 
            mu2 * (0.5 * sigma2^-0.5)/sqrt(sigma2)^2 * (gamma(mu2/sqrt(sigma2)) * 
                digamma(mu2/sqrt(sigma2))) * gamma(y2 + 1) * 
                (1 + sqrt(sigma2))^(y2 + mu2/sqrt(sigma2))) * 
            (2 * ((gamma(mu2/sqrt(sigma2)) * gamma(y2 + 1) * 
                ((1 + sqrt(sigma2))^((y2 + mu2/sqrt(sigma2)) - 
                  1) * ((y2 + mu2/sqrt(sigma2)) * (0.5 * sigma2^-0.5)) - 
                  (1 + sqrt(sigma2))^(y2 + mu2/sqrt(sigma2)) * 
                    (log((1 + sqrt(sigma2))) * (mu2 * (0.5 * 
                      sigma2^-0.5)/sqrt(sigma2)^2))) - mu2 * 
                (0.5 * sigma2^-0.5)/sqrt(sigma2)^2 * (gamma(mu2/sqrt(sigma2)) * 
                digamma(mu2/sqrt(sigma2))) * gamma(y2 + 1) * 
                (1 + sqrt(sigma2))^(y2 + mu2/sqrt(sigma2))) * 
                (gamma(mu2/sqrt(sigma2)) * gamma(y2 + 1) * (1 + 
                  sqrt(sigma2))^(y2 + mu2/sqrt(sigma2)))))/((gamma(mu2/sqrt(sigma2)) * 
            gamma(y2 + 1) * (1 + sqrt(sigma2))^(y2 + mu2/sqrt(sigma2)))^2)^2) 

der2pdf2.dersigma22     <- der2pdf2.dersigma22FUNC(y2, mu2, sigma2) 
   
   
der2pdf2.mu2dersigma2FUNC <- function(y2, mu2, sigma2) (1/sqrt(sigma2) * (gamma(y2 + mu2/sqrt(sigma2)) * digamma(y2 + 
    mu2/sqrt(sigma2))) * (sqrt(sigma2)^(y2 - 1) * (y2 * (0.5 * 
    sigma2^-0.5))) - (1/sqrt(sigma2) * (gamma(y2 + mu2/sqrt(sigma2)) * 
    (mu2 * (0.5 * sigma2^-0.5)/sqrt(sigma2)^2 * trigamma(y2 + 
        mu2/sqrt(sigma2))) + mu2 * (0.5 * sigma2^-0.5)/sqrt(sigma2)^2 * 
    (gamma(y2 + mu2/sqrt(sigma2)) * digamma(y2 + mu2/sqrt(sigma2))) * 
    digamma(y2 + mu2/sqrt(sigma2))) + 0.5 * sigma2^-0.5/sqrt(sigma2)^2 * 
    (gamma(y2 + mu2/sqrt(sigma2)) * digamma(y2 + mu2/sqrt(sigma2)))) * 
    sqrt(sigma2)^y2)/(gamma(mu2/sqrt(sigma2)) * gamma(y2 + 1) * 
    (1 + sqrt(sigma2))^(y2 + mu2/sqrt(sigma2))) - 1/sqrt(sigma2) * 
    (gamma(y2 + mu2/sqrt(sigma2)) * digamma(y2 + mu2/sqrt(sigma2))) * 
    sqrt(sigma2)^y2 * (gamma(mu2/sqrt(sigma2)) * gamma(y2 + 1) * 
    ((1 + sqrt(sigma2))^((y2 + mu2/sqrt(sigma2)) - 1) * ((y2 + 
        mu2/sqrt(sigma2)) * (0.5 * sigma2^-0.5)) - (1 + sqrt(sigma2))^(y2 + 
        mu2/sqrt(sigma2)) * (log((1 + sqrt(sigma2))) * (mu2 * 
        (0.5 * sigma2^-0.5)/sqrt(sigma2)^2))) - mu2 * (0.5 * 
    sigma2^-0.5)/sqrt(sigma2)^2 * (gamma(mu2/sqrt(sigma2)) * 
    digamma(mu2/sqrt(sigma2))) * gamma(y2 + 1) * (1 + sqrt(sigma2))^(y2 + 
    mu2/sqrt(sigma2)))/(gamma(mu2/sqrt(sigma2)) * gamma(y2 + 
    1) * (1 + sqrt(sigma2))^(y2 + mu2/sqrt(sigma2)))^2 - (((gamma(y2 + 
    mu2/sqrt(sigma2)) * (sqrt(sigma2)^(y2 - 1) * (y2 * (0.5 * 
    sigma2^-0.5))) - mu2 * (0.5 * sigma2^-0.5)/sqrt(sigma2)^2 * 
    (gamma(y2 + mu2/sqrt(sigma2)) * digamma(y2 + mu2/sqrt(sigma2))) * 
    sqrt(sigma2)^y2) * (1/sqrt(sigma2) * (gamma(mu2/sqrt(sigma2)) * 
    digamma(mu2/sqrt(sigma2))) * gamma(y2 + 1) * (1 + sqrt(sigma2))^(y2 + 
    mu2/sqrt(sigma2)) + gamma(mu2/sqrt(sigma2)) * gamma(y2 + 
    1) * ((1 + sqrt(sigma2))^(y2 + mu2/sqrt(sigma2)) * (log((1 + 
    sqrt(sigma2))) * (1/sqrt(sigma2))))) + (gamma(y2 + mu2/sqrt(sigma2)) * 
    sqrt(sigma2)^y2) * (1/sqrt(sigma2) * (gamma(mu2/sqrt(sigma2)) * 
    digamma(mu2/sqrt(sigma2))) * gamma(y2 + 1) * ((1 + sqrt(sigma2))^((y2 + 
    mu2/sqrt(sigma2)) - 1) * ((y2 + mu2/sqrt(sigma2)) * (0.5 * 
    sigma2^-0.5)) - (1 + sqrt(sigma2))^(y2 + mu2/sqrt(sigma2)) * 
    (log((1 + sqrt(sigma2))) * (mu2 * (0.5 * sigma2^-0.5)/sqrt(sigma2)^2))) - 
    (1/sqrt(sigma2) * (gamma(mu2/sqrt(sigma2)) * (mu2 * (0.5 * 
        sigma2^-0.5)/sqrt(sigma2)^2 * trigamma(mu2/sqrt(sigma2))) + 
        mu2 * (0.5 * sigma2^-0.5)/sqrt(sigma2)^2 * (gamma(mu2/sqrt(sigma2)) * 
            digamma(mu2/sqrt(sigma2))) * digamma(mu2/sqrt(sigma2))) + 
        0.5 * sigma2^-0.5/sqrt(sigma2)^2 * (gamma(mu2/sqrt(sigma2)) * 
            digamma(mu2/sqrt(sigma2)))) * gamma(y2 + 1) * (1 + 
        sqrt(sigma2))^(y2 + mu2/sqrt(sigma2)) + (gamma(mu2/sqrt(sigma2)) * 
    gamma(y2 + 1) * (((1 + sqrt(sigma2))^((y2 + mu2/sqrt(sigma2)) - 
    1) * ((y2 + mu2/sqrt(sigma2)) * (0.5 * sigma2^-0.5)) - (1 + 
    sqrt(sigma2))^(y2 + mu2/sqrt(sigma2)) * (log((1 + sqrt(sigma2))) * 
    (mu2 * (0.5 * sigma2^-0.5)/sqrt(sigma2)^2))) * (log((1 + 
    sqrt(sigma2))) * (1/sqrt(sigma2))) + (1 + sqrt(sigma2))^(y2 + 
    mu2/sqrt(sigma2)) * (0.5 * sigma2^-0.5/(1 + sqrt(sigma2)) * 
    (1/sqrt(sigma2)) - log((1 + sqrt(sigma2))) * (0.5 * sigma2^-0.5/sqrt(sigma2)^2))) - 
    mu2 * (0.5 * sigma2^-0.5)/sqrt(sigma2)^2 * (gamma(mu2/sqrt(sigma2)) * 
        digamma(mu2/sqrt(sigma2))) * gamma(y2 + 1) * ((1 + sqrt(sigma2))^(y2 + 
        mu2/sqrt(sigma2)) * (log((1 + sqrt(sigma2))) * (1/sqrt(sigma2)))))))/(gamma(mu2/sqrt(sigma2)) * 
    gamma(y2 + 1) * (1 + sqrt(sigma2))^(y2 + mu2/sqrt(sigma2)))^2 - 
    (gamma(y2 + mu2/sqrt(sigma2)) * sqrt(sigma2)^y2) * (1/sqrt(sigma2) * 
        (gamma(mu2/sqrt(sigma2)) * digamma(mu2/sqrt(sigma2))) * 
        gamma(y2 + 1) * (1 + sqrt(sigma2))^(y2 + mu2/sqrt(sigma2)) + 
        gamma(mu2/sqrt(sigma2)) * gamma(y2 + 1) * ((1 + sqrt(sigma2))^(y2 + 
            mu2/sqrt(sigma2)) * (log((1 + sqrt(sigma2))) * (1/sqrt(sigma2))))) * 
        (2 * ((gamma(mu2/sqrt(sigma2)) * gamma(y2 + 1) * ((1 + 
            sqrt(sigma2))^((y2 + mu2/sqrt(sigma2)) - 1) * ((y2 + 
            mu2/sqrt(sigma2)) * (0.5 * sigma2^-0.5)) - (1 + sqrt(sigma2))^(y2 + 
            mu2/sqrt(sigma2)) * (log((1 + sqrt(sigma2))) * (mu2 * 
            (0.5 * sigma2^-0.5)/sqrt(sigma2)^2))) - mu2 * (0.5 * 
            sigma2^-0.5)/sqrt(sigma2)^2 * (gamma(mu2/sqrt(sigma2)) * 
            digamma(mu2/sqrt(sigma2))) * gamma(y2 + 1) * (1 + 
            sqrt(sigma2))^(y2 + mu2/sqrt(sigma2))) * (gamma(mu2/sqrt(sigma2)) * 
            gamma(y2 + 1) * (1 + sqrt(sigma2))^(y2 + mu2/sqrt(sigma2)))))/((gamma(mu2/sqrt(sigma2)) * 
        gamma(y2 + 1) * (1 + sqrt(sigma2))^(y2 + mu2/sqrt(sigma2)))^2)^2) 

der2pdf2.mu2dersigma2     <- der2pdf2.mu2dersigma2FUNC(y2, mu2, sigma2)




if(naive == FALSE){   
 
p2  <- pNBII(y2, mu = mu2, sigma = sqrt(sigma2)) 
 
ly2 <- length(y2)
if(length(sigma2) == 1) sigma2 <- c(rep(sigma2, ly2))

mu2 <- c(mu2)
sigma2 <- c(sigma2)

derp2.dermu2           <- rowSums( derpdf2.dermu2FUNC(        y2m, mu2, sigma2 ) , na.rm = TRUE )
derp2.dersigma2        <- rowSums( derpdf2.sigma2FUNC(        y2m, mu2, sigma2 ) , na.rm = TRUE )
der2p2.dermu22         <- rowSums( der2pdf2.dermu2FUNC(       y2m, mu2, sigma2 ) , na.rm = TRUE )
der2p2.dersigma22      <- rowSums( der2pdf2.dersigma22FUNC(   y2m, mu2, sigma2 ) , na.rm = TRUE )
der2p2.derdermu2sigma2 <- rowSums( der2pdf2.mu2dersigma2FUNC( y2m, mu2, sigma2 ) , na.rm = TRUE )
                      
                        
                        
                   }

}





#############################

if(margin2 == "PIG"){ # K all numerical as well
                      # sigma <- ifelse(sigma>precision, sigma, precision)
                      # sigma <- ifelse(sigma<10^7, sigma, 10^7)
                      # tolerances for derivatives may be changed to get better
                      # performace, difficult to find a general rule here
    
    

mu2 <- exp(eta2)
dermu2.dereta2 <- exp(eta2)
der2mu2.dereta2eta2 <- exp(eta2) 

dersigma2.dersigma2.st  <- exp(sigma2.st)  
dersigma2.dersigma2.st2 <- exp(sigma2.st)  


pdf2 <- dPIG(y2, mu = mu2, sigma = sqrt(sigma2))    

derpdf2.dermu2F <- derpdf2.dermu2FUNC2p(function(mu2) dPIG(y2, mu = mu2, sigma = sqrt(sigma2)), y2, mu2, sigma2) 
derpdf2.dermu2  <- derpdf2.dermu2F$fi 
der2pdf2.dermu2 <- derpdf2.dermu2F$se     
    
 
derpdf2.sigma2F     <- derpdf2.sigma2FUNC2p(function(sigma2) dPIG(y2, mu = mu2, sigma = sqrt(sigma2)), y2, mu2, sigma2) 
derpdf2.sigma2      <- derpdf2.sigma2F$fi      
der2pdf2.dersigma22 <- derpdf2.sigma2F$se 
   
der2pdf2.mu2dersigma2 <- der2pdf2.mu2dersigma2FUNC2p(function(mu2, sigma2) dPIG(y2, mu = mu2, sigma = sqrt(sigma2)), y2, mu2, sigma2)


if(naive == FALSE){   
 
p2  <- pPIG(y2, mu = mu2, sigma = sqrt(sigma2)) 
 
derp2.dermu2F  <- derpdf2.dermu2FUNC2p(function(mu2) pPIG(y2, mu = mu2, sigma = sqrt(sigma2)), y2, mu2, sigma2)
derp2.dermu2   <- derp2.dermu2F$fi
der2p2.dermu22 <- derp2.dermu2F$se
    
derp2.dersigma2F <- derpdf2.sigma2FUNC2p(function(sigma2) pPIG(y2, mu = mu2, sigma = sqrt(sigma2)), y2, mu2, sigma2) 
derp2.dersigma2  <- derp2.dersigma2F$fi 
der2p2.dersigma22<- derp2.dersigma2F$se 
   
der2p2.derdermu2sigma2 <- der2pdf2.mu2dersigma2FUNC2p(function(mu2, sigma2) pPIG(y2, mu = mu2, sigma = sqrt(sigma2)), y2, mu2, sigma2) 
  

                   }

}



if(margin2 == "DGP"){
    

                   mu2  <- eta2 # exi
dermu2.dereta2          <- 1
der2mu2.dereta2eta2     <- 0 

dersigma2.dersigma2.st  <- exp(sigma2.st)  # mu
dersigma2.dersigma2.st2 <- exp(sigma2.st)  


indx1 <- as.numeric( ((1 + mu2*y2/sqrt(sigma2))     > 0) == FALSE ) 
indx2 <- as.numeric( ((1 + mu2*(y2+1)/sqrt(sigma2)) > 0) == FALSE )
indx  <- rowSums(cbind(indx1, indx2))



pdf2FUNC <- function(y2, mu2, sigma2) suppressWarnings(   (1 + mu2*y2/sqrt(sigma2))^(-1/mu2) - (1 + mu2*(1+y2)/sqrt(sigma2))^(-1/mu2)    )  



derpdf2.dermu2FUNC        <- function(y2, mu2, sigma2) suppressWarnings(   (((1 + y2)/(1 + mu2 * (1 + y2)/sqrt(sigma2))^(1 + 1/mu2) - y2/(1 + 
    mu2 * y2/sqrt(sigma2))^(1 + 1/mu2))/sqrt(sigma2) + (log1p(mu2 * 
    y2/sqrt(sigma2))/(1 + mu2 * y2/sqrt(sigma2))^(1/mu2) - log1p(mu2 * 
    (1 + y2)/sqrt(sigma2))/(1 + mu2 * (1 + y2)/sqrt(sigma2))^(1/mu2))/mu2)/mu2   )


derpdf2.sigma2FUNC        <- function(y2, mu2, sigma2)  suppressWarnings(    -((0.5 * ((1 + y2)/(1 + mu2 * (1 + y2)/sqrt(sigma2))^(1 + 1/mu2)) - 
    0.5 * (y2/(1 + mu2 * y2/sqrt(sigma2))^(1 + 1/mu2)))/(sigma2 * 
    sqrt(sigma2)))     )

       
der2pdf2.dermu2FUNC       <- function(y2, mu2, sigma2) suppressWarnings(   (((log1p(mu2 * y2/sqrt(sigma2)) * (log1p(mu2 * y2/sqrt(sigma2))/(mu2 * 
    (1 + mu2 * y2/sqrt(sigma2))^(1/mu2)) - y2/((1 + mu2 * y2/sqrt(sigma2))^(1 + 
    1/mu2) * sqrt(sigma2))) - log1p(mu2 * (1 + y2)/sqrt(sigma2)) * 
    (log1p(mu2 * (1 + y2)/sqrt(sigma2))/(mu2 * (1 + mu2 * (1 + 
        y2)/sqrt(sigma2))^(1/mu2)) - (1 + y2)/((1 + mu2 * (1 + 
        y2)/sqrt(sigma2))^(1 + 1/mu2) * sqrt(sigma2))))/mu2 + 
    (y2/((1 + mu2 * y2/sqrt(sigma2)) * sqrt(sigma2)) - 2 * (log1p(mu2 * 
        y2/sqrt(sigma2))/mu2))/(1 + mu2 * y2/sqrt(sigma2))^(1/mu2) - 
    ((1 + y2)/((1 + mu2 * (1 + y2)/sqrt(sigma2)) * sqrt(sigma2)) - 
        2 * (log1p(mu2 * (1 + y2)/sqrt(sigma2))/mu2))/(1 + mu2 * 
        (1 + y2)/sqrt(sigma2))^(1/mu2))/mu2 + (y2 * ((1/(1 + 
    mu2 * y2/sqrt(sigma2))^(1 + 1/mu2) - log1p(mu2 * y2/sqrt(sigma2))/(mu2 * 
    (1 + mu2 * y2/sqrt(sigma2))^(1 + 1/mu2)))/mu2 + y2 * (1 + 
    1/mu2)/((1 + mu2 * y2/sqrt(sigma2))^(1/mu2 + 2) * sqrt(sigma2))) - 
    ((1 + 1/mu2) * (1 + y2)/((1 + mu2 * (1 + y2)/sqrt(sigma2))^(1/mu2 + 
        2) * sqrt(sigma2)) + (1/(1 + mu2 * (1 + y2)/sqrt(sigma2))^(1 + 
        1/mu2) - log1p(mu2 * (1 + y2)/sqrt(sigma2))/(mu2 * (1 + 
        mu2 * (1 + y2)/sqrt(sigma2))^(1 + 1/mu2)))/mu2) * (1 + 
        y2))/sqrt(sigma2))/mu2    )

    
der2pdf2.dersigma22FUNC   <- function(y2, mu2, sigma2) suppressWarnings(    -((0.5 * ((1 + y2) * ((1 + mu2 * (1 + y2)/sqrt(sigma2))^((1 + 
    1/mu2) - 1) * ((1 + 1/mu2) * (mu2 * (1 + y2) * (0.5 * sigma2^-0.5)/sqrt(sigma2)^2)))/((1 + 
    mu2 * (1 + y2)/sqrt(sigma2))^(1 + 1/mu2))^2) - 0.5 * (y2 * 
    ((1 + mu2 * y2/sqrt(sigma2))^((1 + 1/mu2) - 1) * ((1 + 1/mu2) * 
        (mu2 * y2 * (0.5 * sigma2^-0.5)/sqrt(sigma2)^2)))/((1 + 
    mu2 * y2/sqrt(sigma2))^(1 + 1/mu2))^2))/(sigma2 * sqrt(sigma2)) - 
    (0.5 * ((1 + y2)/(1 + mu2 * (1 + y2)/sqrt(sigma2))^(1 + 1/mu2)) - 
        0.5 * (y2/(1 + mu2 * y2/sqrt(sigma2))^(1 + 1/mu2))) * 
        (sqrt(sigma2) + sigma2 * (0.5 * sigma2^-0.5))/(sigma2 * 
        sqrt(sigma2))^2)     )

    
der2pdf2.mu2dersigma2FUNC <- function(y2, mu2, sigma2) suppressWarnings(    -((0.5 * (((1 + 1/mu2) * (1 + mu2 * (1 + y2)/sqrt(sigma2))^(1/mu2) * 
    (1 + y2)/sqrt(sigma2) - (1 + mu2 * (1 + y2)/sqrt(sigma2))^(1 + 
    1/mu2) * log1p(mu2 * (1 + y2)/sqrt(sigma2))/mu2^2) * (1 + 
    y2)/(1 + mu2 * (1 + y2)/sqrt(sigma2))^(2 * (1 + 1/mu2))) - 
    0.5 * (y2 * (y2 * (1 + 1/mu2) * (1 + mu2 * y2/sqrt(sigma2))^(1/mu2)/sqrt(sigma2) - 
        (1 + mu2 * y2/sqrt(sigma2))^(1 + 1/mu2) * log1p(mu2 * 
            y2/sqrt(sigma2))/mu2^2)/(1 + mu2 * y2/sqrt(sigma2))^(2 * 
        (1 + 1/mu2))))/(sigma2 * sqrt(sigma2)))      )




pdf2                  <- as.numeric( pdf2FUNC(y2, mu2, sigma2) )  
derpdf2.dermu2        <- derpdf2.dermu2FUNC(y2, mu2, sigma2) 
derpdf2.sigma2        <- derpdf2.sigma2FUNC(y2, mu2, sigma2) 
der2pdf2.dermu2       <- der2pdf2.dermu2FUNC(y2, mu2, sigma2) 
der2pdf2.dersigma22   <- der2pdf2.dersigma22FUNC(y2, mu2, sigma2) 
der2pdf2.mu2dersigma2 <- der2pdf2.mu2dersigma2FUNC(y2, mu2, sigma2)



pdf2                  <- ifelse( indx == 0, pdf2, 1) # 'cause in output log(1) = 0 hence it will not add anything to the lik
derpdf2.dermu2        <- ifelse( indx == 0, derpdf2.dermu2, 0) 
derpdf2.sigma2        <- ifelse( indx == 0, derpdf2.sigma2, 0) 
der2pdf2.dermu2       <- ifelse( indx == 0, der2pdf2.dermu2, 0) 
der2pdf2.dersigma22   <- ifelse( indx == 0, der2pdf2.dersigma22, 0) 
der2pdf2.mu2dersigma2 <- ifelse( indx == 0, der2pdf2.mu2dersigma2, 0) 








if(naive == FALSE){   # needs y2m
 
ly2 <- length(y2)
if(length(sigma2) == 1) sigma2 <- c(rep(sigma2, ly2))

mu2    <- c(mu2)
sigma2 <- c(sigma2)

derp2.dermu2           <- rowSums( derpdf2.dermu2FUNC(        y2m, mu2, sigma2 ) , na.rm = TRUE )
derp2.dersigma2        <- rowSums( derpdf2.sigma2FUNC(        y2m, mu2, sigma2 ) , na.rm = TRUE )
der2p2.dermu22         <- rowSums( der2pdf2.dermu2FUNC(       y2m, mu2, sigma2 ) , na.rm = TRUE )
der2p2.dersigma22      <- rowSums( der2pdf2.dersigma22FUNC(   y2m, mu2, sigma2 ) , na.rm = TRUE )
der2p2.derdermu2sigma2 <- rowSums( der2pdf2.mu2dersigma2FUNC( y2m, mu2, sigma2 ) , na.rm = TRUE )


                   }

}






##########################################################




if(margin2 %in% c(cont2par,cont3par)){ 

derpdf2.dereta2              <- derpdf2.dermu2*dermu2.dereta2       
der2pdf2.dereta2             <- der2pdf2.dermu2* dermu2.dereta2^2 + derpdf2.dermu2*der2mu2.dereta2eta2        


der2pdf2.dereta2dersigma2    <- der2pdf2.mu2dersigma2* dermu2.dereta2
der2pdf2.dereta2dersigma2.st <- der2pdf2.dereta2dersigma2 *  dersigma2.dersigma2.st

  
derpdf2.dersigma2.st         <- derpdf2.sigma2 * dersigma2.dersigma2.st   
der2pdf2.dersigma2.st2       <- der2pdf2.dersigma22 * dersigma2.dersigma2.st^2 + derpdf2.sigma2  * dersigma2.dersigma2.st2     

                  

if(naive == FALSE){  



derp2.dereta2                <- derp2.dermu2*dermu2.dereta2
der2p2.dereta2eta2           <- der2p2.dermu22*dermu2.dereta2^2+derp2.dermu2*der2mu2.dereta2eta2      


der2p2.dereta2dersigma2      <- der2p2.derdermu2sigma2* dermu2.dereta2    
der2p2.dereta2dersigma2.st   <- der2p2.dereta2dersigma2 *  dersigma2.dersigma2.st  

derp2.dersigma.st            <- derp2.dersigma2 *  dersigma2.dersigma2.st 
der2p2.dersigma2.st2         <- der2p2.dersigma22 * dersigma2.dersigma2.st^2 + derp2.dersigma2 * dersigma2.dersigma2.st2





                 }
                 
                 
                 

}###############



epsilon <- 0.0000001 
max.p   <- 0.9999999

  pdf2 <- ifelse(pdf2 < epsilon, epsilon, pdf2 )

  p2   <- mm(p2) 



ifef <- function(dv){

epsilon <- 0.0000001 
dv <- ifelse(is.na(dv), epsilon, dv ) 
dv <- ifelse(dv == Inf ,  8.218407e+20, dv )
dv <- ifelse(dv == -Inf ,  -8.218407e+20, dv )
dv

}


# for safety

pdf2                         = ifef(pdf2)
p2                           = ifef(p2) 
derpdf2.dereta2              = ifef(derpdf2.dereta2) 
derpdf2.dersigma2.st         = ifef(derpdf2.dersigma2.st) 
derp2.dersigma.st            = ifef(derp2.dersigma.st)
derp2.dereta2                = ifef(derp2.dereta2)
der2p2.dereta2eta2           = ifef(der2p2.dereta2eta2) 
der2pdf2.dereta2             = ifef(der2pdf2.dereta2)
der2p2.dersigma2.st2         = ifef(der2p2.dersigma2.st2) 
der2pdf2.dersigma2.st2       = ifef(der2pdf2.dersigma2.st2)
der2p2.dereta2dersigma2.st   = ifef(der2p2.dereta2dersigma2.st)  
der2pdf2.dereta2dersigma2.st = ifef(der2pdf2.dereta2dersigma2.st)
der2pdf2.dereta2dernu.st     = ifef(der2pdf2.dereta2dernu.st)   
der2pdf2.sigma2.st2dernu.st  = ifef(der2pdf2.sigma2.st2dernu.st)
derpdf2.dernu.st             = ifef(derpdf2.dernu.st)           
der2pdf2.dernu.st2           = ifef(der2pdf2.dernu.st2)         
derp2.nu.st                  = ifef(derp2.nu.st)                
der2p2.dernu.st2             = ifef(der2p2.dernu.st2)           
der2p2.dereta2dernu.st       = ifef(der2p2.dereta2dernu.st)     
der2p2.dersigma2.stdernu.st  = ifef(der2p2.dersigma2.stdernu.st)




list(pdf2                         = pdf2,
     p2                           = p2, 
     derpdf2.dereta2              = derpdf2.dereta2, 
     derpdf2.dersigma2.st         = derpdf2.dersigma2.st, 
     derp2.dersigma.st            = derp2.dersigma.st,
     derp2.dereta2                = derp2.dereta2,
     der2p2.dereta2eta2           = der2p2.dereta2eta2, 
     der2pdf2.dereta2             = der2pdf2.dereta2,
     der2p2.dersigma2.st2         = der2p2.dersigma2.st2, 
     der2pdf2.dersigma2.st2       = der2pdf2.dersigma2.st2,
     der2p2.dereta2dersigma2.st   = der2p2.dereta2dersigma2.st,            
     der2pdf2.dereta2dersigma2.st = der2pdf2.dereta2dersigma2.st,
     der2pdf2.dereta2dernu.st     = der2pdf2.dereta2dernu.st,   
     der2pdf2.sigma2.st2dernu.st  = der2pdf2.sigma2.st2dernu.st,
     derpdf2.dernu.st             = derpdf2.dernu.st,           
     der2pdf2.dernu.st2           = der2pdf2.dernu.st2,         
     derp2.nu.st                  = derp2.nu.st,                
     der2p2.dernu.st2             = der2p2.dernu.st2,           
     der2p2.dereta2dernu.st       = der2p2.dereta2dernu.st,     
     der2p2.dersigma2.stdernu.st  = der2p2.dersigma2.stdernu.st, 
     indx = indx == 0)     


}




    