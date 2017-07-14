distrHsATDiscr <- function(y2, eta2, sigma2, nu, margin2, y2m, robust = FALSE){

p2 <- 1

if(margin2 %in% c("NBI","NBIa")){

pdf2 <- dNBI(y2, mu = exp(eta2), sigma = sqrt(sigma2))  # gamma(y2+1/sqrt(sigma2))/(gamma(1/sqrt(sigma2))*gamma(y2+1))*(sqrt(sigma2)*exp(eta2)/(1+sqrt(sigma2)*exp(eta2)))^y2*(1/(1+sqrt(sigma2)*exp(eta2)))^(1/sqrt(sigma2))
p2   <- pNBI(y2, mu = exp(eta2), sigma = sqrt(sigma2))  

}

if(margin2 %in% c("NBII","NBIIa")){

pdf2 <- dNBII(y2, mu = exp(eta2), sigma = sqrt(sigma2)) # (gamma(y2 + exp(eta2)/sqrt(sigma2))*sqrt(sigma2)^y2)/(gamma(exp(eta2)/sqrt(sigma2))*gamma(y2+1)*(1+sqrt(sigma2))^(y2+exp(eta2)/sqrt(sigma2)))  
p2   <- pNBII(y2, mu = exp(eta2), sigma = sqrt(sigma2))  

}

if(margin2 == "PIG"){

pdf2 <- dPIG(y2, mu = exp(eta2), sigma = sqrt(sigma2))  
p2   <- pPIG(y2, mu = exp(eta2), sigma = sqrt(sigma2))  

}

if(margin2 == "DEL"){

pdf2 <- dDEL(y2, mu = exp(eta2), sigma = sqrt(sigma2), nu = nu)  
p2   <- pDEL(y2, mu = exp(eta2), sigma = sqrt(sigma2), nu = nu)  

}

if(margin2 == "SICHEL"){

pdf2 <- dSICHEL(y2, mu = exp(eta2), sigma = sqrt(sigma2))  
p2   <- pSICHEL(y2, mu = exp(eta2), sigma = sqrt(sigma2))  

}

if(margin2 == "PO"){


mu2 <- c(exp(eta2))
     
if(max(y2) > 170){

prec <- pmax(53, getPrec(mu2), getPrec(y2))
        
mu2 <- mpfr(mu2, prec)
y2  <- mpfr( y2, prec)        
        
} 


pdf2 <- as.numeric( (exp(-mu2)*mu2^y2)/factorial(y2) )  

p2  <- pPO(as.numeric(y2), mu = as.numeric(mu2)) 

           
}



if(margin2 == "ZTP"){# we need y2m especially here as there is
                     # no other function I can use



mu2 <- c(exp(eta2))
     
if(max(y2) > 170){

prec <- pmax(53, getPrec(mu2), getPrec(y2))
        
mu2 <- mpfr(mu2, prec)
y2  <- mpfr( y2, prec)        
        
} 


pdf2FUNC <- function(y2, mu2) mu2^y2/(exp(mu2)-1)*1/factorial(y2)  
pdf2     <- as.numeric( pdf2FUNC(y2, mu2) ) 

if(robust == FALSE) p2  <- rowSums( matrix(as.numeric( pdf2FUNC(y2m, mu2)),dim(y2m)[1],dim(y2m)[2]), na.rm = TRUE ) 


}













epsilon <- 0.0000001 
pdf2 <- ifelse(pdf2 < epsilon, epsilon, pdf2 )

p2 <- mm(p2) 


ifef <- function(dv){

epsilon <- 0.0000001 
dv <- ifelse(is.na(dv), epsilon, dv ) 
dv <- ifelse(dv == Inf ,  8.218407e+20, dv )
dv <- ifelse(dv == -Inf ,  -8.218407e+20, dv )
dv

}

# for safety

pdf2 = ifef(pdf2)
p2   = ifef(p2)





list(pdf2 = pdf2, p2 = p2)     


}




     























