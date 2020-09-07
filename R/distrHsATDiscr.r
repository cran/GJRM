distrHsATDiscr <- function(y2, eta2, sigma2, nu, margin2, y2m, robust = FALSE,
                           min.dn, min.pr, max.pr){



p2 <- 1

sigma <- sigma2

if(margin2 %in% c("NBI","NBIa")){

pdf2 <- dNBI(y2, mu = exp(eta2), sigma = sigma)  # gamma(y2+1/sqrt(sigma2))/(gamma(1/sqrt(sigma2))*gamma(y2+1))*(sqrt(sigma2)*exp(eta2)/(1+sqrt(sigma2)*exp(eta2)))^y2*(1/(1+sqrt(sigma2)*exp(eta2)))^(1/sqrt(sigma2))
p2   <- pNBI(y2, mu = exp(eta2), sigma = sigma)  


}

if(margin2 %in% c("NBII","NBIIa")){

pdf2 <- dNBII(y2, mu = exp(eta2), sigma = sigma) # (gamma(y2 + exp(eta2)/sqrt(sigma2))*sqrt(sigma2)^y2)/(gamma(exp(eta2)/sqrt(sigma2))*gamma(y2+1)*(1+sqrt(sigma2))^(y2+exp(eta2)/sqrt(sigma2)))  
p2   <- pNBII(y2, mu = exp(eta2), sigma = sigma)  



}



if( margin2 %in% c("DGP","DGPII") ){




if( margin2 == "DGP" )   mu2 <- c(eta2)
if( margin2 == "DGPII" ) mu2 <- c(eta2^2)




sigma <- sigma2 <- c(sigma2)

     
pdf2FUNC <- function(y2, mu2, sigma) suppressWarnings(   (1 + mu2*y2/sigma)^(-1/mu2) - (1 + mu2*(1+y2)/sigma)^(-1/mu2)    )
pdf2     <-   as.numeric( pdf2FUNC(y2, mu2, sigma) )  

if(robust == FALSE) p2  <- suppressWarnings(     rowSums( matrix(as.numeric( pdf2FUNC(y2m, mu2, sigma)),dim(y2m)[1],dim(y2m)[2]), na.rm = TRUE )      )


indx1 <- as.numeric( ((1 + mu2*y2/sigma)     > 0) == FALSE ) 
indx2 <- as.numeric( ((1 + mu2*(y2+1)/sigma) > 0) == FALSE ) # not needed
indx  <- rowSums(cbind(indx1, indx2))



pdf2 <- ifelse( indx == 0, pdf2, 0)                 #
if(robust == FALSE) p2 <- ifelse( indx == 0, p2, 0) #



}










if(margin2 == "PIG"){

pdf2 <- dPIG(y2, mu = exp(eta2), sigma = sigma)  

if(length(y2) == 1){ # pPIG can not use duplicate rule automatically 

   ml <- max( c( length(y2), length(eta2), length(sigma) ) )

   if( length(eta2) == 1 )   eta2 <- rep(eta2,  length = ml)
   if( length(sigma) == 1 ) sigma <- rep(sigma, length = ml)

   p2 <- NA

   for(i in 1:ml) p2[i] <- pPIG(y2, mu = exp(eta2[i]), sigma = sigma[i])  

}else p2   <- pPIG(y2, mu = exp(eta2), sigma = sigma)  







}

if(margin2 == "DEL"){

pdf2 <- dDEL(y2, mu = exp(eta2), sigma = sigma, nu = nu)  
p2   <- pDEL(y2, mu = exp(eta2), sigma = sigma, nu = nu)  




}

if(margin2 == "SICHEL"){

pdf2 <- dSICHEL(y2, mu = exp(eta2), sigma = sigma)  
p2   <- pSICHEL(y2, mu = exp(eta2), sigma = sigma)  




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


pdf2FUNC2 <- function(y2, mu2) mu2^y2/(exp(mu2)-1)*1/factorial(y2)  
pdf2     <- as.numeric( pdf2FUNC2(y2, mu2) ) 

if(robust == FALSE) p2  <- rowSums( matrix(as.numeric( pdf2FUNC2(y2m, mu2)),dim(y2m)[1],dim(y2m)[2]), na.rm = TRUE ) 




}



pdf2 <- ifelse(pdf2 < min.dn, min.dn, pdf2)
p2   <- mm(p2, min.pr = min.pr, max.pr = max.pr) 




list(pdf2 = ifef(pdf2), p2 = ifef(p2))     


}




     























