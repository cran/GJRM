distrHsATDiscr2 <- function(y2, eta2, sigma2, nu, margin2, min.dn, left.trunc = 0){




sigma <- sigma2


if(margin2 %in% c("NBI")){

pdf2 <- dNBI(y2, mu = exp(eta2), sigma = sigma)   

}

if(margin2 %in% c("NBII")){

pdf2 <- dNBII(y2, mu = exp(eta2), sigma = sigma) 

}

if(margin2 == "PIG"){

pdf2 <- dPIG(y2, mu = exp(eta2), sigma = sigma)  

}



if(margin2 %in% c("tNBI")){

pdf2 <- dNBItr(y2, mu = exp(eta2), sigma = sigma, left.trunc)   

}

if(margin2 %in% c("tNBII")){

pdf2 <- dNBIItr(y2, mu = exp(eta2), sigma = sigma, left.trunc) 

}

if(margin2 == "tPIG"){

pdf2 <- dPIGtr(y2, mu = exp(eta2), sigma = sigma, left.trunc)  

}








if(margin2 == "DEL"){

pdf2 <- dDEL(y2, mu = exp(eta2), sigma = sigma, nu = nu)  

}

if(margin2 == "SICHEL"){

pdf2 <- dSICHEL(y2, mu = exp(eta2), sigma = sigma) 

}





if(margin2 == "P"){



mu2 <- exp(eta2)

pdf2 <- dPO(y2, mu = mu2)   


           
}



if(margin2 == "tP"){


mu2 <- exp(eta2)

pdf2 <- dPOtr(y2, mu = mu2, left.trunc) 

}



if(margin2 == "DGP0"){

mu2 <- c(exp(eta2))
     
if(max(y2) > 170){

prec <- pmax(53, getPrec(mu2), getPrec(y2))
        
mu2 <- mpfr(mu2, prec)
y2  <- mpfr( y2, prec)        
        
} 


pdf2FUNC <- function(y2, mu2) exp(-y2/mu2) - exp(-(y2+1)/mu2)  
pdf2     <- as.numeric( pdf2FUNC(y2, mu2) ) 



}




if(margin2 %in% c("DGP","DGPII")){

if(margin2 == "DGP")   mu2 <- c(eta2)
if(margin2 == "DGPII") mu2 <- c( exp(eta2) )  # mu2 <- c(eta2^2)

sigma <- c(sigma)

     
pdf2FUNC2 <- function(y2, mu2, sigma) suppressWarnings(    (1 + mu2*y2/sigma)^(-1/mu2) - (1 + mu2*(1+y2)/sigma)^(-1/mu2)     )
pdf2     <-     as.numeric( pdf2FUNC2(y2, mu2, sigma) )     


# done but exclusions not implemented/checked


}


pdf2 <- ifelse(pdf2 < min.dn, min.dn, pdf2 )


list(pdf2 = ifef(pdf2))     


}




     























