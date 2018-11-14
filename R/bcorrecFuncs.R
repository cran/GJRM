intB <- function(y, eta, sigma2, sigma2.st, nu, nu.st, margin, rc, discr = FALSE, ym = NULL){ 
                                   
 if(discr == FALSE) pdf <- distrHsAT(y, eta, sigma2, nu, margin)$pdf2
 if(discr == TRUE)  pdf <- distrHsDiscr(y, eta, sigma2, 1, 1, 1, margin, naive = TRUE, ym)$pdf2

 log( 1 + exp( log( pdf ) + rc ) )
   
}

gradBbit1 <- function(y, eta, sigma2, sigma2.st, nu, nu.st, margin, rc, discr = FALSE, ym = NULL){ 

       if(discr == FALSE) dHs <- distrHs(y, eta, sigma2, sigma2.st, nu, nu.st, margin, naive = TRUE)
       if(discr == TRUE)  dHs <- distrHsDiscr(y, eta, sigma2, sigma2.st, 1, 1, margin, naive = TRUE, ym)

       pdf2            <- dHs$pdf2
       derpdf2.dereta2 <- dHs$derpdf2.dereta2 

       comp1 <- 1 + exp(log( pdf2 ) + rc) 
       comp2 <- pdf2/comp1

       dl.dbe <- derpdf2.dereta2/pdf2

       comp2*dl.dbe
     
}

gradBbit2 <- function(y, eta, sigma2, sigma2.st, nu, nu.st, margin, rc, discr = FALSE, ym = NULL){ 

       if(discr == FALSE) dHs <- distrHs(y, eta, sigma2, sigma2.st, nu, nu.st, margin, naive = TRUE)
       if(discr == TRUE)  dHs <- distrHsDiscr(y, eta, sigma2, sigma2.st, 1, 1, margin, naive = TRUE, ym)

       pdf2                 <- dHs$pdf2
       derpdf2.dersigma2.st <- dHs$derpdf2.dersigma2.st  

       comp1 <- 1 + exp(log( pdf2 ) + rc) 
       comp2 <- pdf2/comp1

       dl.dsigma.st <- derpdf2.dersigma2.st/pdf2

       comp2*dl.dsigma.st
       
}

gradBbit3 <- function(y, eta, sigma2, sigma2.st, nu, nu.st, margin, rc, discr = FALSE, ym = NULL){ 

       dHs <- distrHs(y, eta, sigma2, sigma2.st, nu, nu.st, margin, naive = TRUE)

       pdf2                 <- dHs$pdf2
       derpdf2.dernu.st     <- dHs$derpdf2.dernu.st  

       comp1 <- 1 + exp(log( pdf2 ) + rc) 
       comp2 <- pdf2/comp1

       dl.dnu.st <- derpdf2.dernu.st/pdf2

       comp2*dl.dnu.st
       
       }
       
hessBbit1 <- function(y, eta, sigma2, sigma2.st, nu, nu.st, margin, rc, discr = FALSE, ym = NULL){ 

       if(discr == FALSE) dHs <- distrHs(y, eta, sigma2, sigma2.st, nu, nu.st, margin, naive = TRUE)
       if(discr == TRUE)  dHs <- distrHsDiscr(y, eta, sigma2, sigma2.st, 1, 1, margin, naive = TRUE, ym)

        pdf2                         <- dHs$pdf2
        derpdf2.dereta2              <- dHs$derpdf2.dereta2
        der2pdf2.dereta2             <- dHs$der2pdf2.dereta2
        
        comp1 <- 1 + exp(log(pdf2) + rc)
        comp2 <- pdf2/comp1
        comp3 <- pdf2/comp1^2
        
        d2l.be.be <- (der2pdf2.dereta2 * pdf2 - (derpdf2.dereta2)^2)/pdf2^2
        dl.dbe    <- derpdf2.dereta2/pdf2
        
        comp2*d2l.be.be + dl.dbe^2*comp3

} 

hessBbit2 <- function(y, eta, sigma2, sigma2.st, nu, nu.st, margin, rc, discr = FALSE, ym = NULL){ 

       if(discr == FALSE) dHs <- distrHs(y, eta, sigma2, sigma2.st, nu, nu.st, margin, naive = TRUE)
       if(discr == TRUE)  dHs <- distrHsDiscr(y, eta, sigma2, sigma2.st, 1, 1, margin, naive = TRUE, ym)

        pdf2                         <- dHs$pdf2
        derpdf2.dersigma2.st         <- dHs$derpdf2.dersigma2.st
        der2pdf2.dersigma2.st2       <- dHs$der2pdf2.dersigma2.st2
        
        comp1 <- 1 + exp(log(pdf2) + rc)
        comp2 <- pdf2/comp1
        comp3 <- pdf2/comp1^2
        
        d2l.sigma.sigma <- (der2pdf2.dersigma2.st2 * pdf2 - (derpdf2.dersigma2.st)^2)/pdf2^2
        dl.dsigma.st    <- derpdf2.dersigma2.st/pdf2

        comp2*d2l.sigma.sigma + dl.dsigma.st^2*comp3
   
}

hessBbit3 <- function(y, eta, sigma2, sigma2.st, nu, nu.st, margin, rc, discr = FALSE, ym = NULL){ 

       if(discr == FALSE) dHs <- distrHs(y, eta, sigma2, sigma2.st, nu, nu.st, margin, naive = TRUE)
       if(discr == TRUE)  dHs <- distrHsDiscr(y, eta, sigma2, sigma2.st, 1, 1, margin, naive = TRUE, ym)

        pdf2                         <- dHs$pdf2
        derpdf2.dereta2              <- dHs$derpdf2.dereta2
        derpdf2.dersigma2.st         <- dHs$derpdf2.dersigma2.st
        der2pdf2.dereta2dersigma2.st <- dHs$der2pdf2.dereta2dersigma2.st
        
        comp1 <- 1 + exp(log(pdf2) + rc)
        comp2 <- pdf2/comp1
        comp3 <- pdf2/comp1^2
        
        d2l.be.sigma <- (der2pdf2.dereta2dersigma2.st * pdf2 - derpdf2.dereta2 * derpdf2.dersigma2.st)/pdf2^2
        
        dl.dbe       <- derpdf2.dereta2/pdf2
        dl.dsigma.st <- derpdf2.dersigma2.st/pdf2

        comp2*d2l.be.sigma + dl.dbe*dl.dsigma.st*comp3
    
}

hessBbit4 <- function(y, eta, sigma2, sigma2.st, nu, nu.st, margin, rc, discr = FALSE, ym = NULL){ 

        dHs <- distrHs(y, eta, sigma2, sigma2.st, nu, nu.st, margin, naive = TRUE)

        pdf2               <- dHs$pdf2
        derpdf2.dernu.st   <- dHs$derpdf2.dernu.st           
        der2pdf2.dernu.st2 <- dHs$der2pdf2.dernu.st2           
        
        comp1 <- 1 + exp(log(pdf2) + rc)
        comp2 <- pdf2/comp1
        comp3 <- pdf2/comp1^2
        
        d2l.nu.nu <- (der2pdf2.dernu.st2*pdf2-(derpdf2.dernu.st)^2)/pdf2^2
        dl.dnu.st <- derpdf2.dernu.st/pdf2

        comp2*d2l.nu.nu + dl.dnu.st^2*comp3
 
}

hessBbit5 <- function(y, eta, sigma2, sigma2.st, nu, nu.st, margin, rc, discr = FALSE, ym = NULL){ 

        dHs <- distrHs(y, eta, sigma2, sigma2.st, nu, nu.st, margin, naive = TRUE)

        pdf2                     <- dHs$pdf2
        derpdf2.dereta2          <- dHs$derpdf2.dereta2
        der2pdf2.dereta2dernu.st <- dHs$der2pdf2.dereta2dernu.st   
        derpdf2.dernu.st         <- dHs$derpdf2.dernu.st           
        
        comp1 <- 1 + exp(log(pdf2) + rc)
        comp2 <- pdf2/comp1
        comp3 <- pdf2/comp1^2
        
        d2l.be.nu  <- (der2pdf2.dereta2dernu.st*pdf2 - derpdf2.dereta2*derpdf2.dernu.st)/pdf2^2 
        dl.dbe     <- derpdf2.dereta2/pdf2
        dl.dnu.st  <- derpdf2.dernu.st/pdf2

        comp2*d2l.be.nu + dl.dbe*dl.dnu.st*comp3

}

hessBbit6 <- function(y, eta, sigma2, sigma2.st, nu, nu.st, margin, rc, discr = FALSE, ym = NULL){ 

        dHs <- distrHs(y, eta, sigma2, sigma2.st, nu, nu.st, margin, naive = TRUE)

        pdf2                          <- dHs$pdf2
        derpdf2.dersigma2.st          <- dHs$derpdf2.dersigma2.st
        der2pdf2.dersigma2.stdernu.st <- dHs$der2pdf2.sigma2.st2dernu.st
        derpdf2.dernu.st              <- dHs$derpdf2.dernu.st           
        
        comp1 <- 1 + exp(log(pdf2) + rc)
        comp2 <- pdf2/comp1
        comp3 <- pdf2/comp1^2
        
        d2l.sigma.nu  <- (der2pdf2.dersigma2.stdernu.st*pdf2-(derpdf2.dersigma2.st*derpdf2.dernu.st))/pdf2^2 

        dl.dsigma.st  <- derpdf2.dersigma2.st/pdf2
        dl.dnu.st     <- derpdf2.dernu.st/pdf2

        comp2*d2l.sigma.nu + dl.dsigma.st*dl.dnu.st*comp3
       
}


###################################################################################################
###################################################################################################

int1f <- function(y, eta, sigma2, nu, margin, rc){ 
   pdf <- distrHsAT(y, eta, sigma2, nu, margin)$pdf2
   log( 1 + exp( log( pdf ) + rc ) )
}




d.bpsi <- function(y, X1, X2, X3, eta, sigma2, sigma2.st, nu, nu.st, margin, rc, j){ 

       dHs <- distrHs(y, eta, sigma2, sigma2.st, nu, nu.st, margin, naive = TRUE)

       pdf2                 <- dHs$pdf2
       derpdf2.dereta2      <- dHs$derpdf2.dereta2 
       derpdf2.dersigma2.st <- dHs$derpdf2.dersigma2.st  
       derpdf2.dernu.st     <- dHs$derpdf2.dernu.st  

       comp1 <- 1 + exp(log( pdf2 ) + rc) 
       comp2 <- pdf2/comp1

       dl.dbe       <- derpdf2.dereta2/pdf2
       dl.dsigma.st <- derpdf2.dersigma2.st/pdf2
       dl.dnu.st    <- derpdf2.dernu.st/pdf2


       if( margin %in% c("DAGUM","SM") ) res <- cbind( comp2*as.numeric(dl.dbe)%*%t(X1), comp2*as.numeric(dl.dsigma.st)%*%t(X2), comp2*as.numeric(dl.dnu.st)%*%t(X3) ) else
                                         res <- cbind( comp2*as.numeric(dl.dbe)%*%t(X1), comp2*as.numeric(dl.dsigma.st)%*%t(X2) )
      
       res[, j]
}




gradF <- function(params, n, VC, margin, lB, uB, rc){

  G <- matrix(NA, n, length(params))

for(i in 1:n){

  X1 <- VC$X1[i,]
  X2 <- VC$X2[i,]
  X3 <- VC$X3[i,]
  nu <- nu.st <- 1
  
  eta       <- X1%*%params[1:VC$X1.d2]
  sigma2.st <- X2%*%params[(1+VC$X1.d2):(VC$X1.d2+VC$X2.d2)]
  
     if( margin %in% c("DAGUM","SM") ){ 
  
       nu.st <- X3%*%params[(1+VC$X1.d2+VC$X2.d2):(VC$X1.d2+VC$X2.d2+VC$X3.d2)] 
 
       ss <- enu.tr(nu.st, margin)  
       nu.st <- ss$vrb.st 
       nu    <- ss$vrb   
                                      }
 
  ss        <- esp.tr(sigma2.st, margin)
  sigma2.st <- ss$vrb.st
  sigma2    <- ss$vrb
                                    
  for(j in 1:length(params)) G[i, j] <- integrate(d.bpsi, lB, uB, X1, X2, X3, eta, sigma2, sigma2.st, nu, nu.st, margin, rc, j)$value

             }
             
  colSums(G)
  
}


