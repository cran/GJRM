bcontSurv <- function(params, respvec, VC, ps, AT = FALSE){

    eta1 <- VC$X1%*%params[1:VC$X1.d2]
    eta2 <- VC$X2%*%params[(VC$X1.d2 + 1):(VC$X1.d2 + VC$X2.d2)]
    etad <- etas1 <- etas2 <- l.ln <- NULL 
  
    epsilon <- 0.0000001 
    max.p   <- 0.9999999
    
    
  if(is.null(VC$X3)){  
    sigma21.st <- etas1 <- params[(VC$X1.d2 + VC$X2.d2 + 1)]
    sigma22.st <- etas2 <- params[(VC$X1.d2 + VC$X2.d2 + 2)]
    teta.st    <- etad  <- params[(VC$X1.d2 + VC$X2.d2 + 3)]
  } 
  
  
  if(!is.null(VC$X3)){  
    sigma21.st <- etas1 <- VC$X3%*%params[(VC$X1.d2 + VC$X2.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2)]
    sigma22.st <- etas2 <- VC$X4%*%params[(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2)]
    teta.st    <- etad  <- VC$X5%*%params[(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2 + VC$X5.d2)]
  }  
  
  
##################
## Transformations
##################
  
sstr1 <- esp.tr(sigma21.st, VC$margins[1])  
sstr2 <- esp.tr(sigma22.st, VC$margins[2])  
  
sigma21.st <- sstr1$vrb.st 
sigma22.st <- sstr2$vrb.st 
    
sigma21    <- sstr1$vrb 
sigma22    <- sstr2$vrb 
    
eta1 <- eta.tr(eta1, VC$margins[1])
eta2 <- eta.tr(eta2, VC$margins[2])
    
resT  <- teta.tr(VC, teta.st)

teta.st1 <- teta.st2 <- teta.st <- resT$teta.st
teta1 <- teta2 <- teta <- resT$teta 

##################

Cop1 <- Cop2 <- VC$BivD 
nC1 <- nC2 <- VC$nC

teta.ind1 <- as.logical(c(1,0,round(runif(VC$n-2))) ) 
teta.ind2 <- teta.ind1 == FALSE 


if(!(VC$BivD %in% VC$BivD2) && length(teta.st) > 1){

teta.st1 <- teta.st[teta.ind1]
teta.st2 <- teta.st[teta.ind2]

teta1 <- teta[teta.ind1]
teta2 <- teta[teta.ind2]

}


 
if(VC$BivD %in% VC$BivD2){

if(VC$BivD %in% VC$BivD2[1:4])  teta.ind1 <- ifelse(VC$my.env$signind*teta > exp(VC$zerov), TRUE, FALSE)
if(VC$BivD %in% VC$BivD2[5:12]) teta.ind1 <- ifelse(VC$my.env$signind*teta > exp(VC$zerov) + 1, TRUE, FALSE) 
teta.ind2 <- teta.ind1 == FALSE 

VC$my.env$signind <- ifelse(teta.ind1 == TRUE,  1, -1) 

teta1 <-  teta[teta.ind1]
teta2 <- -teta[teta.ind2]

teta.st1 <- teta.st[teta.ind1]
teta.st2 <- teta.st[teta.ind2]

if(length(teta) == 1) teta.ind2 <- teta.ind1 <- rep(TRUE, VC$n)  

Cop1Cop2R <- Cop1Cop2(VC$BivD)
Cop1 <- Cop1Cop2R$Cop1
Cop2 <- Cop1Cop2R$Cop2

nC1 <- VC$ct[which(VC$ct[,1] == Cop1),2] 
nC2 <- VC$ct[which(VC$ct[,1] == Cop2),2]

} 


##################

  dHs1 <- distrHs(respvec$y1, eta1, sigma21, sigma21.st, nu = 1, nu.st = 1, margin2=VC$margins[1], naive = FALSE)
  dHs2 <- distrHs(respvec$y2, eta2, sigma22, sigma22.st, nu = 1, nu.st = 1, margin2=VC$margins[2], naive = FALSE)

  pdf1 <- dHs1$pdf2
  pdf2 <- dHs2$pdf2

  p1 <- dHs1$p2
  p2 <- dHs2$p2
  

  if( length(teta1) != 0) dH1 <- copgHs(1-p1[teta.ind1], 1-p2[teta.ind1], eta1=NULL, eta2=NULL, teta1, teta.st1, Cop1, VC$dof, VC$dof)
  if( length(teta2) != 0) dH2 <- copgHs(1-p1[teta.ind2], 1-p2[teta.ind2], eta1=NULL, eta2=NULL, teta2, teta.st2, Cop2, VC$dof, VC$dof)
  c.copula2.be1be2 <- c.copula.be1 <- c.copula.be2 <- p00 <- c.copula.theta <- NA
  
  if( length(teta1) != 0){ c.copula2.be1be2[teta.ind1] <- dH1$c.copula2.be1be2
                           c.copula.be1[teta.ind1]     <- dH1$c.copula.be1
                           c.copula.be2[teta.ind1]     <- dH1$c.copula.be2
                           c.copula.theta[teta.ind1]   <- dH1$c.copula.theta
                           p00[teta.ind1] <- BiCDF(1-p1[teta.ind1], 1-p2[teta.ind1], nC1, teta1, VC$dof) 
                         }
  
  if( length(teta2) != 0){ c.copula2.be1be2[teta.ind2] <- dH2$c.copula2.be1be2  
                           c.copula.be1[teta.ind2]     <- dH2$c.copula.be1
                           c.copula.be2[teta.ind2]     <- dH2$c.copula.be2
                           c.copula.theta[teta.ind2]   <- dH2$c.copula.theta
                           p00[teta.ind2] <- BiCDF(1-p1[teta.ind2], 1-p2[teta.ind2], nC2, teta2, VC$dof) 
                         }
 


l.par <- VC$weights*( VC$c11*( log(c.copula2.be1be2) + log(pdf1) + log(pdf2) ) + 
                      VC$c10*( log(c.copula.be1) + log(pdf1) ) + 
                      VC$c01*( log(c.copula.be2) + log(pdf2) ) + 
                      VC$c00*log(p00) 
                     )
 
 
##################
##################

# cont-cont

 derpdf1.dereta1              <- dHs1$derpdf2.dereta2 
 derpdf1.dersigma21.st        <- dHs1$derpdf2.dersigma2.st 
 
 derpdf2.dereta2              <- dHs2$derpdf2.dereta2 
 derpdf2.dersigma22.st        <- dHs2$derpdf2.dersigma2.st 
 
 derp1.dereta1                <- dHs1$derp2.dereta2
 derp1.dersigma21.st          <- dHs1$derp2.dersigma.st
 
 derp2.dereta2                <- dHs2$derp2.dereta2
 derp2.dersigma22.st          <- dHs2$derp2.dersigma.st 
 
 
 if( length(teta1) != 0) BITS1 <- copgHsCont(1-p1[teta.ind1], 1-p2[teta.ind1], teta1, teta.st1, Cop1, Cont = TRUE, VC$dof, VC$dof)
 if( length(teta2) != 0) BITS2 <- copgHsCont(1-p1[teta.ind2], 1-p2[teta.ind2], teta2, teta.st2, Cop2, Cont = TRUE, VC$dof, VC$dof) 
 

   der2h.derp1p1 <- NA
   if( length(teta1) != 0) der2h.derp1p1[teta.ind1]   <- BITS1$der2h.derp1p1
   if( length(teta2) != 0) der2h.derp1p1[teta.ind2]   <- BITS2$der2h.derp1p1 
   
   derc.dereta1               <- der2h.derp1p1 * derp1.dereta1 
   derc.dersigma21.st         <- der2h.derp1p1 * derp1.dersigma21.st
   
   der2h.derp1p2 <- NA
   if( length(teta1) != 0) der2h.derp1p2[teta.ind1] <- BITS1$der2h.derp1p2 
   if( length(teta2) != 0) der2h.derp1p2[teta.ind2] <- BITS2$der2h.derp1p2 
   
   
   derc.dereta2               <- der2h.derp1p2 * derp2.dereta2    
   derc.dersigma22.st         <- der2h.derp1p2 * derp2.dersigma22.st



   der2h.derp1teta            <- NA
   derteta.derteta.st         <- NA
   if( length(teta1) != 0) der2h.derp1teta[teta.ind1]    <- BITS1$der2h.derp1teta
   if( length(teta2) != 0) der2h.derp1teta[teta.ind2]    <- BITS2$der2h.derp1teta   
   if( length(teta1) != 0) derteta.derteta.st[teta.ind1] <- BITS1$derteta.derteta.st
   if( length(teta2) != 0) derteta.derteta.st[teta.ind2] <- BITS2$derteta.derteta.st 


   der2h.derp1teta.st         <- der2h.derp1teta * derteta.derteta.st # new bit
   
#################
#################
   
c.copula2.be1 <- c.copula2.be2 <- c.copula2.be1th <- c.copula2.be2th <- bit1.th2 <- c.copula2.be1t <- c.copula2.be2t <- NA

if( length(teta1) != 0){
 
c.copula2.be1[teta.ind1]    <- dH1$c.copula2.be1   
c.copula2.be2[teta.ind1]    <- dH1$c.copula2.be2 
c.copula2.be1th[teta.ind1]  <- dH1$c.copula2.be1th 
c.copula2.be2th[teta.ind1]  <- dH1$c.copula2.be2th
c.copula2.be1t[teta.ind1]   <- dH1$c.copula2.be1t 
c.copula2.be2t[teta.ind1]   <- dH1$c.copula2.be2t
bit1.th2[teta.ind1]         <- dH1$bit1.th2

}



if( length(teta2) != 0){

c.copula2.be1[teta.ind2]    <- dH2$c.copula2.be1   
c.copula2.be2[teta.ind2]    <- dH2$c.copula2.be2 
c.copula2.be1th[teta.ind2]  <- dH2$c.copula2.be1th 
c.copula2.be2th[teta.ind2]  <- dH2$c.copula2.be2th
c.copula2.be1t[teta.ind2]   <- dH2$c.copula2.be1t 
c.copula2.be2t[teta.ind2]   <- dH2$c.copula2.be2t
bit1.th2[teta.ind2]         <- dH2$bit1.th2

}

#################
#################
   

 
   dl.dbe1 <- VC$weights*(   
   
   VC$c11*(  derpdf1.dereta1/pdf1 + -derc.dereta1/c.copula2.be1be2  ) +
   
   VC$c10*( c.copula.be1^-1 * -c.copula2.be1 * derp1.dereta1 + derpdf1.dereta1 * pdf1^-1) + 
   
   VC$c01*( c.copula.be2^-1 * -c.copula2.be1be2 * derp1.dereta1   ) +
   
   VC$c00*(  -derp1.dereta1*c.copula.be1/p00  )
   
   )
   
   
   
   
   
   dl.dbe2 <- VC$weights*(   
   
   VC$c11*(  derpdf2.dereta2/pdf2 + -derc.dereta2/c.copula2.be1be2  ) +
   
   VC$c10*( c.copula.be1^-1 * -c.copula2.be1be2 * derp2.dereta2   ) +
   
   VC$c01*( c.copula.be2^-1 * -c.copula2.be2 * derp2.dereta2 + derpdf2.dereta2 * pdf2^-1   ) +
   
   VC$c00*(  -derp2.dereta2*c.copula.be2/p00  )           
   
   )
   
   
   
   dl.dsigma21.st <- VC$weights*(   
   
   VC$c11*(  derpdf1.dersigma21.st/pdf1 + -derc.dersigma21.st/c.copula2.be1be2  ) +
   
   VC$c10*( c.copula.be1^-1 * -c.copula2.be1 * derp1.dersigma21.st + derpdf1.dersigma21.st * pdf1^-1   ) +
   
   VC$c01*( c.copula.be2^-1 * -c.copula2.be1be2 * derp1.dersigma21.st   ) +
   
   VC$c00*(  -derp1.dersigma21.st*c.copula.be1/p00  ) 
   
   )
   
   
   
   dl.dsigma22.st <- VC$weights*(   
   
   VC$c11*(  derpdf2.dersigma22.st/pdf2 + -derc.dersigma22.st/c.copula2.be1be2  ) +
   
   VC$c10*( c.copula.be1^-1 * -c.copula2.be1be2 * derp2.dersigma22.st   ) +
   
   VC$c01*( c.copula.be2^-1 * -c.copula2.be2 * derp2.dersigma22.st + derpdf2.dersigma22.st * pdf2^-1    ) +
   
   VC$c00*(  -derp2.dersigma22.st*c.copula.be2/p00  ) 
   
   )  
   

   dl.dteta.st    <- VC$weights*(   
   
   VC$c11*(  der2h.derp1teta.st/c.copula2.be1be2  ) +
   
   VC$c10*( c.copula.be1^-1 * c.copula2.be1th   ) + 
   
   VC$c01*( c.copula.be2^-1 * c.copula2.be2th   ) +
   
   VC$c00*(  c.copula.theta/p00  )
   
   )
  
             
#################################################################################################
# cont-cont


der2c.derrho.derrho    <- NA
der2c.derp1.derp1      <- NA
der2c.derp2.derp2      <- NA
der2c.derp1.derp2      <- NA 
der2c.derp1.derrho     <- NA
der2c.derp2.derrho     <- NA
der2teta.derteta.stteta.st <- NA 


if( length(teta1) != 0){ der2c.derrho.derrho[teta.ind1] <- BITS1$der2c.derrho.derrho
 der2c.derp1.derp1[teta.ind1]      <- BITS1$der2c.derp1.derp1  
 der2c.derp2.derp2[teta.ind1]      <- BITS1$der2c.derp2.derp2  
 der2c.derp1.derp2[teta.ind1]      <- BITS1$der2c.derp1.derp2  
 der2c.derp1.derrho[teta.ind1]     <- BITS1$der2c.derp1.derrho 
 der2c.derp2.derrho[teta.ind1]     <- BITS1$der2c.derp2.derrho }


if( length(teta2) != 0){ der2c.derrho.derrho[teta.ind2] <- BITS2$der2c.derrho.derrho
 der2c.derp1.derp1[teta.ind2]      <- BITS2$der2c.derp1.derp1  
 der2c.derp2.derp2[teta.ind2]      <- BITS2$der2c.derp2.derp2  
 der2c.derp1.derp2[teta.ind2]      <- BITS2$der2c.derp1.derp2  
 der2c.derp1.derrho[teta.ind2]     <- BITS2$der2c.derp1.derrho 
 der2c.derp2.derrho[teta.ind2]     <- BITS2$der2c.derp2.derrho }

if( length(teta1) != 0) der2teta.derteta.stteta.st[teta.ind1] <- BITS1$der2teta.derteta.stteta.st
if( length(teta2) != 0) der2teta.derteta.stteta.st[teta.ind2] <- BITS2$der2teta.derteta.stteta.st 



der2pdf1.dereta1 <- dHs1$der2pdf2.dereta2
der2pdf2.dereta2 <- dHs2$der2pdf2.dereta2

der2pdf1.dersigma21.st2  <- dHs1$der2pdf2.dersigma2.st2
der2pdf2.dersigma22.st2  <- dHs2$der2pdf2.dersigma2.st2

der2p1.dereta1eta1 <- dHs1$der2p2.dereta2eta2
der2p2.dereta2eta2 <- dHs2$der2p2.dereta2eta2
 
der2p1.dersigma21.st2 <-  dHs1$der2p2.dersigma2.st2
der2p2.dersigma22.st2 <-  dHs2$der2p2.dersigma2.st2
 
der2pdf1.dereta1dersigma21.st <-  dHs1$der2pdf2.dereta2dersigma2.st
der2pdf2.dereta2dersigma22.st <-  dHs2$der2pdf2.dereta2dersigma2.st
 
der2p1.dereta1dersigma21.st <-  dHs1$der2p2.dereta2dersigma2.st
der2p2.dereta2dersigma22.st <-  dHs2$der2p2.dereta2dersigma2.st


########################
########################


bit1.b1b1 <- c.copula2.be1*derp1.dereta1^2 + -c.copula.be1*der2p1.dereta1eta1   
bit1.s1s1 <- c.copula2.be1*derp1.dersigma21.st^2 + -c.copula.be1*der2p1.dersigma21.st2
bit1.b1s1 <- c.copula2.be1*derp1.dereta1*derp1.dersigma21.st + -c.copula.be1*der2p1.dereta1dersigma21.st                                                                                                                                

bit1.b2s1 <- c.copula2.be1be2*derp2.dereta2*derp1.dersigma21.st
bit1.b1s2 <- c.copula2.be1be2*derp1.dereta1*derp2.dersigma22.st

bit1.b2b2 <- c.copula2.be2*derp2.dereta2^2 + -c.copula.be2*der2p2.dereta2eta2 
bit1.s2s2 <- c.copula2.be2*derp2.dersigma22.st^2 + -c.copula.be2*der2p2.dersigma22.st2
bit1.b2s2 <- c.copula2.be2*derp2.dereta2*derp2.dersigma22.st + -c.copula.be2*der2p2.dereta2dersigma22.st

bit1.b1b2 <- c.copula2.be1be2*derp1.dereta1*derp2.dereta2

bit1.b1th <- -c.copula2.be1th*derp1.dereta1

bit1.b2th <- -c.copula2.be2th*derp2.dereta2

bit1.s1th <- -c.copula2.be1th*derp1.dersigma21.st
bit1.s2th <- -c.copula2.be2th*derp2.dersigma22.st

bit1.s1s2 <- c.copula2.be1be2*derp1.dersigma21.st*derp2.dersigma22.st

      
########################
########################


der3C.derp1p1p1 <- der3C.derp1tetateta <- der2h.derteta.teta.st <- der3C.p1p1teta <- der2h.derp2teta <- der2h.derp2p2 <- NA


if( length(teta1) != 0){der3C.derp1p1p1[teta.ind1]       <- BITS1$der3C.derp1p1p1 
                        der2h.derteta.teta.st[teta.ind1] <- BITS1$der2h.derteta.teta.st
                        der3C.derp1tetateta[teta.ind1]   <- BITS1$der3C.derp1tetateta
                        der3C.p1p1teta[teta.ind1]        <- BITS1$der3C.p1p1teta  
                        der2h.derp2teta[teta.ind1]       <- BITS1$der2h.derp2teta
                        der2h.derp2p2[teta.ind1]         <- BITS1$der2h.derp2p2
                       }


if( length(teta2) != 0){der3C.derp1p1p1[teta.ind2]       <- BITS2$der3C.derp1p1p1
                        der2h.derteta.teta.st[teta.ind2] <- BITS2$der2h.derteta.teta.st
                        der3C.derp1tetateta[teta.ind2]   <- BITS2$der3C.derp1tetateta
                        der3C.p1p1teta[teta.ind2]        <- BITS2$der3C.p1p1teta
                        der2h.derp2teta[teta.ind2]       <- BITS2$der2h.derp2teta
                        der2h.derp2p2[teta.ind2]         <- BITS2$der2h.derp2p2
                       }


########################
########################                                                                                                 der2h.derp1p1 * derp1.dereta1
    
          
  d2l.be1.be1      <-  -VC$weights*( 
  
   VC$c11*(  (der2pdf1.dereta1 * pdf1 - derpdf1.dereta1^2) / pdf1^2  + 
                       ((der2c.derp1.derp1 * derp1.dereta1^2 + der2h.derp1p1 * -der2p1.dereta1eta1) * c.copula2.be1be2 + -derc.dereta1^2) /c.copula2.be1be2^2   ) +
   
   VC$c10*( -(c.copula.be1)^-2 * c.copula2.be1^2 * derp1.dereta1^2 + c.copula.be1^-1 * der3C.derp1p1p1 * derp1.dereta1^2 + c.copula.be1^-1 * c.copula2.be1 * - der2p1.dereta1eta1 +
            der2pdf1.dereta1 * pdf1^-1 +  derpdf1.dereta1^2*-pdf1^-2 ) +
   
   VC$c01*( -(c.copula.be2)^-2 * c.copula2.be1be2^2 * derp1.dereta1^2 + c.copula.be2^-1 * der2h.derp1p1  * derp1.dereta1^2 +  c.copula.be2^-1 * c.copula2.be1be2 * -der2p1.dereta1eta1   ) +
   
   VC$c00*(  (bit1.b1b1*p00 - (-c.copula.be1*derp1.dereta1)^2)/p00^2  ) 
  
  ) 
  
                   
                       
  d2l.be2.be2  <-  -VC$weights*( 
 
    VC$c11*(  (der2pdf2.dereta2 * pdf2 - derpdf2.dereta2^2) / pdf2^2   + 
                       ((der2c.derp2.derp2 * derp2.dereta2^2 + der2h.derp1p2 * -der2p2.dereta2eta2) * c.copula2.be1be2 - derc.dereta2^2) /c.copula2.be1be2^2  ) +
    
    VC$c10*( -(c.copula.be1)^-2 * c.copula2.be1be2^2 * derp2.dereta2^2 + c.copula.be1^-1 * der2h.derp1p2  * derp2.dereta2^2 +  c.copula.be1^-1 * c.copula2.be1be2 * -der2p2.dereta2eta2) +
    
    VC$c01*( -(c.copula.be2)^-2 * c.copula2.be2^2 * derp2.dereta2^2 + c.copula.be2^-1 * der2h.derp2p2 * derp2.dereta2^2 + c.copula.be2^-1 * c.copula2.be2 * - der2p2.dereta2eta2 +
            der2pdf2.dereta2 * pdf2^-1 +  derpdf2.dereta2^2*-pdf2^-2   ) +
    
    VC$c00*(  (bit1.b2b2*p00 - (-c.copula.be2*derp2.dereta2)^2)/p00^2  )
  
  )
       
       
       
 
 
 
####################################################
####################################################

if( !(VC$BivD %in% VC$BivD2) ){ 

if(VC$BivD %in% c("C90","C180","C270","J90","J180","J270","G90","G180","G270")){ 

    bit2rho <- VC$c10*( -(c.copula.be1)^-2 * c.copula2.be1th^2 + c.copula.be1^-1 * der3C.derp1tetateta *  derteta.derteta.st^2    - c.copula.be1^-1 * c.copula2.be1t * der2teta.derteta.stteta.st)    
    bit3rho <- VC$c01*( -(c.copula.be2)^-2 * c.copula2.be2th^2 + c.copula.be2^-1 * der2h.derteta.teta.st  *  derteta.derteta.st^2 - c.copula.be2^-1 * c.copula2.be2t * der2teta.derteta.stteta.st)
 
                                                                               }else{
                                                                               
    bit2rho <- VC$c10*( -(c.copula.be1)^-2 * c.copula2.be1th^2 + c.copula.be1^-1 * der3C.derp1tetateta *  derteta.derteta.st^2    + c.copula.be1^-1 * c.copula2.be1t * der2teta.derteta.stteta.st)    
    bit3rho <- VC$c01*( -(c.copula.be2)^-2 * c.copula2.be2th^2 + c.copula.be2^-1 * der2h.derteta.teta.st  *  derteta.derteta.st^2 + c.copula.be2^-1 * c.copula2.be2t * der2teta.derteta.stteta.st)
                                                                                
                                                                               }                                                                               
}  


if( VC$BivD %in% VC$BivD2 ){ 

bit2rho <- bit3rho <- NA

if(Cop1 %in% c("C0","J0","G0")){ 
                                                           
    bit2rho[teta.ind1] <- VC$c10[teta.ind1]*( -(c.copula.be1[teta.ind1])^-2 * c.copula2.be1th[teta.ind1]^2 + c.copula.be1[teta.ind1]^-1 * der3C.derp1tetateta[teta.ind1] *  derteta.derteta.st[teta.ind1]^2    + c.copula.be1[teta.ind1]^-1 * c.copula2.be1t[teta.ind1] * der2teta.derteta.stteta.st[teta.ind1])    
    bit3rho[teta.ind1] <- VC$c01[teta.ind1]*( -(c.copula.be2[teta.ind1])^-2 * c.copula2.be2th[teta.ind1]^2 + c.copula.be2[teta.ind1]^-1 * der2h.derteta.teta.st[teta.ind1]  *  derteta.derteta.st[teta.ind1]^2 + c.copula.be2[teta.ind1]^-1 * c.copula2.be2t[teta.ind1] * der2teta.derteta.stteta.st[teta.ind1])
                                                                                
                               } 

if(Cop1 %in% c("C180","J180","G180")){ 
                                                                        
    bit2rho[teta.ind1] <- VC$c10[teta.ind1]*( -(c.copula.be1[teta.ind1])^-2 * c.copula2.be1th[teta.ind1]^2 + c.copula.be1[teta.ind1]^-1 * der3C.derp1tetateta[teta.ind1] *  derteta.derteta.st[teta.ind1]^2    - c.copula.be1[teta.ind1]^-1 * c.copula2.be1t[teta.ind1] * der2teta.derteta.stteta.st[teta.ind1])    
    bit3rho[teta.ind1] <- VC$c01[teta.ind1]*( -(c.copula.be2[teta.ind1])^-2 * c.copula2.be2th[teta.ind1]^2 + c.copula.be2[teta.ind1]^-1 * der2h.derteta.teta.st[teta.ind1]  *  derteta.derteta.st[teta.ind1]^2 - c.copula.be2[teta.ind1]^-1 * c.copula2.be2t[teta.ind1] * der2teta.derteta.stteta.st[teta.ind1])
                                                                                
                                     } 


if(Cop2 %in% c("C90","J90","G90","C270","J270","G270")){ 
                                                           
    bit2rho[teta.ind2] <- VC$c10[teta.ind2]*( -(c.copula.be1[teta.ind2])^-2 * c.copula2.be1th[teta.ind2]^2 + c.copula.be1[teta.ind2]^-1 * der3C.derp1tetateta[teta.ind2] *  derteta.derteta.st[teta.ind2]^2    - c.copula.be1[teta.ind2]^-1 * c.copula2.be1t[teta.ind2] * der2teta.derteta.stteta.st[teta.ind2])    
    bit3rho[teta.ind2] <- VC$c01[teta.ind2]*( -(c.copula.be2[teta.ind2])^-2 * c.copula2.be2th[teta.ind2]^2 + c.copula.be2[teta.ind2]^-1 * der2h.derteta.teta.st[teta.ind2]  *  derteta.derteta.st[teta.ind2]^2 - c.copula.be2[teta.ind2]^-1 * c.copula2.be2t[teta.ind2] * der2teta.derteta.stteta.st[teta.ind2])
                                                                                
                                                       } 
   
}                                                                               
   
                                                                               
                                                                               
                                                                               
############################ 
            
d2l.rho.rho <- -VC$weights*( 
         
 VC$c11*(  ((der2c.derrho.derrho*derteta.derteta.st^2 + der2h.derp1teta *der2teta.derteta.stteta.st)*c.copula2.be1be2 - der2h.derp1teta.st^2) /c.copula2.be1be2^2  ) +
            
 bit2rho +
            
 bit3rho + 
            
 VC$c00*(  (bit1.th2*p00  - (c.copula.theta)^2)/p00^2  )
          
          ) 

###################################################          
###################################################  

                                    
  d2l.sigma21.sigma21  <- -VC$weights*( 
  
   VC$c11*(  (der2pdf1.dersigma21.st2 * pdf1 - derpdf1.dersigma21.st^2) / pdf1^2   + 
                          ((der2c.derp1.derp1*derp1.dersigma21.st^2+der2h.derp1p1*-der2p1.dersigma21.st2)*c.copula2.be1be2 - derc.dersigma21.st^2) /c.copula2.be1be2^2  ) +
   
   VC$c10*( -(c.copula.be1)^-2 * c.copula2.be1^2 * derp1.dersigma21.st^2 + c.copula.be1^-1 * der3C.derp1p1p1 * derp1.dersigma21.st^2 + c.copula.be1^-1 * c.copula2.be1 * - der2p1.dersigma21.st2 +
            der2pdf1.dersigma21.st2 * pdf1^-1 +  derpdf1.dersigma21.st^2*-pdf1^-2    ) +
   
   VC$c01*( -(c.copula.be2)^-2 * c.copula2.be1be2^2 * derp1.dersigma21.st^2 + c.copula.be2^-1 * der2h.derp1p1  * derp1.dersigma21.st^2 +  c.copula.be2^-1 * c.copula2.be1be2 * -der2p1.dersigma21.st2   ) +
   
   VC$c00*(  (bit1.s1s1*p00 - (-c.copula.be1*derp1.dersigma21.st)^2)/p00^2  )  
  
  )  
  
  
  
  
  
         
  
   
   
   
  d2l.sigma22.sigma22  <- -VC$weights*( 

   VC$c11*(  (der2pdf2.dersigma22.st2 * pdf2 - derpdf2.dersigma22.st^2) / pdf2^2   + 
                          ((der2c.derp2.derp2*derp2.dersigma22.st^2+der2h.derp1p2*-der2p2.dersigma22.st2)*c.copula2.be1be2 - derc.dersigma22.st^2) /c.copula2.be1be2^2  ) +
   
   VC$c10*( -(c.copula.be1)^-2 * c.copula2.be1be2^2 * derp2.dersigma22.st^2 + c.copula.be1^-1 * der2h.derp1p2  * derp2.dersigma22.st^2 +  c.copula.be1^-1 * c.copula2.be1be2 * 
             -der2p2.dersigma22.st2   ) +
   
   VC$c01*( -(c.copula.be2)^-2 * c.copula2.be2^2 * derp2.dersigma22.st^2 + c.copula.be2^-1 * der2h.derp2p2 * derp2.dersigma22.st^2 + c.copula.be2^-1 * c.copula2.be2 * - der2p2.dersigma22.st2 +
            der2pdf2.dersigma22.st2 * pdf2^-1 +  derpdf2.dersigma22.st^2*-pdf2^-2     ) +
   
   VC$c00*(  (bit1.s2s2*p00 - (-c.copula.be2*derp2.dersigma22.st)^2)/p00^2  )  
  
  )

















  d2l.be1.be2    <- -VC$weights*( 


   VC$c11*(  (der2c.derp1.derp2*derp1.dereta1*derp2.dereta2*c.copula2.be1be2 - derc.dereta1*derc.dereta2) /c.copula2.be1be2^2  ) +
   
   VC$c10*( -(c.copula.be1)^-2 * c.copula2.be1be2 * c.copula2.be1 * derp1.dereta1 * derp2.dereta2 +  c.copula.be1^-1 * der2h.derp1p1 *  derp1.dereta1 * derp2.dereta2 ) +
   
   VC$c01*( -(c.copula.be2)^-2 * c.copula2.be1be2 * c.copula2.be2 * derp1.dereta1 * derp2.dereta2 +  c.copula.be2^-1 * der2h.derp1p2 *  derp1.dereta1 * derp2.dereta2  ) +
   
   VC$c00*( (bit1.b1b2*p00 - (c.copula.be1*derp1.dereta1*c.copula.be2*derp2.dereta2))/p00^2   )
  
  
  )
  
  
  
  
  
  
  

  
  d2l.be1.rho   <- -VC$weights*(


   VC$c11*(  (der2c.derp1.derrho *-derp1.dereta1*derteta.derteta.st* c.copula2.be1be2 - -derc.dereta1*der2h.derp1teta*derteta.derteta.st)  /c.copula2.be1be2^2  ) +
   
   VC$c10*(  -(c.copula.be1)^-2 * c.copula2.be1th *  c.copula2.be1 * - derp1.dereta1 + c.copula.be1^-1 * der3C.p1p1teta  *derteta.derteta.st*- derp1.dereta1) +
   
   VC$c01*(  -(c.copula.be2)^-2 * c.copula2.be2th *  c.copula2.be1be2 * - derp1.dereta1 + c.copula.be2^-1 * der2h.derp1teta *derteta.derteta.st*- derp1.dereta1  ) +
   
   VC$c00*(  (bit1.b1th*p00 - (-c.copula.be1*derp1.dereta1*c.copula.theta))/p00^2  )  
  
  )
  
  
  
  
  
  
  
  

    
  
  d2l.be2.rho   <- -VC$weights*( 
  
   VC$c11*(  (der2c.derp2.derrho *-derp2.dereta2*derteta.derteta.st* c.copula2.be1be2 - -derc.dereta2*der2h.derp1teta*derteta.derteta.st) /c.copula2.be1be2^2  ) +
   
   VC$c10*( -(c.copula.be1)^-2 * c.copula2.be1th *  c.copula2.be1be2 * - derp2.dereta2 + c.copula.be1^-1 * der2h.derp1teta *derteta.derteta.st*- derp2.dereta2  ) +
  
   VC$c01*(  -(c.copula.be2)^-2 * c.copula2.be2th *  c.copula2.be2 * - derp2.dereta2 + c.copula.be2^-1 * der2h.derp2teta *derteta.derteta.st*- derp2.dereta2   ) +
   
   VC$c00*(  (bit1.b2th*p00 - (-c.copula.be2*derp2.dereta2*c.copula.theta))/p00^2  ) 
  
  )
  
  
  
  
  
    
  




  
  d2l.be1.be1      <-  -VC$weights*( 
  
   VC$c11*(  (der2pdf1.dereta1 * pdf1 - derpdf1.dereta1^2) / pdf1^2  + 
                       ((der2c.derp1.derp1 * derp1.dereta1^2 + der2h.derp1p1 * -der2p1.dereta1eta1) * c.copula2.be1be2 + -derc.dereta1^2) /c.copula2.be1be2^2   ) +
   
   VC$c10*( -(c.copula.be1)^-2 * c.copula2.be1^2 * derp1.dereta1^2 + c.copula.be1^-1 * der3C.derp1p1p1 * derp1.dereta1^2 + c.copula.be1^-1 * c.copula2.be1 * - der2p1.dereta1eta1 +
            der2pdf1.dereta1 * pdf1^-1 +  derpdf1.dereta1^2*-pdf1^-2 ) +
   
   VC$c01*( -(c.copula.be2)^-2 * c.copula2.be1be2^2 * derp1.dereta1^2 + c.copula.be2^-1 * der2h.derp1p1  * derp1.dereta1^2 +  c.copula.be2^-1 * c.copula2.be1be2 * -der2p1.dereta1eta1   ) +
   
   VC$c00*(  (bit1.b1b1*p00 - (-c.copula.be1*derp1.dereta1)^2)/p00^2  ) 
  
  ) 
  
  


  d2l.be1.sigma21  <- -VC$weights*( 
 
   VC$c11*(  (der2pdf1.dereta1dersigma21.st * pdf1 - derpdf1.dereta1*derpdf1.dersigma21.st) / pdf1^2   + 
                       ((der2c.derp1.derp1 * derp1.dereta1* derp1.dersigma21.st+der2h.derp1p1*-der2p1.dereta1dersigma21.st) * c.copula2.be1be2 - derc.dereta1*derc.dersigma21.st) /c.copula2.be1be2^2  ) +
   
   VC$c10*( -(c.copula.be1)^-2 * c.copula2.be1^2 * derp1.dereta1*derp1.dersigma21.st + c.copula.be1^-1 * der3C.derp1p1p1 * derp1.dereta1*derp1.dersigma21.st + c.copula.be1^-1 * c.copula2.be1 * - der2p1.dereta1dersigma21.st +
            der2pdf1.dereta1dersigma21.st * pdf1^-1 +  derpdf1.dereta1*derpdf1.dersigma21.st*-pdf1^-2    ) +
   
   VC$c01*( -(c.copula.be2)^-2 * c.copula2.be1be2^2 * derp1.dereta1*derp1.dersigma21.st + c.copula.be2^-1 * der2h.derp1p1  * derp1.dereta1*derp1.dersigma21.st +  c.copula.be2^-1 * c.copula2.be1be2 * -der2p1.dereta1dersigma21.st   ) +
    
   VC$c00*(  (bit1.b1s1*p00 - c.copula.be1^2*derp1.dereta1*derp1.dersigma21.st)/p00^2  )
 
 
  )
  
  
  
  
  
  
    
   
    
                   
  d2l.be2.sigma22 <- -VC$weights*( 
  
 
    VC$c11*(  (der2pdf2.dereta2dersigma22.st * pdf2 - derpdf2.dereta2*derpdf2.dersigma22.st) / pdf2^2   + 
                      ((der2c.derp2.derp2 * derp2.dereta2* derp2.dersigma22.st+der2h.derp1p2*-der2p2.dereta2dersigma22.st) * c.copula2.be1be2 - derc.dereta2*derc.dersigma22.st) /c.copula2.be1be2^2  ) +
    
    VC$c10*(-(c.copula.be1)^-2 * c.copula2.be1be2^2 * derp2.dereta2*derp2.dersigma22.st + c.copula.be1^-1 * der2h.derp1p2  * derp2.dereta2*derp2.dersigma22.st +  c.copula.be1^-1 * c.copula2.be1be2 * -der2p2.dereta2dersigma22.st    ) +
    
    VC$c01*(  -(c.copula.be2)^-2 * c.copula2.be2^2 * derp2.dereta2*derp2.dersigma22.st + c.copula.be2^-1 *  der2h.derp2p2 * derp2.dereta2*derp2.dersigma22.st + c.copula.be2^-1 * c.copula2.be2 * - der2p2.dereta2dersigma22.st +
            der2pdf2.dereta2dersigma22.st * pdf2^-1 +  derpdf2.dereta2*derpdf2.dersigma22.st*-pdf2^-2   ) +
    
    VC$c00*(  (bit1.b2s2*p00 - c.copula.be2^2*derp2.dereta2*derp2.dersigma22.st)/p00^2  )
 
  
  )
                        
       
    
    
 
         
       
  d2l.be2.sigma21 <- -VC$weights*( 

   VC$c11*(  (der2c.derp1.derp2*derp2.dereta2*derp1.dersigma21.st* c.copula2.be1be2 - derc.dereta2*derc.dersigma21.st)/c.copula2.be1be2^2  ) +
   
   VC$c10*( -(c.copula.be1)^-2 * c.copula2.be1be2 * c.copula2.be1 * derp1.dersigma21.st * derp2.dereta2 +  c.copula.be1^-1 * der2h.derp1p1 *  derp1.dersigma21.st * derp2.dereta2     ) +
   
   VC$c01*( -(c.copula.be2)^-2 * c.copula2.be1be2 * c.copula2.be2 * derp1.dersigma21.st * derp2.dereta2 +  c.copula.be2^-1 * der2h.derp1p2 *  derp1.dersigma21.st * derp2.dereta2   ) +
   
   VC$c00*(  (bit1.b2s1*p00 - (c.copula.be2*derp2.dereta2*c.copula.be1*derp1.dersigma21.st))/p00^2  )  
  
  )
  
  
  
  
  
  d2l.be1.sigma22 <- -VC$weights*( 

   VC$c11*(  (der2c.derp1.derp2*derp1.dereta1*derp2.dersigma22.st* c.copula2.be1be2 - derc.dereta1*derc.dersigma22.st)/c.copula2.be1be2^2  ) +
   
   VC$c10*( -(c.copula.be1)^-2 * c.copula2.be1be2 * c.copula2.be1 * derp1.dereta1 * derp2.dersigma22.st +  c.copula.be1^-1 * der2h.derp1p1 *  derp1.dereta1 * derp2.dersigma22.st   ) +
   
   VC$c01*( -(c.copula.be2)^-2 * c.copula2.be1be2 * c.copula2.be2 * derp1.dereta1 * derp2.dersigma22.st +  c.copula.be2^-1 * der2h.derp1p2 *  derp1.dereta1 * derp2.dersigma22.st   ) +
   
   VC$c00*(  (bit1.b1s2*p00 - (c.copula.be1*derp1.dereta1*c.copula.be2*derp2.dersigma22.st))/p00^2  )

  
  )
  
  
  
  

  
  d2l.rho.sigma21 <-  -VC$weights*( 

   VC$c11*(  (der2c.derp1.derrho*-derp1.dersigma21.st*derteta.derteta.st*c.copula2.be1be2 - -derc.dersigma21.st*der2h.derp1teta *derteta.derteta.st)/c.copula2.be1be2^2  ) +
   
   VC$c10*( -(c.copula.be1)^-2 * c.copula2.be1th *  c.copula2.be1 * - derp1.dersigma21.st + c.copula.be1^-1 * der3C.p1p1teta *derteta.derteta.st*- derp1.dersigma21.st   ) +
   
   VC$c01*(  -(c.copula.be2)^-2 * c.copula2.be2th *  c.copula2.be1be2 * - derp1.dersigma21.st + c.copula.be2^-1 * der2h.derp1teta *derteta.derteta.st*- derp1.dersigma21.st   ) +
   
   VC$c00*(  (bit1.s1th*p00 - (-c.copula.be1*derp1.dersigma21.st*c.copula.theta))/p00^2  )  
  
  )
  
  
  
  
  
  
  
  d2l.rho.sigma22 <-  -VC$weights*( 

   VC$c11*(  (der2c.derp2.derrho*-derp2.dersigma22.st*derteta.derteta.st*c.copula2.be1be2 - -derc.dersigma22.st*der2h.derp1teta *derteta.derteta.st)/c.copula2.be1be2^2  ) +
   
   VC$c10*(  -(c.copula.be1)^-2 * c.copula2.be1th  *  c.copula2.be1be2 * - derp2.dersigma22.st + c.copula.be1^-1 * der2h.derp1teta *derteta.derteta.st*- derp2.dersigma22.st   ) +
   
   VC$c01*( -(c.copula.be2)^-2 * c.copula2.be2th  *  c.copula2.be2 * - derp2.dersigma22.st + c.copula.be2^-1 * der2h.derp2teta *derteta.derteta.st*- derp2.dersigma22.st    ) +
   
   VC$c00*(  (bit1.s2th*p00 - (-c.copula.be2*derp2.dersigma22.st*c.copula.theta))/p00^2  )
  
  
  )
  
  
  
  
    
  
  
  
  
  
  d2l.sigma21.sigma22 <- -VC$weights*( 


   VC$c11*(  (der2c.derp1.derp2*derp1.dersigma21.st*derp2.dersigma22.st* c.copula2.be1be2 - derc.dersigma21.st*derc.dersigma22.st ) /c.copula2.be1be2^2  ) +
   
   VC$c10*( -(c.copula.be1)^-2 * c.copula2.be1be2 * c.copula2.be1 * derp1.dersigma21.st * derp2.dersigma22.st +  c.copula.be1^-1 * der2h.derp1p1 *  derp1.dersigma21.st * derp2.dersigma22.st    ) +
   
   VC$c01*( -(c.copula.be2)^-2 * c.copula2.be1be2 * c.copula2.be2 * derp1.dersigma21.st * derp2.dersigma22.st +  c.copula.be2^-1 * der2h.derp1p2 *  derp1.dersigma21.st * derp2.dersigma22.st    ) +
   
   VC$c00*(  (bit1.s1s2*p00 - (c.copula.be1*derp1.dersigma21.st*c.copula.be2*derp2.dersigma22.st))/p00^2  )
   
  
  )
         
         
         
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
         
         
         
         
 

if( is.null(VC$X3) ){

  G   <- -c(colSums( c(dl.dbe1)*VC$X1 ),
            colSums( c(dl.dbe2)*VC$X2 ),
            sum( dl.dsigma21.st ),
            sum( dl.dsigma22.st ),
            sum( dl.dteta.st )            )



  be1.be1 <- crossprod(VC$X1*c(d2l.be1.be1),VC$X1)
  be2.be2 <- crossprod(VC$X2*c(d2l.be2.be2),VC$X2)
  be1.be2 <- crossprod(VC$X1*c(d2l.be1.be2),VC$X2)
  be1.rho <-   t(t(rowSums(t(VC$X1*c(d2l.be1.rho)))))
  be1.sigma21 <- t(t(rowSums(t(VC$X1*c(d2l.be1.sigma21))))) 
  be1.sigma22 <- t(t(rowSums(t(VC$X1*c(d2l.be1.sigma22))))) 
  be2.rho <-   t(t(rowSums(t(VC$X2*c(d2l.be2.rho)))))
  be2.sigma21 <- t(t(rowSums(t(VC$X2*c(d2l.be2.sigma21))))) 
  be2.sigma22 <- t(t(rowSums(t(VC$X2*c(d2l.be2.sigma22)))))

  H <- rbind( cbind( be1.be1    ,   be1.be2    ,   be1.sigma21,    be1.sigma22,       be1.rho  ), 
              cbind( t(be1.be2) ,   be2.be2    ,   be2.sigma21,    be2.sigma22,        be2.rho  ), 
              cbind( t(be1.sigma21) , t(be2.sigma21) , sum(d2l.sigma21.sigma21), sum(d2l.sigma21.sigma22), sum(d2l.rho.sigma21)),
              cbind( t(be1.sigma22) , t(be2.sigma22) , sum(d2l.sigma21.sigma22), sum(d2l.sigma22.sigma22), sum(d2l.rho.sigma22)),
              cbind( t(be1.rho) ,   t(be2.rho) ,   sum(d2l.rho.sigma21), sum(d2l.rho.sigma22), sum(d2l.rho.rho) ) 
              
              ) 


}



if( !is.null(VC$X3) ){



G   <- -c( colSums(       c(dl.dbe1)*VC$X1 ) ,
           colSums(       c(dl.dbe2)*VC$X2 ) ,
           colSums(c(dl.dsigma21.st)*VC$X3 ) ,
           colSums(c(dl.dsigma22.st)*VC$X4 ) ,
           colSums(   c(dl.dteta.st)*VC$X5 )  )  

                 
    be1.be1         <- crossprod(VC$X1*c(d2l.be1.be1),VC$X1)
    be2.be2         <- crossprod(VC$X2*c(d2l.be2.be2),VC$X2)
    be1.be2         <- crossprod(VC$X1*c(d2l.be1.be2),VC$X2)
    be1.rho         <- crossprod(VC$X1*c(d2l.be1.rho),VC$X5)    
    be2.rho         <- crossprod(VC$X2*c(d2l.be2.rho),VC$X5)   
    be1.sigma21     <- crossprod(VC$X1*c(d2l.be1.sigma21),VC$X3)      
    be1.sigma22     <- crossprod(VC$X1*c(d2l.be1.sigma22),VC$X4)      
    be2.sigma21     <- crossprod(VC$X2*c(d2l.be2.sigma21),VC$X3)   
    be2.sigma22     <- crossprod(VC$X2*c(d2l.be2.sigma22),VC$X4)   
    sigma21.sigma21 <- crossprod(VC$X3*c(d2l.sigma21.sigma21),VC$X3)     
    sigma21.sigma22 <- crossprod(VC$X3*c(d2l.sigma21.sigma22),VC$X4)
    rho.sigma21     <- crossprod(VC$X3*c(d2l.rho.sigma21),VC$X5)
    sigma22.sigma22 <- crossprod(VC$X4*c(d2l.sigma22.sigma22),VC$X4)
    rho.sigma22     <- crossprod(VC$X4*c(d2l.rho.sigma22),VC$X5)   
    rho.rho         <- crossprod(VC$X5*c(d2l.rho.rho),VC$X5)    
    
    
    H <- rbind( cbind( be1.be1        ,   be1.be2      ,   be1.sigma21,      be1.sigma22,       be1.rho    ), 
                cbind( t(be1.be2)     ,   be2.be2      ,   be2.sigma21,      be2.sigma22,       be2.rho    ), 
                cbind( t(be1.sigma21) , t(be2.sigma21) ,   sigma21.sigma21,  sigma21.sigma22,   rho.sigma21),
                cbind( t(be1.sigma22) , t(be2.sigma22) , t(sigma21.sigma22), sigma22.sigma22,   rho.sigma22),
                cbind( t(be1.rho)     ,   t(be2.rho)   ,   t(rho.sigma21),   t(rho.sigma22),    rho.rho    ) ) 
                
                 

}



res <- -sum(l.par)

 
if(VC$extra.regI == "pC") H <- regH(H, type = 1)
  
  
  S.h  <- ps$S.h  
  
  if( length(S.h) != 1){
  
  S.h1 <- 0.5*crossprod(params,S.h)%*%params # as(S.h, "sparseMatrix") not sure about advantage...
  S.h2 <- S.h%*%params
  
  } else S.h <- S.h1 <- S.h2 <- 0   
  
  S.res <- res
  res   <- S.res + S.h1
  G     <- G + S.h2
  H     <- H + S.h  
        
      
if(VC$extra.regI == "sED") H <- regH(H, type = 2)   
  
  
  
if( VC$margins[1] == "LN" || VC$margins[2] == "LN"){
  
  
  if(VC$margins[1] == "LN") dHs1 <- distrHsAT(exp(respvec$y1), eta1, sigma21, 1, margin2=VC$margins[1])
  if(VC$margins[2] == "LN") dHs2 <- distrHsAT(exp(respvec$y2), eta2, sigma22, 1, margin2=VC$margins[2])

  pdf1 <- dHs1$pdf2
  pdf2 <- dHs2$pdf2

  p1 <- dHs1$p2
  p2 <- dHs2$p2
 
 
  if( length(teta1) != 0) dH1 <- copgHsAT(p1[teta.ind1], p2[teta.ind1], teta1, Cop1, Ln = TRUE)
  if( length(teta2) != 0) dH2 <- copgHsAT(p1[teta.ind2], p2[teta.ind2], teta2, Cop2, Ln = TRUE)
  c.copula2.be1be2 <- NA
  if( length(teta1) != 0) c.copula2.be1be2[teta.ind1] <- dH1$c.copula2.be1be2
  if( length(teta2) != 0) c.copula2.be1be2[teta.ind2] <- dH2$c.copula2.be1be2    
  

  l.ln <- -sum( VC$weights*( log(pdf1) + log(pdf2) + log(c.copula2.be1be2) ) )

  }
  
  
  
  

         list(value=res, gradient=G, hessian=H, S.h=S.h, S.h1=S.h1, S.h2=S.h2, 
              l=S.res, l.ln = l.ln, l.par=l.par, ps = ps, 
              eta1=eta1, eta2=eta2, etad=etad, etas1 = etas1, etas2 = etas2, 
              BivD=VC$BivD, p1 = p1, p2 = p2,
              dl.dbe1          =dl.dbe1,       
              dl.dbe2          =dl.dbe2,       
              dl.dsigma21.st   =dl.dsigma21.st,
              dl.dsigma22.st   =dl.dsigma22.st,
              dl.dteta.st      =dl.dteta.st, 
              teta.ind2 = teta.ind2, teta.ind1 = teta.ind1,
              Cop1 = Cop1, Cop2 = Cop2, teta1 = teta1, teta2 = teta2) 
              
              

  }
  

  
