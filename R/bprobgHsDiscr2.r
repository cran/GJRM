bprobgHsDiscr2 <- function(params, respvec, VC, ps, AT = FALSE){


  eta1 <- VC$X1%*%params[1:VC$X1.d2]
  eta2 <- VC$X2%*%params[(VC$X1.d2+1):(VC$X1.d2+VC$X2.d2)]
  etad <- etas <- l.ln <- NULL 

  epsilon <- sqrt(.Machine$double.eps)
  max.p   <- 0.9999999
  
  
if(is.null(VC$X3)){  
  sigma2.st <- etas <- params[(VC$X1.d2 + VC$X2.d2 + 1)]
  teta.st   <- etad <- params[(VC$X1.d2 + VC$X2.d2 + 2)]
} 

if(!is.null(VC$X3)){  
  sigma2.st <- etas <- VC$X3%*%params[(VC$X1.d2+VC$X2.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2)]
  teta.st   <- etad <- VC$X4%*%params[(VC$X1.d2+VC$X2.d2+VC$X3.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2+VC$X4.d2)]
}  
  
    sstr1 <- esp.tr(sigma2.st, VC$margins[2])  
    sigma2.st <- sstr1$vrb.st 
    sigma2    <- sstr1$vrb 

    eta2 <- eta.tr(eta2, VC$margins[2])
    
  
 dHs <- distrHsDiscr(respvec$y2, eta2, sigma2, sigma2.st, nu = 1, nu.st = 1, margin2=VC$margins[2], naive = FALSE, y2m = VC$y2m)
  
  
 pdf2                         <- dHs$pdf2
 p2                           <- dHs$p2 
 derpdf2.dereta2              <- dHs$derpdf2.dereta2 
 derpdf2.dersigma2.st         <- dHs$derpdf2.dersigma2.st 
 derp2.dersigma.st            <- dHs$derp2.dersigma.st
 derp2.dereta2                <- dHs$derp2.dereta2
 der2p2.dereta2eta2           <- dHs$der2p2.dereta2eta2 
 der2pdf2.dereta2             <- dHs$der2pdf2.dereta2
 der2p2.dersigma2.st2         <- dHs$der2p2.dersigma2.st2
 der2pdf2.dersigma2.st2       <- dHs$der2pdf2.dersigma2.st2
 der2p2.dereta2dersigma2.st   <- dHs$der2p2.dereta2dersigma2.st            
 der2pdf2.dereta2dersigma2.st <- dHs$der2pdf2.dereta2dersigma2.st  
  
  
  pd1 <- probm(eta1, VC$margins[1], bc = TRUE) 
  p1  <- 1 - pd1$pr                           #   pnorm(-eta1), p(y1=0)


########################################################################################################  
  
resT    <- teta.tr(VC, teta.st)

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
    

########################################################################################################


  C1 <- C2 <- A <- B <- NA
  
  
if( length(teta1) != 0){  
  
  C1[teta.ind1] <- BiCDF(p1[teta.ind1], p2[teta.ind1],          nC1, teta1, VC$dof)
  C2[teta.ind1] <- BiCDF(p1[teta.ind1], mm(p2[teta.ind1]-pdf2[teta.ind1]), nC1, teta1, VC$dof)
  
  A[teta.ind1] <- ifelse(C1[teta.ind1] - C2[teta.ind1] < epsilon, epsilon, C1[teta.ind1] - C2[teta.ind1])
  B[teta.ind1] <- ifelse( pdf2[teta.ind1] - A[teta.ind1] < epsilon, epsilon, pdf2[teta.ind1] - A[teta.ind1])
  
}  


if( length(teta2) != 0){  
  
  C1[teta.ind2] <- BiCDF(p1[teta.ind2], p2[teta.ind2],          nC2, teta2, VC$dof)
  C2[teta.ind2] <- BiCDF(p1[teta.ind2], mm(p2[teta.ind2]-pdf2[teta.ind2]), nC2, teta2, VC$dof)
  
  A[teta.ind2] <- ifelse(C1[teta.ind2] - C2[teta.ind2] < epsilon, epsilon, C1[teta.ind2] - C2[teta.ind2])
  B[teta.ind2] <- ifelse( pdf2[teta.ind2] - A[teta.ind2] < epsilon, epsilon, pdf2[teta.ind2] - A[teta.ind2])
  
}

  

l.par <- VC$weights*( respvec$cy*log( A ) + respvec$y1*log( B ) )
  
########################################################################################################
 
   c.copula.be1.C1 <- c.copula.be1.C2 <- c.copula.be2.C1 <- c.copula.be2.C2 <- c.copula.theta.C1 <- c.copula.theta.C2 <- NA
   
   
   
   
   if( length(teta1) != 0) dH1F <- copgHs(p1[teta.ind1], p2[teta.ind1], eta1=NULL, eta2=NULL, teta1, teta.st1, Cop1, VC$dof)
   if( length(teta2) != 0) dH1S <- copgHs(p1[teta.ind2], p2[teta.ind2], eta1=NULL, eta2=NULL, teta2, teta.st2, Cop2, VC$dof)
   
   if( length(teta1) != 0) dH2F <- copgHs(p1[teta.ind1], mm(p2[teta.ind1]-pdf2[teta.ind1]), eta1=NULL, eta2=NULL, teta1, teta.st1, Cop1, VC$dof) 
   if( length(teta2) != 0) dH2S <- copgHs(p1[teta.ind2], mm(p2[teta.ind2]-pdf2[teta.ind2]), eta1=NULL, eta2=NULL, teta2, teta.st2, Cop2, VC$dof) 
   
   
 if( length(teta1) != 0){  
    
   c.copula.be1.C1[teta.ind1] <- dH1F$c.copula.be1 
   c.copula.be1.C2[teta.ind1] <- dH2F$c.copula.be1 
   
   c.copula.be2.C1[teta.ind1] <- dH1F$c.copula.be2 
   c.copula.be2.C2[teta.ind1] <- dH2F$c.copula.be2 
   
   c.copula.theta.C1[teta.ind1] <- dH1F$c.copula.theta # here there is theta star already
   c.copula.theta.C2[teta.ind1] <- dH2F$c.copula.theta
   
 }  
   
   
   if( length(teta2) != 0){  
      
     c.copula.be1.C1[teta.ind2] <- dH1S$c.copula.be1 
     c.copula.be1.C2[teta.ind2] <- dH2S$c.copula.be1 
     
     c.copula.be2.C1[teta.ind2] <- dH1S$c.copula.be2 
     c.copula.be2.C2[teta.ind2] <- dH2S$c.copula.be2 
     
     c.copula.theta.C1[teta.ind2] <- dH1S$c.copula.theta # here there is theta star already
     c.copula.theta.C2[teta.ind2] <- dH2S$c.copula.theta
     
   }  
     
  

derp2m1.dereta2     <- derp2.dereta2 - derpdf2.dereta2
derp2m1.dersigma.st <- derp2.dersigma.st - derpdf2.dersigma2.st 



derp1.dereta1   <- pd1$derp1.dereta1    # -dnorm(-eta1) 


Cc <- mm(c.copula.be1.C1 - c.copula.be1.C2) 
C  <- Cc*derp1.dereta1

Cs    <- c.copula.theta.C1  - c.copula.theta.C2
Cssb2 <- c.copula.be2.C1*derp2.dereta2 - c.copula.be2.C2*derp2m1.dereta2  
CssSI <- c.copula.be2.C1*derp2.dersigma.st - c.copula.be2.C2*derp2m1.dersigma.st  



  dl.dbe1      <- VC$weights*( ( respvec$cy/A - respvec$y1/B )*C  ) 
  dl.dbe2      <- VC$weights*( respvec$cy/A*Cssb2 + respvec$y1/B*(derpdf2.dereta2 - Cssb2)    )
  dl.dsigma.st <- VC$weights*( respvec$cy/A*CssSI + respvec$y1/B*(derpdf2.dersigma2.st - CssSI)  )
  dl.dteta.st  <- VC$weights*( ( respvec$cy/A - respvec$y1/B )*Cs   )                     
 
 
######################################################################################################## 
 
  c.copula2.be1.C1 <- c.copula2.be1.C2 <- c.copula2.be2.C1 <- c.copula2.be2.C2 <- c.copula2.be1be2.C1 <- c.copula2.be1be2.C2 <- c.copula2.be2th.C1 <- c.copula2.be2th.C2 <- c.copula2.theta.C1 <- c.copula2.theta.C2 <- c.copula.thet.C1 <- c.copula.thet.C2 <- derteta.derteta.st <- der2teta.derteta.stteta.st <- c.copula2.be1th.C1 <- c.copula2.be1th.C2 <- NA
    
    
    
 if( length(teta1) != 0){     
    
     
     c.copula2.be1.C1[teta.ind1]           <- dH1F$c.copula2.be1
     c.copula2.be1.C2[teta.ind1]           <- dH2F$c.copula2.be1
     
     c.copula2.be2.C1[teta.ind1]           <- dH1F$c.copula2.be2
     c.copula2.be2.C2[teta.ind1]           <- dH2F$c.copula2.be2
     
     c.copula2.be1be2.C1[teta.ind1]        <- dH1F$c.copula2.be1be2
     c.copula2.be1be2.C2[teta.ind1]        <- dH2F$c.copula2.be1be2
     
     c.copula2.be2th.C1[teta.ind1]         <- dH1F$c.copula2.be2th
     c.copula2.be2th.C2[teta.ind1]         <- dH2F$c.copula2.be2th  
     
     c.copula2.theta.C1[teta.ind1]         <- dH1F$bit1.th2ATE 
     c.copula2.theta.C2[teta.ind1]         <- dH2F$bit1.th2ATE 
     
     c.copula.thet.C1[teta.ind1]           <- dH1F$c.copula.thet # NO star
     c.copula.thet.C2[teta.ind1]           <- dH2F$c.copula.thet    
    
     derteta.derteta.st[teta.ind1]         <- dH1F$derteta.derteta.st         # does not matter dH1 or dH2 
     der2teta.derteta.stteta.st[teta.ind1] <- dH1F$der2teta.derteta.stteta.st   
    
     c.copula2.be1th.C1[teta.ind1]         <- dH1F$c.copula2.be1th 
     c.copula2.be1th.C2[teta.ind1]         <- dH2F$c.copula2.be1th    
    
 }   
 
 
 
 if( length(teta2) != 0){     
    
     c.copula2.be1.C1[teta.ind2]           <- dH1S$c.copula2.be1
     c.copula2.be1.C2[teta.ind2]           <- dH2S$c.copula2.be1
     
     c.copula2.be2.C1[teta.ind2]           <- dH1S$c.copula2.be2
     c.copula2.be2.C2[teta.ind2]           <- dH2S$c.copula2.be2
     
     c.copula2.be1be2.C1[teta.ind2]        <- dH1S$c.copula2.be1be2
     c.copula2.be1be2.C2[teta.ind2]        <- dH2S$c.copula2.be1be2
     
     c.copula2.be2th.C1[teta.ind2]         <- dH1S$c.copula2.be2th
     c.copula2.be2th.C2[teta.ind2]         <- dH2S$c.copula2.be2th  
     
     c.copula2.theta.C1[teta.ind2]         <- dH1S$bit1.th2ATE 
     c.copula2.theta.C2[teta.ind2]         <- dH2S$bit1.th2ATE 
     
     c.copula.thet.C1[teta.ind2]           <- dH1S$c.copula.thet # NO star
     c.copula.thet.C2[teta.ind2]           <- dH2S$c.copula.thet    
    
     derteta.derteta.st[teta.ind2]         <- dH1S$derteta.derteta.st         # does not matter dH1 or dH2 
     der2teta.derteta.stteta.st[teta.ind2] <- dH1S$der2teta.derteta.stteta.st   
    
     c.copula2.be1th.C1[teta.ind2]         <- dH1S$c.copula2.be1th 
     c.copula2.be1th.C2[teta.ind2]         <- dH2S$c.copula2.be1th    
    
 }   
 
 

  
  
 
  der2p1.dereta1eta1 <- pd1$der2p1.dereta1eta1
  
  derC.dereta1 <- (c.copula2.be1.C1 - c.copula2.be1.C2)*derp1.dereta1^2 + Cc*der2p1.dereta1eta1
  
 
  
  derCs.dertheta.st <- (c.copula2.theta.C1 - c.copula2.theta.C2)*derteta.derteta.st^2 + (c.copula.thet.C1 - c.copula.thet.C2)*der2teta.derteta.stteta.st
  
  derA.dereta2 <- c.copula.be2.C1*derp2.dereta2   - c.copula.be2.C2*derp2m1.dereta2 
  derB.dereta2 <- derpdf2.dereta2 - derA.dereta2  
  der2p2m1.dereta2eta2 <- der2p2.dereta2eta2 - der2pdf2.dereta2
  derCssb2.dereta2 <- c.copula2.be2.C1*derp2.dereta2^2 + c.copula.be2.C1*der2p2.dereta2eta2 - (c.copula2.be2.C2*derp2m1.dereta2^2 + c.copula.be2.C2*der2p2m1.dereta2eta2)                                                                     
      
  derA.dersigma2.st <- c.copula.be2.C1*derp2.dersigma.st - c.copula.be2.C2*derp2m1.dersigma.st 
  derB.dersigma2.st <- derpdf2.dersigma2.st - derA.dersigma2.st 
  der2p2m1.dersigma2.st2 <- der2p2.dersigma2.st2 - der2pdf2.dersigma2.st2
  derCssSI.dersigma2.st <- c.copula2.be2.C1*derp2.dersigma.st^2 + c.copula.be2.C1*der2p2.dersigma2.st2 - (c.copula2.be2.C2*derp2m1.dersigma.st^2 + c.copula.be2.C2*der2p2m1.dersigma2.st2)                                                                     
      
  derC.dereta2 <- (c.copula2.be1be2.C1*derp2.dereta2 - c.copula2.be1be2.C2*derp2m1.dereta2)*derp1.dereta1    
  derC.dersigma2.st <- (c.copula2.be1be2.C1*derp2.dersigma.st - c.copula2.be1be2.C2*derp2m1.dersigma.st)*derp1.dereta1 
 
 
  derC.dertheta.st <- (c.copula2.be1th.C1 - c.copula2.be1th.C2)*derp1.dereta1 
  
  der2p2m1.dereta2dersigma2.st <- der2p2.dereta2dersigma2.st - der2pdf2.dereta2dersigma2.st 
  
  derCssb2.dersigma2.st <- c.copula2.be2.C1*derp2.dereta2*derp2.dersigma.st + c.copula.be2.C1*der2p2.dereta2dersigma2.st - (c.copula2.be2.C2*derp2m1.dereta2*derp2m1.dersigma.st + c.copula.be2.C2*der2p2m1.dereta2dersigma2.st)   
  
  
  derCs.dereta2      <- c.copula2.be2th.C1*derp2.dereta2 - c.copula2.be2th.C2*derp2m1.dereta2
  derCs.dersigma2.st <- c.copula2.be2th.C1*derp2.dersigma.st - c.copula2.be2th.C2*derp2m1.dersigma.st
  

  d2l.be1.be1      <- -VC$weights*( (-respvec$cy/A^2 - respvec$y1/B^2)*C^2 + ( respvec$cy/A - respvec$y1/B )*derC.dereta1 )
  d2l.rho.rho      <- -VC$weights*( (-respvec$cy/A^2 - respvec$y1/B^2)*Cs^2 + ( respvec$cy/A - respvec$y1/B )*derCs.dertheta.st )
  d2l.be2.be2      <- -VC$weights*( -respvec$cy*derA.dereta2*Cssb2/A^2 + respvec$cy*derCssb2.dereta2/A - respvec$y1*derB.dereta2/B^2*(derpdf2.dereta2-Cssb2) + respvec$y1/B*( der2pdf2.dereta2 - derCssb2.dereta2)  )
  d2l.sigma.sigma  <- -VC$weights*( -respvec$cy*derA.dersigma2.st*CssSI/A^2 + respvec$cy*derCssSI.dersigma2.st/A - respvec$y1*derB.dersigma2.st/B^2*(derpdf2.dersigma2.st-CssSI) + respvec$y1/B*( der2pdf2.dersigma2.st2 - derCssSI.dersigma2.st)  )
  d2l.be1.be2      <- -VC$weights*( (-respvec$cy*derA.dereta2/A^2 + respvec$y1*derB.dereta2/B^2)*C + (respvec$cy/A - respvec$y1/B)*derC.dereta2 )
  d2l.be1.sigma    <- -VC$weights*( (-respvec$cy*derA.dersigma2.st/A^2 + respvec$y1*derB.dersigma2.st/B^2)*C + (respvec$cy/A - respvec$y1/B)*derC.dersigma2.st )
  d2l.be1.rho      <- -VC$weights*( (-respvec$cy/A^2 - respvec$y1/B^2)*C*Cs + ( respvec$cy/A - respvec$y1/B )*derC.dertheta.st )
  d2l.be2.sigma    <- -VC$weights*( -respvec$cy*derA.dersigma2.st*Cssb2/A^2 + respvec$cy*derCssb2.dersigma2.st/A - respvec$y1*derB.dersigma2.st/B^2*(derpdf2.dereta2-Cssb2) + respvec$y1/B*( der2pdf2.dereta2dersigma2.st - derCssb2.dersigma2.st)    )

  d2l.be2.rho      <- -VC$weights*( (-respvec$cy*derA.dereta2/A^2 + respvec$y1*derB.dereta2/B^2)*Cs + ( respvec$cy/A - respvec$y1/B )*derCs.dereta2 )
  

  d2l.rho.sigma    <- -VC$weights*( (-respvec$cy*derA.dersigma2.st/A^2 + respvec$y1*derB.dersigma2.st/B^2)*Cs + ( respvec$cy/A - respvec$y1/B )*derCs.dersigma2.st )
  
  
  

if( is.null(VC$X3) ){


  be1.be1   <- crossprod(VC$X1*c(d2l.be1.be1),VC$X1)
  be2.be2   <- crossprod(VC$X2*c(d2l.be2.be2),VC$X2)
  be1.be2   <- crossprod(VC$X1*c(d2l.be1.be2),VC$X2)
  be1.rho   <- t(t(rowSums(t(VC$X1*c(d2l.be1.rho)))))
  be1.sigma <- t(t(rowSums(t(VC$X1*c(d2l.be1.sigma))))) 
  be2.rho   <- t(t(rowSums(t(VC$X2*c(d2l.be2.rho)))))
  be2.sigma <- t(t(rowSums(t(VC$X2*c(d2l.be2.sigma))))) 

  H <- rbind( cbind( be1.be1    ,   be1.be2    ,   be1.sigma,            be1.rho  ), 
              cbind( t(be1.be2) ,   be2.be2    ,   be2.sigma,            be2.rho  ), 
              cbind( t(be1.sigma),  t(be2.sigma),  sum(d2l.sigma.sigma), sum(d2l.rho.sigma) ),
              cbind( t(be1.rho) ,   t(be2.rho),    sum(d2l.rho.sigma),   sum(d2l.rho.rho)   ) 
              ) 
         
  G   <- -c( colSums( c(dl.dbe1)*VC$X1 ) ,
             colSums( c(dl.dbe2)*VC$X2 ) ,
             sum( dl.dsigma.st ),
             sum( dl.dteta.st ) )
    
}




if( !is.null(VC$X3) ){

  be1.be1   <- crossprod(VC$X1*c(d2l.be1.be1),VC$X1)
  be2.be2   <- crossprod(VC$X2*c(d2l.be2.be2),VC$X2)
  be1.be2   <- crossprod(VC$X1*c(d2l.be1.be2),VC$X2)
  
  be1.rho   <- crossprod(VC$X1*c(d2l.be1.rho),  VC$X4)                                     
  be1.sigma <- crossprod(VC$X1*c(d2l.be1.sigma),VC$X3)                                   
  be2.rho   <- crossprod(VC$X2*c(d2l.be2.rho),  VC$X4)                                     
  be2.sigma <- crossprod(VC$X2*c(d2l.be2.sigma),VC$X3)  
  
  sigma.sigma <- crossprod(VC$X3*c(d2l.sigma.sigma),VC$X3)   
  sigma.rho   <- crossprod(VC$X3*c(d2l.rho.sigma),VC$X4)  
  rho.rho     <- crossprod(VC$X4*c(d2l.rho.rho),    VC$X4)    
  
  

  H <- rbind( cbind( be1.be1     ,  be1.be2     ,  be1.sigma   , be1.rho   ), 
              cbind( t(be1.be2)  ,  be2.be2     ,  be2.sigma   , be2.rho   ), 
              cbind( t(be1.sigma),  t(be2.sigma),  sigma.sigma , sigma.rho ),
              cbind( t(be1.rho)  ,  t(be2.rho)  ,  t(sigma.rho), rho.rho   ) 
             )  
            
   
  G   <- -c( colSums(      c(dl.dbe1)*VC$X1 ) ,
             colSums(      c(dl.dbe2)*VC$X2 ) ,
             colSums( c(dl.dsigma.st)*VC$X3 ) ,
             colSums(  c(dl.dteta.st)*VC$X4 ) )   
   
    
}




      res <- -sum(l.par)





 
if(VC$extra.regI == "pC") H <- regH(H, type = 1)
  
  S.h  <- ps$S.h  


  if( length(S.h) != 1){
  
  S.h1 <- 0.5*crossprod(params,S.h)%*%params
  S.h2 <- S.h%*%params
  
  } else S.h <- S.h1 <- S.h2 <- 0   
  
  S.res <- res
  res   <- S.res + S.h1
  G     <- G + S.h2
  H     <- H + S.h  
        
if(VC$extra.regI == "sED") H <- regH(H, type = 2)   
  
 
         list(value=res, gradient=G, hessian=H, S.h=S.h, S.h1=S.h1, S.h2=S.h2, l=S.res, l.par=l.par, ps = ps, etas = etas,
              eta1=eta1, eta2=eta2, etad=etad,
              dl.dbe1=dl.dbe1, dl.dbe2=dl.dbe2, dl.dsigma.st = dl.dsigma.st, dl.dteta.st = dl.dteta.st,
              BivD=VC$BivD, p1=1-p1, p2=p2, theta.star = teta.st,
              teta.ind2 = teta.ind2, teta.ind1 = teta.ind1,
              Cop1 = Cop1, Cop2 = Cop2, teta1 = teta1, teta2 = teta2)      

}




     























