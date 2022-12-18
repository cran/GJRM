bprobgHsCont3binTW <- function(params, respvec, VC, ps, AT = FALSE){

###################################################################################

p1 <- p2 <- pdf1 <- pdf2 <- c.copula.be2 <- c.copula.be1 <- c.copula2.be1be2 <- NA

y00 <- respvec$y00
y10 <- respvec$y10
y0p <- respvec$y0p
y1p <- respvec$y1p

YP <- y0p + y1p  
Y0 <- y00 + y10 



# note: code below will be perhaps a bit inefficient but
# easily allows to combine two files reducing the chances of mistakes
###################################################################################


  eta1 <- VC$X1%*%params[1:VC$X1.d2]
  eta2 <- VC$X2%*%params[(VC$X1.d2+1):(VC$X1.d2+VC$X2.d2)]
  etad <- etas <- etan <- NULL 
 
if(is.null(VC$X3)){  
  sigma2.st <- etas <- params[(VC$X1.d2 + VC$X2.d2 + 1)]
  nu.st     <- etan <- params[(VC$X1.d2 + VC$X2.d2 + 2)]
  teta.st   <- etad <- params[(VC$X1.d2 + VC$X2.d2 + 3)]
} 

if(!is.null(VC$X3)){  
  sigma2.st <- etas <- VC$X3%*%params[(VC$X1.d2+VC$X2.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2)]
  nu.st     <- etan <- VC$X4%*%params[(VC$X1.d2+VC$X2.d2+VC$X3.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2+VC$X4.d2)]
  teta.st   <- etad <- VC$X5%*%params[(VC$X1.d2+VC$X2.d2+VC$X3.d2+VC$X4.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2+VC$X4.d2+VC$X5.d2)]
}  
  
  
    sstr1 <- esp.tr(sigma2.st, VC$margins[2])  
    sigma2.st <- sstr1$vrb.st 
    sigma2    <- sstr1$vrb 

    
sstr1 <- enu.tr(nu.st, VC$margins[2])  
nu.st <- sstr1$vrb.st 
nu    <- sstr1$vrb 

    eta2 <- eta.tr(eta2, VC$margins[2])
    

 dHs  <- distrHs(respvec$y2, eta2, sigma2, sigma2.st, nu, nu.st, margin2=VC$margins[2], naive = FALSE, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
  
 pdf2                         <- dHs$pdf2
 p2                           <- dHs$p2 
 derpdf2.dereta2              <- dHs$derpdf2.dereta2 
 derpdf2.dersigma2.st         <- dHs$derpdf2.dersigma2.st 
 derp2.dersigma.st            <- dHs$derp2.dersigma.st
 derpdf2.dernu.st             <- dHs$derpdf2.dernu.st 
 derp2.dernu.st               <- dHs$derp2.nu.st
 derp2.dereta2                <- dHs$derp2.dereta2
 der2p2.dereta2eta2           <- dHs$der2p2.dereta2eta2 
 der2pdf2.dereta2             <- dHs$der2pdf2.dereta2
 der2p2.dersigma2.st2         <- dHs$der2p2.dersigma2.st2
 der2pdf2.dersigma2.st2       <- dHs$der2pdf2.dersigma2.st2
 
 der2p2.dernu.st2         <- dHs$der2p2.dernu.st2
 der2pdf2.dernu.st2       <- dHs$der2pdf2.dernu.st2
 
 der2p2.dereta2dersigma2.st   <- dHs$der2p2.dereta2dersigma2.st            
 der2pdf2.dereta2dersigma2.st <- dHs$der2pdf2.dereta2dersigma2.st  
 
 
 der2p2.dereta2dernu.st   <- dHs$der2p2.dereta2dernu.st            
 der2pdf2.dereta2dernu.st <- dHs$der2pdf2.dereta2dernu.st 
  
 der2p2.dersigma2.stdernu.st   <- dHs$der2p2.dersigma2.stdernu.st            
 der2pdf2.dersigma2.stdernu.st <- dHs$der2pdf2.sigma2.st2dernu.st 
  
  pd1 <- probm(eta1, VC$margins[1], bc = TRUE, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr) 
  
  p1 <- 1 - pd1$pr #   pnorm(-eta1)

########################################################################################################  
   
 resT    <- teta.tr(VC, teta.st)
 
 teta.st1 <- teta.st2 <- teta.st <- resT$teta.st
 teta1 <- teta2 <- teta <- resT$teta 
     
 ########################################################################################################  
 
 nC1 <- nC2 <- VC$nC 
 
 
 Cop1 <- Cop2 <- VC$BivD 
 
 teta.ind1 <- as.logical(c(1,0,round(runif(VC$n-2))) ) 
 teta.ind2 <- teta.ind1 == FALSE 
 
 
 if(!(VC$BivD %in% VC$BivD2) && length(teta.st) > 1){
 
 teta.st1 <- teta.st[teta.ind1]
 teta.st2 <- teta.st[teta.ind2]
 
 teta1 <- teta[teta.ind1]
 teta2 <- teta[teta.ind2]
 
 }
 
 
  
 if(VC$BivD %in% VC$BivD2){
 
if(VC$BivD %in% VC$BivD2[c(1:4,13:16)])  teta.ind1 <- ifelse(VC$my.env$signind*teta > exp(VC$zerov), TRUE, FALSE)
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
 
 
   if( length(teta1) != 0) dH1 <- copgHs(p1[teta.ind1], p2[teta.ind1], eta1=NULL, eta2=NULL, teta1, teta.st1, Cop1, VC$dof, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
   if( length(teta2) != 0) dH2 <- copgHs(p1[teta.ind2], p2[teta.ind2], eta1=NULL, eta2=NULL, teta2, teta.st2, Cop2, VC$dof, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
   
   h <- NA
   if( length(teta1) != 0) h[teta.ind1] <- dH1$c.copula.be2 
   if( length(teta2) != 0) h[teta.ind2] <- dH2$c.copula.be2 
  
 

# LIK HERE

  
########################################################################################################
 

bit1.th2 <- c.copula.be1 <- c.copula.be2 <- c.copula.theta <- c.copula2.be1 <- c.copula2.be1th <- c.copula2.be2 <- c.copula2.be1be2 <- c.copula2.be2th <- NA
 

if( length(teta1) != 0){
 
  c.copula2.be2[teta.ind1]    <- dH1$c.copula2.be2 
  c.copula2.be1be2[teta.ind1] <- dH1$c.copula2.be1be2   
  c.copula2.be2th[teta.ind1]  <- dH1$c.copula2.be2th  
  c.copula.be1[teta.ind1]     <- dH1$c.copula.be1
  c.copula.be2[teta.ind1]     <- dH1$c.copula.be2
  c.copula.theta[teta.ind1]   <- dH1$c.copula.theta 
  c.copula2.be1[teta.ind1]    <- dH1$c.copula2.be1   
  c.copula2.be1th[teta.ind1]  <- dH1$c.copula2.be1th 
  bit1.th2[teta.ind1]         <- dH1$bit1.th2

}  
  

if( length(teta2) != 0){
 
  c.copula2.be2[teta.ind2]    <- dH2$c.copula2.be2 
  c.copula2.be1be2[teta.ind2] <- dH2$c.copula2.be1be2   
  c.copula2.be2th[teta.ind2]  <- dH2$c.copula2.be2th  
  c.copula.be1[teta.ind2]     <- dH2$c.copula.be1
  c.copula.be2[teta.ind2]     <- dH2$c.copula.be2
  c.copula.theta[teta.ind2]   <- dH2$c.copula.theta 
  c.copula2.be1[teta.ind2]    <- dH2$c.copula2.be1   
  c.copula2.be1th[teta.ind2]  <- dH2$c.copula2.be1th 
  bit1.th2[teta.ind2]         <- dH2$bit1.th2

    
}  
   
  
  
  derp1.dereta1    <- pd1$derp1.dereta1 
  
  derh.dereta1     <- c.copula2.be1be2 * derp1.dereta1
  derh.dereta2     <- c.copula2.be2 * derp2.dereta2





  # GRAD HERE


 
######################################################################################################## 
 

 
 der2h.derp2p2 <- der2h.derteta.teta.st <- derteta.derteta.st <- der2teta.derteta.stteta.st <- der2h.derp1p2 <- der2h.derp1teta <- der2h.derp2teta <- der2h.derp1p1 <- NA
 
  if( length(teta1) != 0) BITS1 <- copgHsCont(p1[teta.ind1], p2[teta.ind1], teta1, teta.st1, Cop1, par2 = VC$dof, nu.st = log(VC$dof-2))
  if( length(teta2) != 0) BITS2 <- copgHsCont(p1[teta.ind2], p2[teta.ind2], teta2, teta.st2, Cop2, par2 = VC$dof, nu.st = log(VC$dof-2)) 
 
 
if( length(teta1) != 0){ 
   
  der2h.derp2p2[teta.ind1]              <- BITS1$der2h.derp2p2 
  der2h.derteta.teta.st[teta.ind1]      <- BITS1$der2h.derteta.teta.st  
  derteta.derteta.st[teta.ind1]         <- BITS1$derteta.derteta.st 
  der2teta.derteta.stteta.st[teta.ind1] <- BITS1$der2teta.derteta.stteta.st  
  der2h.derp1p2[teta.ind1]              <- BITS1$der2h.derp1p2  
  der2h.derp1teta[teta.ind1]            <- BITS1$der2h.derp1teta                                     
  der2h.derp2teta[teta.ind1]            <- BITS1$der2h.derp2teta  
  der2h.derp1p1[teta.ind1]              <- BITS1$der2h.derp1p1
  
}


if( length(teta2) != 0){ 
   
  der2h.derp2p2[teta.ind2]              <- BITS2$der2h.derp2p2 
  der2h.derteta.teta.st[teta.ind2]      <- BITS2$der2h.derteta.teta.st  
  derteta.derteta.st[teta.ind2]         <- BITS2$derteta.derteta.st 
  der2teta.derteta.stteta.st[teta.ind2] <- BITS2$der2teta.derteta.stteta.st  
  der2h.derp1p2[teta.ind2]              <- BITS2$der2h.derp1p2  
  der2h.derp1teta[teta.ind2]            <- BITS2$der2h.derp1teta                                     
  der2h.derp2teta[teta.ind2]            <- BITS2$der2h.derp2teta  
  der2h.derp1p1[teta.ind2]              <- BITS2$der2h.derp1p1
  
}
  
  
 
  der2p1.dereta1eta1 <- pd1$der2p1.dereta1eta1     
                    
  der2h.dereta2.dereta2         <- der2h.derp2p2*derp2.dereta2^2 + c.copula2.be2*der2p2.dereta2eta2                                        
  der2h.derteta.st2             <- der2h.derteta.teta.st*derteta.derteta.st^2 + c.copula2.be2th/derteta.derteta.st* der2teta.derteta.stteta.st                                                                                        
  der2h.derp2dersigma2.st       <- der2h.derp2p2*derp2.dersigma.st 
  
  der2h.derp2dernu.st           <- der2h.derp2p2*derp2.dernu.st 

  der2h.dersigma2.st2           <- der2h.derp2dersigma2.st*derp2.dersigma.st + c.copula2.be2*der2p2.dersigma2.st2
  der2h.dernu.st2               <- der2h.derp2dernu.st*derp2.dernu.st        + c.copula2.be2*der2p2.dernu.st2
  
  derh.dersigma2.st             <- c.copula2.be2 * derp2.dersigma.st
  derh.dernu.st                 <- c.copula2.be2 * derp2.dernu.st
  
  der2h.dereta1.dereta2         <- der2h.derp1p2*derp1.dereta1*derp2.dereta2                                                                   
  der2h.dereta1.derteta.st      <- der2h.derp1teta*derp1.dereta1*derteta.derteta.st  
  der2h.dereta1.dersigma2.st    <- der2h.derp1p2 * derp2.dersigma.st*derp1.dereta1  
  
  der2h.dereta1.dernu.st        <- der2h.derp1p2 * derp2.dernu.st*derp1.dereta1 
  
  der2h.dereta2.derteta.st      <- der2h.derp2teta*derp2.dereta2*derteta.derteta.st 
  
  der2h.derteta.st.dersigma2.st <- der2h.derp2teta* derteta.derteta.st*derp2.dersigma.st 
  der2h.derteta.st.dernu.st     <- der2h.derp2teta* derteta.derteta.st*derp2.dernu.st 
  
  der2h.dersigma2.st.dernu.st   <- der2h.derp2dernu.st*derp2.dersigma.st + c.copula2.be2*der2p2.dersigma2.stdernu.st
  
  der2h.dereta2.dersigma2.st    <- der2h.derp2dersigma2.st*derp2.dereta2 + c.copula2.be2*der2p2.dereta2dersigma2.st  
  der2h.dereta2.dernu.st        <- der2h.derp2dernu.st*derp2.dereta2 + c.copula2.be2*der2p2.dereta2dernu.st 
  
  der2h.dereta1.dereta1         <- der2h.derp1p1*derp1.dereta1^2 + c.copula2.be1be2*der2p1.dereta1eta1      
  

# HESS here



############################################################################################################################################################
## 0 0 and 1 0 cases now
##

  
  p1   <- pd1$pr  # before p1 was 1 - p1, now I make this explicit


  p00 <- p10 <- NA


# no need to transform p2 as that is alredy evaluated at 0 when y2 = 0

if( length(teta1) != 0) p00[teta.ind1] <- mm(BiCDF(1-p1[teta.ind1], p2[teta.ind1], nC1, teta1, VC$dof), min.pr = VC$min.pr, max.pr = VC$max.pr  )
if( length(teta2) != 0) p00[teta.ind2] <- mm(BiCDF(1-p1[teta.ind2], p2[teta.ind2], nC2, teta2, VC$dof), min.pr = VC$min.pr, max.pr = VC$max.pr  )
                        p10            <- mm(p2 - p00, min.pr = VC$min.pr, max.pr = VC$max.pr)
  

l.par <- Y0*VC$weights*( y00*log(p00) + y10*log(p10) ) + YP*VC$weights*( y0p*log(h) + y1p*log(1 - h) + log(pdf2) )





dl.dbe1      <- Y0*VC$weights*( c.copula.be1*derp1.dereta1*(y00/p00 - y10/p10) )                          + YP*VC$weights*( derh.dereta1 *(y0p/h - y1p/(1-h)) ) 
dl.dbe2      <- Y0*VC$weights*( derp2.dereta2*( y00*c.copula.be2/p00 + y10*(1 - c.copula.be2)/p10 ) )     + YP*VC$weights*( derh.dereta2 * (y0p/h - y1p/(1-h)) + derpdf2.dereta2/pdf2 )
dl.dsigma.st <- Y0*VC$weights*( derp2.dersigma.st*( y00*c.copula.be2/p00 + y10*(1 - c.copula.be2)/p10 ) ) + YP*VC$weights*( c.copula2.be2 * derp2.dersigma.st *(y0p/h - y1p/(1-h)) + derpdf2.dersigma2.st/pdf2 )
dl.dnu.st    <- Y0*VC$weights*( derp2.dernu.st*( y00*c.copula.be2/p00 + y10*(1 - c.copula.be2)/p10 ) )    + YP*VC$weights*( c.copula2.be2 * derp2.dernu.st *(y0p/h - y1p/(1-h)) + derpdf2.dernu.st/pdf2 )
dl.dteta.st  <- Y0*VC$weights*( c.copula.theta*(y00/p00 - y10/p10) )                                      + YP*VC$weights*( c.copula2.be2th*(y0p/h - y1p/(1-h)) )                     
 



  d2l.be1.be1      <- -Y0*VC$weights*(  c.copula2.be1*derp1.dereta1^2*(y00/p00 - y10/p10) + c.copula.be1*der2p1.dereta1eta1*(y00/p00 - y10/p10) + c.copula.be1*derp1.dereta1*( -y00/p00^2*c.copula.be1*derp1.dereta1 - y10/p10^2*c.copula.be1*derp1.dereta1)     )    -YP*VC$weights*(der2h.dereta1.dereta1 *(y0p/h - y1p/(1-h)) - derh.dereta1^2 * (y0p/h^2 + y1p/(1-h)^2) )
  d2l.be2.be2      <- -Y0*VC$weights*( der2p2.dereta2eta2*(y00/p00*c.copula.be2 + y10/p10*(1-c.copula.be2)) + derp2.dereta2*(-y00/p00^2*derp2.dereta2*c.copula.be2^2 + y00/p00*c.copula2.be2*derp2.dereta2 - y10/p10^2*(derp2.dereta2-c.copula.be2*derp2.dereta2)*(1-c.copula.be2) - y10/p10*c.copula2.be2*derp2.dereta2   )  )  -YP*VC$weights*(der2h.dereta2.dereta2 *(y0p/h - y1p/(1-h)) - derh.dereta2^2 * (y0p/h^2 + y1p/(1-h)^2) + (der2pdf2.dereta2*pdf2-(derpdf2.dereta2)^2)/(pdf2)^2 )
  d2l.be1.be2      <- -Y0*VC$weights*( c.copula2.be1be2*derp1.dereta1*derp2.dereta2*(y00/p00 - y10/p10) + c.copula.be1*derp1.dereta1*(  -y00/p00^2*c.copula.be2*derp2.dereta2 + y10/p10^2*(derp2.dereta2-c.copula.be2*derp2.dereta2)   )   )      -YP*VC$weights*(der2h.dereta1.dereta2 *(y0p/h - y1p/(1-h)) - derh.dereta1*derh.dereta2 * (y0p/h^2 + y1p/(1-h)^2)  )
  d2l.sigma.sigma  <- -Y0*VC$weights*( der2p2.dersigma2.st2*(y00/p00*c.copula.be2 + y10/p10*(1-c.copula.be2)) + derp2.dersigma.st*(-y00/p00^2*derp2.dersigma.st*c.copula.be2^2 + y00/p00*c.copula2.be2*derp2.dersigma.st - y10/p10^2*(derp2.dersigma.st-c.copula.be2*derp2.dersigma.st)*(1-c.copula.be2) - y10/p10*c.copula2.be2*derp2.dersigma.st  )  )    -YP*VC$weights*(der2h.dersigma2.st2 *(y0p/h - y1p/(1-h)) - derh.dersigma2.st^2 * (y0p/h^2 + y1p/(1-h)^2) + (der2pdf2.dersigma2.st2*pdf2-(derpdf2.dersigma2.st)^2)/(pdf2)^2 )
  d2l.nu.nu        <- -Y0*VC$weights*( der2p2.dernu.st2*(y00/p00*c.copula.be2 + y10/p10*(1-c.copula.be2)) + derp2.dernu.st*(-y00/p00^2*derp2.dernu.st*c.copula.be2^2 + y00/p00*c.copula2.be2*derp2.dernu.st - y10/p10^2*(derp2.dernu.st-c.copula.be2*derp2.dernu.st)*(1-c.copula.be2) - y10/p10*c.copula2.be2*derp2.dernu.st  )  )    -YP*VC$weights*(der2h.dernu.st2 *(y0p/h - y1p/(1-h)) - derh.dernu.st^2 * (y0p/h^2 + y1p/(1-h)^2) + (der2pdf2.dernu.st2*pdf2-(derpdf2.dernu.st)^2)/(pdf2)^2 )
  d2l.be1.sigma    <- -Y0*VC$weights*( c.copula2.be1be2*derp1.dereta1*derp2.dersigma.st*(y00/p00 - y10/p10) + c.copula.be1*derp1.dereta1*(  -y00/p00^2*c.copula.be2*derp2.dersigma.st + y10/p10^2*(derp2.dersigma.st-c.copula.be2*derp2.dersigma.st)   )   )   -YP*VC$weights*(der2h.dereta1.dersigma2.st *(y0p/h - y1p/(1-h)) - derh.dereta1 * derh.dersigma2.st* (y0p/h^2 + y1p/(1-h)^2) )
  d2l.be1.nu       <- -Y0*VC$weights*( c.copula2.be1be2*derp1.dereta1*derp2.dernu.st*(y00/p00 - y10/p10) + c.copula.be1*derp1.dereta1*(  -y00/p00^2*c.copula.be2*derp2.dernu.st + y10/p10^2*(derp2.dernu.st-c.copula.be2*derp2.dernu.st)   )   )   -YP*VC$weights*(der2h.dereta1.dernu.st *(y0p/h - y1p/(1-h)) - derh.dereta1 * derh.dernu.st* (y0p/h^2 + y1p/(1-h)^2) )  
  d2l.be2.sigma    <- -Y0*VC$weights*(  der2p2.dereta2dersigma2.st*(y00/p00*c.copula.be2 + y10/p10*(1-c.copula.be2)) + derp2.dereta2*(-y00/p00^2*derp2.dersigma.st*c.copula.be2^2 + y00/p00*c.copula2.be2*derp2.dersigma.st - y10/p10^2*(derp2.dersigma.st-c.copula.be2*derp2.dersigma.st)*(1-c.copula.be2) - y10/p10*c.copula2.be2*derp2.dersigma.st   )  )   -YP*VC$weights*(der2h.dereta2.dersigma2.st *(y0p/h - y1p/(1-h)) - derh.dereta2*derh.dersigma2.st * (y0p/h^2 + y1p/(1-h)^2) + (der2pdf2.dereta2dersigma2.st*pdf2-(derpdf2.dereta2*derpdf2.dersigma2.st))/(pdf2)^2 ) 
  d2l.be2.nu       <- -Y0*VC$weights*(  der2p2.dereta2dernu.st*(y00/p00*c.copula.be2 + y10/p10*(1-c.copula.be2)) + derp2.dereta2*(-y00/p00^2*derp2.dernu.st*c.copula.be2^2 + y00/p00*c.copula2.be2*derp2.dernu.st - y10/p10^2*(derp2.dernu.st-c.copula.be2*derp2.dernu.st)*(1-c.copula.be2) - y10/p10*c.copula2.be2*derp2.dernu.st   )  )  -YP*VC$weights*(der2h.dereta2.dernu.st *(y0p/h - y1p/(1-h)) - derh.dereta2*derh.dernu.st * (y0p/h^2 + y1p/(1-h)^2) + (der2pdf2.dereta2dernu.st*pdf2-(derpdf2.dereta2*derpdf2.dernu.st))/(pdf2)^2 ) 
  d2l.rho.rho      <- -Y0*VC$weights*(  bit1.th2*(y00/p00 - y10/p10) + c.copula.theta^2*(-y00/p00^2 - y10/p10^2) ) -YP*VC$weights*(der2h.derteta.st2*(y0p/h - y1p/(1-h)) - c.copula2.be2th^2 * (y0p/h^2 + y1p/(1-h)^2) )
  d2l.be1.rho      <- -Y0*VC$weights*( c.copula2.be1th*derp1.dereta1*(y00/p00 - y10/p10) + c.copula.theta*( -y00/p00^2*c.copula.be1*derp1.dereta1 - y10/p10^2*c.copula.be1*derp1.dereta1)     ) -YP*VC$weights*(der2h.dereta1.derteta.st *(y0p/h - y1p/(1-h)) - derh.dereta1 * c.copula2.be2th* (y0p/h^2 + y1p/(1-h)^2) )
  d2l.be2.rho      <- -Y0*VC$weights*( c.copula2.be2th*derp2.dereta2*(y00/p00 - y10/p10) + c.copula.theta*( -y00/p00^2*c.copula.be2*derp2.dereta2 + y10/p10^2*(derp2.dereta2 - c.copula.be2*derp2.dereta2))     )  -YP*VC$weights*(der2h.dereta2.derteta.st *(y0p/h - y1p/(1-h)) - derh.dereta2*c.copula2.be2th * (y0p/h^2 + y1p/(1-h)^2)  )
  d2l.rho.sigma    <- -Y0*VC$weights*( c.copula2.be2th*derp2.dersigma.st*(y00/p00 - y10/p10) + c.copula.theta*( -y00/p00^2*c.copula.be2*derp2.dersigma.st + y10/p10^2*(derp2.dersigma.st - c.copula.be2*derp2.dersigma.st))     ) -YP*VC$weights*(der2h.derteta.st.dersigma2.st *(y0p/h - y1p/(1-h)) - derh.dersigma2.st * c.copula2.be2th* (y0p/h^2 + y1p/(1-h)^2) )
  d2l.rho.nu       <- -Y0*VC$weights*( c.copula2.be2th*derp2.dernu.st*(y00/p00 - y10/p10) + c.copula.theta*( -y00/p00^2*c.copula.be2*derp2.dernu.st + y10/p10^2*(derp2.dernu.st - c.copula.be2*derp2.dernu.st))     )    -YP*VC$weights*(der2h.derteta.st.dernu.st *(y0p/h - y1p/(1-h)) - derh.dernu.st * c.copula2.be2th* (y0p/h^2 + y1p/(1-h)^2) )
  d2l.sigma.nu     <- -Y0*VC$weights*(  der2p2.dersigma2.stdernu.st*(y00/p00*c.copula.be2 + y10/p10*(1-c.copula.be2)) + derp2.dersigma.st*(-y00/p00^2*derp2.dernu.st*c.copula.be2^2 + y00/p00*c.copula2.be2*derp2.dernu.st - y10/p10^2*(derp2.dernu.st-c.copula.be2*derp2.dernu.st)*(1-c.copula.be2) - y10/p10*c.copula2.be2*derp2.dernu.st   )  )  -YP*VC$weights*(der2h.dersigma2.st.dernu.st *(y0p/h - y1p/(1-h)) - derh.dersigma2.st * derh.dernu.st * (y0p/h^2 + y1p/(1-h)^2) + (der2pdf2.dersigma2.stdernu.st*pdf2-(derpdf2.dersigma2.st*derpdf2.dernu.st))/(pdf2)^2 )
     
         
         
         
         
         
         
         
         
         
         
         
if( is.null(VC$X3) ){


  be1.be1   <- crossprod(VC$X1*c(d2l.be1.be1),VC$X1)
  be2.be2   <- crossprod(VC$X2*c(d2l.be2.be2),VC$X2)
  be1.be2   <- crossprod(VC$X1*c(d2l.be1.be2),VC$X2)
  be1.rho   <- t(t(rowSums(t(VC$X1*c(d2l.be1.rho)))))
  be1.sigma <- t(t(rowSums(t(VC$X1*c(d2l.be1.sigma))))) 
  be2.rho   <- t(t(rowSums(t(VC$X2*c(d2l.be2.rho)))))
  be2.sigma <- t(t(rowSums(t(VC$X2*c(d2l.be2.sigma))))) 
  be2.nu <- t(t(rowSums(t(VC$X2*c(d2l.be2.nu)))))
  be1.nu <- t(t(rowSums(t(VC$X1*c(d2l.be1.nu)))))   
  


         
 
   H <- rbind( cbind( be1.be1     ,  be1.be2     ,  be1.sigma   , be1.nu,    be1.rho   ), 
               cbind( t(be1.be2)  ,  be2.be2     ,  be2.sigma   , be2.nu,    be2.rho   ), 
               cbind( t(be1.sigma),  t(be2.sigma),  sum(d2l.sigma.sigma) , sum(d2l.sigma.nu),  sum(d2l.rho.sigma) ),
               cbind( t(be1.nu),     t(be2.nu),     sum(d2l.sigma.nu) ,    sum(d2l.nu.nu), sum(d2l.rho.nu)  ),
               cbind( t(be1.rho)  ,  t(be2.rho)  ,  sum(d2l.rho.sigma), sum(d2l.rho.nu), sum(d2l.rho.rho)   ) 
               
               
             )  
         
         
         
         
         
  G   <- -c( colSums( c(dl.dbe1)*VC$X1 ) ,
             colSums( c(dl.dbe2)*VC$X2 ) ,
             sum( dl.dsigma.st ),
             sum( dl.dnu.st ),
             sum( dl.dteta.st )    )
    
}




if( !is.null(VC$X3) ){

  be1.be1   <- crossprod(VC$X1*c(d2l.be1.be1),VC$X1)
  be2.be2   <- crossprod(VC$X2*c(d2l.be2.be2),VC$X2)
  be1.be2   <- crossprod(VC$X1*c(d2l.be1.be2),VC$X2)
  
  be1.rho   <- crossprod(VC$X1*c(d2l.be1.rho),  VC$X5)                                     
  be1.sigma <- crossprod(VC$X1*c(d2l.be1.sigma),VC$X3)  
  be1.nu    <- crossprod(VC$X1*c(d2l.be1.nu),VC$X4)                                   

  be2.rho   <- crossprod(VC$X2*c(d2l.be2.rho),  VC$X5)                                     
  be2.sigma <- crossprod(VC$X2*c(d2l.be2.sigma),VC$X3) 
  be2.nu    <- crossprod(VC$X2*c(d2l.be2.nu),VC$X4)  
  
  
  sigma.sigma <- crossprod(VC$X3*c(d2l.sigma.sigma),VC$X3)
  sigma.nu    <- crossprod(VC$X3*c(d2l.sigma.nu),VC$X4)   
  
  sigma.rho   <- crossprod(VC$X3*c(d2l.rho.sigma),VC$X5) 
 
  rho.rho     <- crossprod(VC$X5*c(d2l.rho.rho),    VC$X5) 
  rho.nu      <- crossprod(VC$X4*c(d2l.rho.nu),    VC$X5) 
  nu.nu       <- crossprod(VC$X4*c(d2l.nu.nu),    VC$X4)    
  

  
  

  H <- rbind( cbind( be1.be1     ,  be1.be2     ,  be1.sigma   , be1.nu      , be1.rho   ), 
              cbind( t(be1.be2)  ,  be2.be2     ,  be2.sigma   , be2.nu      , be2.rho   ), 
              cbind( t(be1.sigma),  t(be2.sigma),  sigma.sigma , sigma.nu    , sigma.rho ),
              cbind( t(be1.nu)   ,  t(be2.nu)   ,  t(sigma.nu) , nu.nu       , rho.nu    ),
              cbind( t(be1.rho)  ,  t(be2.rho)  ,  t(sigma.rho), t(rho.nu)   , rho.rho   ) 
              
              
             )  
            
   
  G   <- -c( colSums(      c(dl.dbe1)*VC$X1 ) ,
             colSums(      c(dl.dbe2)*VC$X2 ) ,
             colSums( c(dl.dsigma.st)*VC$X3 ) ,
             colSums(    c(dl.dnu.st)*VC$X4 ) ,
             colSums(  c(dl.dteta.st)*VC$X5 ) )   
   
    
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
     eta1 = eta1, eta2 = eta2, etad = etad, etan = etan,
     dl.dbe1=dl.dbe1, dl.dbe2=dl.dbe2, dl.dsigma.st = dl.dsigma.st, 
     dl.dnu.st = dl.dnu.st, dl.dteta.st = dl.dteta.st,           
     BivD=VC$BivD, theta.star = teta.st,
     p1 = p1, p2 = p2, pdf1 = pdf1, pdf2 = pdf2,
     c.copula.be2 = c.copula.be2,
     c.copula.be1 = c.copula.be1,
     c.copula2.be1be2 = c.copula2.be1be2,
     teta.ind2 = teta.ind2, teta.ind1 = teta.ind1,
     Cop1 = Cop1, Cop2 = Cop2, teta1 = teta1, teta2 = teta2)      



}




     























