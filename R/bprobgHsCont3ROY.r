bprobgHsCont3ROY <- function(params, respvec, VC, ps, AT = FALSE){

 p1 <- p2 <- pdf1 <- pdf2 <- c.copula.be2 <- c.copula.be1 <- c.copula2.be1be2 <- l.par <- dl.dbe1 <- d2l.be1.be1 <- l.parL <- NA
 etad1 <- etad2 <- etas1 <- etas2 <- etan1 <- etan2 <- l.ln <- NULL 

 ########

   eta1 <- VC$X1%*%params[1:VC$X1.d2]
   eta2 <- eta.tr( VC$X2%*%params[(VC$X1.d2+1):(VC$X1.d2+VC$X2.d2)],                   VC$margins[2])
   eta3 <- eta.tr( VC$X3%*%params[(VC$X1.d2+VC$X2.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2)], VC$margins[3]) 

   sigma1.st <- VC$X4%*%params[(VC$X1.d2+VC$X2.d2+VC$X3.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2+VC$X4.d2)]
   sigma2.st <- VC$X5%*%params[(VC$X1.d2+VC$X2.d2+VC$X3.d2+VC$X4.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2+VC$X4.d2+VC$X5.d2)]   

   nu1.st <- VC$X6%*%params[(VC$X1.d2+VC$X2.d2+VC$X3.d2+VC$X4.d2+VC$X5.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2+VC$X4.d2+VC$X5.d2 + VC$X6.d2)]
   nu2.st <- VC$X7%*%params[(VC$X1.d2+VC$X2.d2+VC$X3.d2+VC$X4.d2+VC$X5.d2 + VC$X6.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2+VC$X4.d2+VC$X5.d2 + VC$X6.d2+VC$X7.d2)]  

   teta.st1 <- VC$X8%*%params[(VC$X1.d2+VC$X2.d2+VC$X3.d2+VC$X4.d2+VC$X5.d2 + VC$X6.d2+VC$X7.d2 + 1):(VC$X1.d2+VC$X2.d2+VC$X3.d2+VC$X4.d2+VC$X5.d2 + VC$X6.d2 + VC$X7.d2 + VC$X8.d2)]
   teta.st2 <- VC$X9%*%params[(VC$X1.d2+VC$X2.d2+VC$X3.d2+VC$X4.d2+VC$X5.d2 + VC$X6.d2+VC$X7.d2 + VC$X8.d2 + 1):(VC$X1.d2+VC$X2.d2+VC$X3.d2+VC$X4.d2+VC$X5.d2 + VC$X6.d2 + VC$X7.d2 + VC$X8.d2 + VC$X9.d2)]
 
 ########
 
   pd1 <- probm(eta1, VC$margins[1], bc = TRUE, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr) 
        
   p1 <- pd1$pr            
   p0 <- 1 - p1           
    
   derp1.dereta1      <- pd1$derp1.dereta1      # recall that these derivs are wrt to 1 - p1
   der2p1.dereta1eta1 <- pd1$der2p1.dereta1eta1


   s1tr      <- esp.tr(sigma1.st, VC$margins[2])  
   sigma1.st <- s1tr$vrb.st 
   sigma1    <- s1tr$vrb 

   s2tr      <- esp.tr(sigma2.st, VC$margins[3])  
   sigma2.st <- s2tr$vrb.st 
   sigma2    <- s2tr$vrb 
   
   n1tr   <- enu.tr(nu1.st, VC$margins[2])  
   nu1.st <- n1tr$vrb.st 
   nu1    <- n1tr$vrb    

   n2tr   <- enu.tr(nu2.st, VC$margins[3])  
   nu2.st <- n2tr$vrb.st 
   nu2    <- n2tr$vrb

   VC1   <- list(BivD = VC$BivD1, BivD2 = NULL)  
   resT1 <- teta.tr(VC1, teta.st1)
   
   VC2   <- list(BivD = VC$BivD2, BivD2 = NULL)
   resT2 <- teta.tr(VC2, teta.st2) 

   teta.st1 <- resT1$teta.st
   teta.st2 <- resT2$teta.st
   teta1    <- resT1$teta
   teta2    <- resT2$teta 
    
   Cop1 <- VC$BivD1
   Cop2 <- VC$BivD2 
   nC1  <- VC$nC1
   nC2  <- VC$nC2    
  
########  


 dHs  <- distrHs(respvec$y2, eta2, sigma1, sigma1.st, nu1, nu1.st, margin2=VC$margins[2], naive = FALSE, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
  
 pdf2.M2                          <- dHs$pdf2
 p2.M2                            <- dHs$p2 
 derpdf2.dereta2.M2               <- dHs$derpdf2.dereta2 
 derpdf2.dersigma2.st.M2          <- dHs$derpdf2.dersigma2.st 
 derp2.dersigma.st.M2             <- dHs$derp2.dersigma.st
 derpdf2.dernu.st.M2              <- dHs$derpdf2.dernu.st 
 derp2.dernu.st.M2                <- dHs$derp2.nu.st
 derp2.dereta2.M2                 <- dHs$derp2.dereta2
 der2p2.dereta2eta2.M2            <- dHs$der2p2.dereta2eta2 
 der2pdf2.dereta2.M2              <- dHs$der2pdf2.dereta2
 der2p2.dersigma2.st2.M2          <- dHs$der2p2.dersigma2.st2
 der2pdf2.dersigma2.st2.M2        <- dHs$der2pdf2.dersigma2.st2
 der2p2.dernu.st2.M2              <- dHs$der2p2.dernu.st2
 der2pdf2.dernu.st2.M2            <- dHs$der2pdf2.dernu.st2
 der2p2.dereta2dersigma2.st.M2    <- dHs$der2p2.dereta2dersigma2.st            
 der2pdf2.dereta2dersigma2.st.M2  <- dHs$der2pdf2.dereta2dersigma2.st  
 der2p2.dereta2dernu.st.M2        <- dHs$der2p2.dereta2dernu.st            
 der2pdf2.dereta2dernu.st.M2      <- dHs$der2pdf2.dereta2dernu.st 
 der2p2.dersigma2.stdernu.st.M2   <- dHs$der2p2.dersigma2.stdernu.st            
 der2pdf2.dersigma2.stdernu.st.M2 <- dHs$der2pdf2.sigma2.st2dernu.st 
 
 
 dHs  <- distrHs(respvec$y3, eta3, sigma2, sigma2.st, nu2, nu2.st, margin2=VC$margins[2], naive = FALSE, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
  
 pdf2.M3                          <- dHs$pdf2
 p2.M3                            <- dHs$p2 
 derpdf2.dereta2.M3               <- dHs$derpdf2.dereta2 
 derpdf2.dersigma2.st.M3          <- dHs$derpdf2.dersigma2.st 
 derp2.dersigma.st.M3             <- dHs$derp2.dersigma.st
 derpdf2.dernu.st.M3              <- dHs$derpdf2.dernu.st 
 derp2.dernu.st.M3                <- dHs$derp2.nu.st
 derp2.dereta2.M3                 <- dHs$derp2.dereta2
 der2p2.dereta2eta2.M3            <- dHs$der2p2.dereta2eta2 
 der2pdf2.dereta2.M3              <- dHs$der2pdf2.dereta2
 der2p2.dersigma2.st2.M3          <- dHs$der2p2.dersigma2.st2
 der2pdf2.dersigma2.st2.M3        <- dHs$der2pdf2.dersigma2.st2
 der2p2.dernu.st2.M3              <- dHs$der2p2.dernu.st2
 der2pdf2.dernu.st2.M3            <- dHs$der2pdf2.dernu.st2
 der2p2.dereta2dersigma2.st.M3    <- dHs$der2p2.dereta2dersigma2.st            
 der2pdf2.dereta2dersigma2.st.M3  <- dHs$der2pdf2.dereta2dersigma2.st  
 der2p2.dereta2dernu.st.M3        <- dHs$der2p2.dereta2dernu.st            
 der2pdf2.dereta2dernu.st.M3      <- dHs$der2pdf2.dereta2dernu.st 
 der2p2.dersigma2.stdernu.st.M3   <- dHs$der2p2.dersigma2.stdernu.st            
 der2pdf2.dersigma2.stdernu.st.M3 <- dHs$der2pdf2.sigma2.st2dernu.st  
 
 
########

  dH.M2 <- copgHs(p0[VC$inde0], p2.M2, eta1=NULL, eta2=NULL, teta1, teta.st1, Cop1, VC$dof1, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
  dH.M3 <- copgHs(p0[VC$inde1], p2.M3, eta1=NULL, eta2=NULL, teta2, teta.st2, Cop2, VC$dof2, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
  
  h.M2 <- dH.M2$c.copula.be2 
  h.M3 <- dH.M3$c.copula.be2 
  
  l.par[VC$inde0] <- log( h.M2 )     + log(pdf2.M2)
  l.par[VC$inde1] <- log( 1 - h.M3 ) + log(pdf2.M3)  
      
  l.par <- VC$weights*l.par 
 
########

  c.copula2.be2.M2    <- dH.M2$c.copula2.be2 
  c.copula2.be1be2.M2 <- dH.M2$c.copula2.be1be2   
  c.copula2.be2th.M2  <- dH.M2$c.copula2.be2th  
  
  c.copula2.be2.M3    <- dH.M3$c.copula2.be2 
  c.copula2.be1be2.M3 <- dH.M3$c.copula2.be1be2   
  c.copula2.be2th.M3  <- dH.M3$c.copula2.be2th     # th star is in
  

## GRADIENT

  dl.dbe1[VC$inde0] <- 1/h.M2*c.copula2.be1be2.M2*derp1.dereta1[VC$inde0] 
  dl.dbe1[VC$inde1] <- 1/(1 - h.M3)*-c.copula2.be1be2.M3*derp1.dereta1[VC$inde1] 
  dl.dbe1           <- VC$weights*dl.dbe1
  
  dl.dbe2       <- VC$weights[VC$inde0]*(      1/h.M2*c.copula2.be2.M2*derp2.dereta2.M2     + 1/pdf2.M2*derpdf2.dereta2.M2)
  dl.dsigma1.st <- VC$weights[VC$inde0]*(      1/h.M2*c.copula2.be2.M2*derp2.dersigma.st.M2 + 1/pdf2.M2*derpdf2.dersigma2.st.M2)
  dl.dnu1.st    <- VC$weights[VC$inde0]*(      1/h.M2*c.copula2.be2.M2*derp2.dernu.st.M2    + 1/pdf2.M2*derpdf2.dernu.st.M2)

  dl.dbe3       <- VC$weights[VC$inde1]*( 1/(1 - h.M3)*-c.copula2.be2.M3*derp2.dereta2.M3     + 1/pdf2.M3*derpdf2.dereta2.M3 )
  dl.dsigma2.st <- VC$weights[VC$inde1]*( 1/(1 - h.M3)*-c.copula2.be2.M3*derp2.dersigma.st.M3 + 1/pdf2.M3*derpdf2.dersigma2.st.M3 )
  dl.dnu2.st    <- VC$weights[VC$inde1]*( 1/(1 - h.M3)*-c.copula2.be2.M3*derp2.dernu.st.M3    + 1/pdf2.M3*derpdf2.dernu.st.M3 )
  
  dl.dteta1.st  <- VC$weights[VC$inde0]*( 1/h.M2*c.copula2.be2th.M2 )                     
  dl.dteta2.st  <- VC$weights[VC$inde1]*( 1/(1 - h.M3)*-c.copula2.be2th.M3 )                    


  ########################################################################################################


  BITS.M2 <- copgHsCont(p0[VC$inde0], p2.M2, teta1, teta.st1, Cop1, par2 = VC$dof1, nu.st = log(VC$dof1-2))
  BITS.M3 <- copgHsCont(p0[VC$inde1], p2.M3, teta2, teta.st2, Cop2, par2 = VC$dof2, nu.st = log(VC$dof2-2))
  
   der2h.derp2p2.M2              <- BITS.M2$der2h.derp2p2 
   der2h.derteta.teta.st.M2      <- BITS.M2$der2h.derteta.teta.st  
   derteta.derteta.st.M2         <- BITS.M2$derteta.derteta.st 
   der2teta.derteta.stteta.st.M2 <- BITS.M2$der2teta.derteta.stteta.st  
   der2h.derp1p2.M2              <- BITS.M2$der2h.derp1p2  
   der2h.derp1teta.M2            <- BITS.M2$der2h.derp1teta                                     
   der2h.derp2teta.M2            <- BITS.M2$der2h.derp2teta  
   der2h.derp1p1.M2              <- BITS.M2$der2h.derp1p1
   der2h.derteta.st2.M2          <- der2h.derteta.teta.st.M2*derteta.derteta.st.M2^2 + c.copula2.be2th.M2/derteta.derteta.st.M2*der2teta.derteta.stteta.st.M2                                                                                        
   
   der2h.derp2p2.M3              <- BITS.M3$der2h.derp2p2 
   der2h.derteta.teta.st.M3      <- BITS.M3$der2h.derteta.teta.st  
   derteta.derteta.st.M3         <- BITS.M3$derteta.derteta.st 
   der2teta.derteta.stteta.st.M3 <- BITS.M3$der2teta.derteta.stteta.st  
   der2h.derp1p2.M3              <- BITS.M3$der2h.derp1p2  
   der2h.derp1teta.M3            <- BITS.M3$der2h.derp1teta                                     
   der2h.derp2teta.M3            <- BITS.M3$der2h.derp2teta  
   der2h.derp1p1.M3              <- BITS.M3$der2h.derp1p1
   der2h.derteta.st2.M3          <- der2h.derteta.teta.st.M3*derteta.derteta.st.M3^2 + c.copula2.be2th.M3/derteta.derteta.st.M3*der2teta.derteta.stteta.st.M3                                                                                        
   
 
        d2l.be1.be1[VC$inde0] <- ( der2h.derp1p1.M2*derp1.dereta1[VC$inde0]^2 +  c.copula2.be1be2.M2*der2p1.dereta1eta1[VC$inde0])/h.M2       - ( c.copula2.be1be2.M2*derp1.dereta1[VC$inde0])^2/h.M2^2
        d2l.be1.be1[VC$inde1] <- (-der2h.derp1p1.M3*derp1.dereta1[VC$inde1]^2 + -c.copula2.be1be2.M3*der2p1.dereta1eta1[VC$inde1])/(1 - h.M3) - (-c.copula2.be1be2.M3*derp1.dereta1[VC$inde1])^2/(1 - h.M3)^2 
        d2l.be1.be1           <- -VC$weights*d2l.be1.be1         
      
	d2l.be2.be2 <- -VC$weights[VC$inde0]*( (der2h.derp2p2.M2*derp2.dereta2.M2^2     + c.copula2.be2.M2*der2p2.dereta2eta2.M2  )/h.M2 - (c.copula2.be2.M2*derp2.dereta2.M2    )^2/h.M2^2 +       der2pdf2.dereta2.M2/pdf2.M2 -      derpdf2.dereta2.M2^2/pdf2.M2^2 )  
	d2l.si1.si1 <- -VC$weights[VC$inde0]*( (der2h.derp2p2.M2*derp2.dersigma.st.M2^2 + c.copula2.be2.M2*der2p2.dersigma2.st2.M2)/h.M2 - (c.copula2.be2.M2*derp2.dersigma.st.M2)^2/h.M2^2 + der2pdf2.dersigma2.st2.M2/pdf2.M2 - derpdf2.dersigma2.st.M2^2/pdf2.M2^2 )  
	d2l.nu1.nu1 <- -VC$weights[VC$inde0]*( (der2h.derp2p2.M2*derp2.dernu.st.M2^2    + c.copula2.be2.M2*der2p2.dernu.st2.M2    )/h.M2 - (c.copula2.be2.M2*derp2.dernu.st.M2   )^2/h.M2^2 +     der2pdf2.dernu.st2.M2/pdf2.M2 -     derpdf2.dernu.st.M2^2/pdf2.M2^2 )

	d2l.be3.be3 <- -VC$weights[VC$inde1]*(  (-der2h.derp2p2.M3*derp2.dereta2.M3^2     + -c.copula2.be2.M3*der2p2.dereta2eta2.M3  )/(1 - h.M3) - (-c.copula2.be2.M3*derp2.dereta2.M3    )^2/(1 - h.M3)^2 + der2pdf2.dereta2.M3/pdf2.M3       - derpdf2.dereta2.M3^2/pdf2.M3^2    )    
	d2l.si2.si2 <- -VC$weights[VC$inde1]*(  (-der2h.derp2p2.M3*derp2.dersigma.st.M3^2 + -c.copula2.be2.M3*der2p2.dersigma2.st2.M3)/(1 - h.M3) - (-c.copula2.be2.M3*derp2.dersigma.st.M3)^2/(1 - h.M3)^2 + der2pdf2.dersigma2.st2.M3/pdf2.M3 - derpdf2.dersigma2.st.M3^2/pdf2.M3^2    )    
	d2l.nu2.nu2 <- -VC$weights[VC$inde1]*(  (-der2h.derp2p2.M3*derp2.dernu.st.M3^2    + -c.copula2.be2.M3*der2p2.dernu.st2.M3    )/(1 - h.M3) - (-c.copula2.be2.M3*derp2.dernu.st.M3   )^2/(1 - h.M3)^2 + der2pdf2.dernu.st2.M3/pdf2.M3     - derpdf2.dernu.st.M3^2/pdf2.M3^2    )  
	
	d2l.th1.th1 <- -VC$weights[VC$inde0]*(  der2h.derteta.st2.M2/h.M2       - c.copula2.be2th.M2^2/h.M2^2   )  
	d2l.th2.th2 <- -VC$weights[VC$inde1]*( -der2h.derteta.st2.M3/(1 - h.M3) - (-c.copula2.be2th.M3)^2/(1 - h.M3)^2  ) 
	
	d2l.be2.si1 <- -VC$weights[VC$inde0]*(   (der2h.derp2p2.M2*derp2.dereta2.M2*derp2.dersigma.st.M2 +  c.copula2.be2.M2*der2p2.dereta2dersigma2.st.M2  )/h.M2       -    c.copula2.be2.M2^2*derp2.dereta2.M2*derp2.dersigma.st.M2/h.M2^2    + der2pdf2.dereta2dersigma2.st.M2/pdf2.M2 - derpdf2.dereta2.M2*derpdf2.dersigma2.st.M2/pdf2.M2^2 )    
	d2l.si1.nu1 <- -VC$weights[VC$inde0]*(   (der2h.derp2p2.M2*derp2.dernu.st.M2*derp2.dersigma.st.M2 +  c.copula2.be2.M2*der2p2.dersigma2.stdernu.st.M2  )/h.M2       -    c.copula2.be2.M2^2*derp2.dernu.st.M2*derp2.dersigma.st.M2/h.M2^2 + der2pdf2.dersigma2.stdernu.st.M2/pdf2.M2 - derpdf2.dernu.st.M2*derpdf2.dersigma2.st.M2/pdf2.M2^2 ) 

	d2l.be2.nu1 <- -VC$weights[VC$inde0]*(   (der2h.derp2p2.M2*derp2.dereta2.M2*derp2.dernu.st.M2    +  c.copula2.be2.M2*der2p2.dereta2dernu.st.M2      )/h.M2       -    c.copula2.be2.M2^2*derp2.dereta2.M2*derp2.dernu.st.M2/h.M2^2           + der2pdf2.dereta2dernu.st.M2/pdf2.M2 -     derpdf2.dereta2.M2*derpdf2.dernu.st.M2/pdf2.M2^2 ) 
	
	d2l.be3.si2 <- -VC$weights[VC$inde1]*(  (-der2h.derp2p2.M3*derp2.dereta2.M3*derp2.dersigma.st.M3 + -c.copula2.be2.M3*der2p2.dereta2dersigma2.st.M3  )/(1 - h.M3) - (-c.copula2.be2.M3)^2*derp2.dereta2.M3*derp2.dersigma.st.M3/(1 - h.M3)^2 + der2pdf2.dereta2dersigma2.st.M3/pdf2.M3 - derpdf2.dereta2.M3*derpdf2.dersigma2.st.M3/pdf2.M3^2 )      

	d2l.be3.nu2 <- -VC$weights[VC$inde1]*(  (-der2h.derp2p2.M3*derp2.dereta2.M3*derp2.dernu.st.M3     + -c.copula2.be2.M3*der2p2.dereta2dernu.st.M3      )/(1 - h.M3) - (-c.copula2.be2.M3)^2*derp2.dereta2.M3*derp2.dernu.st.M3/(1 - h.M3)^2    + der2pdf2.dereta2dernu.st.M3/pdf2.M3     - derpdf2.dereta2.M3*derpdf2.dernu.st.M3/pdf2.M3^2 )   
        d2l.si2.nu2 <- -VC$weights[VC$inde1]*(  (-der2h.derp2p2.M3*derp2.dersigma.st.M3*derp2.dernu.st.M3 + -c.copula2.be2.M3*der2p2.dersigma2.stdernu.st.M3      )/(1 - h.M3) - (-c.copula2.be2.M3)^2*derp2.dersigma.st.M3*derp2.dernu.st.M3/(1 - h.M3)^2    + der2pdf2.dersigma2.stdernu.st.M3/pdf2.M3     - derpdf2.dersigma2.st.M3*derpdf2.dernu.st.M3/pdf2.M3^2 )   
        
	d2l.be1.be2 <- -VC$weights[VC$inde0]*( der2h.derp1p2.M2*derp2.dereta2.M2*derp1.dereta1[VC$inde0]/h.M2     - c.copula2.be2.M2*derp2.dereta2.M2*c.copula2.be1be2.M2*derp1.dereta1[VC$inde0]/h.M2^2 )
	d2l.be1.si1 <- -VC$weights[VC$inde0]*( der2h.derp1p2.M2*derp2.dersigma.st.M2*derp1.dereta1[VC$inde0]/h.M2 - c.copula2.be2.M2*derp2.dersigma.st.M2*c.copula2.be1be2.M2*derp1.dereta1[VC$inde0]/h.M2^2  )
	d2l.be1.nu1 <- -VC$weights[VC$inde0]*( der2h.derp1p2.M2*derp2.dernu.st.M2*derp1.dereta1[VC$inde0]/h.M2    - c.copula2.be2.M2*derp2.dernu.st.M2*c.copula2.be1be2.M2*derp1.dereta1[VC$inde0]/h.M2^2  )

	d2l.be1.be3 <- -VC$weights[VC$inde1]*( -der2h.derp1p2.M3*derp1.dereta1[VC$inde1]*derp2.dereta2.M3/(1 - h.M3)     - -c.copula2.be2.M3*derp2.dereta2.M3*-c.copula2.be1be2.M3*derp1.dereta1[VC$inde1]/(1 - h.M3)^2 )  
	d2l.be1.si2 <- -VC$weights[VC$inde1]*( -der2h.derp1p2.M3*derp1.dereta1[VC$inde1]*derp2.dersigma.st.M3/(1 - h.M3) - -c.copula2.be2.M3*derp2.dersigma.st.M3*-c.copula2.be1be2.M3*derp1.dereta1[VC$inde1]/(1 - h.M3)^2 )	
	d2l.be1.nu2 <- -VC$weights[VC$inde1]*( -der2h.derp1p2.M3*derp1.dereta1[VC$inde1]*derp2.dernu.st.M3/(1 - h.M3)    - -c.copula2.be2.M3*derp2.dernu.st.M3*-c.copula2.be1be2.M3*derp1.dereta1[VC$inde1]/(1 - h.M3)^2 )	
        
	d2l.be3.th2 <- -VC$weights[VC$inde1]*( -der2h.derp2teta.M3*derteta.derteta.st.M3*derp2.dereta2.M3/(1 - h.M3)     - -c.copula2.be2.M3*derp2.dereta2.M3*-c.copula2.be2th.M3/(1 - h.M3)^2 )
	d2l.si2.th2 <- -VC$weights[VC$inde1]*( -der2h.derp2teta.M3*derteta.derteta.st.M3*derp2.dersigma.st.M3/(1 - h.M3) - -c.copula2.be2.M3*derp2.dersigma.st.M3*-c.copula2.be2th.M3/(1 - h.M3)^2 )    
	d2l.nu2.th2 <- -VC$weights[VC$inde1]*( -der2h.derp2teta.M3*derteta.derteta.st.M3*derp2.dernu.st.M3/(1 - h.M3)    - -c.copula2.be2.M3*derp2.dernu.st.M3*-c.copula2.be2th.M3/(1 - h.M3)^2 ) 

	d2l.be2.th1 <- -VC$weights[VC$inde0]*( der2h.derp2teta.M2*derteta.derteta.st.M2*derp2.dereta2.M2/h.M2     -     c.copula2.be2.M2*derp2.dereta2.M2*c.copula2.be2th.M2/h.M2^2  ) 
	d2l.si1.th1 <- -VC$weights[VC$inde0]*( der2h.derp2teta.M2*derteta.derteta.st.M2*derp2.dersigma.st.M2/h.M2 - c.copula2.be2.M2*derp2.dersigma.st.M2*c.copula2.be2th.M2/h.M2^2  )      
        d2l.nu1.th1 <- -VC$weights[VC$inde0]*( der2h.derp2teta.M2*derteta.derteta.st.M2*derp2.dernu.st.M2/h.M2    -    c.copula2.be2.M2*derp2.dernu.st.M2*c.copula2.be2th.M2/h.M2^2  ) 

	d2l.be1.th1 <- -VC$weights[VC$inde0]*( der2h.derp1teta.M2*derteta.derteta.st.M2*derp1.dereta1[VC$inde0]/h.M2 - c.copula2.be2th.M2*c.copula2.be1be2.M2*derp1.dereta1[VC$inde0]/h.M2^2  )
	d2l.be1.th2 <- -VC$weights[VC$inde1]*( -der2h.derp1teta.M3*derteta.derteta.st.M3*derp1.dereta1[VC$inde1]/(1 - h.M3) - -c.copula2.be2th.M3*-c.copula2.be1be2.M3*derp1.dereta1[VC$inde1]/(1 - h.M3)^2 )    
	
	
	
	
	
	#d2l.be2.th2 <- 0      
	#d2l.be3.si1 <- 0    
	#d2l.be3.th1 <- 0   
	#d2l.si1.si2 <- 0    
	#d2l.si1.th2 <- 0    
	#d2l.si2.th1 <- 0    
	#d2l.th1.th2 <- 0   
 	#d2l.be2.be3 <- 0 
	#d2l.be2.si2 <- 0
     
  be1.be1 <- crossprod(VC$X1*c(d2l.be1.be1),           VC$X1)
  be1.be2 <- crossprod(VC$X1[VC$inde0,]*c(d2l.be1.be2),VC$X2)
  be1.be3 <- crossprod(VC$X1[VC$inde1,]*c(d2l.be1.be3),VC$X3)
  be1.si1 <- crossprod(VC$X1[VC$inde0,]*c(d2l.be1.si1),VC$X4)
  be1.si2 <- crossprod(VC$X1[VC$inde1,]*c(d2l.be1.si2),VC$X5)
  be1.nu1 <- crossprod(VC$X1[VC$inde0,]*c(d2l.be1.nu1),VC$X6)
  be1.nu2 <- crossprod(VC$X1[VC$inde1,]*c(d2l.be1.nu2),VC$X7)
  be1.th1 <- crossprod(VC$X1[VC$inde0,]*c(d2l.be1.th1),VC$X8)
  be1.th2 <- crossprod(VC$X1[VC$inde1,]*c(d2l.be1.th2),VC$X9)
  
  be2.be2 <-    crossprod(VC$X2*c(d2l.be2.be2),VC$X2)
  be2.be3 <- matrix(0,dim(VC$X2)[2],       dim(VC$X3)[2])
  be2.si1 <-    crossprod(VC$X2*c(d2l.be2.si1),VC$X4)
  be2.si2 <- matrix(0,dim(VC$X2)[2],       dim(VC$X5)[2]) 
  be2.nu1 <-    crossprod(VC$X2*c(d2l.be2.nu1),VC$X6)
  be2.nu2 <- matrix(0,dim(VC$X2)[2],       dim(VC$X7)[2])  
  be2.th1 <-    crossprod(VC$X2*c(d2l.be2.th1),VC$X8)  
  be2.th2 <- matrix(0,dim(VC$X2)[2],       dim(VC$X9)[2])   
  
  be3.be3 <-    crossprod(VC$X3*c(d2l.be3.be3),VC$X3)
  be3.si1 <- matrix(0,dim(VC$X3)[2],       dim(VC$X4)[2])
  be3.si2 <-    crossprod(VC$X3*c(d2l.be3.si2),VC$X5)
  be3.nu1 <- matrix(0,dim(VC$X3)[2],       dim(VC$X6)[2])   
  be3.nu2 <-    crossprod(VC$X3*c(d2l.be3.nu2),VC$X7)   
  be3.th1 <- matrix(0,dim(VC$X3)[2],       dim(VC$X8)[2]) 
  be3.th2 <-    crossprod(VC$X3*c(d2l.be3.th2),VC$X9)   
  
  si1.si1 <-    crossprod(VC$X4*c(d2l.si1.si1),VC$X4)
  si1.si2 <- matrix(0,dim(VC$X4)[2],       dim(VC$X5)[2])
  si1.nu1 <-    crossprod(VC$X4*c(d2l.si1.nu1),VC$X6)
  si1.nu2 <- matrix(0,dim(VC$X4)[2],       dim(VC$X7)[2])  
  si1.th1 <-    crossprod(VC$X4*c(d2l.si1.th1),VC$X8)
  si1.th2 <- matrix(0,dim(VC$X4)[2],       dim(VC$X9)[2])   
  
  si2.si2 <-    crossprod(VC$X5*c(d2l.si2.si2),VC$X5)
  si2.nu1 <- matrix(0,dim(VC$X5)[2],       dim(VC$X6)[2])
  si2.nu2 <-    crossprod(VC$X5*c(d2l.si2.nu2),VC$X7)
  si2.th1 <- matrix(0,dim(VC$X5)[2],       dim(VC$X8)[2])
  si2.th2 <-    crossprod(VC$X5*c(d2l.si2.th2),VC$X9)    
  
  nu1.nu1 <-    crossprod(VC$X6*c(d2l.nu1.nu1),VC$X6)
  nu1.nu2 <- matrix(0,dim(VC$X6)[2],       dim(VC$X7)[2]) 
  nu1.th1 <-    crossprod(VC$X6*c(d2l.nu1.th1),VC$X8)
  nu1.th2 <- matrix(0,dim(VC$X6)[2],       dim(VC$X9)[2])  
  
  nu2.nu2 <-    crossprod(VC$X7*c(d2l.nu2.nu2),VC$X7)
  nu2.th1 <- matrix(0,dim(VC$X7)[2],       dim(VC$X8)[2]) 
  nu2.th2 <-    crossprod(VC$X7*c(d2l.nu2.th2),VC$X9)  
  
  th1.th1 <-    crossprod(VC$X8*c(d2l.th1.th1),VC$X8)  
  th1.th2 <- matrix(0,dim(VC$X8)[2],       dim(VC$X9)[2])   

  th2.th2 <- crossprod(VC$X9*c(d2l.th2.th2),VC$X9) 
  
  

  
  H <- rbind( cbind(   be1.be1 ,   be1.be2 ,   be1.be3 ,   be1.si1 ,   be1.si2 ,   be1.nu1 ,   be1.nu2 ,   be1.th1 ,  be1.th2 ),
              cbind( t(be1.be2),   be2.be2 ,   be2.be3 ,   be2.si1 ,   be2.si2 ,   be2.nu1 ,   be2.nu2 ,   be2.th1 ,  be2.th2 ), 
              cbind( t(be1.be3), t(be2.be3),   be3.be3 ,   be3.si1 ,   be3.si2 ,   be3.nu1 ,   be3.nu2 ,   be3.th1 ,  be3.th2 ),
              cbind( t(be1.si1), t(be2.si1), t(be3.si1),   si1.si1 ,   si1.si2 ,   si1.nu1 ,   si1.nu2 ,   si1.th1 ,  si1.th2 ),
              cbind( t(be1.si2), t(be2.si2), t(be3.si2), t(si1.si2),   si2.si2 ,   si2.nu1 ,   si2.nu2 ,   si2.th1 ,  si2.th2 ),
              cbind( t(be1.nu1), t(be2.nu1), t(be3.nu1), t(si1.nu1), t(si2.nu1),   nu1.nu1 ,   nu1.nu2 ,   nu1.th1 ,  nu1.th2 ),
              cbind( t(be1.nu2), t(be2.nu2), t(be3.nu2), t(si1.nu2), t(si2.nu2), t(nu1.nu2),   nu2.nu2 ,   nu2.th1 ,  nu2.th2 ),
              cbind( t(be1.th1), t(be2.th1), t(be3.th1), t(si1.th1), t(si2.th1), t(nu1.th1), t(nu2.th1),   th1.th1 ,  th1.th2 ),
              cbind( t(be1.th2), t(be2.th2), t(be3.th2), t(si1.th2), t(si2.th2), t(nu1.th2), t(nu2.th2), t(th1.th2),  th2.th2 )
            ) 
            
            
            
           
  G   <- -c( colSums(c(dl.dbe1      )*VC$X1),
             colSums(c(dl.dbe2      )*VC$X2),
             colSums(c(dl.dbe3      )*VC$X3),
             colSums(c(dl.dsigma1.st)*VC$X4),
             colSums(c(dl.dsigma2.st)*VC$X5),
             colSums(c(dl.dnu1.st   )*VC$X6),
             colSums(c(dl.dnu2.st   )*VC$X7),             
             colSums(c(dl.dteta1.st )*VC$X8), 
             colSums(c(dl.dteta2.st )*VC$X9) 
             
             ) 
             
             
         res <- -sum(l.par)

if(VC$extra.regI == "pC" && VC$hess==FALSE) H <- regH(H, type = 1)
  
  S.h <- ps$S.h
  
  if( length(S.h) != 1){
  
    S.h1 <- 0.5*crossprod(params,S.h)%*%params
    S.h2 <- S.h%*%params
  
  } else S.h <- S.h1 <- S.h2 <- 0   
  
  S.res <- res
  res   <- S.res + S.h1
  G     <- G + S.h2
  H     <- H + S.h 
      
if(VC$extra.regI == "sED") H <- regH(H, type = 2) 



  
         list(value = res, gradient = G, hessian = H, S.h = S.h, S.h1 = S.h1, S.h2 = S.h2, l = S.res, l.par = l.par, ps = ps,
              eta1 = eta1, eta2 = eta2, eta3 = eta3, sigma1.st = sigma1.st, sigma2.st = sigma2.st, nu1.st = nu1.st, nu2.st = nu2.st,
              teta.st1 = teta.st1, teta.st2 = teta.st2, teta1 = teta1, teta2 = teta2,
              dl.dbe1 = dl.dbe1, dl.dbe2 = dl.dbe2, dl.dbe3 = dl.dbe3, dl.dsigma1.st = dl.dsigma1.st, dl.dsigma2.st = dl.dsigma2.st,
              dl.dteta1.st = dl.dteta1.st, dl.dteta2.st = dl.dteta2.st, dl.dnu1.st = dl.dnu1.st, dl.dnu2.st = dl.dnu2.st,
              BivD1 = VC$BivD1, BivD2 = VC$BivD2,                              
              p1 = p1, p0 = p0,           
              Cop1 = Cop1, Cop2 = Cop2)      
              
         
}









