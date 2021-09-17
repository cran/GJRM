bprobgHsDiscr1ROY <- function(params, respvec, VC, ps, AT = FALSE){

 p1 <- p2 <- pdf1 <- pdf2 <- c.copula.be2 <- c.copula.be1 <- c.copula2.be1be2 <- l.par <- dl.dbe1 <- d2l.be1.be1 <- NA

   eta1 <-         VC$X1%*%params[1:VC$X1.d2]
   eta2 <- eta.tr( VC$X2%*%params[(VC$X1.d2+1):(VC$X1.d2+VC$X2.d2)],                   VC$margins[2])
   eta3 <- eta.tr( VC$X3%*%params[(VC$X1.d2+VC$X2.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2)], VC$margins[3]) 
   
   
   etad1 <- etad2 <- etas1 <- etas2 <- etan1 <- etan2 <- l.ln <- NULL 
    
   
    teta.st1 <- etad1 <- VC$X4%*%params[(VC$X1.d2+VC$X2.d2+VC$X3.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2+VC$X4.d2)]
    teta.st2 <- etad2 <- VC$X5%*%params[(VC$X1.d2+VC$X2.d2+VC$X3.d2+VC$X4.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2+VC$X4.d2+VC$X5.d2)]
 
 
    pd1 <- probm(eta1, VC$margins[1], bc = TRUE, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr) 
        
    p1 <- pd1$pr    #        
    p0 <- 1 - p1    #       
    
    derp1.dereta1      <- pd1$derp1.dereta1  # recall that these derivs are wrt to 1 - p1
    der2p1.dereta1eta1 <- pd1$der2p1.dereta1eta1


  dHs <- distrHsDiscr(respvec$y2, eta2, 1, 1, nu = 1, nu.st = 1, margin2=VC$margins[2], naive = FALSE, y2m = VC$y2m, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
   
   
   pdf2.M2               <- dHs$pdf2
   p2.M2                 <- dHs$p2 
   derpdf2.dereta2.M2    <- dHs$derpdf2.dereta2 
   derp2.dereta2.M2      <- dHs$derp2.dereta2
   der2p2.dereta2eta2.M2 <- dHs$der2p2.dereta2eta2 
   der2pdf2.dereta2.M2   <- dHs$der2pdf2.dereta2
   
   
  #******** SET UP y3m  ************ DONE
   
  dHs <- distrHsDiscr(respvec$y3, eta3, 1, 1, nu = 1, nu.st = 1, margin2=VC$margins[3], naive = FALSE, y2m = VC$y3m, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
   
   
   pdf2.M3               <- dHs$pdf2
   p2.M3                 <- dHs$p2 
   derpdf2.dereta2.M3    <- dHs$derpdf2.dereta2 
   derp2.dereta2.M3      <- dHs$derp2.dereta2
   der2p2.dereta2eta2.M3 <- dHs$der2p2.dereta2eta2 
   der2pdf2.dereta2.M3   <- dHs$der2pdf2.dereta2   
    
    
  ########################################################################################################  
    
    
  VC1   <- list(BivD = VC$BivD1, BivD2 = NULL)  
  resT1 <- teta.tr(VC1, teta.st1)
  
  VC2   <- list(BivD = VC$BivD2, BivD2 = NULL)
  resT2 <- teta.tr(VC2, teta.st2) # ********** careful here the correct copula is selected DONE

  teta.st1 <- resT1$teta.st
  teta.st2 <- resT2$teta.st

  teta1 <- resT1$teta
  teta2 <- resT2$teta 
    
  ##################

  Cop1 <- VC$BivD1
  Cop2 <- VC$BivD2 

  nC1 <- VC$nC1
  nC2 <- VC$nC2    
    
  ########################################################################################################
  
    C1.M2 <- mm(BiCDF(p0[VC$inde0], p2.M2,                                                       nC1, teta1, VC$dof1), min.pr = VC$min.pr, max.pr = VC$max.pr  )
    C2.M2 <- mm(BiCDF(p0[VC$inde0], mm(p2.M2 - pdf2.M2, min.pr = VC$min.pr, max.pr = VC$max.pr), nC1, teta1, VC$dof1), min.pr = VC$min.pr, max.pr = VC$max.pr  )
    
    A.M2  <- mm(C1.M2 - C2.M2, min.pr = VC$min.pr, max.pr = VC$max.pr)
    
    
    C1.M3 <- mm(BiCDF(p0[VC$inde1], p2.M3,                                                       nC2, teta2, VC$dof2), min.pr = VC$min.pr, max.pr = VC$max.pr  )
    C2.M3 <- mm(BiCDF(p0[VC$inde1], mm(p2.M3 - pdf2.M3, min.pr = VC$min.pr, max.pr = VC$max.pr), nC2, teta2, VC$dof2), min.pr = VC$min.pr, max.pr = VC$max.pr  )
    
    A.M3  <- mm(C1.M3 - C2.M3, min.pr = VC$min.pr, max.pr = VC$max.pr)    
    
    B.M3  <- mm( pdf2.M3 - A.M3, min.pr = VC$min.pr, max.pr = VC$max.pr)
    

    l.par[VC$inde0] <- log( A.M2 )
    l.par[VC$inde1] <- log( B.M3 )  
    
    l.par <- VC$weights*l.par 
    
      
  ########################################################################################################
     

    dH1F.M2 <- copgHs(p0[VC$inde0], p2.M2,                                                       eta1=NULL, eta2=NULL, teta1, teta.st1, Cop1, VC$dof1, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)  
    dH2F.M2 <- copgHs(p0[VC$inde0], mm(p2.M2 - pdf2.M2, min.pr = VC$min.pr, max.pr = VC$max.pr), eta1=NULL, eta2=NULL, teta1, teta.st1, Cop1, VC$dof1, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr) 
    
    c.copula.be1.C1.M2   <- dH1F.M2$c.copula.be1 
    c.copula.be1.C2.M2   <- dH2F.M2$c.copula.be1 
    
    c.copula.be2.C1.M2   <- dH1F.M2$c.copula.be2 
    c.copula.be2.C2.M2   <- dH2F.M2$c.copula.be2 
    
    c.copula.theta.C1.M2 <- dH1F.M2$c.copula.theta 
    c.copula.theta.C2.M2 <- dH2F.M2$c.copula.theta
    
    derp2m1.dereta2.M2   <- derp2.dereta2.M2 - derpdf2.dereta2.M2
    
    Cc.M2    <- c.copula.be1.C1.M2 - c.copula.be1.C2.M2 # mm(c.copula.be1.C1.M2 - c.copula.be1.C2.M2, min.pr = VC$min.pr, max.pr = VC$max.pr) 
    C.M2     <- Cc.M2*derp1.dereta1[VC$inde0] # contains already -
    
    Cs.M2    <- c.copula.theta.C1.M2  - c.copula.theta.C2.M2

    Ceta2 <- (c.copula.be2.C1.M2 - c.copula.be2.C2.M2)*derp2.dereta2.M2 + c.copula.be2.C2.M2*derpdf2.dereta2.M2  # mm(c.copula.be2.C1.M2 - c.copula.be2.C2.M2, min.pr = VC$min.pr, max.pr = VC$max.pr)




    dH1F.M3 <- copgHs(p0[VC$inde1], p2.M3,                                                       eta1=NULL, eta2=NULL, teta2, teta.st2, Cop2, VC$dof2, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)  
    dH2F.M3 <- copgHs(p0[VC$inde1], mm(p2.M3 - pdf2.M3, min.pr = VC$min.pr, max.pr = VC$max.pr), eta1=NULL, eta2=NULL, teta2, teta.st2, Cop2, VC$dof2, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr) 
    
    c.copula.be1.C1.M3   <- dH1F.M3$c.copula.be1 
    c.copula.be1.C2.M3   <- dH2F.M3$c.copula.be1 
    
    c.copula.be2.C1.M3   <- dH1F.M3$c.copula.be2 
    c.copula.be2.C2.M3   <- dH2F.M3$c.copula.be2 
    
    c.copula.theta.C1.M3 <- dH1F.M3$c.copula.theta 
    c.copula.theta.C2.M3 <- dH2F.M3$c.copula.theta
    
    derp2m1.dereta2.M3   <- derp2.dereta2.M3 - derpdf2.dereta2.M3
    
    Cc.M3    <- c.copula.be1.C1.M3 - c.copula.be1.C2.M3  # mm(c.copula.be1.C1.M3 - c.copula.be1.C2.M3, min.pr = VC$min.pr, max.pr = VC$max.pr) 
    C.M3     <- Cc.M3*derp1.dereta1[VC$inde1] # contains already -
    
    Cs.M3    <- c.copula.theta.C1.M3  - c.copula.theta.C2.M3
    
    
    Ceta3 <- derpdf2.dereta2.M3 - ((c.copula.be2.C1.M3 -  c.copula.be2.C2.M3)*derp2.dereta2.M3 + c.copula.be2.C2.M3*derpdf2.dereta2.M3)
    
    
    # ************ check from above what is really needed and not
    
    
    
      dl.dbe1[VC$inde0] <-  C.M2/A.M2
      dl.dbe1[VC$inde1] <- -C.M3/B.M3
      dl.dbe1           <-  VC$weights*dl.dbe1
 

      dl.dbe2      <- VC$weights[VC$inde0]*( Ceta2/A.M2)      
      dl.dbe3      <- VC$weights[VC$inde1]*( Ceta3/B.M3)
      
      
      
      dl.dteta1.st <- VC$weights[VC$inde0]*( Cs.M2/A.M2)  
      dl.dteta2.st <- VC$weights[VC$inde1]*(-Cs.M3/B.M3)                    




    ######################################################################################################## 
     
      
      c.copula2.be1.C1.M2           <- dH1F.M2$c.copula2.be1
      c.copula2.be1.C2.M2           <- dH2F.M2$c.copula2.be1
      
      c.copula2.be2.C1.M2           <- dH1F.M2$c.copula2.be2
      c.copula2.be2.C2.M2           <- dH2F.M2$c.copula2.be2
      
      c.copula2.be1be2.C1.M2        <- dH1F.M2$c.copula2.be1be2
      c.copula2.be1be2.C2.M2        <- dH2F.M2$c.copula2.be1be2
      
      c.copula2.be2th.C1.M2         <- dH1F.M2$c.copula2.be2th # with star
      c.copula2.be2th.C2.M2         <- dH2F.M2$c.copula2.be2th  
      
      c.copula2.theta.C1.M2         <- dH1F.M2$bit1.th2ATE # no start
      c.copula2.theta.C2.M2         <- dH2F.M2$bit1.th2ATE 
      
      c.copula.thet.C1.M2           <- dH1F.M2$c.copula.thet # NO star
      c.copula.thet.C2.M2           <- dH2F.M2$c.copula.thet    
     
      derteta.derteta.st.M2         <- dH1F.M2$derteta.derteta.st         # does not matter dH1 or dH2 
      der2teta.derteta.stteta.st.M2 <- dH1F.M2$der2teta.derteta.stteta.st   
     
      c.copula2.be1th.C1.M2         <- dH1F.M2$c.copula2.be1th 
      c.copula2.be1th.C2.M2         <- dH2F.M2$c.copula2.be1th    
     
      c.copula2.be1.C1.M3           <- dH1F.M3$c.copula2.be1
      c.copula2.be1.C2.M3           <- dH2F.M3$c.copula2.be1
      
      c.copula2.be2.C1.M3           <- dH1F.M3$c.copula2.be2
      c.copula2.be2.C2.M3           <- dH2F.M3$c.copula2.be2
      
      c.copula2.be1be2.C1.M3        <- dH1F.M3$c.copula2.be1be2
      c.copula2.be1be2.C2.M3        <- dH2F.M3$c.copula2.be1be2
      
      c.copula2.be2th.C1.M3         <- dH1F.M3$c.copula2.be2th
      c.copula2.be2th.C2.M3         <- dH2F.M3$c.copula2.be2th  
      
      c.copula2.theta.C1.M3         <- dH1F.M3$bit1.th2ATE 
      c.copula2.theta.C2.M3         <- dH2F.M3$bit1.th2ATE 
      
      c.copula.thet.C1.M3           <- dH1F.M3$c.copula.thet # NO star
      c.copula.thet.C2.M3           <- dH2F.M3$c.copula.thet    
     
      derteta.derteta.st.M3         <- dH1F.M3$derteta.derteta.st         # does not matter dH1 or dH2 
      der2teta.derteta.stteta.st.M3 <- dH1F.M3$der2teta.derteta.stteta.st   
     
      c.copula2.be1th.C1.M3         <- dH1F.M3$c.copula2.be1th 
      c.copula2.be1th.C2.M3         <- dH2F.M3$c.copula2.be1th   
  
      der2p2m1.dereta2eta2.M2 <- der2p2.dereta2eta2.M2 - der2pdf2.dereta2.M2
      der2p2m1.dereta2eta2.M3 <- der2p2.dereta2eta2.M3 - der2pdf2.dereta2.M3

##########################

      b1b1CY <- (c.copula2.be1.C1.M2 - c.copula2.be1.C2.M2)*derp1.dereta1[VC$inde0]^2 + Cc.M2*der2p1.dereta1eta1[VC$inde0]
      b1b1Y  <- (c.copula2.be1.C1.M3 - c.copula2.be1.C2.M3)*derp1.dereta1[VC$inde1]^2 + Cc.M3*der2p1.dereta1eta1[VC$inde1]  
      d2l.be1.be1[VC$inde0] <-  b1b1CY/A.M2 - C.M2^2/A.M2^2
      d2l.be1.be1[VC$inde1] <- -b1b1Y/B.M3  - C.M3^2/B.M3^2
      d2l.be1.be1           <- -VC$weights*d2l.be1.be1         # ok verified 
      

     
      b3b3Y  <- der2pdf2.dereta2.M3 - ( c.copula2.be2.C1.M3*derp2.dereta2.M3^2 + c.copula.be2.C1.M3*der2p2.dereta2eta2.M3 - (c.copula2.be2.C2.M3*derp2m1.dereta2.M3^2 + c.copula.be2.C2.M3*der2p2m1.dereta2eta2.M3)  ) 
      d2l.be3.be3 <- -VC$weights[VC$inde1]*( b3b3Y/B.M3  - Ceta3^2/B.M3^2    ) # ok verified
                 

      b2b2CY <- c.copula2.be2.C1.M2*derp2.dereta2.M2^2 + c.copula.be2.C1.M2*der2p2.dereta2eta2.M2 - (c.copula2.be2.C2.M2*derp2m1.dereta2.M2^2 + c.copula.be2.C2.M2*der2p2m1.dereta2eta2.M2)   
      d2l.be2.be2 <- -VC$weights[VC$inde0]*( b2b2CY/A.M2 - Ceta2^2/A.M2^2    ) # ok verified
      
            
      
      b1b3 <- -(c.copula2.be1be2.C1.M3*derp2.dereta2.M3 - c.copula2.be1be2.C2.M3*derp2m1.dereta2.M3)*derp1.dereta1[VC$inde1]
      d2l.be1.be3 <- -VC$weights[VC$inde1]*(b1b3/B.M3 - -C.M3*Ceta3/B.M3^2  )  # ok verified
      
      
      b1t2 <- -(c.copula2.be1th.C1.M3 - c.copula2.be1th.C2.M3)*derp1.dereta1[VC$inde1] 
      d2l.be1.th2 <- -VC$weights[VC$inde1]*( b1t2/B.M3 - -C.M3*-Cs.M3/B.M3^2  ) # ok verified
     
     
      t1t1CY <- c.copula2.theta.C1.M2*derteta.derteta.st.M2^2 + c.copula.thet.C1.M2*der2teta.derteta.stteta.st.M2 - (c.copula2.theta.C2.M2*derteta.derteta.st.M2^2 + c.copula.thet.C2.M2*der2teta.derteta.stteta.st.M2 ) 
      d2l.th1.th1 <- -VC$weights[VC$inde0]*( t1t1CY/A.M2 - Cs.M2^2/A.M2^2    )  # ok verified
     
           
      

     
     
     
     
      t2t2Y  <- -(c.copula2.theta.C1.M3*derteta.derteta.st.M3^2 + c.copula.thet.C1.M3*der2teta.derteta.stteta.st.M3 - (c.copula2.theta.C2.M3*derteta.derteta.st.M3^2 + c.copula.thet.C2.M3*der2teta.derteta.stteta.st.M3)  ) 
      d2l.th2.th2 <- -VC$weights[VC$inde1]*(t2t2Y/B.M3 - Cs.M3^2/B.M3^2    ) # ok but then results are not good


      
      
      
      b1t1 <-  (c.copula2.be1th.C1.M2 - c.copula2.be1th.C2.M2)*derp1.dereta1[VC$inde0] 
      d2l.be1.th1 <- -VC$weights[VC$inde0]*( b1t1/A.M2 - C.M2*Cs.M2/A.M2^2  ) # looks ok





     


      b1b2 <-  (c.copula2.be1be2.C1.M2*derp2.dereta2.M2 - c.copula2.be1be2.C2.M2*derp2m1.dereta2.M2)*derp1.dereta1[VC$inde0]
      d2l.be1.be2 <- -VC$weights[VC$inde0]*( b1b2/A.M2 - C.M2*Ceta2/A.M2^2 ) # looks ok
      
     
     
     
     
      
      b2t1 <-   c.copula2.be2th.C1.M2*derp2.dereta2.M2 - c.copula2.be2th.C2.M2*derp2m1.dereta2.M2
      d2l.be2.th1 <- -VC$weights[VC$inde0]*( b2t1/A.M2 -  Ceta2*Cs.M2/A.M2^2 ) # looks ok
      
      
      
      
      
      b3t2 <- -(c.copula2.be2th.C1.M3*derp2.dereta2.M3 - c.copula2.be2th.C2.M3*derp2m1.dereta2.M3)
      d2l.be3.th2 <- -VC$weights[VC$inde1]*( b3t2/B.M3 - -Ceta3*Cs.M3/B.M3^2 ) # looks ok





      d2l.th1.th2 <- 0
      d2l.be2.be3 <- 0 
      d2l.be2.th2 <- 0
      d2l.be3.th1 <- 0
      
      
      

  be1.be1 <- crossprod(VC$X1*c(d2l.be1.be1),VC$X1)
  be2.be2 <- crossprod(VC$X2*c(d2l.be2.be2),VC$X2)
  be3.be3 <- crossprod(VC$X3*c(d2l.be3.be3),VC$X3)
  th1.th1 <- crossprod(VC$X4*c(d2l.th1.th1),VC$X4)
  th2.th2 <- crossprod(VC$X5*c(d2l.th2.th2),VC$X5)
  be1.be2 <- crossprod(VC$X1[VC$inde0,]*c(d2l.be1.be2),VC$X2)
  be1.be3 <- crossprod(VC$X1[VC$inde1,]*c(d2l.be1.be3),VC$X3)
  be1.th1 <- crossprod(VC$X1[VC$inde0,]*c(d2l.be1.th1),VC$X4)
  be1.th2 <- crossprod(VC$X1[VC$inde1,]*c(d2l.be1.th2),VC$X5)
  be2.th1 <- crossprod(VC$X2*c(d2l.be2.th1),VC$X4)
  be3.th2 <- crossprod(VC$X3*c(d2l.be3.th2),VC$X5)
  
  
 
  th1.th2 <- matrix(0, dim(VC$X4)[2], dim(VC$X5)[2])
  be2.be3 <- matrix(0, dim(VC$X2)[2], dim(VC$X3)[2])
  be2.th2 <- matrix(0, dim(VC$X2)[2], dim(VC$X5)[2])
  be3.th1 <- matrix(0, dim(VC$X3)[2], dim(VC$X4)[2])



  
  H <- rbind( cbind(   be1.be1 ,   be1.be2 ,   be1.be3 ,   be1.th1,  be1.th2 ), 
              cbind( t(be1.be2),   be2.be2 ,   be2.be3 ,   be2.th1,  be2.th2 ), 
              cbind( t(be1.be3), t(be2.be3),   be3.be3 ,   be3.th1,  be3.th2 ),
              cbind( t(be1.th1), t(be2.th1), t(be3.th1),   th1.th1,  th1.th2 ),
              cbind( t(be1.th2), t(be2.th2), t(be3.th2), t(th1.th2), th2.th2 )               
            ) 
            
           
  G   <- -c( colSums(     c(dl.dbe1)*VC$X1),
             colSums(     c(dl.dbe2)*VC$X2),
             colSums(     c(dl.dbe3)*VC$X3),
             colSums(c(dl.dteta1.st)*VC$X4), 
             colSums(c(dl.dteta2.st)*VC$X5) 
             
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
  
  
  
         list(value = res, gradient = G, hessian = H, S.h = S.h, S.h1 = S.h1, S.h2 = S.h2, l = S.res, l.par = l.par, ps = ps, l.ln = NULL,
              eta1 = eta1, eta2 = eta2, eta3 = eta3, etad1 = etad1, 
              etad2 = etad2, etas1 = etas1, etas2 = etas2, etan1 = etan1, etan2 = etan2,
              dl.dbe1 = dl.dbe1, dl.dbe2 = dl.dbe2, dl.dbe3 = dl.dbe3, dl.dteta1.st = dl.dteta1.st, dl.dteta2.st = dl.dteta2.st,
              BivD1 = VC$BivD1, BivD2 = VC$BivD2,                              
              p1 = p1, p0 = p0,           
              teta.st1 = teta.st1, teta.st2 = teta.st2,
              Cop1 = Cop1, Cop2 = Cop2, teta1 = teta1, teta2 = teta2)      

}

