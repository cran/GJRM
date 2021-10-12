bprobgHsBinROY <- function(params, respvec, VC, ps, AT = FALSE){

 p1 <- p2 <- pdf1 <- pdf2 <- c.copula.be2 <- c.copula.be1 <- c.copula2.be1be2 <- l.par <- dl.dbe1 <- d2l.be1.be1 <- NA
 etad1 <- etad2 <- etas1 <- etas2 <- etan1 <- etan2 <- l.ln <- NULL 


   eta1 <- VC$X1%*%params[1:VC$X1.d2]
   eta2 <- VC$X2%*%params[(VC$X1.d2+1):(VC$X1.d2+VC$X2.d2)]
   eta3 <- VC$X3%*%params[(VC$X1.d2+VC$X2.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2)] 
   
   
   teta.st1 <- VC$X4%*%params[(VC$X1.d2+VC$X2.d2+VC$X3.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2+VC$X4.d2)]
   teta.st2 <- VC$X5%*%params[(VC$X1.d2+VC$X2.d2+VC$X3.d2+VC$X4.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2+VC$X4.d2+VC$X5.d2)]
 
   VC1   <- list(BivD = VC$BivD1, BivD2 = NULL)  
   resT1 <- teta.tr(VC1, teta.st1)
  
   VC2   <- list(BivD = VC$BivD2, BivD2 = NULL)
   resT2 <- teta.tr(VC2, teta.st2) 

   teta.st1 <- resT1$teta.st
   teta.st2 <- resT2$teta.st

   teta1 <- resT1$teta
   teta2 <- resT2$teta 
    
   Cop1 <- VC$BivD1
   Cop2 <- VC$BivD2 

   nC1 <- VC$nC1
   nC2 <- VC$nC2    
 
 
  pd1 <- probm(eta1, VC$margins[1], bc = TRUE, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
  pd2 <- probm(eta2, VC$margins[2], bc = TRUE, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
  pd3 <- probm(eta3, VC$margins[3], bc = TRUE, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)

  p1 <- pd1$pr
  p2 <- pd2$pr
  p3 <- pd3$pr
       
  derp1.dereta1.M1      <- pd1$derp1.dereta1      # recall that these derivs are wrt to 1 - p1
  der2p1.dereta1eta1.M1 <- pd1$der2p1.dereta1eta1
  derp1.dereta1.M2      <- pd2$derp1.dereta1      
  der2p1.dereta1eta1.M2 <- pd2$der2p1.dereta1eta1    
  derp1.dereta1.M3      <- pd3$derp1.dereta1     
  der2p1.dereta1eta1.M3 <- pd3$der2p1.dereta1eta1    
    
    
  p00 <- mm(BiCDF(1 - p1[VC$inde0], 1 - p2, nC1, teta1, VC$dof1), min.pr = VC$min.pr, max.pr = VC$max.pr )
  p01 <- (1 - p1[VC$inde0]) - p00
  
  p10 <- (1 - p3) - mm(BiCDF(1 - p1[VC$inde1], 1 - p3, nC2, teta2, VC$dof2), min.pr = VC$min.pr, max.pr = VC$max.pr )
  p11 <- p1[VC$inde1] - p10

  l.par[VC$y10.y20] <- log(p00[VC$y10.y20R]) 
  l.par[VC$y10.y21] <- log(p01[VC$y10.y21R]) 
  l.par[VC$y11.y30] <- log(p10[VC$y11.y30R]) 
  l.par[VC$y11.y31] <- log(p11[VC$y11.y31R]) 
  
  l.par <- VC$weights*l.par 

  ########################################################################################################
  
 dH <- copgHs(1 - p1[VC$inde0], 1 - p2, eta1 = NULL, eta2 = NULL, teta1, teta.st1, Cop1, VC$dof1, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
 
 c.copula.be1.M2     <- dH$c.copula.be1
 c.copula.be2.M2     <- dH$c.copula.be2
 c.copula.theta.M2   <- dH$c.copula.theta # has star in it
 c.copula2.be1.M2    <- dH$c.copula2.be1   
 c.copula2.be2.M2    <- dH$c.copula2.be2 
 c.copula2.be1be2.M2 <- dH$c.copula2.be1be2
 c.copula2.be1th.M2  <- dH$c.copula2.be1th  # has star in it
 c.copula2.be2th.M2  <- dH$c.copula2.be2th
 bit1.th2.M2         <- dH$bit1.th2       # second deriv with star in it
 
 dH <- copgHs(1 - p1[VC$inde1], 1 - p3, eta1 = NULL, eta2 = NULL, teta2, teta.st2, Cop2, VC$dof2, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
 
 c.copula.be1.M3     <- dH$c.copula.be1
 c.copula.be2.M3     <- dH$c.copula.be2
 c.copula.theta.M3   <- dH$c.copula.theta 
 c.copula2.be1.M3    <- dH$c.copula2.be1   
 c.copula2.be2.M3    <- dH$c.copula2.be2 
 c.copula2.be1be2.M3 <- dH$c.copula2.be1be2
 c.copula2.be1th.M3  <- dH$c.copula2.be1th 
 c.copula2.be2th.M3  <- dH$c.copula2.be2th
 bit1.th2.M3         <- dH$bit1.th2 
 
 
    
      dl.dbe1 <- dl.dbe2 <- dl.dbe3 <- dl.dteta1.st <- dl.dteta2.st <- NA

      dl.dbe1[VC$y10.y20] <- (1/p00*c.copula.be1.M2*derp1.dereta1.M1[VC$inde0])[VC$y10.y20R] 
      dl.dbe1[VC$y10.y21] <- (1/p01*(derp1.dereta1.M1[VC$inde0] - c.copula.be1.M2*derp1.dereta1.M1[VC$inde0]))[VC$y10.y21R]
      dl.dbe1[VC$y11.y30] <- (1/p10*-c.copula.be1.M3*derp1.dereta1.M1[VC$inde1])[VC$y11.y30R]
      dl.dbe1[VC$y11.y31] <- (1/p11*(-derp1.dereta1.M1[VC$inde1] + c.copula.be1.M3*derp1.dereta1.M1[VC$inde1]))[VC$y11.y31R] 
      dl.dbe1             <-  VC$weights*dl.dbe1
 
      dl.dbe2[VC$y10.y20R] <- (1/p00*c.copula.be2.M2*derp1.dereta1.M2)[VC$y10.y20R] 
      dl.dbe2[VC$y10.y21R] <- (1/p01*-c.copula.be2.M2*derp1.dereta1.M2)[VC$y10.y21R]
      dl.dbe2              <- VC$weights[VC$inde0]*dl.dbe2
      
      dl.dbe3[VC$y11.y30R] <- (1/p10*(derp1.dereta1.M3 - c.copula.be2.M3*derp1.dereta1.M3))[VC$y11.y30R] 
      dl.dbe3[VC$y11.y31R] <- (1/p11*(-derp1.dereta1.M3 + c.copula.be2.M3*derp1.dereta1.M3))[VC$y11.y31R]
      dl.dbe3              <- VC$weights[VC$inde1]*dl.dbe3      
      
      dl.dteta1.st[VC$y10.y20R] <-  (1/p00*c.copula.theta.M2)[VC$y10.y20R]
      dl.dteta1.st[VC$y10.y21R] <- (1/p01*-c.copula.theta.M2)[VC$y10.y21R]
      dl.dteta1.st              <- VC$weights[VC$inde0]*dl.dteta1.st

      dl.dteta2.st[VC$y11.y30R] <- (1/p10*-c.copula.theta.M3)[VC$y11.y30R]
      dl.dteta2.st[VC$y11.y31R] <- (1/p11*c.copula.theta.M3)[VC$y11.y31R]
      dl.dteta2.st              <- VC$weights[VC$inde1]*dl.dteta2.st                    



      d2l.be1.be1 <- d2l.be2.be2 <- d2l.be3.be3 <- d2l.be1.th1 <- NA  
      d2l.th1.th1 <- d2l.th2.th2 <- d2l.be1.be2 <- d2l.be1.be3 <- d2l.be1.th2 <- d2l.be3.th2 <- d2l.be2.th1 <-  NA    

      d2l.be1.be1[VC$y10.y20] <- ( (c.copula2.be1.M2*derp1.dereta1.M1[VC$inde0]^2 + c.copula.be1.M2*der2p1.dereta1eta1.M1[VC$inde0])/p00 - (c.copula.be1.M2*derp1.dereta1.M1[VC$inde0])^2/p00^2 )[VC$y10.y20R] 
      d2l.be1.be1[VC$y10.y21] <- ( (der2p1.dereta1eta1.M1[VC$inde0] - c.copula2.be1.M2*derp1.dereta1.M1[VC$inde0]^2 + -c.copula.be1.M2*der2p1.dereta1eta1.M1[VC$inde0])/p01 - (derp1.dereta1.M1[VC$inde0] - c.copula.be1.M2*derp1.dereta1.M1[VC$inde0])^2/p01^2 )[VC$y10.y21R]
      d2l.be1.be1[VC$y11.y30] <- ( (-c.copula2.be1.M3*derp1.dereta1.M1[VC$inde1]^2 + -c.copula.be1.M3*der2p1.dereta1eta1.M1[VC$inde1])/p10 -  (-c.copula.be1.M3*derp1.dereta1.M1[VC$inde1])^2/p10^2 )[VC$y11.y30R]
      d2l.be1.be1[VC$y11.y31] <- ( (-der2p1.dereta1eta1.M1[VC$inde1] + c.copula2.be1.M3*derp1.dereta1.M1[VC$inde1]^2 + c.copula.be1.M3*der2p1.dereta1eta1.M1[VC$inde1])/p11 - (-derp1.dereta1.M1[VC$inde1] +c.copula.be1.M3*derp1.dereta1.M1[VC$inde1])^2/p11^2 )[VC$y11.y31R] 
      d2l.be1.be1             <- -VC$weights*d2l.be1.be1         
      
      d2l.be2.be2[VC$y10.y20R] <- ( (c.copula2.be2.M2*derp1.dereta1.M2^2 + c.copula.be2.M2*der2p1.dereta1eta1.M2)/p00  - (c.copula.be2.M2*derp1.dereta1.M2)^2/p00^2     )[VC$y10.y20R] 
      d2l.be2.be2[VC$y10.y21R] <- ( (-c.copula2.be2.M2*derp1.dereta1.M2^2 + -c.copula.be2.M2*der2p1.dereta1eta1.M2)/p01  - (-c.copula.be2.M2*derp1.dereta1.M2)^2/p01^2  )[VC$y10.y21R]
      d2l.be2.be2              <- -VC$weights[VC$inde0]*d2l.be2.be2     
     
      d2l.be3.be3[VC$y11.y30R] <- (  (der2p1.dereta1eta1.M3 - c.copula2.be2.M3*derp1.dereta1.M3^2 + -c.copula.be2.M3*der2p1.dereta1eta1.M3)/p10  - (derp1.dereta1.M3 - c.copula.be2.M3*derp1.dereta1.M3)^2/p10^2   )[VC$y11.y30R] 
      d2l.be3.be3[VC$y11.y31R] <- ( (-der2p1.dereta1eta1.M3 + c.copula2.be2.M3*derp1.dereta1.M3^2 + c.copula.be2.M3*der2p1.dereta1eta1.M3)/p11  - (-derp1.dereta1.M3 + c.copula.be2.M3*derp1.dereta1.M3)^2/p11^2   )[VC$y11.y31R]
      d2l.be3.be3 <- -VC$weights[VC$inde1]*d2l.be3.be3     

      d2l.th1.th1[VC$y10.y20R] <- ( bit1.th2.M2/p00  - c.copula.theta.M2^2/p00^2      )[VC$y10.y20R]
      d2l.th1.th1[VC$y10.y21R] <- ( -bit1.th2.M2/p01  - (-c.copula.theta.M2)^2/p01^2  )[VC$y10.y21R]
      d2l.th1.th1              <- -VC$weights[VC$inde0]*d2l.th1.th1
     
      d2l.th2.th2[VC$y11.y30R] <- (  -bit1.th2.M3/p10  - (-c.copula.theta.M3)^2/p10^2  )[VC$y11.y30R] 
      d2l.th2.th2[VC$y11.y31R] <- (   bit1.th2.M3/p11  - (c.copula.theta.M3)^2/p11^2   )[VC$y11.y31R]
      d2l.th2.th2              <- -VC$weights[VC$inde1]*d2l.th2.th2

   

      d2l.be1.be2[VC$y10.y20R] <- ( (c.copula2.be1be2.M2*derp1.dereta1.M1[VC$inde0]*derp1.dereta1.M2)/p00  - (c.copula.be1.M2*c.copula.be2.M2*derp1.dereta1.M1[VC$inde0]*derp1.dereta1.M2)/p00^2 )[VC$y10.y20R] 
      d2l.be1.be2[VC$y10.y21R] <- ( (-c.copula2.be1be2.M2*derp1.dereta1.M1[VC$inde0]*derp1.dereta1.M2)/p01  - ((derp1.dereta1.M1[VC$inde0] - c.copula.be1.M2*derp1.dereta1.M1[VC$inde0])*(-c.copula.be2.M2*derp1.dereta1.M2))/p01^2  )[VC$y10.y21R]
      d2l.be1.be2              <- -VC$weights[VC$inde0]*d2l.be1.be2

      d2l.be1.th1[VC$y10.y20R] <- ( (c.copula2.be1th.M2*derp1.dereta1.M1[VC$inde0])/p00  - (c.copula.be1.M2*c.copula.theta.M2*derp1.dereta1.M1[VC$inde0])/p00^2 )[VC$y10.y20R] 
      d2l.be1.th1[VC$y10.y21R] <- ( (-c.copula2.be1th.M2*derp1.dereta1.M1[VC$inde0])/p01  - (derp1.dereta1.M1[VC$inde0] - c.copula.be1.M2*derp1.dereta1.M1[VC$inde0])*-c.copula.theta.M2/p01^2  )[VC$y10.y21R]
      d2l.be1.th1              <- -VC$weights[VC$inde0]*d2l.be1.th1

   
      d2l.be2.th1[VC$y10.y20R] <- ( (c.copula2.be2th.M2*derp1.dereta1.M2)/p00  - (c.copula.be2.M2*c.copula.theta.M2*derp1.dereta1.M2)/p00^2 )[VC$y10.y20R] 
      d2l.be2.th1[VC$y10.y21R] <- ( (-c.copula2.be2th.M2*derp1.dereta1.M2)/p01  - (-c.copula.be2.M2*-c.copula.theta.M2*derp1.dereta1.M2)/p01^2 )[VC$y10.y21R]
      d2l.be2.th1              <- -VC$weights[VC$inde0]*d2l.be2.th1
 
      d2l.be1.be3[VC$y11.y30R] <- ( (-c.copula2.be1be2.M3*derp1.dereta1.M3*derp1.dereta1.M1[VC$inde1])/p10 - ( (-c.copula.be1.M3*derp1.dereta1.M1[VC$inde1])*(derp1.dereta1.M3 - c.copula.be2.M3*derp1.dereta1.M3) )/p10^2 )[VC$y11.y30R] 
      d2l.be1.be3[VC$y11.y31R] <- (  (c.copula2.be1be2.M3*derp1.dereta1.M3*derp1.dereta1.M1[VC$inde1])/p11 - ( (-derp1.dereta1.M1[VC$inde1] + c.copula.be1.M3*derp1.dereta1.M1[VC$inde1])*(-derp1.dereta1.M3 + c.copula.be2.M3*derp1.dereta1.M3) )/p11^2 )[VC$y11.y31R]
      d2l.be1.be3              <- -VC$weights[VC$inde1]*d2l.be1.be3 
 
 

      d2l.be1.th2[VC$y11.y30R] <- ( (-c.copula2.be1th.M3*derp1.dereta1.M1[VC$inde1])/p10 - ( (-c.copula.be1.M3*derp1.dereta1.M1[VC$inde1])*(-c.copula.theta.M3) )/p10^2 )[VC$y11.y30R] 
      d2l.be1.th2[VC$y11.y31R] <- ( (c.copula2.be1th.M3*derp1.dereta1.M1[VC$inde1])/p11 - ( (-derp1.dereta1.M1[VC$inde1]+c.copula.be1.M3*derp1.dereta1.M1[VC$inde1])*(c.copula.theta.M3) )/p11^2 )[VC$y11.y31R]
      d2l.be1.th2              <- -VC$weights[VC$inde1]*d2l.be1.th2

 
       d2l.be3.th2[VC$y11.y30R] <- ( (-c.copula2.be2th.M3*derp1.dereta1.M3)/p10 - ( (derp1.dereta1.M3-c.copula.be2.M3*derp1.dereta1.M3)*(-c.copula.theta.M3) )/p10^2 )[VC$y11.y30R] 
       d2l.be3.th2[VC$y11.y31R] <- ( (c.copula2.be2th.M3*derp1.dereta1.M3)/p11 - ( (-derp1.dereta1.M3+c.copula.be2.M3*derp1.dereta1.M3)*(c.copula.theta.M3) )/p11^2 )[VC$y11.y31R]
       d2l.be3.th2              <- -VC$weights[VC$inde1]*d2l.be3.th2
  
 
 

    ######################################################################################################## 
     
      #d2l.th1.th2 <- 0
      #d2l.be2.be3 <- 0 
      #d2l.be2.th2 <- 0
      #d2l.be3.th1 <- 0
      
  be1.be1 <- crossprod(VC$X1*c(d2l.be1.be1),           VC$X1)
  be1.be2 <- crossprod(VC$X1[VC$inde0,]*c(d2l.be1.be2),VC$X2)
  be1.be3 <- crossprod(VC$X1[VC$inde1,]*c(d2l.be1.be3),VC$X3)
  be1.th1 <- crossprod(VC$X1[VC$inde0,]*c(d2l.be1.th1),VC$X4)
  be1.th2 <- crossprod(VC$X1[VC$inde1,]*c(d2l.be1.th2),VC$X5)  
  
  
  be2.be2 <-     crossprod(VC$X2*c(d2l.be2.be2),VC$X2)
  be2.be3 <- matrix(0, dim(VC$X2)[2],       dim(VC$X3)[2])
  be2.th1 <-     crossprod(VC$X2*c(d2l.be2.th1),VC$X4)
  be2.th2 <- matrix(0, dim(VC$X2)[2],       dim(VC$X5)[2])  
  
  
  be3.be3 <-     crossprod(VC$X3*c(d2l.be3.be3),VC$X3)
  be3.th1 <- matrix(0, dim(VC$X3)[2],       dim(VC$X4)[2])
  be3.th2 <-     crossprod(VC$X3*c(d2l.be3.th2),VC$X5)
  
  th1.th1 <-     crossprod(VC$X4*c(d2l.th1.th1),VC$X4)
  th1.th2 <- matrix(0, dim(VC$X4)[2],       dim(VC$X5)[2])
  
  th2.th2 <-     crossprod(VC$X5*c(d2l.th2.th2),VC$X5)

  

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
              p1 = p1, p2 = p2, p3 = p3, c.copula2.be1be2 = c(c.copula2.be1be2.M2, c.copula2.be1be2.M3),          
              teta.st1 = teta.st1, teta.st2 = teta.st2,
              Cop1 = Cop1, Cop2 = Cop2, teta1 = teta1, teta2 = teta2)    

}




     























