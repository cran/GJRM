triprobgHsESS <- function (params, respvec, VC, ps, AT = FALSE){
  
  eta1 <- VC$X1%*%params[1:VC$X1.d2]
  eta2 <- VC$X2%*%params[(VC$X1.d2 + 1):(VC$X1.d2 + VC$X2.d2)]
  eta3 <- VC$X3%*%params[(VC$X1.d2 + VC$X2.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2)]
  
  etad <- p111 <- A <- NULL
  
  theta12.st <- params[(VC$X1.d2 + VC$X2.d2 + VC$X3.d2+1)]    
  theta13.st <- params[(VC$X1.d2 + VC$X2.d2 + VC$X3.d2+2)]    
  theta23.st <- params[(VC$X1.d2 + VC$X2.d2 + VC$X3.d2+3)]    
  
  p1 <- probm(eta1, VC$margins[1], min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)$pr
  p2 <- probm(eta2, VC$margins[2], min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)$pr
  p3 <- probm(eta3, VC$margins[3], min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)$pr 
  
  mar1 <- qnorm(p1)
  mar2 <- qnorm(p2)
  mar3 <- qnorm(p3)
  
  ###########################
  
  if(VC$Chol == FALSE){
    
    theta12    <- tanh(theta12.st)
    theta13    <- tanh(theta13.st)
    theta23    <- tanh(theta23.st)
    
    Sigma <-  matrix( c( 1,        theta12, theta13,
                         theta12,        1, theta23,
                         theta13,  theta23,        1), 3 , 3) 
    
   Sigma <- PosDefCor(Sigma, Chol = FALSE, theta12.st, theta13.st, theta23.st)
    
    theta12 <- Sigma[1,2]
    theta13 <- Sigma[1,3]
    theta23 <- Sigma[2,3]   
    
    theta12.st <- atanh(theta12)
    theta13.st <- atanh(theta13)
    theta23.st <- atanh(theta23)
    
  }
  
  if(VC$Chol == TRUE){
    
      Sigma <- PosDefCor(Sigma = 1, Chol = TRUE, theta12.st, theta13.st, theta23.st)
    
    theta12 <- Sigma[1,2]
    theta13 <- Sigma[1,3]
    theta23 <- Sigma[2,3]   
    
  }
  
  
  theta12 <- mmf(theta12, max.pr = VC$max.pr)
  theta13 <- mmf(theta13, max.pr = VC$max.pr)
  theta23 <- mmf(theta23, max.pr = VC$max.pr)
  
  p11 <- mm(pbinorm(mar1[VC$inde], mar2, cov12 = theta12), min.pr = VC$min.pr, max.pr = VC$max.pr)
  p13 <- mm(pbinorm(mar1[VC$inde], mar3, cov12 = theta13), min.pr = VC$min.pr, max.pr = VC$max.pr)
  p23 <- mm(pbinorm(mar2,          mar3, cov12 = theta23), min.pr = VC$min.pr, max.pr = VC$max.pr)
  
  #params <- c(params[1:(VC$X1.d2 + VC$X2.d2 + VC$X3.d2)],theta12.st,theta13.st,theta23.st)
  
  
  
  if(VC$approx == FALSE){ 
  
  eta123 <- cbind(mar1[VC$inde],mar2,mar3)
  
  #for(i in 1:length(mar3)) p111[i] <- mm( pmnorm(x = eta123[i, ], varcov = Sigma)[1] , min.pr = VC$min.pr, max.pr = VC$max.pr)
  
                            p111    <- mm( pmnorm(x = eta123,      varcov = Sigma),     min.pr = VC$min.pr, max.pr = VC$max.pr)

  }
  
  if(VC$approx == TRUE) p111 <- mm(  TRIapprox(mar1[VC$inde], mar2, mar3, Sigma) , min.pr = VC$min.pr, max.pr = VC$max.pr)
  
  ############################################################################################
  
  
  p110 <- mm(p11 - p111, min.pr = VC$min.pr, max.pr = VC$max.pr)
  p101 <- mm(p13 - p111, min.pr = VC$min.pr, max.pr = VC$max.pr)
  p100 <- mm(p1[VC$inde] - p11 - p101, min.pr = VC$min.pr, max.pr = VC$max.pr)
  
  p0   <- mm(1 - p1, min.pr = VC$min.pr, max.pr = VC$max.pr) 
  
  ##########################################################################################
  l.par1            <- respvec$cy1*log(p0) 
  l.par1[VC$inde]   <- respvec$y1.y2.y3  * log(p111) + respvec$y1.cy2.y3  * log(p101) + respvec$y1.y2.cy3 * log(p110) + respvec$y1.cy2.cy3 * log(p100)
  l.par              <- VC$weights * l.par1 
  
  res <- -sum(l.par) 
  ##########################################################################################
  
  
  TIn <- list(eta1 = eta1, eta2 = eta2, eta3 = eta3, 
              theta12 = theta12, theta13 = theta13, theta23 = theta23, 
              theta12.st = theta12.st, theta13.st = theta13.st, theta23.st = theta23.st, 
              mar1 = mar1, mar2 = mar2, mar3 = mar3,
              p111 = p111, p110 = p110, p100 = p100, p101 = p101, p0 = p0)
  
  gTRI <- g.triESS(respvec = respvec, VC = VC, TIn = TIn)
  
  G <- -c( colSums(c(gTRI$dl.de1) * VC$X1), 
           colSums(c(gTRI$dl.de2) * VC$X2), 
           colSums(c(gTRI$dl.de3) * VC$X3), 
           sum(gTRI$dl.dtheta12.st), 
           sum(gTRI$dl.dtheta13.st), 
           sum(gTRI$dl.dtheta23.st) )
  
  
  LgTRI <- list(p12.g = gTRI$p12.g, p13.g = gTRI$p13.g, p23.g = gTRI$p23.g, 
                p12.g.c = gTRI$p12.g.c, p13.g.c = gTRI$p13.g.c, p23.g.c = gTRI$p23.g.c, 
                d.1 = gTRI$d.1, d.2 = gTRI$d.2, d.3 = gTRI$d.3,
                dmar1 = gTRI$dmar1, dmar2 = gTRI$dmar2, dmar3 = gTRI$dmar3,
                d11.12 = gTRI$d11.12, d11.13 = gTRI$d11.13, d11.23 = gTRI$d11.23,
                p.1.11 = gTRI$p.1.11, p.1.10 = gTRI$p.1.10, p.1.00 = gTRI$p.1.00, 
                p.1.01 = gTRI$p.1.01, p.2.11 = gTRI$p.2.11, p.2.10 = gTRI$p.2.10, 
                p.2.00 = gTRI$p.2.00, p.2.01 = gTRI$p.2.01, p.3.11 = gTRI$p.3.11, 
                p.3.10 = gTRI$p.3.10, p.3.00 = gTRI$p.3.00, p.3.01 = gTRI$p.3.01,
                d11.12 = gTRI$d11.12,
                mean.12 = gTRI$mean.12,
                mean.13 = gTRI$mean.13,
                mean.23 = gTRI$mean.23, sd.12 = gTRI$sd.12,
                sd.13 = gTRI$sd.13,
                sd.23 = gTRI$sd.23,
                upst.1 = gTRI$upst.1,
                upst.2 = gTRI$upst.2,
                dF1.de1 = gTRI$dF1.de1, dF2.de2 = gTRI$dF2.de2, dF3.de3 = gTRI$dF3.de3,
                dl.dF1 = gTRI$dl.dF1, dl.dF2 = gTRI$dl.dF2, dl.dF3 = gTRI$dl.dF3,
                dl.dtheta12 = gTRI$dl.dtheta12, dl.dtheta13 = gTRI$dl.dtheta13, dl.dtheta23 = gTRI$dl.dtheta23) 
  
  HTRI <- H.triESS(respvec = respvec, VC = VC, TIn = TIn, LgTRI = LgTRI)
  H    <- - rbind(HTRI$h1, HTRI$h2, HTRI$h3, HTRI$h4, HTRI$h5, HTRI$h6)
  
  ##########################################################################################
  ##########################################################################################
  
  S.h  <- ps$S.h
  
  
  if( VC$penCor %in% c("lasso", "alasso") ){
    
    if(VC$penCor %in% c("lasso")) A <- diag(1/(sqrt(params[(length(params)-2):length(params)]^2 + 1e-08)))  
    
    if(VC$penCor %in% c("alasso")){
      
      wc <- 1/abs(VC$wc)^VC$gamma
      A  <- diag(wc * 1/(sqrt(params[(length(params)-2):length(params)]^2 + 1e-08)))
      
    }
    
    if( VC$l.sp1==0 && VC$l.sp2==0 && VC$l.sp3==0) VC$qu.mag$Ss[[1]]                      <- A 
    if( VC$l.sp1!=0 || VC$l.sp2!=0 || VC$l.sp3!=0) VC$qu.mag$Ss[[length(VC$qu.mag$Ss)+1]] <- A
    
    S.h <- adiag( S.h, VC$sp[length(VC$sp)]*VC$qu.mag$Ss[[length(VC$qu.mag$Ss)]])
    
  }
  
  
  
  if (VC$extra.regI == "pC") H <- regH(H, type = 1)
  
  
  if( length(S.h) != 1){
    
    S.h1 <- 0.5*crossprod(params,S.h)%*%params
    S.h2 <- S.h%*%params
    
  } else S.h <- S.h1 <- S.h2 <- 0   
  
  
  S.res <- res
  res   <- S.res + S.h1
  G     <- G + S.h2
  H     <- H + S.h  
  
  
  if (VC$extra.regI == "sED") H <- regH(H, type = 2)
  
  dl.de2 <- dl.de3 <- dl.dtheta12.st <- dl.dtheta13.st <- dl.dtheta23.st <- rep(0, length(eta1))
  
  dl.de1 <- gTRI$dl.de1
  
  dl.de2[VC$inde] <- gTRI$dl.de2
  dl.de3[VC$inde] <- gTRI$dl.de3
  
  #dl.dtheta12.st[VC$inde2] <- gTRI$dl.dtheta12.st 
  #dl.dtheta13.st[VC$inde2] <- gTRI$dl.dtheta13.st 
  #dl.dtheta23.st[VC$inde2] <- gTRI$dl.dtheta23.st 
  
  dl.dtheta12.st[VC$inde] <- gTRI$dl.dtheta12.st 
  dl.dtheta13.st[VC$inde] <- gTRI$dl.dtheta13.st 
  dl.dtheta23.st[VC$inde] <- gTRI$dl.dtheta23.st 
  
  list(value = res, gradient = G, hessian = H, S.h=S.h, S.h1=S.h1, S.h2=S.h2, qu.mag = VC$qu.mag,
       l = S.res, l.par = l.par, ps = ps, 
       eta1 = eta1, eta2 = eta2, eta3 = eta3, 
       p111 = p111,
       p110 = p110,
       p10  = NULL,
       p0   = p0,
       p011 = NULL,
       p101 = NULL,
       p100 = NULL,
       p010 = NULL,
       p001 = NULL,
       p000 = NULL,
       theta12 = theta12,
       theta13 = theta13,
       theta23 = theta23,
       dl.de1 = dl.de1,
       dl.de2 = dl.de2,
       dl.de3 = dl.de3,
       dl.dtheta12.st = dl.dtheta12.st, 
       dl.dtheta13.st = dl.dtheta13.st, 
       dl.dtheta23.st = dl.dtheta23.st,
       p1 = p1, p2 = p2, p3 = p3)
  
}
