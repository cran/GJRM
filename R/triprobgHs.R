triprobgHs <- function(params, respvec, VC, ps, AT = FALSE){
  
  eta1 <- VC$X1%*%params[1:VC$X1.d2]
  eta2 <- VC$X2%*%params[(VC$X1.d2 + 1):(VC$X1.d2 + VC$X2.d2)]
  eta3 <- VC$X3%*%params[(VC$X1.d2 + VC$X2.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2)]
  
  theta12 <- theta13 <- theta23 <- NA
  
  
  etad <- p111 <- A <- NULL
  
  if(is.null(VC$X4)){
    
    theta12.st <- params[(VC$X1.d2 + VC$X2.d2 + VC$X3.d2+1)]    
    theta13.st <- params[(VC$X1.d2 + VC$X2.d2 + VC$X3.d2+2)]    
    theta23.st <- params[(VC$X1.d2 + VC$X2.d2 + VC$X3.d2+3)]  
    
  }
  
  if(!is.null(VC$X4)){
    
    theta12.st <- VC$X4%*%params[(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2)]    
    theta13.st <- VC$X5%*%params[(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2 + VC$X5.d2)]    
    theta23.st <- VC$X6%*%params[(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2 + VC$X5.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2 + VC$X5.d2 + VC$X6.d2)]  
    
  }  
  
  
  #####
  # ! #
  #####################################################################
  ## I replaced "probit" with  margins[1], margins[2] and margins[3] ##                   
  #####################################################################
  
  p1 <- probm(eta1, VC$margins[1], min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)$pr
  p2 <- probm(eta2, VC$margins[2], min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)$pr
  p3 <- probm(eta3, VC$margins[3], min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)$pr 
  
  #####################################################################
  
  #####
  # ! #
  ###########################
  ## The following is new: ##                   
  ###########################
  
  mar1 <- qnorm(p1)
  mar2 <- qnorm(p2)
  mar3 <- qnorm(p3)
  
  ###########################
  
  if(VC$Chol == FALSE){
    
    theta12 <- tanh(theta12.st)
    theta13 <- tanh(theta13.st)
    theta23 <- tanh(theta23.st)
    
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
    
    if(is.null(VC$X4)){
  
      Sigma <- PosDefCor(Sigma = 1, Chol = TRUE, theta12.st, theta13.st, theta23.st)
 
      theta12 <- Sigma[1,2]
      theta13 <- Sigma[1,3]
      theta23 <- Sigma[2,3]   
    } 
    
    
    if(!is.null(VC$X4)){
      
      Sigma <- list()
      
      for(i in 1:VC$n){
        
        Sigma[[i]] <- PosDefCor(Sigma = NULL, Chol = TRUE, theta12.st = theta12.st[i], theta13.st = theta13.st[i], theta23.st = theta23.st[i]) 
       
        theta12[i] <- Sigma[[i]][1,2]
        theta13[i] <- Sigma[[i]][1,3]
        theta23[i] <- Sigma[[i]][2,3]   
        
      }
      
    }   
    
    
  }
  
  
  theta12 <- mmf(theta12, max.pr = VC$max.pr)
  theta13 <- mmf(theta13, max.pr = VC$max.pr)
  theta23 <- mmf(theta23, max.pr = VC$max.pr)
  
  p11 <- mm( pbinorm( mar1, mar2, cov12 = theta12) , min.pr = VC$min.pr, max.pr = VC$max.pr)
  p13 <- mm( pbinorm( mar1, mar3, cov12 = theta13) , min.pr = VC$min.pr, max.pr = VC$max.pr)
  p23 <- mm( pbinorm( mar2, mar3, cov12 = theta23) , min.pr = VC$min.pr, max.pr = VC$max.pr)
  
  #if(VC$Chol == FALSE) params <- c(params[1:(VC$X1.d2 + VC$X2.d2 + VC$X3.d2)],theta12.st,theta13.st,theta23.st)
  
  if(is.null(VC$X4)){
  
    #if(VC$approx == FALSE){ for(i in 1:VC$n) p111[i] <- mm( pmnorm(x = c(mar1[i], mar2[i], mar3[i]), varcov = Sigma)[1] , min.pr = VC$min.pr, max.pr = VC$max.pr)  }
    
    if(VC$approx == FALSE){                     p111 <- mm( pmnorm(x = cbind(mar1, mar2, mar3),      varcov = Sigma),     min.pr = VC$min.pr, max.pr = VC$max.pr)  }

     
    if(VC$approx == TRUE) p111 <- mm( TRIapprox(mar1, mar2, mar3, Sigma), min.pr = VC$min.pr, max.pr = VC$max.pr )
  }
  
  if(!is.null(VC$X4)){
    for(i in 1:VC$n) p111[i] <- mm( pmnorm(x = c(mar1[i], mar2[i], mar3[i]), varcov = Sigma[[i]])[1] , min.pr = VC$min.pr, max.pr = VC$max.pr) 
  }  
  
  
  ############################################################################################
  
  p011 <- mm(p23 - p111, min.pr = VC$min.pr, max.pr = VC$max.pr)
  p101 <- mm(p13 - p111, min.pr = VC$min.pr, max.pr = VC$max.pr)
  p110 <- mm(p11 - p111, min.pr = VC$min.pr, max.pr = VC$max.pr)
  p100 <- mm(p1 - p11 - p101, min.pr = VC$min.pr, max.pr = VC$max.pr)
  p010 <- mm(p2 - p11 - p011, min.pr = VC$min.pr, max.pr = VC$max.pr)
  p001 <- mm(p3 - p23 - p101, min.pr = VC$min.pr, max.pr = VC$max.pr)
  p000 <- mm(1 - p111 - p011 - p101 - p110 - p001 - p010 - p100, min.pr = VC$min.pr, max.pr = VC$max.pr)
  
  sall.p <- rowSums(cbind(p111, p011, p101, p110, p100, p010, p001, p000))
  
  p111 <- p111/sall.p 
  p011 <- p011/sall.p 
  p101 <- p101/sall.p 
  p110 <- p110/sall.p 
  p100 <- p100/sall.p 
  p010 <- p010/sall.p 
  p001 <- p001/sall.p 
  p000 <- p000/sall.p 

  ##########################################################################################
  
  l.par <- VC$weights * (respvec$y1.y2.y3    * log(p111) + respvec$y1.y2.cy3  * log(p110) +
                           respvec$cy1.y2.y3   * log(p011) + respvec$cy1.y2.cy3 * log(p010) +
                           respvec$cy1.cy2.cy3 * log(p000) + respvec$cy1.cy2.y3 * log(p001) + 
                           respvec$y1.cy2.cy3  * log(p100) + respvec$y1.cy2.y3  * log(p101)   ) 
  
  res <- -sum(l.par) 
  
  ##########################################################################################
  
  
  #####
  # ! #
  #######################################################
  ## In TIn I added the 4th line because we need these ##
  ## quantities for the derivatives                    ##
  #######################################################
  
  TIn <- list(eta1 = eta1, eta2 = eta2, eta3 = eta3, 
              theta12 = theta12, theta13 = theta13, theta23 = theta23, 
              theta12.st = theta12.st, theta13.st = theta13.st, theta23.st = theta23.st, 
              mar1 = mar1, mar2 = mar2, mar3 = mar3,
              p111 = p111, p110 = p110, p011 = p011, p010 = p010, 
              p000 = p000, p001 = p001, p100 = p100, p101 = p101)
  
  #####################################################################
  
  gTRI <- g.tri(respvec = respvec, VC = VC, TIn = TIn)
  
  if(is.null(VC$X4)){  
    
    G <- -c( colSums(c(gTRI$dl.de1) * VC$X1), 
             colSums(c(gTRI$dl.de2) * VC$X2), 
             colSums(c(gTRI$dl.de3) * VC$X3), 
             sum(gTRI$dl.dtheta12.st), 
             sum(gTRI$dl.dtheta13.st), 
             sum(gTRI$dl.dtheta23.st) )
    
  }
  
  
  if(!is.null(VC$X4)){  
    
    G <- -c( colSums(c(gTRI$dl.de1) * VC$X1), 
             colSums(c(gTRI$dl.de2) * VC$X2), 
             colSums(c(gTRI$dl.de3) * VC$X3), 
             colSums(c(gTRI$dl.dtheta12.st) * rbind(VC$X4, VC$X4, VC$X4)), 
             colSums(c(gTRI$dl.dtheta13.st) * rbind(VC$X5, VC$X5, VC$X5)), 
             colSums(c(gTRI$dl.dtheta23.st) * rbind(VC$X6, VC$X6, VC$X6)) )
    
  }  
  
  
  
  #####
  # ! #
  ###############################################################
  ## In LgTRI I added the 3rd, 4th, 16th and 17th line because ## 
  ## we need these quantities for the Hessian                  ##                 
  ###############################################################
  
  LgTRI <- list(p12.g = gTRI$p12.g, p13.g = gTRI$p13.g, p23.g = gTRI$p23.g, 
                p12.g.c = gTRI$p12.g.c, p13.g.c = gTRI$p13.g.c, p23.g.c = gTRI$p23.g.c, 
                d.1 = gTRI$d.1, d.2 = gTRI$d.2, d.3 = gTRI$d.3,
                dmar1 = gTRI$dmar1, dmar2 = gTRI$dmar2, dmar3 = gTRI$dmar3,
                d11.12 = gTRI$d11.12, d11.13 = gTRI$d11.13, d11.23 = gTRI$d11.23,
                p.1.11 = gTRI$p.1.11, p.1.10 = gTRI$p.1.10, p.1.00 = gTRI$p.1.00, 
                p.1.01 = gTRI$p.1.01, p.2.11 = gTRI$p.2.11, p.2.10 = gTRI$p.2.10, 
                p.2.00 = gTRI$p.2.00, p.2.01 = gTRI$p.2.01, p.3.11 = gTRI$p.3.11, 
                p.3.10 = gTRI$p.3.10, p.3.00 = gTRI$p.3.00, p.3.01 = gTRI$p.3.01,
                mean.12 = gTRI$mean.12,
                mean.13 = gTRI$mean.13,
                mean.23 = gTRI$mean.23, 
                sd.12 = gTRI$sd.12,
                sd.13 = gTRI$sd.13,
                sd.23 = gTRI$sd.23,
                dF1.de1 = gTRI$dF1.de1, dF2.de2 = gTRI$dF2.de2, dF3.de3 = gTRI$dF3.de3,
                dl.dF1 = gTRI$dl.dF1, dl.dF2 = gTRI$dl.dF2, dl.dF3 = gTRI$dl.dF3,
                dl.dtheta12 = gTRI$dl.dtheta12, dl.dtheta13 = gTRI$dl.dtheta13, dl.dtheta23 = gTRI$dl.dtheta23) 
  
  #####################################################################
  
  HTRI <- H.tri(respvec = respvec, VC = VC, TIn = TIn, LgTRI = LgTRI)
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
  
  
  list(value = res, gradient = G, hessian = H, S.h=S.h, S.h1=S.h1, S.h2=S.h2, qu.mag = VC$qu.mag,
       l = S.res, l.par = l.par, ps = ps, 
       eta1 = eta1, eta2 = eta2, eta3 = eta3, 
       p111 = p111,
       p011 = p011,
       p101 = p101,
       p110 = p110,
       p100 = p100,
       p010 = p010,
       p001 = p001,
       p000 = p000,
       theta12 = theta12,
       theta13 = theta13,
       theta23 = theta23,
       dl.de1 = gTRI$dl.de1,
       dl.de2 = gTRI$dl.de2,
       dl.de3 = gTRI$dl.de3,
       dl.dtheta12.st = gTRI$dl.dtheta12.st, 
       dl.dtheta13.st = gTRI$dl.dtheta13.st, 
       dl.dtheta23.st = gTRI$dl.dtheta23.st,
       p1 = p1, p2 = p2, p3 = p3)
  
}
