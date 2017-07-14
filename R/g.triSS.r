g.triSS <- function(params, respvec, VC, TIn){
  
  #####
  # ! #
  ########################################################################
  ## I replaced TIn$eta1 with TIn$mar1. Same for TIn$eta2 and TIn$eta3  ##                 
  ########################################################################
  
  mean1 <- TIn$theta12 * TIn$mar1
  mean2 <- TIn$theta13 * TIn$mar1
  mean3 <- TIn$theta12 * TIn$mar2
  mean4 <- TIn$theta23 * TIn$mar2
  mean5 <- TIn$theta13 * TIn$mar3
  mean6 <- TIn$theta23 * TIn$mar3
  
  ################################################################
  
  var1 <- 1 - TIn$theta12^2   
  var2 <- 1 - TIn$theta13^2
  var3 <- 1 - TIn$theta23^2
  
  cov1 <- TIn$theta23 - TIn$theta12 * TIn$theta13
  cov2 <- TIn$theta13 - TIn$theta12 * TIn$theta23
  cov3 <- TIn$theta12 - TIn$theta13 * TIn$theta23
  
  cov1 <- mmf(cov1)
  cov2 <- mmf(cov2)
  cov3 <- mmf(cov3)  
  
  #####
  # ! #
  ########################################################################
  ## I replaced TIn$eta1 with TIn$mar1. Same for TIn$eta2 and TIn$eta3  ##                 
  ########################################################################
  
  d.1 <- dnorm(TIn$mar1)
  d.2 <- dnorm(TIn$mar2)
  d.3 <- dnorm(TIn$mar3)
  
  p.1.11  <- mm(pbinorm( TIn$mar2[VC$inde2.1],  TIn$mar3,             mean1 =   mean1[VC$inde2],   mean2 =   mean2[VC$inde2], var1 = var1, var2 = var2, cov12 =   cov1) )            
  p.1.10  <- mm(pbinorm( TIn$mar2[VC$inde2.1], -TIn$mar3,             mean1 =   mean1[VC$inde2],   mean2 =  -mean2[VC$inde2], var1 = var1, var2 = var2, cov12 =  -cov1) ) 
  
  p.2.11  <- mm(pbinorm( TIn$mar1[VC$inde2],    TIn$mar3,             mean1 =   mean3[VC$inde2.1], mean2 =   mean4[VC$inde2.1], var1 = var1, var2 = var3, cov12 =  cov2) )
  p.2.10  <- mm(pbinorm( TIn$mar1[VC$inde2],   -TIn$mar3,             mean1 =   mean3[VC$inde2.1], mean2 =  -mean4[VC$inde2.1], var1 = var1, var2 = var3, cov12 = -cov2) )
  
  p.3.11  <- mm(pbinorm( TIn$mar1[VC$inde2],    TIn$mar2[VC$inde2.1], mean1 =   mean5,             mean2 =   mean6, var1 = var2, var2 = var3, cov12 =   cov3) )
  p.3.10  <- mm(pbinorm( TIn$mar1[VC$inde2],   -TIn$mar2[VC$inde2.1], mean1 =   mean5,             mean2 =  -mean6, var1 = var2, var2 = var3, cov12 =  -cov3) ) 
  
  upst.1   <- mm( pnorm( (TIn$mar2            - TIn$theta12 * TIn$mar1[VC$inde1])/sqrt(1 - TIn$theta12^2)) )
  upst.2   <- mm( pnorm( (TIn$mar1[VC$inde1]  - TIn$theta12 * TIn$mar2          )/sqrt(1 - TIn$theta12^2)) )
  
  ################################################################################################
  
  #####
  # ! #
  ###############################
  ## The next 6 lines are new  ##                 
  ###############################
  
  dmar1 <- probm(TIn$eta1, VC$margins[1], only.pr = FALSE)$d.n
  dmar2 <- probm(TIn$eta2, VC$margins[2], only.pr = FALSE)$d.n
  dmar3 <- probm(TIn$eta3, VC$margins[3], only.pr = FALSE)$d.n
  
  dF1.de1 <- (1/d.1) * dmar1
  dF2.de2 <- (1/d.2) * dmar2
  dF3.de3 <- (1/d.3) * dmar3
  
  ###################################################################
  
  dl.dF1.1           <- - respvec$cy1/TIn$p0 * d.1
  dl.dF1.1[VC$inde1] <-   respvec$y1.cy2/TIn$p10 * (d.1[VC$inde1] - d.1[VC$inde1] * upst.1) 
  dl.dF1.1[VC$inde2] <-   respvec$y1.y2.cy3/TIn$p110 * d.1[VC$inde2] * p.1.10 + 
    respvec$y1.y2.y3/TIn$p111  * d.1[VC$inde2] * p.1.11 
  dl.dF1 <- dl.dF1.1
  
  #####
  # ! #
  ###########################
  ## The following is new: ##                 
  ###########################
  
  dl.de1 <- dl.dF1 * dF1.de1
  
  ################################
  
  dl.dF2.1             <- - respvec$y1.cy2/TIn$p10 * d.2 * upst.2
  dl.dF2.1[VC$inde2.1] <-   respvec$y1.y2.cy3/TIn$p110   * d.2[VC$inde2.1] * p.2.10 +
    respvec$y1.y2.y3/TIn$p111    * d.2[VC$inde2.1] * p.2.11 
  dl.dF2 <- dl.dF2.1
  
  #####
  # ! #
  ###########################
  ## The following is new: ##                 
  ###########################
  
  dl.de2 <- dl.dF2 * dF2.de2
  
  ################################
  
  dl.dF3 <- - respvec$y1.y2.cy3/TIn$p110 * d.3 * p.3.11 +
    respvec$y1.y2.y3/TIn$p111 * d.3 * p.3.11 
  
  #####
  # ! #
  ###########################
  ## The following is new: ##                 
  ###########################
  
  dl.de3 <- dl.dF3 * dF3.de3
  
  ################################
  
  #####
  # ! #
  ########################################################################
  ## I replaced TIn$eta1 with TIn$mar1. Same for TIn$eta2 and TIn$eta3  ##                 
  ########################################################################
  
  mean.12 <- ( TIn$mar1[VC$inde1]   * (TIn$theta13 - TIn$theta12 * TIn$theta23) + TIn$mar2 * (TIn$theta23 - TIn$theta12 * TIn$theta13) )/( 1 - TIn$theta12^2 )
  mean.13 <- ( TIn$mar1[VC$inde2]   * (TIn$theta12 - TIn$theta13 * TIn$theta23) + TIn$mar3 * (TIn$theta23 - TIn$theta12 * TIn$theta13) )/( 1 - TIn$theta13^2 )
  mean.23 <- ( TIn$mar2[VC$inde2.1] * (TIn$theta12 - TIn$theta13 * TIn$theta23) + TIn$mar3 * (TIn$theta13 - TIn$theta12 * TIn$theta23) )/( 1 - TIn$theta23^2 )
  
  ############################################################################################################
  
  deno <- 1 - TIn$theta12^2 - TIn$theta13^2 - TIn$theta23^2 + 2 * TIn$theta12 * TIn$theta13 * TIn$theta23
  sd.12   <- sqrt( deno / ( 1 - TIn$theta12^2 ) )
  sd.13   <- sqrt( deno / ( 1 - TIn$theta13^2 ) )
  sd.23   <- sqrt( deno / ( 1 - TIn$theta23^2 ) )
  
  #####
  # ! #
  ########################################################################
  ## I replaced TIn$eta1 with TIn$mar1. Same for TIn$eta2 and TIn$eta3  ##                 
  ########################################################################
  
  p12.g <- mm( pnorm( (TIn$mar3             - mean.12[VC$inde2.1])/sd.12) )
  p13.g <- mm( pnorm( (TIn$mar2[VC$inde2.1] - mean.13            )/sd.13) )
  p23.g <- mm( pnorm( (TIn$mar1[VC$inde2]   - mean.23            )/sd.23) )
  
  ########################################################################
  
  p12.g.c <- mm(1 - p12.g)
  p13.g.c <- mm(1 - p13.g)
  p23.g.c <- mm(1 - p23.g)
  
  #####
  # ! #
  ########################################################################
  ## I replaced TIn$eta1 with TIn$mar1. Same for TIn$eta2 and TIn$eta3  ##                 
  ########################################################################
  
  d11.12 <- dbinorm( TIn$mar1[VC$inde1]  , TIn$mar2, cov12 =  TIn$theta12)
  d11.13 <- dbinorm( TIn$mar1[VC$inde2]  , TIn$mar3, cov12 =  TIn$theta13)
  d11.23 <- dbinorm( TIn$mar2[VC$inde2.1], TIn$mar3, cov12 =  TIn$theta23)
  
  ########################################################################
  
  dl.dtheta12.1             <- - respvec$y1.cy2/TIn$p10   * d11.12
  dl.dtheta12.1[VC$inde2.1] <-   respvec$y1.y2.cy3/TIn$p110 * d11.12[VC$inde2.1] * p12.g.c +
    respvec$y1.y2.y3/TIn$p111   * d11.12[VC$inde2.1] * p12.g 
  
  dl.dtheta12 <- dl.dtheta12.1
  
  
  dl.dtheta13 <- - respvec$y1.y2.cy3/TIn$p110 * d11.13 * p13.g +
    respvec$y1.y2.y3/TIn$p111 * d11.13 * p13.g
  
  dl.dtheta23 <- - respvec$y1.y2.cy3/TIn$p110 * d11.23 * p23.g +
    respvec$y1.y2.y3/TIn$p111 * d11.23 * p23.g
  
  if(VC$Chol == FALSE){
    dtheta12.dtheta12.st <- 4 * exp( 2 * TIn$theta12.st )/( exp(2 * TIn$theta12.st) + 1 )^2
    dtheta13.dtheta13.st <- 4 * exp( 2 * TIn$theta13.st )/( exp(2 * TIn$theta13.st) + 1 )^2
    dtheta23.dtheta23.st <- 4 * exp( 2 * TIn$theta23.st )/( exp(2 * TIn$theta23.st) + 1 )^2
    
    dl.dtheta12.st <- dl.dtheta12 * dtheta12.dtheta12.st
    dl.dtheta13.st <- dl.dtheta13 * dtheta13.dtheta13.st
    dl.dtheta23.st <- dl.dtheta23 * dtheta23.dtheta23.st
  }
  
  if(VC$Chol == TRUE){
    
    dl.dtheta <- matrix(0,length(VC$inde1),3)
    dl.dtheta[VC$inde1, 1] <- dl.dtheta12
    dl.dtheta[VC$inde2, 2]   <- dl.dtheta13
    dl.dtheta[VC$inde2, 3]   <- dl.dtheta23
    
    dth12.dth12.st <- 1/(1 + TIn$theta12.st^2)^(3/2) 
    dth12.dth13.st <- 0
    dth12.dth23.st <- 0
    
    dth13.dth12.st <- 0
    dth13.dth13.st <- (1 + TIn$theta23.st^2)/(1 + TIn$theta13.st^2 + TIn$theta23.st^2)^(3/2)
    dth13.dth23.st <- - (TIn$theta13.st * TIn$theta23.st)/(1 + TIn$theta13.st^2 + TIn$theta23.st^2)^(3/2)
    
    
    dth23.dth12.st <- TIn$theta13.st/sqrt((1 + TIn$theta12.st^2) * (1 + TIn$theta13.st^2 + TIn$theta23.st^2)) - (TIn$theta12.st * (TIn$theta12.st * TIn$theta13.st + TIn$theta23.st))/((1 + TIn$theta12.st^2)^(3/2) * sqrt(1 + TIn$theta13.st^2 + TIn$theta23.st^2))
    dth23.dth13.st <- TIn$theta12.st/sqrt((1 + TIn$theta12.st^2) * (1 + TIn$theta13.st^2 + TIn$theta23.st^2)) - (TIn$theta13.st * (TIn$theta12.st * TIn$theta13.st + TIn$theta23.st))/(sqrt(1 + TIn$theta12.st^2) * (1 + TIn$theta13.st^2 + TIn$theta23.st^2)^(3/2))
    dth23.dth23.st <- 1/sqrt((1 + TIn$theta12.st^2) * (1 + TIn$theta13.st^2 + TIn$theta23.st^2)) - (TIn$theta23.st * (TIn$theta12.st * TIn$theta13.st + TIn$theta23.st))/(sqrt(1 + TIn$theta12.st^2) * (1 + TIn$theta13.st^2 + TIn$theta23.st^2)^(3/2))
    
    dtheta.theta.st <- matrix( c(  dth12.dth12.st,  dth13.dth12.st, dth23.dth12.st,
                                   dth12.dth13.st,  dth13.dth13.st, dth23.dth13.st,
                                   dth12.dth23.st,  dth13.dth23.st, dth23.dth23.st ), 3 , 3)
    
    
    
    dl.dtheta.st   <- dl.dtheta %*% dtheta.theta.st
    
    dl.dtheta12.st <- dl.dtheta.st[, 1]
    dl.dtheta12.st <- dl.dtheta12.st[VC$inde1]
    dl.dtheta13.st <- dl.dtheta.st[, 2]
    dl.dtheta13.st <- dl.dtheta13.st[VC$inde2]
    dl.dtheta23.st <- dl.dtheta.st[, 3]
    dl.dtheta23.st <- dl.dtheta23.st[VC$inde2]
    
    
  }
  
  
  #####
  # ! #
  #######################################################
  ## In GTRIVec: 3rd, 4th, 9th and 10th lines are new  ##
  #######################################################
  
  GTRIVec <- list(p12.g = p12.g, p13.g = p13.g, p23.g = p23.g, 
                  p12.g.c = p12.g.c, p13.g.c = p13.g.c, p23.g.c = p23.g.c, 
                  d.1 = d.1, d.2 = d.2, d.3 = d.3,
                  dmar1 = dmar1, dmar2 = dmar2, dmar3 = dmar3,
                  d11.12 = d11.12, d11.13 = d11.13, d11.23 = d11.23, 
                  p.1.11 = p.1.11, p.1.10 = p.1.10,  
                  p.2.11 = p.2.11, p.2.10 = p.2.10,  
                  p.3.11 = p.3.11, p.3.10 = p.3.10, 
                  dF1.de1 = dF1.de1, dF2.de2 = dF2.de2, dF3.de3 = dF3.de3,
                  dl.dF1 = dl.dF1, dl.dF2 = dl.dF2, dl.dF3 = dl.dF3, 
                  dl.de1 = VC$weights*dl.de1, dl.de2 = VC$weights[VC$inde1]*dl.de2, dl.de3 = VC$weights[VC$inde2]*dl.de3, 
                  dl.dtheta12.st = VC$weights[VC$inde1]*dl.dtheta12.st, dl.dtheta13.st = VC$weights[VC$inde2]*dl.dtheta13.st,
                  dl.dtheta23.st = VC$weights[VC$inde2]*dl.dtheta23.st,
                  mean.12 = mean.12,
                  mean.13 = mean.13,
                  mean.23 = mean.23, 
                  sd.12 = sd.12,
                  sd.13 = sd.13,
                  sd.23 = sd.23,
                  upst.1 = upst.1,
                  upst.2 = upst.2,
                  dl.dtheta12 =dl.dtheta12, dl.dtheta13 = dl.dtheta13, dl.dtheta23 = dl.dtheta23)
  
  GTRIVec
  
}   
