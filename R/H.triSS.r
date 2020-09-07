H.triSS <- function(params, respvec, VC, TIn, LgTRI){
  
  #####
  # ! #
  ########################################################################
  ## I replaced TIn$eta1 with TIn$mar1. Same for TIn$eta2 and TIn$eta3  ##                 
  ########################################################################
  
  dst.1 <- dnorm( (TIn$mar2  - TIn$theta12 * TIn$mar1[VC$inde1] )/sqrt(1 - TIn$theta12^2) )  
  pst.1 <- pnorm( ( ((TIn$mar3 - TIn$theta13 * TIn$mar1[VC$inde2])/sqrt(1 - TIn$theta13^2)) - ((TIn$theta23 - TIn$theta12 * TIn$theta13)/sqrt((1 - TIn$theta12^2) * (1 - TIn$theta13^2))) * ((TIn$mar2[VC$inde2.1]  - TIn$theta12 * TIn$mar1[VC$inde2])/sqrt(1 - TIn$theta12^2)) )/sqrt(1 - ((TIn$theta23 - TIn$theta12 * TIn$theta13)/sqrt((1 - TIn$theta12^2) * (1 - TIn$theta13^2)))^2))
  pst.1 <- mm(pst.1, min.pr = VC$min.pr, max.pr = VC$max.pr)
  st.1 <- -TIn$theta12/sqrt(1 - TIn$theta12^2)
  
  dst.2 <- dnorm((TIn$mar3 - TIn$theta13 * TIn$mar1[VC$inde2])/sqrt(1 - TIn$theta13^2))
  pst.2 <- pnorm( ( ((TIn$mar2[VC$inde2.1]  - TIn$theta12 * TIn$mar1[VC$inde2])/sqrt(1 - TIn$theta12^2)) - ((TIn$theta23 - TIn$theta12 * TIn$theta13)/sqrt((1 - TIn$theta12^2) * (1 - TIn$theta13^2))) * ((TIn$mar3 - TIn$theta13 * TIn$mar1[VC$inde2])/sqrt(1 - TIn$theta13^2)) )/sqrt(1 - ((TIn$theta23 - TIn$theta12 * TIn$theta13)/sqrt((1 - TIn$theta12^2) * (1 - TIn$theta13^2)))^2))
  pst.2 <- mm(pst.2, min.pr = VC$min.pr, max.pr = VC$max.pr)
  st.2 <- -TIn$theta13/sqrt(1 - TIn$theta13^2)
  
  dst.3 <- dnorm((TIn$mar1[VC$inde1] - TIn$theta12 * TIn$mar2)/sqrt(1 - TIn$theta12^2))
  pst.3 <- pnorm( ( ((TIn$mar3 - TIn$theta23 * TIn$mar2[VC$inde2.1])/sqrt(1 - TIn$theta23^2)) - ((TIn$theta13 - TIn$theta12 * TIn$theta23)/sqrt((1 - TIn$theta12^2) * (1 - TIn$theta23^2))) * ((TIn$mar1[VC$inde2] - TIn$theta12 * TIn$mar2[VC$inde2.1])/sqrt(1 - TIn$theta12^2)) )/sqrt(1 - ((TIn$theta13 - TIn$theta12 * TIn$theta23)/sqrt((1 - TIn$theta12^2) * (1 - TIn$theta23^2)))^2))
  pst.3 <- mm(pst.3, min.pr = VC$min.pr, max.pr = VC$max.pr)
  st.3 <- -TIn$theta12/sqrt(1 - TIn$theta12^2)
  
  dst.4 <- dnorm((TIn$mar3 - TIn$theta23 * TIn$mar2[VC$inde2.1])/sqrt(1 - TIn$theta23^2))
  pst.4 <- pnorm( ( ((TIn$mar1[VC$inde2] - TIn$theta12 * TIn$mar2[VC$inde2.1])/sqrt(1 - TIn$theta12^2)) - ((TIn$theta13 - TIn$theta12 * TIn$theta23)/sqrt((1 - TIn$theta12^2) * (1 - TIn$theta23^2))) * ((TIn$mar3 - TIn$theta23 * TIn$mar2[VC$inde2.1])/sqrt(1 - TIn$theta23^2)) )/sqrt(1 - ((TIn$theta13 - TIn$theta12 * TIn$theta23)/sqrt((1 - TIn$theta12^2) * (1 - TIn$theta23^2)))^2))
  pst.4 <- mm(pst.4, min.pr = VC$min.pr, max.pr = VC$max.pr)
  st.4  <- -TIn$theta23/sqrt(1 - TIn$theta23^2)
  
  dst.5 <- dnorm((TIn$mar1[VC$inde2] - TIn$theta13 * TIn$mar3)/sqrt(1 - TIn$theta13^2))
  pst.5 <- pnorm( ( ((TIn$mar2[VC$inde2.1]  - TIn$theta23 * TIn$mar3)/sqrt(1 - TIn$theta23^2)) - ((TIn$theta12 - TIn$theta13 * TIn$theta23)/sqrt((1 - TIn$theta13^2) * (1 - TIn$theta23^2))) * ((TIn$mar1[VC$inde2] - TIn$theta13 * TIn$mar3)/sqrt(1 - TIn$theta13^2)) )/sqrt(1 - ((TIn$theta12 - TIn$theta13 * TIn$theta23)/sqrt((1 - TIn$theta13^2) * (1 - TIn$theta23^2)))^2))
  pst.5 <- mm(pst.5, min.pr = VC$min.pr, max.pr = VC$max.pr)
  st.5  <- -TIn$theta13/sqrt(1 - TIn$theta13^2)
  
  dst.6 <- dnorm((TIn$mar2[VC$inde2.1] - TIn$theta23 * TIn$mar3)/sqrt(1 - TIn$theta23^2))
  pst.6 <- pnorm( ( ((TIn$mar1[VC$inde2] - TIn$theta13 * TIn$mar3)/sqrt(1 - TIn$theta13^2)) - ((TIn$theta12 - TIn$theta13 * TIn$theta23)/sqrt((1 - TIn$theta13^2) * (1 - TIn$theta23^2))) * ((TIn$mar2[VC$inde2.1] - TIn$theta23 * TIn$mar3)/sqrt(1 - TIn$theta23^2)) )/sqrt(1 - ((TIn$theta12 - TIn$theta13 * TIn$theta23)/sqrt((1 - TIn$theta13^2) * (1 - TIn$theta23^2)))^2))
  pst.6 <- mm(pst.6, min.pr = VC$min.pr, max.pr = VC$max.pr)
  st.6  <- -TIn$theta23/sqrt(1 - TIn$theta23^2)
  
  ########################################################################
  
  dp.1.11.de1 <- dst.1[VC$inde2.1] * pst.1        * st.1      + dst.2 * pst.2       * st.2
  dp.1.10.de1 <- dst.1[VC$inde2.1] * (1 - pst.1)  * st.1      + dst.2 * pst.2       * ( -st.2 )
  
  dp.2.11.de2 <- dst.3[VC$inde2.1] * pst.3        * st.3      + dst.4 * pst.4       * st.4
  dp.2.10.de2 <- dst.3[VC$inde2.1] * (1 - pst.3)  * st.3      + dst.4 * pst.4       * ( -st.4 )
  
  dp.3.11.de3 <- dst.5 * pst.5        * st.5      + dst.6 * pst.6       * st.6
  
  #####
  # ! #
  ###############################
  ## The next 6 lines are new  ##                 
  ###############################
  
  der2p.dereta1 <- probm(TIn$eta1, VC$margins[1], only.pr = FALSE, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)$der2p.dereta
  der2p.dereta2 <- probm(TIn$eta2, VC$margins[2], only.pr = FALSE, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)$der2p.dereta
  der2p.dereta3 <- probm(TIn$eta3, VC$margins[3], only.pr = FALSE, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)$der2p.dereta
  
  d2F1.de1 <- (TIn$mar1 * LgTRI$dmar1^2)/LgTRI$d.1^2 + der2p.dereta1/LgTRI$d.1
  d2F2.de2 <- (TIn$mar2 * LgTRI$dmar2^2)/LgTRI$d.2^2 + der2p.dereta2/LgTRI$d.2
  d2F3.de3 <- (TIn$mar3 * LgTRI$dmar3^2)/LgTRI$d.3^2 + der2p.dereta3/LgTRI$d.3
  
  ###################################################################
  
  d2l.dF1.F1.1           <- respvec$cy1       * ( -1/TIn$p0^2   * (LgTRI$d.1)^2 + 1/TIn$p0 * (TIn$mar1 * LgTRI$d.1) )
  d2l.dF1.F1.1[VC$inde1] <- respvec$y1.cy2    * ( -1/TIn$p10^2  * (LgTRI$d.1[VC$inde1] - LgTRI$d.1[VC$inde1] * LgTRI$upst.1)^2 + 1/TIn$p10  * ( - TIn$mar1[VC$inde1] * LgTRI$d.1[VC$inde1] + TIn$mar1[VC$inde1] * LgTRI$d.1[VC$inde1] * LgTRI$upst.1 - LgTRI$d.1[VC$inde1] * dst.1 * st.1) ) 
  d2l.dF1.F1.1[VC$inde2] <- respvec$y1.y2.cy3 * ( -1/TIn$p110^2 * (LgTRI$d.1[VC$inde2] * LgTRI$p.1.10)^2   + 1/TIn$p110 * ( -TIn$mar1[VC$inde2] * LgTRI$d.1[VC$inde2] * LgTRI$p.1.10   + LgTRI$d.1[VC$inde2] * dp.1.10.de1)  ) +
    respvec$y1.y2.y3  * ( -1/TIn$p111^2 * (LgTRI$d.1[VC$inde2] * LgTRI$p.1.11)^2   + 1/TIn$p111 * ( -TIn$mar1[VC$inde2] * LgTRI$d.1[VC$inde2] * LgTRI$p.1.11   + LgTRI$d.1[VC$inde2] * dp.1.11.de1)  ) 
  d2l.dF1.F1 <- d2l.dF1.F1.1
  
  #####
  # ! #
  ###########################
  ## The following is new: ##                 
  ###########################
  
  d2l.de1.e1 <- d2l.dF1.F1 * LgTRI$dF1.de1^2 + LgTRI$dl.dF1 * d2F1.de1
  
  ################################################################
  
  d2l.dF2.F2.1        <- respvec$y1.cy2 * ( -1/TIn$p10^2  * ( - LgTRI$d.2 * LgTRI$upst.2)^2 + 1/TIn$p10  * ( TIn$mar2 * LgTRI$d.2 * LgTRI$upst.2 - LgTRI$d.2 * dst.3 * st.3) ) 
  d2l.dF2.F2.1[VC$inde2.1] <- respvec$y1.y2.cy3 * ( -1/TIn$p110^2 * (LgTRI$d.2[VC$inde2.1] * LgTRI$p.2.10)^2 + 1/TIn$p110 * ( -TIn$mar2[VC$inde2.1]  * LgTRI$d.2[VC$inde2.1] * LgTRI$p.2.10 + LgTRI$d.2[VC$inde2.1] * dp.2.10.de2) ) +
    respvec$y1.y2.y3  * ( -1/TIn$p111^2 * (LgTRI$d.2[VC$inde2.1] * LgTRI$p.2.11)^2 + 1/TIn$p111 * ( -TIn$mar2[VC$inde2.1]  * LgTRI$d.2[VC$inde2.1] * LgTRI$p.2.11 + LgTRI$d.2[VC$inde2.1] * dp.2.11.de2) ) 
  d2l.dF2.F2 <- d2l.dF2.F2.1
  
  #####
  # ! #
  ###########################
  ## The following is new: ##                 
  ###########################
  
  d2l.de2.e2 <- d2l.dF2.F2 * LgTRI$dF2.de2^2 + LgTRI$dl.dF2 * d2F2.de2
  
  ################################################################
  
  d2l.dF3.F3 <-  - respvec$y1.y2.cy3 * (  1/TIn$p110^2 * (LgTRI$d.3 * LgTRI$p.3.11)^2 + 1/TIn$p110 * ( -TIn$mar3 * LgTRI$d.3 * LgTRI$p.3.11 + LgTRI$d.3 * dp.3.11.de3) ) +
    respvec$y1.y2.y3 * ( -1/TIn$p111^2 * (LgTRI$d.3 * LgTRI$p.3.11)^2 + 1/TIn$p111 * ( -TIn$mar3 * LgTRI$d.3 * LgTRI$p.3.11 + LgTRI$d.3 * dp.3.11.de3) ) 
  
  #####
  # ! #
  ###########################
  ## The following is new: ##                 
  ###########################
  
  d2l.de3.e3 <- d2l.dF3.F3 * LgTRI$dF3.de3^2 + LgTRI$dl.dF3 * d2F3.de3
  
  ################################################################
  
  
  #####
  # ! #
  ########################################################################
  ## I replaced TIn$eta1 with TIn$mar1. Same for TIn$eta2 and TIn$eta3  ##                 
  ########################################################################
  
  mean11 <- TIn$theta23 * TIn$mar2 + ((TIn$theta13 - TIn$theta12 * TIn$theta23) * (TIn$mar1[VC$inde1]   - TIn$theta12 * TIn$mar2))/(1 - TIn$theta12^2)
  mean22 <- TIn$theta23 * TIn$mar3 + ((TIn$theta12 - TIn$theta13 * TIn$theta23) * (TIn$mar1[VC$inde2]   - TIn$theta13 * TIn$mar3))/(1 - TIn$theta13^2)
  mean33 <- TIn$theta13 * TIn$mar3 + ((TIn$theta12 - TIn$theta13 * TIn$theta23) * (TIn$mar2[VC$inde2.1] - TIn$theta23 * TIn$mar3))/(1 - TIn$theta23^2)
  
  ###################################################################################
  
  deno <- 1 - TIn$theta12^2 - TIn$theta13^2 - TIn$theta23^2 + 2 * TIn$theta12 * TIn$theta13 * TIn$theta23
  
  sd11 <- sqrt( deno / ( 1 - TIn$theta12^2 ) )
  sd22 <- sqrt( deno / ( 1 - TIn$theta13^2 ) )
  sd33 <- sqrt( deno / ( 1 - TIn$theta23^2 ) )
  
  #####
  # ! #
  ########################################################################
  ## I replaced TIn$eta1 with TIn$mar1. Same for TIn$eta2 and TIn$eta3  ##                 
  ########################################################################
  
  p1.1 <- mm( pnorm((TIn$mar3             - mean11[VC$inde2.1])/sd11), min.pr = VC$min.pr, max.pr = VC$max.pr )
  p2.2 <- mm( pnorm((TIn$mar2[VC$inde2.1] - mean22)/sd22), min.pr = VC$min.pr, max.pr = VC$max.pr )
  p3.3 <- mm( pnorm((TIn$mar1[VC$inde2]   - mean33)/sd33) , min.pr = VC$min.pr, max.pr = VC$max.pr)
  
  ###########################################################################
  
  p1.1.c <- mm(1 - p1.1, min.pr = VC$min.pr, max.pr = VC$max.pr)
  
  #####
  # ! #
  ########################################################################
  ## I replaced TIn$eta1 with TIn$mar1. Same for TIn$eta2 and TIn$eta3  ##                 
  ########################################################################
  
  d.1.1  <- dnorm(TIn$mar1[VC$inde1]  , mean = TIn$theta12 * TIn$mar2, sd = sqrt(1 - TIn$theta12^2))
  d.1.2  <- dnorm(TIn$mar1[VC$inde2]  , mean = TIn$theta13 * TIn$mar3, sd = sqrt(1 - TIn$theta13^2))
  d.1.3  <- dnorm(TIn$mar2[VC$inde2.1], mean = TIn$theta23 * TIn$mar3, sd = sqrt(1 - TIn$theta23^2))
  
  ###########################################################################
  
  d2l.dF1.F2.1          <- respvec$y1.cy2 * ( -1/TIn$p10^2  * (LgTRI$d.1[VC$inde1] - LgTRI$d.1[VC$inde1] * LgTRI$upst.1) * ( - LgTRI$d.2 * LgTRI$upst.2) - 1/TIn$p10 * LgTRI$d.2 * dst.3 * 1/sqrt(1 - TIn$theta12^2)) 
  d2l.dF1.F2.1[VC$inde2.1] <- respvec$y1.y2.cy3 * ( -1/TIn$p110^2 * LgTRI$d.2[VC$inde2.1] * LgTRI$p.2.10 * LgTRI$d.1[VC$inde2] * LgTRI$p.1.10 + 1/TIn$p110 * LgTRI$d.2[VC$inde2.1] * d.1.1[VC$inde2.1] * p1.1.c ) + 
    respvec$y1.y2.y3  * ( -1/TIn$p111^2 * LgTRI$d.2[VC$inde2.1] * LgTRI$p.2.11 * LgTRI$d.1[VC$inde2] * LgTRI$p.1.11 + 1/TIn$p111 * LgTRI$d.2[VC$inde2.1] * d.1.1[VC$inde2.1] * p1.1   )
  d2l.dF1.F2 <- d2l.dF1.F2.1
  
  #####
  # ! #
  ###########################
  ## The following is new: ##                 
  ###########################
  
  d2l.de1.e2 <- d2l.dF1.F2 * LgTRI$dF1.de1[VC$inde1] * LgTRI$dF2.de2  
  
  ###########################################################
  
  d2l.dF1.F3 <- - respvec$y1.y2.cy3 * ( -1/TIn$p110^2 * LgTRI$d.3 * LgTRI$p.3.11 * LgTRI$d.1[VC$inde2] * LgTRI$p.1.10 + 1/TIn$p110 * LgTRI$d.3 * d.1.2 * p2.2   ) + 
    respvec$y1.y2.y3 * ( -1/TIn$p111^2 * LgTRI$d.3 * LgTRI$p.3.11 * LgTRI$d.1[VC$inde2] * LgTRI$p.1.11 + 1/TIn$p111 * LgTRI$d.3 * d.1.2 * p2.2   )
  
  #####
  # ! #
  ###########################
  ## The following is new: ##                 
  ###########################
  
  d2l.de1.e3 <- d2l.dF1.F3 * LgTRI$dF1.de1[VC$inde2] * LgTRI$dF3.de3  
  
  ###########################################################
  
  d2l.dF2.F3 <- - respvec$y1.y2.cy3 * ( -1/TIn$p110^2 * LgTRI$d.3 * LgTRI$p.3.11 * LgTRI$d.2[VC$inde2.1] * LgTRI$p.2.10 + 1/TIn$p110 * LgTRI$d.3 * d.1.3 * p3.3   ) + 
    respvec$y1.y2.y3 * ( -1/TIn$p111^2 * LgTRI$d.3 * LgTRI$p.3.11 * LgTRI$d.2[VC$inde2.1] * LgTRI$p.2.11 + 1/TIn$p111 * LgTRI$d.3 * d.1.3 * p3.3   ) 
  
  #####
  # ! #
  ###########################
  ## The following is new: ##                 
  ###########################
  
  d2l.de2.e3 <- d2l.dF2.F3 * LgTRI$dF2.de2[VC$inde2.1] * LgTRI$dF3.de3  
  
  ###########################################################
  
  #####
  # ! #
  ########################################################################
  ## I replaced TIn$eta1 with TIn$mar1. Same for TIn$eta2 and TIn$eta3  ##                 
  ########################################################################
  
  d12 <- dnorm( (TIn$mar3             - LgTRI$mean.12[VC$inde2.1])/LgTRI$sd.12 )
  d13 <- dnorm( (TIn$mar2[VC$inde2.1] - LgTRI$mean.13)/LgTRI$sd.13 )
  d23 <- dnorm( (TIn$mar1[VC$inde2]   - LgTRI$mean.23)/LgTRI$sd.23 )
  
  #############################################################################
  
  d2l.dF1.theta12.1          <- respvec$y1.cy2 * ( -1/TIn$p10^2 * (LgTRI$d.1[VC$inde1] - LgTRI$d.1[VC$inde1] * LgTRI$upst.1) * ( - LgTRI$d11.12) + 1/TIn$p10 * LgTRI$d11.12 * (TIn$mar1[VC$inde1] - TIn$theta12 * TIn$mar2)/(1 - TIn$theta12^2) )
  d2l.dF1.theta12.1[VC$inde2.1] <- respvec$y1.y2.cy3 * ( -1/TIn$p110^2 * LgTRI$d.1[VC$inde2] * LgTRI$p.1.10 * LgTRI$d11.12[VC$inde2.1] * LgTRI$p12.g.c + 1/TIn$p110 * ( LgTRI$d11.12[VC$inde2.1] * (TIn$theta12*TIn$mar2[VC$inde2.1]  - TIn$mar1[VC$inde2])/(1 - TIn$theta12^2) * LgTRI$p12.g.c + LgTRI$d11.12[VC$inde2.1] * d12/LgTRI$sd.12 * ((TIn$theta13 - TIn$theta12 * TIn$theta23)/(1 - TIn$theta12^2)) )) + 
    respvec$y1.y2.y3  * ( -1/TIn$p111^2 * LgTRI$d.1[VC$inde2] * LgTRI$p.1.11 * LgTRI$d11.12[VC$inde2.1] * LgTRI$p12.g   + 1/TIn$p111 * ( LgTRI$d11.12[VC$inde2.1] * (TIn$theta12*TIn$mar2[VC$inde2.1]  - TIn$mar1[VC$inde2])/(1 - TIn$theta12^2) * LgTRI$p12.g   - LgTRI$d11.12[VC$inde2.1] * d12/LgTRI$sd.12 * ((TIn$theta13 - TIn$theta12 * TIn$theta23)/(1 - TIn$theta12^2)) )) 
  d2l.dF1.theta12 <- d2l.dF1.theta12.1
  
  #####
  # ! #
  ###########################
  ## The following is new: ##                 
  ###########################
  
  d2l.de1.theta12 <- d2l.dF1.theta12 * LgTRI$dF1.de1[VC$inde1]
  
  ###########################################################
  
  
  d2l.dF1.theta13 <- - respvec$y1.y2.cy3   * ( -1/TIn$p110^2 * LgTRI$d.1[VC$inde2] * LgTRI$p.1.10 * LgTRI$d11.13 * LgTRI$p13.g   + 1/TIn$p110 * ( LgTRI$d11.13 * (TIn$theta13 * TIn$mar3 - TIn$mar1[VC$inde2])/(1 - TIn$theta13^2) * LgTRI$p13.g   - LgTRI$d11.13 * d13/LgTRI$sd.13 * ((TIn$theta12 - TIn$theta13 * TIn$theta23)/(1 - TIn$theta13^2)) )) +
    respvec$y1.y2.y3 * ( -1/TIn$p111^2 * LgTRI$d.1[VC$inde2] * LgTRI$p.1.11 * LgTRI$d11.13 * LgTRI$p13.g   + 1/TIn$p111 * ( LgTRI$d11.13 * (TIn$theta13 * TIn$mar3 - TIn$mar1[VC$inde2])/(1 - TIn$theta13^2) * LgTRI$p13.g   - LgTRI$d11.13 * d13/LgTRI$sd.13 * ((TIn$theta12 - TIn$theta13 * TIn$theta23)/(1 - TIn$theta13^2)) )) 
  
  #####
  # ! #
  ###########################
  ## The following is new: ##                 
  ###########################
  
  d2l.de1.theta13 <- d2l.dF1.theta13 * LgTRI$dF1.de1[VC$inde2]
  
  ###########################################################
  
  d2l.dF1.theta23 <- - respvec$y1.y2.cy3  * ( -1/TIn$p110^2 * LgTRI$d.1[VC$inde2] * LgTRI$p.1.10 * LgTRI$d11.23 * LgTRI$p23.g   + 1/TIn$p110 * LgTRI$d11.23 * d23/LgTRI$sd.23 ) +
    respvec$y1.y2.y3 * ( -1/TIn$p111^2 * LgTRI$d.1[VC$inde2] * LgTRI$p.1.11 * LgTRI$d11.23 * LgTRI$p23.g   + 1/TIn$p111 * LgTRI$d11.23 * d23/LgTRI$sd.23) 
  
  #####
  # ! #
  ###########################
  ## The following is new: ##                 
  ###########################
  
  d2l.de1.theta23 <- d2l.dF1.theta23 * LgTRI$dF1.de1[VC$inde2]
  
  ###########################################################
  
  d2l.dF2.theta12.1          <- respvec$y1.cy2 * ( -1/TIn$p10^2 * LgTRI$d.2 *  LgTRI$upst.2 * LgTRI$d11.12 - 1/TIn$p10  * LgTRI$d11.12  * (TIn$theta12 * TIn$mar1[VC$inde1] - TIn$mar2)/(1 - TIn$theta12^2) )
  d2l.dF2.theta12.1[VC$inde2.1] <- respvec$y1.y2.cy3 * ( -1/TIn$p110^2 * LgTRI$d.2[VC$inde2.1] * LgTRI$p.2.10 * LgTRI$d11.12[VC$inde2.1] * LgTRI$p12.g.c + 1/TIn$p110 * (LgTRI$d11.12[VC$inde2.1] * (TIn$theta12 * TIn$mar1[VC$inde2] - TIn$mar2[VC$inde2.1])/(1 - TIn$theta12^2) * LgTRI$p12.g.c + LgTRI$d11.12[VC$inde2.1] * (d12/LgTRI$sd.12) * (TIn$theta23 - TIn$theta12 * TIn$theta13)/(1 - TIn$theta12^2) ) ) + 
    respvec$y1.y2.y3  * ( -1/TIn$p111^2 * LgTRI$d.2[VC$inde2.1] * LgTRI$p.2.11 * LgTRI$d11.12[VC$inde2.1] * LgTRI$p12.g   + 1/TIn$p111 * (LgTRI$d11.12[VC$inde2.1] * (TIn$theta12 * TIn$mar1[VC$inde2] - TIn$mar2[VC$inde2.1])/(1 - TIn$theta12^2) * LgTRI$p12.g   - LgTRI$d11.12[VC$inde2.1] * (d12/LgTRI$sd.12) * (TIn$theta23 - TIn$theta12 * TIn$theta13)/(1 - TIn$theta12^2) ) ) 
  d2l.dF2.theta12 <- d2l.dF2.theta12.1
  
  #####
  # ! #
  ###########################
  ## The following is new: ##                 
  ###########################
  
  d2l.de2.theta12 <- d2l.dF2.theta12 * LgTRI$dF2.de2
  
  ###########################################################
  
  
  d2l.dF2.theta13 <- - respvec$y1.y2.cy3 * ( -1/TIn$p110^2 * LgTRI$d.2[VC$inde2.1] * LgTRI$p.2.10 * LgTRI$d11.13 * LgTRI$p13.g   + 1/TIn$p110 * LgTRI$d11.13 * d13/LgTRI$sd.13 ) + 
    respvec$y1.y2.y3 * ( -1/TIn$p111^2 * LgTRI$d.2[VC$inde2.1] * LgTRI$p.2.11 * LgTRI$d11.13 * LgTRI$p13.g   + 1/TIn$p111 * LgTRI$d11.13 * d13/LgTRI$sd.13 )
  
  #####
  # ! #
  ###########################
  ## The following is new: ##                 
  ###########################
  
  d2l.de2.theta13 <- d2l.dF2.theta13 * LgTRI$dF2.de2[VC$inde2.1]
  
  ###########################################################
  
  d2l.dF2.theta23 <- - respvec$y1.y2.cy3 * ( -1/TIn$p110^2 * LgTRI$d.2[VC$inde2.1] * LgTRI$p.2.10 * LgTRI$d11.23 * LgTRI$p23.g   + 1/TIn$p110 * (LgTRI$d11.23 * (TIn$theta23 * TIn$mar3 - TIn$mar2[VC$inde2.1])/(1 - TIn$theta23^2) * LgTRI$p23.g   - LgTRI$d11.23 * d23/LgTRI$sd.23 * (TIn$theta12 - TIn$theta13 * TIn$theta23)/(1 - TIn$theta23^2) ) ) + 
    respvec$y1.y2.y3 * ( -1/TIn$p111^2 * LgTRI$d.2[VC$inde2.1] * LgTRI$p.2.11 * LgTRI$d11.23 * LgTRI$p23.g   + 1/TIn$p111 * (LgTRI$d11.23 * (TIn$theta23 * TIn$mar3 - TIn$mar2[VC$inde2.1])/(1 - TIn$theta23^2) * LgTRI$p23.g   - LgTRI$d11.23 * d23/LgTRI$sd.23 * (TIn$theta12 - TIn$theta13 * TIn$theta23)/(1 - TIn$theta23^2) ) )
  
  #####
  # ! #
  ###########################
  ## The following is new: ##                 
  ###########################
  
  d2l.de2.theta23 <- d2l.dF2.theta23 * LgTRI$dF2.de2[VC$inde2.1]
  
  ###########################################################
  
  d2l.dF3.theta12 <- - respvec$y1.y2.cy3 * ( -1/TIn$p110^2 * LgTRI$d.3 * LgTRI$p.3.11 * LgTRI$d11.12[VC$inde2.1] * LgTRI$p12.g.c + 1/TIn$p110 * LgTRI$d11.12[VC$inde2.1] * d12/LgTRI$sd.12 ) +
    respvec$y1.y2.y3 * ( -1/TIn$p111^2 * LgTRI$d.3 * LgTRI$p.3.11 * LgTRI$d11.12[VC$inde2.1] * LgTRI$p12.g   + 1/TIn$p111 * LgTRI$d11.12[VC$inde2.1] * d12/LgTRI$sd.12 ) 
  
  
  #####
  # ! #
  ###########################
  ## The following is new: ##                 
  ###########################
  
  d2l.de3.theta12 <- d2l.dF3.theta12 * LgTRI$dF3.de3
  
  ###########################################################
  
  
  d2l.dF3.theta13 <- - respvec$y1.y2.cy3 * (  1/TIn$p110^2 * LgTRI$d.3 * LgTRI$p.3.11 * LgTRI$d11.13 * LgTRI$p13.g   + 1/TIn$p110 * (LgTRI$d11.13 * ((TIn$theta13 * TIn$mar1[VC$inde2] - TIn$mar3)/(1 - TIn$theta13^2)) * LgTRI$p13.g   - LgTRI$d11.13 * (d13/LgTRI$sd.13) * (TIn$theta23 - TIn$theta12 * TIn$theta13)/(1 - TIn$theta13^2) ) ) +
    respvec$y1.y2.y3 * ( -1/TIn$p111^2 * LgTRI$d.3 * LgTRI$p.3.11 * LgTRI$d11.13 * LgTRI$p13.g   + 1/TIn$p111 * (LgTRI$d11.13 * ((TIn$theta13*TIn$mar1[VC$inde2] - TIn$mar3)/(1 - TIn$theta13^2)) * LgTRI$p13.g   - LgTRI$d11.13 * (d13/LgTRI$sd.13) * (TIn$theta23 - TIn$theta12 * TIn$theta13)/(1 - TIn$theta13^2) ) ) 
  
  #####
  # ! #
  ###########################
  ## The following is new: ##                 
  ###########################
  
  d2l.de3.theta13 <- d2l.dF3.theta13 * LgTRI$dF3.de3
  
  ###########################################################
  
  
  d2l.dF3.theta23 <- - respvec$y1.y2.cy3 * (  1/TIn$p110^2 * LgTRI$d.3 * LgTRI$p.3.11 * LgTRI$d11.23 * LgTRI$p23.g   + 1/TIn$p110 * (LgTRI$d11.23 * ((TIn$theta23*TIn$mar2[VC$inde2.1]  - TIn$mar3)/(1 - TIn$theta23^2)) * LgTRI$p23.g   - LgTRI$d11.23 * d23/LgTRI$sd.23 * (TIn$theta13 - TIn$theta12*TIn$theta23)/(1 - TIn$theta23^2) ) ) +
    respvec$y1.y2.y3 * ( -1/TIn$p111^2 * LgTRI$d.3 * LgTRI$p.3.11 * LgTRI$d11.23 * LgTRI$p23.g   + 1/TIn$p111 * (LgTRI$d11.23 * ((TIn$theta23*TIn$mar2[VC$inde2.1]  - TIn$mar3)/(1 - TIn$theta23^2)) * LgTRI$p23.g   - LgTRI$d11.23 * d23/LgTRI$sd.23 * (TIn$theta13 - TIn$theta12 * TIn$theta23)/(1 - TIn$theta23^2) ) ) 
  
  #####
  # ! #
  ###########################
  ## The following is new: ##                 
  ###########################
  
  d2l.de3.theta23 <- d2l.dF3.theta23 * LgTRI$dF3.de3
  
  ###########################################################
  
  
  #####
  # ! #
  ########################################################################
  ## I replaced TIn$eta1 with TIn$mar1. Same for TIn$eta2 and TIn$eta3  ##                 
  ########################################################################
  
  dd11.12.dtheta12 <- (LgTRI$d11.12/(1 - TIn$theta12^2)) * ( ( (TIn$theta12 * TIn$mar2 - TIn$mar1[VC$inde1]  ) * (TIn$theta12 * TIn$mar1[VC$inde1]   - TIn$mar2)/(1 - TIn$theta12^2) ) + TIn$theta12)
  dd11.13.dtheta13 <- (LgTRI$d11.13/(1 - TIn$theta13^2)) * ( ( (TIn$theta13 * TIn$mar3 - TIn$mar1[VC$inde2]  ) * (TIn$theta13 * TIn$mar1[VC$inde2]   - TIn$mar3)/(1 - TIn$theta13^2) ) + TIn$theta13)
  dd11.23.dtheta23 <- (LgTRI$d11.23/(1 - TIn$theta23^2)) * ( ( (TIn$theta23 * TIn$mar3 - TIn$mar2[VC$inde2.1]) * (TIn$theta23 * TIn$mar2[VC$inde2.1] - TIn$mar3)/(1 - TIn$theta23^2) ) + TIn$theta23)
  
  
  dmean12.dtheta12 <- ( (-TIn$mar1[VC$inde1]   * TIn$theta23 - TIn$mar2 * TIn$theta13) * (1 - TIn$theta12^2) + 2 * TIn$theta12 * (TIn$mar1[VC$inde1]   * (TIn$theta13 - TIn$theta12 * TIn$theta23) + TIn$mar2 * (TIn$theta23 - TIn$theta12 * TIn$theta13)))/(1 - TIn$theta12^2)^2
  dmean13.dtheta13 <- ( (-TIn$mar1[VC$inde2]   * TIn$theta23 - TIn$mar3 * TIn$theta12) * (1 - TIn$theta13^2) + 2 * TIn$theta13 * (TIn$mar1[VC$inde2]   * (TIn$theta12 - TIn$theta13 * TIn$theta23) + TIn$mar3 * (TIn$theta23 - TIn$theta12 * TIn$theta13)))/(1 - TIn$theta13^2)^2
  dmean23.dtheta23 <- ( (-TIn$mar2[VC$inde2.1] * TIn$theta13 - TIn$mar3 * TIn$theta12) * (1 - TIn$theta23^2) + 2 * TIn$theta23 * (TIn$mar2[VC$inde2.1] * (TIn$theta12 - TIn$theta13 * TIn$theta23) + TIn$mar3 * (TIn$theta13 - TIn$theta12 * TIn$theta23)))/(1 - TIn$theta23^2)^2
  
  dmean13.dtheta12 <- ( TIn$mar1[VC$inde2]   - TIn$mar3              * TIn$theta13 )/( 1 - TIn$theta13^2 )
  dmean23.dtheta12 <- ( TIn$mar2[VC$inde2.1] - TIn$mar3              * TIn$theta23 )/( 1 - TIn$theta23^2 )
  dmean23.dtheta13 <- ( TIn$mar3             - TIn$mar2[VC$inde2.1]  * TIn$theta23 )/( 1 - TIn$theta23^2 )
  
  #############################################################################################################
  
  dvar12.dtheta12 <- ( 2 * (TIn$theta13 * TIn$theta23 - TIn$theta12) * (1 - TIn$theta12^2) + 2 * TIn$theta12 * deno )/(1 - TIn$theta12^2)^2
  dvar13.dtheta13 <- ( 2 * (TIn$theta12 * TIn$theta23 - TIn$theta13) * (1 - TIn$theta13^2) + 2 * TIn$theta13 * deno )/(1 - TIn$theta13^2)^2
  dvar23.dtheta23 <- ( 2 * (TIn$theta12 * TIn$theta13 - TIn$theta23) * (1 - TIn$theta23^2) + 2 * TIn$theta23 * deno )/(1 - TIn$theta23^2)^2
  
  dvar13.dtheta12 <- ( 2 * (TIn$theta13 * TIn$theta23 - TIn$theta12 ) )/(1 - TIn$theta13^2)
  dvar23.dtheta12 <- ( 2 * (TIn$theta13 * TIn$theta23 - TIn$theta12 ) )/(1 - TIn$theta23^2)
  dvar23.dtheta13 <- ( 2 * (TIn$theta12 * TIn$theta23 - TIn$theta13 ) )/(1 - TIn$theta23^2)
  
  
  d2l.dtheta12.theta12.1          <- respvec$y1.cy2 * ( -1/TIn$p10^2  * LgTRI$d11.12^2 - 1/TIn$p10 * dd11.12.dtheta12)
  d2l.dtheta12.theta12.1[VC$inde2.1] <- respvec$y1.y2.cy3 * ( -1/TIn$p110^2 * (LgTRI$d11.12[VC$inde2.1] * LgTRI$p12.g.c )^2 + 1/TIn$p110 * (( dd11.12.dtheta12[VC$inde2.1] * LgTRI$p12.g.c) + ((d12/LgTRI$sd.12) * (dmean12.dtheta12[VC$inde2.1]  + (((TIn$mar3 - LgTRI$mean.12[VC$inde2.1])/(2 * LgTRI$sd.12^2)) * dvar12.dtheta12)) * LgTRI$d11.12[VC$inde2.1]) )) + 
    respvec$y1.y2.y3  * ( -1/TIn$p111^2 * (LgTRI$d11.12[VC$inde2.1] * LgTRI$p12.g   )^2 + 1/TIn$p111 * (( dd11.12.dtheta12[VC$inde2.1] * LgTRI$p12.g  ) - ((d12/LgTRI$sd.12) * (dmean12.dtheta12[VC$inde2.1]  + (((TIn$mar3 - LgTRI$mean.12[VC$inde2.1])/(2 * LgTRI$sd.12^2)) * dvar12.dtheta12)) * LgTRI$d11.12[VC$inde2.1]) ))
  
  d2l.dtheta12.theta12 <- d2l.dtheta12.theta12.1
  
  d2l.dtheta13.theta13 <- - respvec$y1.y2.cy3 * (  1/TIn$p110^2 * (LgTRI$d11.13 * LgTRI$p13.g)^2 + 1/TIn$p110 * ((dd11.13.dtheta13 * LgTRI$p13.g  ) - ((d13/LgTRI$sd.13) * (dmean13.dtheta13 + (((TIn$mar2[VC$inde2.1] - LgTRI$mean.13)/(2 * LgTRI$sd.13^2)) * dvar13.dtheta13)) * LgTRI$d11.13) )) + 
    respvec$y1.y2.y3 * ( -1/TIn$p111^2 * (LgTRI$d11.13 * LgTRI$p13.g)^2 + 1/TIn$p111 * (( dd11.13.dtheta13 * LgTRI$p13.g  ) - ((d13/LgTRI$sd.13) * (dmean13.dtheta13 + (((TIn$mar2[VC$inde2.1] - LgTRI$mean.13)/(2 * LgTRI$sd.13^2)) * dvar13.dtheta13)) * LgTRI$d11.13) )) 
  
  d2l.dtheta23.theta23 <- - respvec$y1.y2.cy3 * (  1/TIn$p110^2 * (LgTRI$d11.23 * LgTRI$p23.g   )^2 + 1/TIn$p110 * ((dd11.23.dtheta23 * LgTRI$p23.g  ) - ((d23/LgTRI$sd.23) * (dmean23.dtheta23 + ((TIn$mar1[VC$inde2] - LgTRI$mean.23)/(2*LgTRI$sd.23^2) * dvar23.dtheta23)) * LgTRI$d11.23) )) + 
    respvec$y1.y2.y3 * ( -1/TIn$p111^2 * (LgTRI$d11.23 * LgTRI$p23.g   )^2 + 1/TIn$p111 * ((dd11.23.dtheta23 * LgTRI$p23.g  ) - ((d23/LgTRI$sd.23) * (dmean23.dtheta23 + ((TIn$mar1[VC$inde2] - LgTRI$mean.23)/(2*LgTRI$sd.23^2) * dvar23.dtheta23)) * LgTRI$d11.23) )) 
  
  d2l.dtheta12.theta13 <- - respvec$y1.y2.cy3 * ( -1/TIn$p110^2 * LgTRI$d11.12[VC$inde2.1] * LgTRI$p12.g.c * LgTRI$d11.13 * LgTRI$p13.g   - 1/TIn$p110 * LgTRI$d11.13 * (d13/LgTRI$sd.13) * (dmean13.dtheta12 + ((TIn$mar2[VC$inde2.1] - LgTRI$mean.13)/(2*LgTRI$sd.13^2) * dvar13.dtheta12) )) +
    respvec$y1.y2.y3 * ( -1/TIn$p111^2 * LgTRI$d11.12[VC$inde2.1] * LgTRI$p12.g   * LgTRI$d11.13 * LgTRI$p13.g   - 1/TIn$p111 * LgTRI$d11.13 * (d13/LgTRI$sd.13) * (dmean13.dtheta12 + ((TIn$mar2[VC$inde2.1] - LgTRI$mean.13)/(2*LgTRI$sd.13^2) * dvar13.dtheta12) )) 
  
  d2l.dtheta12.theta23 <- - respvec$y1.y2.cy3 * ( -1/TIn$p110^2 * LgTRI$d11.12[VC$inde2.1] * LgTRI$p12.g.c * LgTRI$d11.23 * LgTRI$p23.g   - 1/TIn$p110 * LgTRI$d11.23 * (d23/LgTRI$sd.23) * (dmean23.dtheta12 + ((TIn$mar1[VC$inde2] - LgTRI$mean.23)/(2*LgTRI$sd.23^2) * dvar23.dtheta12) )) +
    respvec$y1.y2.y3 * ( -1/TIn$p111^2 * LgTRI$d11.12[VC$inde2.1] * LgTRI$p12.g   * LgTRI$d11.23 * LgTRI$p23.g   - 1/TIn$p111 * LgTRI$d11.23 * (d23/LgTRI$sd.23) * (dmean23.dtheta12 + ((TIn$mar1[VC$inde2] - LgTRI$mean.23)/(2*LgTRI$sd.23^2) * dvar23.dtheta12) )) 
  
  d2l.dtheta13.theta23 <- - respvec$y1.y2.cy3 * (  1/TIn$p110^2 * LgTRI$d11.13 * LgTRI$p13.g   * LgTRI$d11.23 * LgTRI$p23.g   - 1/TIn$p110 * LgTRI$d11.23 * (d23/LgTRI$sd.23) * (dmean23.dtheta13 + ((TIn$mar1[VC$inde2] - LgTRI$mean.23)/(2 * LgTRI$sd.23^2) * dvar23.dtheta13) )) +
    respvec$y1.y2.y3 * ( -1/TIn$p111^2 * LgTRI$d11.13 * LgTRI$p13.g   * LgTRI$d11.23 * LgTRI$p23.g   - 1/TIn$p111 * LgTRI$d11.23 * (d23/LgTRI$sd.23) * (dmean23.dtheta13 + ((TIn$mar1[VC$inde2] - LgTRI$mean.23)/(2 * LgTRI$sd.23^2) * dvar23.dtheta13) )) 
  
  
  d2l.dbe1.be1 <- crossprod( VC$X1 * c(VC$weights*d2l.de1.e1), VC$X1 )
  d2l.dbe2.be2 <- crossprod( VC$X2 * c(VC$weights[VC$inde1]*d2l.de2.e2), VC$X2 )
  d2l.dbe3.be3 <- crossprod( VC$X3 * c(VC$weights[VC$inde2]*d2l.de3.e3), VC$X3 )
  d2l.dbe1.be2 <- crossprod( VC$X1[VC$inde1,   ] * c(VC$weights[VC$inde1]*d2l.de1.e2), VC$X2 )
  d2l.dbe1.be3 <- crossprod( VC$X1[VC$inde2,   ] * c(VC$weights[VC$inde2]*d2l.de1.e3), VC$X3 )
  d2l.dbe2.be3 <- crossprod( VC$X2[VC$inde2.1, ] * c(VC$weights[VC$inde2]*d2l.de2.e3), VC$X3 )
  
  if(VC$Chol == FALSE){
    dtheta12.dtheta12.st <- 4 * exp( 2 * TIn$theta12.st )/( exp(2 * TIn$theta12.st) + 1 )^2
    dtheta13.dtheta13.st <- 4 * exp( 2 * TIn$theta13.st )/( exp(2 * TIn$theta13.st) + 1 )^2
    dtheta23.dtheta23.st <- 4 * exp( 2 * TIn$theta23.st )/( exp(2 * TIn$theta23.st) + 1 )^2
    
    d2theta12.theta12.st <- (8 * exp( 2 * TIn$theta12.st ) - 8 * exp( 4 * TIn$theta12.st ))/( exp(2 * TIn$theta12.st) + 1 )^3
    d2theta13.theta13.st <- (8 * exp( 2 * TIn$theta13.st ) - 8 * exp( 4 * TIn$theta13.st ))/( exp(2 * TIn$theta13.st) + 1 )^3
    d2theta23.theta23.st <- (8 * exp( 2 * TIn$theta23.st ) - 8 * exp( 4 * TIn$theta23.st ))/( exp(2 * TIn$theta23.st) + 1 )^3
    
    d2l.dbe1.theta12.st <- colSums(c(VC$weights[VC$inde1]*d2l.de1.theta12* dtheta12.dtheta12.st) * VC$X1[VC$inde1, ]) 
    d2l.dbe1.theta13.st <- colSums(c(VC$weights[VC$inde2]*d2l.de1.theta13* dtheta13.dtheta13.st) * VC$X1[VC$inde2, ]) 
    d2l.dbe1.theta23.st <- colSums(c(VC$weights[VC$inde2]*d2l.de1.theta23* dtheta23.dtheta23.st) * VC$X1[VC$inde2, ]) 
    
    d2l.dbe2.theta12.st <- colSums(c(VC$weights[VC$inde1]*d2l.de2.theta12* dtheta12.dtheta12.st) * VC$X2) 
    d2l.dbe2.theta13.st <- colSums(c(VC$weights[VC$inde2]*d2l.de2.theta13* dtheta13.dtheta13.st) * VC$X2[VC$inde2.1, ]) 
    d2l.dbe2.theta23.st <- colSums(c(VC$weights[VC$inde2]*d2l.de2.theta23* dtheta23.dtheta23.st) * VC$X2[VC$inde2.1, ]) 
    
    d2l.dbe3.theta12.st <- colSums(c(VC$weights[VC$inde2]*d2l.de3.theta12* dtheta12.dtheta12.st) * VC$X3) 
    d2l.dbe3.theta13.st <- colSums(c(VC$weights[VC$inde2]*d2l.de3.theta13* dtheta13.dtheta13.st) * VC$X3) 
    d2l.dbe3.theta23.st <- colSums(c(VC$weights[VC$inde2]*d2l.de3.theta23* dtheta23.dtheta23.st) * VC$X3) 
    
    d2l.dtheta12.st <- sum( VC$weights[VC$inde1]*(d2l.dtheta12.theta12 * (dtheta12.dtheta12.st)^2  + LgTRI$dl.dtheta12 * d2theta12.theta12.st) )
    d2l.dtheta13.st <- sum( VC$weights[VC$inde2]*(d2l.dtheta13.theta13 * (dtheta13.dtheta13.st)^2  + LgTRI$dl.dtheta13 * d2theta13.theta13.st) )
    d2l.dtheta23.st <- sum( VC$weights[VC$inde2]*(d2l.dtheta23.theta23 * (dtheta23.dtheta23.st)^2  + LgTRI$dl.dtheta23 * d2theta23.theta23.st) )
    
    d2l.dtheta12.st.theta13.st <- sum(VC$weights[VC$inde2]*d2l.dtheta12.theta13 * dtheta12.dtheta12.st * dtheta13.dtheta13.st )
    d2l.dtheta12.st.theta23.st <- sum(VC$weights[VC$inde2]*d2l.dtheta12.theta23 * dtheta12.dtheta12.st * dtheta23.dtheta23.st )
    d2l.dtheta13.st.theta23.st <- sum(VC$weights[VC$inde2]*d2l.dtheta13.theta23 * dtheta13.dtheta13.st * dtheta23.dtheta23.st )
    
  }
  
  if(VC$Chol == TRUE){
    
    d2l.de1.theta <- matrix(0,length(VC$inde1),3)
    d2l.de1.theta[VC$inde1, 1] <- VC$weights[VC$inde1]*d2l.de1.theta12
    d2l.de1.theta[VC$inde2, 2] <- VC$weights[VC$inde2]*d2l.de1.theta13
    d2l.de1.theta[VC$inde2, 3] <- VC$weights[VC$inde2]*d2l.de1.theta23
    
    d2l.de2.theta <-  matrix(0,length(VC$inde1),3)
    d2l.de2.theta[VC$inde1, 1] <- VC$weights[VC$inde1]*d2l.de2.theta12
    d2l.de2.theta[VC$inde2, 2] <- VC$weights[VC$inde2]*d2l.de2.theta13
    d2l.de2.theta[VC$inde2, 3] <- VC$weights[VC$inde2]*d2l.de2.theta23
    
    d2l.de3.theta <-  matrix(0,length(VC$inde1),3)
    d2l.de3.theta[VC$inde2, 1] <- VC$weights[VC$inde2]*d2l.de3.theta12
    d2l.de3.theta[VC$inde2, 2] <- VC$weights[VC$inde2]*d2l.de3.theta13
    d2l.de3.theta[VC$inde2, 3] <- VC$weights[VC$inde2]*d2l.de3.theta23
    
    dth12.dth12.st <- 1/(1 + TIn$theta12.st^2)^(3/2) 
    dth12.dth13.st <- 0
    dth12.dth23.st <- 0
    
    dth13.dth12.st <- 0
    dth13.dth13.st <- (1 + TIn$theta23.st^2)/(1 + TIn$theta13.st^2 + TIn$theta23.st^2)^(3/2)
    dth13.dth23.st <- - (TIn$theta13.st * TIn$theta23.st)/(1 + TIn$theta13.st^2 + TIn$theta23.st^2)^(3/2)
    
    
    dth23.dth12.st <- TIn$theta13.st/sqrt( (1 + TIn$theta12.st^2) * (1 + TIn$theta13.st^2 + TIn$theta23.st^2)) - (TIn$theta12.st * (TIn$theta12.st * TIn$theta13.st + TIn$theta23.st))/((1 + TIn$theta12.st^2)^(3/2) * sqrt(1 + TIn$theta13.st^2 + TIn$theta23.st^2))
    dth23.dth13.st <- TIn$theta12.st/sqrt((1 + TIn$theta12.st^2) * (1 + TIn$theta13.st^2 + TIn$theta23.st^2)) - (TIn$theta13.st * (TIn$theta12.st * TIn$theta13.st + TIn$theta23.st))/(sqrt(1 + TIn$theta12.st^2) * (1 + TIn$theta13.st^2 + TIn$theta23.st^2)^(3/2))
    dth23.dth23.st <- 1/sqrt((1 + TIn$theta12.st^2) * (1 + TIn$theta13.st^2 + TIn$theta23.st^2)) - (TIn$theta23.st * (TIn$theta12.st * TIn$theta13.st + TIn$theta23.st))/(sqrt(1 + TIn$theta12.st^2) * (1 + TIn$theta13.st^2 + TIn$theta23.st^2)^(3/2))
    
    dtheta.theta.st <- matrix( c(  dth12.dth12.st,  dth13.dth12.st, dth23.dth12.st,
                                   dth12.dth13.st,  dth13.dth13.st, dth23.dth13.st,
                                   dth12.dth23.st,  dth13.dth23.st, dth23.dth23.st ), 3 , 3)
    
    
    d2l.dbe1.theta12.st <- (t(VC$X1) %*%  d2l.de1.theta %*% dtheta.theta.st)[, 1]
    d2l.dbe1.theta13.st <- (t(VC$X1) %*%  d2l.de1.theta %*% dtheta.theta.st)[, 2]
    d2l.dbe1.theta23.st <- (t(VC$X1) %*%  d2l.de1.theta %*% dtheta.theta.st)[, 3]
    
    d2l.dbe2.theta12.st <- (t(VC$X2) %*%  d2l.de2.theta[VC$inde1, ] %*% dtheta.theta.st)[, 1]
    d2l.dbe2.theta13.st <- (t(VC$X2) %*%  d2l.de2.theta[VC$inde1, ] %*% dtheta.theta.st)[, 2]
    d2l.dbe2.theta23.st <- (t(VC$X2) %*%  d2l.de2.theta[VC$inde1, ] %*% dtheta.theta.st)[, 3]
    
    d2l.dbe3.theta12.st <- (t(VC$X3) %*%  d2l.de3.theta[VC$inde2, ] %*% dtheta.theta.st)[, 1]
    d2l.dbe3.theta13.st <- (t(VC$X3) %*%  d2l.de3.theta[VC$inde2, ] %*% dtheta.theta.st)[, 2]
    d2l.dbe3.theta23.st <- (t(VC$X3) %*%  d2l.de3.theta[VC$inde2, ] %*% dtheta.theta.st)[, 3]
    
    ############################
    ## d2l.dtheta.st.theta.st ## 
    ############################
    d2th12.dth12.st.th12.st <- -(3 * TIn$theta12.st)/(1 + TIn$theta12.st^2)^(5/2)
    d2th12.dth13.st.th12.st <- 0
    d2th12.dth23.st.th12.st <- 0
    d2th12.dth12.st.th13.st <- 0
    d2th12.dth13.st.th13.st <- 0
    d2th12.dth23.st.th13.st <- 0
    d2th12.dth12.st.th23.st <- 0
    d2th12.dth13.st.th23.st <- 0
    d2th12.dth23.st.th23.st <- 0
    
    d2th13.dth12.st.th12.st <- 0
    d2th13.dth13.st.th12.st <- 0
    d2th13.dth23.st.th12.st <- 0
    d2th13.dth12.st.th13.st <- 0
    d2th13.dth13.st.th13.st<- -(3 * TIn$theta13.st * (1 + TIn$theta23.st^2))/(1 + TIn$theta13.st^2 + TIn$theta23.st^2)^(5/2)
    d2th13.dth23.st.th13.st <- (2 * TIn$theta23.st)/(1 + TIn$theta13.st^2 + TIn$theta23.st^2)^(3/2) - (3 * TIn$theta23.st * (1 + TIn$theta23.st^2))/(1 + TIn$theta13.st^2 + TIn$theta23.st^2)^(5/2)
    d2th13.dth12.st.th23.st <- 0
    d2th13.dth13.st.th23.st <- -TIn$theta23.st/(1 + TIn$theta13.st^2 + TIn$theta23.st^2)^(3/2) + (3 * TIn$theta13.st^2 * TIn$theta23.st)/(1 + TIn$theta13.st^2 + TIn$theta23.st^2)^(5/2)
    d2th13.dth23.st.th23.st <- -TIn$theta13.st/(1 + TIn$theta13.st^2 + TIn$theta23.st^2)^(3/2) + (3 * TIn$theta13.st * TIn$theta23.st^2)/(1 + TIn$theta13.st^2 + TIn$theta23.st^2)^(5/2)
    
    d2th23.dth12.st.th12.st <- -(3 * TIn$theta12.st * TIn$theta13.st + TIn$theta23.st)/((1 + TIn$theta12.st^2)^(3/2) * sqrt(1 + TIn$theta13.st^2 + TIn$theta23.st^2)) + (3 * TIn$theta12.st^2 * (TIn$theta12.st * TIn$theta13.st + TIn$theta23.st))/((1 + TIn$theta12.st^2)^(5/2) * sqrt(1 + TIn$theta13.st^2 + TIn$theta23.st^2))
    d2th23.dth13.st.th12.st <- (1 + TIn$theta12.st^2 + TIn$theta23.st^2 + TIn$theta12.st^2 * TIn$theta23.st^2 + TIn$theta12.st^2 * TIn$theta13.st^2 + TIn$theta12.st * TIn$theta13.st * TIn$theta23.st)/((1 + TIn$theta12.st^2) * (1 + TIn$theta13.st^2 + TIn$theta23.st^2))^(3/2) - TIn$theta12.st^2/((1 + TIn$theta12.st^2)^(3/2) * sqrt(1 + TIn$theta13.st^2 + TIn$theta23.st^2))
    d2th23.dth23.st.th12.st <- -(TIn$theta13.st * TIn$theta23.st)/(sqrt(1 + TIn$theta12.st^2) * (1 + TIn$theta13.st^2 + TIn$theta23.st^2)^(3/2)) - TIn$theta12.st/((1 + TIn$theta12.st^2)^(3/2) * sqrt(1 + TIn$theta13.st^2 + TIn$theta23.st^2)) + (TIn$theta12.st * TIn$theta23.st * (TIn$theta12.st * TIn$theta13.st + TIn$theta23.st))/((1 + TIn$theta12.st^2) * (1 + TIn$theta13.st^2 + TIn$theta23.st^2))^(3/2)
    d2th23.dth12.st.th13.st <- (1 + TIn$theta13.st^2 + TIn$theta23.st^2 + TIn$theta12.st * TIn$theta13.st * (TIn$theta12.st * TIn$theta13.st + TIn$theta23.st))/((1 + TIn$theta12.st^2) * (1 + TIn$theta13.st^2 + TIn$theta23.st^2))^(3/2) - TIn$theta13.st^2/(sqrt(1 + TIn$theta12.st^2) * (1 + TIn$theta13.st^2 + TIn$theta23.st^2)^(3/2))
    d2th23.dth13.st.th13.st <- -(3 * TIn$theta12.st * TIn$theta13.st + TIn$theta23.st)/(sqrt(1 + TIn$theta12.st^2) * (1 + TIn$theta13.st^2 + TIn$theta23.st^2)^(3/2)) + (3 * TIn$theta13.st^2 * (TIn$theta12.st * TIn$theta13.st + TIn$theta23.st))/(sqrt(1 + TIn$theta12.st^2) * (1 + TIn$theta13.st^2 + TIn$theta23.st^2)^(5/2))
    d2th23.dth23.st.th13.st <- -(TIn$theta12.st * TIn$theta23.st + TIn$theta13.st)/(sqrt(1 + TIn$theta12.st^2) * (1 + TIn$theta13.st^2 + TIn$theta23.st^2)^(3/2)) + (3 * TIn$theta13.st * TIn$theta23.st * (TIn$theta12.st * TIn$theta13.st + TIn$theta23.st))/(sqrt(1 + TIn$theta12.st^2) * (1 + TIn$theta13.st^2 + TIn$theta23.st^2)^(5/2)) 
    d2th23.dth12.st.th23.st <- -TIn$theta12.st/((1 + TIn$theta12.st^2)^(3/2) * sqrt(1 + TIn$theta13.st^2 + TIn$theta23.st^2)) - (TIn$theta13.st * TIn$theta23.st)/(sqrt(1 + TIn$theta12.st^2) * (1 + TIn$theta13.st^2 + TIn$theta23.st^2)^(3/2)) + (TIn$theta12.st * TIn$theta23.st * (TIn$theta12.st * TIn$theta13.st + TIn$theta23.st))/((1 + TIn$theta12.st^2) * (1 + TIn$theta13.st^2 + TIn$theta23.st^2))^(3/2)
    d2th23.dth13.st.th23.st <- -(TIn$theta13.st + TIn$theta12.st * TIn$theta23.st)/(sqrt(1 + TIn$theta12.st^2) * (1 + TIn$theta13.st^2 + TIn$theta23.st^2)^(3/2)) + (3 * TIn$theta13.st * TIn$theta23.st * (TIn$theta12.st * TIn$theta13.st + TIn$theta23.st))/(sqrt(1 + TIn$theta12.st^2) * (1 + TIn$theta13.st^2 + TIn$theta23.st^2)^(5/2))
    d2th23.dth23.st.th23.st <- -(3 * TIn$theta23.st + TIn$theta12.st * TIn$theta13.st)/(sqrt(1 + TIn$theta12.st^2) * (1 + TIn$theta13.st^2 + TIn$theta23.st^2)^(3/2)) + (3 * TIn$theta23.st^2 * (TIn$theta12.st * TIn$theta13.st + TIn$theta23.st))/(sqrt(1 + TIn$theta12.st^2) * (1 + TIn$theta13.st^2 + TIn$theta23.st^2)^(5/2)) 
    
    
    d2theta12.theta12.st <- matrix( c(          d2th12.dth12.st.th12.st,  d2th13.dth12.st.th12.st, d2th23.dth12.st.th12.st,
                                                d2th12.dth12.st.th13.st,  d2th13.dth12.st.th13.st, d2th23.dth12.st.th13.st,
                                                d2th12.dth12.st.th23.st,  d2th13.dth12.st.th23.st, d2th23.dth12.st.th23.st ), 3 , 3)
    
    d2theta13.theta13.st <- matrix( c(          d2th12.dth13.st.th12.st,  d2th13.dth13.st.th12.st, d2th23.dth13.st.th12.st,
                                                d2th12.dth13.st.th13.st,  d2th13.dth13.st.th13.st, d2th23.dth13.st.th13.st,
                                                d2th12.dth13.st.th23.st,  d2th13.dth13.st.th23.st, d2th23.dth13.st.th23.st ), 3 , 3)
    
    d2theta23.theta23.st <- matrix( c(          d2th12.dth23.st.th12.st,  d2th13.dth23.st.th12.st, d2th23.dth23.st.th12.st,
                                                d2th12.dth23.st.th13.st,  d2th13.dth23.st.th13.st, d2th23.dth23.st.th13.st,
                                                d2th12.dth23.st.th23.st,  d2th13.dth23.st.th23.st, d2th23.dth23.st.th23.st ), 3 , 3)
    
    
    
    d2l.theta <- matrix( c(          sum(VC$weights[VC$inde1]*d2l.dtheta12.theta12),  sum(VC$weights[VC$inde2]*d2l.dtheta12.theta13), sum(VC$weights[VC$inde2]*d2l.dtheta12.theta23),
                                     sum(VC$weights[VC$inde2]*d2l.dtheta12.theta13),  sum(VC$weights[VC$inde2]*d2l.dtheta13.theta13), sum(VC$weights[VC$inde2]*d2l.dtheta13.theta23),
                                     sum(VC$weights[VC$inde2]*d2l.dtheta12.theta23),  sum(VC$weights[VC$inde2]*d2l.dtheta13.theta23), sum(VC$weights[VC$inde2]*d2l.dtheta23.theta23) ), 3 , 3)
    
    
    
    dl.dtheta <- matrix(0,length(VC$inde1),3)
    dl.dtheta[VC$inde1, 1] <- LgTRI$dl.dtheta12
    dl.dtheta[VC$inde2, 2]   <- LgTRI$dl.dtheta13
    dl.dtheta[VC$inde2, 3]   <- LgTRI$dl.dtheta23
    
    
    dl.dtheta.d2theta.dtheta.st <- matrix(c(colSums(dl.dtheta) %*% d2theta12.theta12.st,
                                            colSums(dl.dtheta) %*% d2theta13.theta13.st,
                                            colSums(dl.dtheta) %*% d2theta23.theta23.st), 3, 3)
    
    d2l.dtheta12.st <- ( t(dtheta.theta.st) %*% d2l.theta %*% dtheta.theta.st + dl.dtheta.d2theta.dtheta.st )[1, 1]
    
    d2l.dtheta13.st <- ( t(dtheta.theta.st) %*% d2l.theta %*% dtheta.theta.st + dl.dtheta.d2theta.dtheta.st )[2, 2]
    
    d2l.dtheta23.st <- ( t(dtheta.theta.st) %*% d2l.theta %*% dtheta.theta.st + dl.dtheta.d2theta.dtheta.st )[3, 3]
    
    d2l.dtheta12.st.theta13.st <- ( t(dtheta.theta.st) %*% d2l.theta %*% dtheta.theta.st + dl.dtheta.d2theta.dtheta.st )[1, 2]
    
    d2l.dtheta12.st.theta23.st <- ( t(dtheta.theta.st) %*% d2l.theta %*% dtheta.theta.st + dl.dtheta.d2theta.dtheta.st )[1, 3]
    
    d2l.dtheta13.st.theta23.st <- ( t(dtheta.theta.st) %*% d2l.theta %*% dtheta.theta.st + dl.dtheta.d2theta.dtheta.st )[2, 3]
    
  }
  
  
  h1 <- cbind(d2l.dbe1.be1, d2l.dbe1.be2, d2l.dbe1.be3, d2l.dbe1.theta12.st, d2l.dbe1.theta13.st, d2l.dbe1.theta23.st)
  h2 <- cbind(t(d2l.dbe1.be2), d2l.dbe2.be2, d2l.dbe2.be3, d2l.dbe2.theta12.st, d2l.dbe2.theta13.st, d2l.dbe2.theta23.st)
  h3 <- cbind(t(d2l.dbe1.be3), t(d2l.dbe2.be3), d2l.dbe3.be3, d2l.dbe3.theta12.st, d2l.dbe3.theta13.st, d2l.dbe3.theta23.st)
  h4 <- cbind(t(d2l.dbe1.theta12.st), t(d2l.dbe2.theta12.st), t(d2l.dbe3.theta12.st), d2l.dtheta12.st, d2l.dtheta12.st.theta13.st, d2l.dtheta12.st.theta23.st)
  h5 <- cbind(t(d2l.dbe1.theta13.st), t(d2l.dbe2.theta13.st), t(d2l.dbe3.theta13.st), t(d2l.dtheta12.st.theta13.st), d2l.dtheta13.st, d2l.dtheta13.st.theta23.st)
  h6 <- cbind(t(d2l.dbe1.theta23.st), t(d2l.dbe2.theta23.st), t(d2l.dbe3.theta23.st), t(d2l.dtheta12.st.theta23.st), t(d2l.dtheta13.st.theta23.st), d2l.dtheta23.st)
  
  HTRIVec <- list(h1 = h1, h2 = h2, h3 = h3, h4 = h4, h5 = h5, h6 = h6)
  
  HTRIVec
  
}
