H.tri <- function(respvec, VC, TIn, LgTRI){
  
  #####
  # ! #
  ########################################################################
  ## I replaced TIn$eta1 with TIn$mar1. Same for TIn$eta2 and TIn$eta3  ##                 
  ########################################################################
  
  dst.1 <- dnorm( (TIn$mar2  - TIn$theta12 * TIn$mar1 )/sqrt(1 - TIn$theta12^2) )  
  pst.1 <- pnorm( ( ((TIn$mar3 - TIn$theta13 * TIn$mar1)/sqrt(1 - TIn$theta13^2)) - ((TIn$theta23 - TIn$theta12 * TIn$theta13)/sqrt((1 - TIn$theta12^2) * (1 - TIn$theta13^2))) * ((TIn$mar2  - TIn$theta12 * TIn$mar1)/sqrt(1 - TIn$theta12^2)) )/sqrt(1 - ((TIn$theta23 - TIn$theta12 * TIn$theta13)/sqrt((1 - TIn$theta12^2) * (1 - TIn$theta13^2)))^2))
  pst.1 <- mm(pst.1)
  st.1 <- -TIn$theta12/sqrt(1 - TIn$theta12^2)
  
  dst.2 <- dnorm((TIn$mar3 - TIn$theta13 * TIn$mar1)/sqrt(1 - TIn$theta13^2))
  pst.2 <- pnorm( ( ((TIn$mar2  - TIn$theta12 * TIn$mar1)/sqrt(1 - TIn$theta12^2)) - ((TIn$theta23 - TIn$theta12 * TIn$theta13)/sqrt((1 - TIn$theta12^2) * (1 - TIn$theta13^2))) * ((TIn$mar3 - TIn$theta13 * TIn$mar1)/sqrt(1 - TIn$theta13^2)) )/sqrt(1 - ((TIn$theta23 - TIn$theta12 * TIn$theta13)/sqrt((1 - TIn$theta12^2) * (1 - TIn$theta13^2)))^2))
  pst.2 <- mm(pst.2)
  st.2 <- -TIn$theta13/sqrt(1 - TIn$theta13^2)
  
  dst.3 <- dnorm((TIn$mar1 - TIn$theta12 * TIn$mar2)/sqrt(1 - TIn$theta12^2))
  pst.3 <- pnorm( ( ((TIn$mar3 - TIn$theta23 * TIn$mar2)/sqrt(1 - TIn$theta23^2)) - ((TIn$theta13 - TIn$theta12 * TIn$theta23)/sqrt((1 - TIn$theta12^2) * (1 - TIn$theta23^2))) * ((TIn$mar1 - TIn$theta12 * TIn$mar2)/sqrt(1 - TIn$theta12^2)) )/sqrt(1 - ((TIn$theta13 - TIn$theta12 * TIn$theta23)/sqrt((1 - TIn$theta12^2) * (1 - TIn$theta23^2)))^2))
  pst.3 <- mm(pst.3)
  st.3 <- -TIn$theta12/sqrt(1 - TIn$theta12^2)
  
  dst.4 <- dnorm((TIn$mar3 - TIn$theta23 * TIn$mar2)/sqrt(1 - TIn$theta23^2))
  pst.4 <- pnorm( ( ((TIn$mar1 - TIn$theta12 * TIn$mar2)/sqrt(1 - TIn$theta12^2)) - ((TIn$theta13 - TIn$theta12 * TIn$theta23)/sqrt((1 - TIn$theta12^2) * (1 - TIn$theta23^2))) * ((TIn$mar3 - TIn$theta23 * TIn$mar2)/sqrt(1 - TIn$theta23^2)) )/sqrt(1 - ((TIn$theta13 - TIn$theta12 * TIn$theta23)/sqrt((1 - TIn$theta12^2) * (1 - TIn$theta23^2)))^2))
  pst.4 <- mm(pst.4)
  st.4  <- -TIn$theta23/sqrt(1 - TIn$theta23^2)
  
  dst.5 <- dnorm((TIn$mar1 - TIn$theta13 * TIn$mar3)/sqrt(1 - TIn$theta13^2))
  pst.5 <- pnorm( ( ((TIn$mar2  - TIn$theta23 * TIn$mar3)/sqrt(1 - TIn$theta23^2)) - ((TIn$theta12 - TIn$theta13 * TIn$theta23)/sqrt((1 - TIn$theta13^2) * (1 - TIn$theta23^2))) * ((TIn$mar1 - TIn$theta13 * TIn$mar3)/sqrt(1 - TIn$theta13^2)) )/sqrt(1 - ((TIn$theta12 - TIn$theta13 * TIn$theta23)/sqrt((1 - TIn$theta13^2) * (1 - TIn$theta23^2)))^2))
  pst.5 <- mm(pst.5)
  st.5  <- -TIn$theta13/sqrt(1 - TIn$theta13^2)
  
  dst.6 <- dnorm((TIn$mar2 - TIn$theta23 * TIn$mar3)/sqrt(1 - TIn$theta23^2))
  pst.6 <- pnorm( ( ((TIn$mar1 - TIn$theta13 * TIn$mar3)/sqrt(1 - TIn$theta13^2)) - ((TIn$theta12 - TIn$theta13 * TIn$theta23)/sqrt((1 - TIn$theta13^2) * (1 - TIn$theta23^2))) * ((TIn$mar2 - TIn$theta23 * TIn$mar3)/sqrt(1 - TIn$theta23^2)) )/sqrt(1 - ((TIn$theta12 - TIn$theta13 * TIn$theta23)/sqrt((1 - TIn$theta13^2) * (1 - TIn$theta23^2)))^2))
  pst.6 <- mm(pst.6)
  st.6  <- -TIn$theta23/sqrt(1 - TIn$theta23^2)
  
  ########################################################################
  
  dp.1.11.de1 <- dst.1 * pst.1        * st.1      + dst.2 * pst.2       * st.2
  dp.1.10.de1 <- dst.1 * (1 - pst.1)  * st.1      + dst.2 * pst.2       * ( -st.2 )
  dp.1.00.de1 <- dst.1 * (1 - pst.1)  * ( -st.1 ) + dst.2 * (1 - pst.2) * ( -st.2 )
  dp.1.01.de1 <- dst.1 * pst.1        * ( -st.1 ) + dst.2 * (1 - pst.2) * st.2 
  
  dp.2.11.de2 <- dst.3 * pst.3        * st.3      + dst.4 * pst.4       * st.4
  dp.2.10.de2 <- dst.3 * (1 - pst.3)  * st.3      + dst.4 * pst.4       * ( -st.4 )
  dp.2.00.de2 <- dst.3 * (1 - pst.3)  * ( -st.3 ) + dst.4 * (1 - pst.4) * ( -st.4 )
  dp.2.01.de2 <- dst.3 * pst.3        * ( -st.3 ) + dst.4 * (1 - pst.4) * st.4 
  
  dp.3.11.de3 <- dst.5 * pst.5        * st.5      + dst.6 * pst.6       * st.6
  dp.3.10.de3 <- dst.5 * (1 - pst.5)  * st.5      + dst.6 * pst.6       * ( -st.6 )
  dp.3.00.de3 <- dst.5 * (1 - pst.5)  * ( -st.5 ) + dst.6 * (1 - pst.6) * ( -st.6 )
  dp.3.01.de3 <- dst.5 * pst.5        * ( -st.5 ) + dst.6 * (1 - pst.6) * st.6 
  
  #####
  # ! #
  ###############################
  ## The next 6 lines are new  ##                 
  ###############################
  
  der2p.dereta1 <- probm(TIn$eta1, VC$margins[1], only.pr = FALSE)$der2p.dereta
  der2p.dereta2 <- probm(TIn$eta2, VC$margins[2], only.pr = FALSE)$der2p.dereta
  der2p.dereta3 <- probm(TIn$eta3, VC$margins[3], only.pr = FALSE)$der2p.dereta
  
  
  d2F1.de1 <- (TIn$mar1 * LgTRI$dmar1^2)/LgTRI$d.1^2 + der2p.dereta1/LgTRI$d.1
  d2F2.de2 <- (TIn$mar2 * LgTRI$dmar2^2)/LgTRI$d.2^2 + der2p.dereta2/LgTRI$d.2
  d2F3.de3 <- (TIn$mar3 * LgTRI$dmar3^2)/LgTRI$d.3^2 + der2p.dereta3/LgTRI$d.3
  
  
  ###################################################################
  
  d2l.dF1.F1 <- respvec$y1.y2.y3  * ( -1/TIn$p111^2 * (LgTRI$d.1 * LgTRI$p.1.11)^2 + 1/TIn$p111 * ( -TIn$mar1 * LgTRI$d.1 * LgTRI$p.1.11 + LgTRI$d.1 * dp.1.11.de1) ) + 
    respvec$y1.y2.cy3   * ( -1/TIn$p110^2 * (LgTRI$d.1 * LgTRI$p.1.10)^2 + 1/TIn$p110 * ( -TIn$mar1 * LgTRI$d.1 * LgTRI$p.1.10 + LgTRI$d.1 * dp.1.10.de1) ) -
    respvec$cy1.y2.y3   * (  1/TIn$p011^2 * (LgTRI$d.1 * LgTRI$p.1.11)^2 + 1/TIn$p011 * ( -TIn$mar1 * LgTRI$d.1 * LgTRI$p.1.11 + LgTRI$d.1 * dp.1.11.de1) ) -
    respvec$cy1.y2.cy3  * (  1/TIn$p010^2 * (LgTRI$d.1 * LgTRI$p.1.10)^2 + 1/TIn$p010 * ( -TIn$mar1 * LgTRI$d.1 * LgTRI$p.1.10 + LgTRI$d.1 * dp.1.10.de1) ) -
    respvec$cy1.cy2.cy3 * (  1/TIn$p000^2 * (LgTRI$d.1 * LgTRI$p.1.00)^2 + 1/TIn$p000 * ( -TIn$mar1 * LgTRI$d.1 * LgTRI$p.1.00 + LgTRI$d.1 * dp.1.00.de1) ) -
    respvec$cy1.cy2.y3  * (  1/TIn$p001^2 * (LgTRI$d.1 * LgTRI$p.1.01)^2 + 1/TIn$p001 * ( -TIn$mar1 * LgTRI$d.1 * LgTRI$p.1.01 + LgTRI$d.1 * dp.1.01.de1) ) +
    respvec$y1.cy2.cy3  * ( -1/TIn$p100^2 * (LgTRI$d.1 * LgTRI$p.1.00)^2 + 1/TIn$p100 * ( -TIn$mar1 * LgTRI$d.1 * LgTRI$p.1.00 + LgTRI$d.1 * dp.1.00.de1) ) +
    respvec$y1.cy2.y3   * ( -1/TIn$p101^2 * (LgTRI$d.1 * LgTRI$p.1.01)^2 + 1/TIn$p101 * ( -TIn$mar1 * LgTRI$d.1 * LgTRI$p.1.01 + LgTRI$d.1 * dp.1.01.de1) )
  
  #####
  # ! #
  ######################
  ## This bit is new  ##                 
  ######################
  
  d2l.de1.e1 <- d2l.dF1.F1 * LgTRI$dF1.de1^2 + LgTRI$dl.dF1 * d2F1.de1
  
  ################################################################
  
  d2l.dF2.F2 <-  respvec$y1.y2.y3 * ( -1/TIn$p111^2 * (LgTRI$d.2 * LgTRI$p.2.11)^2 + 1/TIn$p111 * ( -TIn$mar2   * LgTRI$d.2 * LgTRI$p.2.11 + LgTRI$d.2 * dp.2.11.de2) ) + 
    respvec$y1.y2.cy3 * ( -1/TIn$p110^2 * (LgTRI$d.2 * LgTRI$p.2.10)^2 + 1/TIn$p110 * ( -TIn$mar2  * LgTRI$d.2 * LgTRI$p.2.10 + LgTRI$d.2 * dp.2.10.de2) ) +
    respvec$cy1.y2.y3 * ( -1/TIn$p011^2 * (LgTRI$d.2 * LgTRI$p.2.01)^2 + 1/TIn$p011 * ( -TIn$mar2  * LgTRI$d.2 * LgTRI$p.2.01 + LgTRI$d.2 * dp.2.01.de2) ) +
    respvec$cy1.y2.cy3 * ( -1/TIn$p010^2 * (LgTRI$d.2 * LgTRI$p.2.00)^2 + 1/TIn$p010 * ( -TIn$mar2  * LgTRI$d.2 * LgTRI$p.2.00 + LgTRI$d.2 * dp.2.00.de2) ) -
    respvec$cy1.cy2.cy3 * (  1/TIn$p000^2 * (LgTRI$d.2 * LgTRI$p.2.00)^2 + 1/TIn$p000 * ( -TIn$mar2  * LgTRI$d.2 * LgTRI$p.2.00 + LgTRI$d.2 * dp.2.00.de2) ) -
    respvec$cy1.cy2.y3 * (  1/TIn$p001^2 * (LgTRI$d.2 * LgTRI$p.2.01)^2 + 1/TIn$p001 * ( -TIn$mar2  * LgTRI$d.2 * LgTRI$p.2.01 + LgTRI$d.2 * dp.2.01.de2) ) -
    respvec$y1.cy2.cy3 * (  1/TIn$p100^2 * (LgTRI$d.2 * LgTRI$p.2.10)^2 + 1/TIn$p100 * ( -TIn$mar2  * LgTRI$d.2 * LgTRI$p.2.10 + LgTRI$d.2 * dp.2.10.de2) ) -
    respvec$y1.cy2.y3 * (  1/TIn$p101^2 * (LgTRI$d.2 * LgTRI$p.2.11)^2 + 1/TIn$p101 * ( -TIn$mar2  * LgTRI$d.2 * LgTRI$p.2.11 + LgTRI$d.2 * dp.2.11.de2) )
  
  #####
  # ! #
  ######################
  ## This bit is new  ##                 
  ######################
  
  d2l.de2.e2 <- d2l.dF2.F2 * LgTRI$dF2.de2^2 + LgTRI$dl.dF2 * d2F2.de2
  
  ################################################################
  
  d2l.dF3.F3 <-  respvec$y1.y2.y3 * ( -1/TIn$p111^2 * (LgTRI$d.3 * LgTRI$p.3.11)^2 + 1/TIn$p111 * ( -TIn$mar3 * LgTRI$d.3 * LgTRI$p.3.11 + LgTRI$d.3 * dp.3.11.de3) ) - 
    respvec$y1.y2.cy3 * (  1/TIn$p110^2 * (LgTRI$d.3 * LgTRI$p.3.11)^2 + 1/TIn$p110 * ( -TIn$mar3 * LgTRI$d.3 * LgTRI$p.3.11 + LgTRI$d.3 * dp.3.11.de3) ) +
    respvec$cy1.y2.y3 * ( -1/TIn$p011^2 * (LgTRI$d.3 * LgTRI$p.3.01)^2 + 1/TIn$p011 * ( -TIn$mar3 * LgTRI$d.3 * LgTRI$p.3.01 + LgTRI$d.3 * dp.3.01.de3) ) -
    respvec$cy1.y2.cy3 * (  1/TIn$p010^2 * (LgTRI$d.3 * LgTRI$p.3.01)^2 + 1/TIn$p010 * ( -TIn$mar3 * LgTRI$d.3 * LgTRI$p.3.01 + LgTRI$d.3 * dp.3.01.de3) ) -
    respvec$cy1.cy2.cy3 * (  1/TIn$p000^2 * (LgTRI$d.3 * LgTRI$p.3.00)^2 + 1/TIn$p000 * ( -TIn$mar3 * LgTRI$d.3 * LgTRI$p.3.00 + LgTRI$d.3 * dp.3.00.de3) ) +
    respvec$cy1.cy2.y3 * ( -1/TIn$p001^2 * (LgTRI$d.3 * LgTRI$p.3.00)^2 + 1/TIn$p001 * ( -TIn$mar3 * LgTRI$d.3 * LgTRI$p.3.00 + LgTRI$d.3 * dp.3.00.de3) ) -
    respvec$y1.cy2.cy3 * (  1/TIn$p100^2 * (LgTRI$d.3 * LgTRI$p.3.10)^2 + 1/TIn$p100 * ( -TIn$mar3 * LgTRI$d.3 * LgTRI$p.3.10 + LgTRI$d.3 * dp.3.10.de3) ) +
    respvec$y1.cy2.y3 * ( -1/TIn$p101^2 * (LgTRI$d.3 * LgTRI$p.3.10)^2 + 1/TIn$p101 * ( -TIn$mar3 * LgTRI$d.3 * LgTRI$p.3.10 + LgTRI$d.3 * dp.3.10.de3) )
  
  #####
  # ! #
  ######################
  ## This bit is new  ##                 
  ######################
  
  d2l.de3.e3 <- d2l.dF3.F3 * LgTRI$dF3.de3^2 + LgTRI$dl.dF3 * d2F3.de3
  
  ################################################################
  
  
  #####
  # ! #
  ########################################################################
  ## I replaced TIn$eta1 with TIn$mar1. Same for TIn$eta2 and TIn$eta3  ##                 
  ########################################################################
  
  mean11 <- TIn$theta23 * TIn$mar2 + ((TIn$theta13 - TIn$theta12 * TIn$theta23) * (TIn$mar1 - TIn$theta12 * TIn$mar2))/(1 - TIn$theta12^2)
  mean22 <- TIn$theta23 * TIn$mar3 + ((TIn$theta12 - TIn$theta13 * TIn$theta23) * (TIn$mar1 - TIn$theta13 * TIn$mar3))/(1 - TIn$theta13^2)
  mean33 <- TIn$theta13 * TIn$mar3 + ((TIn$theta12 - TIn$theta13 * TIn$theta23) * (TIn$mar2 - TIn$theta23 * TIn$mar3))/(1 - TIn$theta23^2)
  
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
  
  p1.1 <- mm( pnorm((TIn$mar3 - mean11)/sd11) )
  p2.2 <- mm( pnorm((TIn$mar2 - mean22)/sd22) )
  p3.3 <- mm( pnorm((TIn$mar1 - mean33)/sd33) )
  
  ###########################################################################
  
  p1.1.c <- mm(1 - p1.1)
  p2.2.c <- mm(1 - p2.2)
  p3.3.c <- mm(1 - p3.3)
  
  #####
  # ! #
  ########################################################################
  ## I replaced TIn$eta1 with TIn$mar1. Same for TIn$eta2 and TIn$eta3  ##                 
  ########################################################################
  
  d.1.1  <- dnorm(TIn$mar1, mean = TIn$theta12 * TIn$mar2, sd = sqrt(1 - TIn$theta12^2))
  d.1.2  <- dnorm(TIn$mar1, mean = TIn$theta13 * TIn$mar3, sd = sqrt(1 - TIn$theta13^2))
  d.1.3  <- dnorm(TIn$mar2, mean = TIn$theta23 * TIn$mar3, sd = sqrt(1 - TIn$theta23^2))
  
  ###########################################################################
  
  d2l.dF1.F2 <-  respvec$y1.y2.y3 * ( -1/TIn$p111^2 * LgTRI$d.2 * LgTRI$p.2.11 * LgTRI$d.1 * LgTRI$p.1.11 + 1/TIn$p111 * LgTRI$d.2 * d.1.1 * p1.1   ) +
    respvec$y1.y2.cy3 * ( -1/TIn$p110^2 * LgTRI$d.2 * LgTRI$p.2.10 * LgTRI$d.1 * LgTRI$p.1.10 + 1/TIn$p110 * LgTRI$d.2 * d.1.1 * p1.1.c ) -
    respvec$cy1.y2.y3 * ( -1/TIn$p011^2 * LgTRI$d.2 * LgTRI$p.2.01 * LgTRI$d.1 * LgTRI$p.1.11 + 1/TIn$p011 * LgTRI$d.2 * d.1.1 * p1.1   ) -
    respvec$cy1.y2.cy3 * ( -1/TIn$p010^2 * LgTRI$d.2 * LgTRI$p.2.00 * LgTRI$d.1 * LgTRI$p.1.10 + 1/TIn$p010 * LgTRI$d.2 * d.1.1 * p1.1.c ) +
    respvec$cy1.cy2.cy3 * ( -1/TIn$p000^2 * LgTRI$d.2 * LgTRI$p.2.00 * LgTRI$d.1 * LgTRI$p.1.00 + 1/TIn$p000 * LgTRI$d.2 * d.1.1 * p1.1.c ) +
    respvec$cy1.cy2.y3 * ( -1/TIn$p001^2 * LgTRI$d.2 * LgTRI$p.2.01 * LgTRI$d.1 * LgTRI$p.1.01 + 1/TIn$p001 * LgTRI$d.2 * d.1.1 * p1.1   ) -
    respvec$y1.cy2.cy3 * ( -1/TIn$p100^2 * LgTRI$d.2 * LgTRI$p.2.10 * LgTRI$d.1 * LgTRI$p.1.00 + 1/TIn$p100 * LgTRI$d.2 * d.1.1 * p1.1.c ) -
    respvec$y1.cy2.y3 * ( -1/TIn$p101^2 * LgTRI$d.2 * LgTRI$p.2.11 * LgTRI$d.1 * LgTRI$p.1.01 + 1/TIn$p101 * LgTRI$d.2 * d.1.1 * p1.1   ) 
  
  #####
  # ! #
  ######################
  ## This bit is new  ##                 
  ######################
  
  d2l.de1.e2 <- d2l.dF1.F2 * LgTRI$dF1.de1 * LgTRI$dF2.de2  
  
  ###########################################################
  
  d2l.dF1.F3 <-  respvec$y1.y2.y3 * ( -1/TIn$p111^2 * LgTRI$d.3 * LgTRI$p.3.11 * LgTRI$d.1 * LgTRI$p.1.11 + 1/TIn$p111 * LgTRI$d.3 * d.1.2 * p2.2   ) -
    respvec$y1.y2.cy3 * ( -1/TIn$p110^2 * LgTRI$d.3 * LgTRI$p.3.11 * LgTRI$d.1 * LgTRI$p.1.10 + 1/TIn$p110 * LgTRI$d.3 * d.1.2 * p2.2   ) -
    respvec$cy1.y2.y3 * ( -1/TIn$p011^2 * LgTRI$d.3 * LgTRI$p.3.01 * LgTRI$d.1 * LgTRI$p.1.11 + 1/TIn$p011 * LgTRI$d.3 * d.1.2 * p2.2   ) +
    respvec$cy1.y2.cy3 * ( -1/TIn$p010^2 * LgTRI$d.3 * LgTRI$p.3.01 * LgTRI$d.1 * LgTRI$p.1.10 + 1/TIn$p010 * LgTRI$d.3 * d.1.2 * p2.2   ) +
    respvec$cy1.cy2.cy3 * ( -1/TIn$p000^2 * LgTRI$d.3 * LgTRI$p.3.00 * LgTRI$d.1 * LgTRI$p.1.00 + 1/TIn$p000 * LgTRI$d.3 * d.1.2 * p2.2.c ) -
    respvec$cy1.cy2.y3 * ( -1/TIn$p001^2 * LgTRI$d.3 * LgTRI$p.3.00 * LgTRI$d.1 * LgTRI$p.1.01 + 1/TIn$p001 * LgTRI$d.3 * d.1.2 * p2.2.c ) -
    respvec$y1.cy2.cy3 * ( -1/TIn$p100^2 * LgTRI$d.3 * LgTRI$p.3.10 * LgTRI$d.1 * LgTRI$p.1.00 + 1/TIn$p100 * LgTRI$d.3 * d.1.2 * p2.2.c ) +
    respvec$y1.cy2.y3 * ( -1/TIn$p101^2 * LgTRI$d.3 * LgTRI$p.3.10 * LgTRI$d.1 * LgTRI$p.1.01 + 1/TIn$p101 * LgTRI$d.3 * d.1.2 * p2.2.c ) 
  
  #####
  # ! #
  ######################
  ## This bit is new  ##                 
  ######################
  
  d2l.de1.e3 <- d2l.dF1.F3 * LgTRI$dF1.de1 * LgTRI$dF3.de3  
  
  ###########################################################
  
  d2l.dF2.F3 <-  respvec$y1.y2.y3 * ( -1/TIn$p111^2 * LgTRI$d.3 * LgTRI$p.3.11 * LgTRI$d.2 * LgTRI$p.2.11 + 1/TIn$p111 * LgTRI$d.3 * d.1.3 * p3.3  ) -
    respvec$y1.y2.cy3 * ( -1/TIn$p110^2 * LgTRI$d.3 * LgTRI$p.3.11 * LgTRI$d.2 * LgTRI$p.2.10 + 1/TIn$p110 * LgTRI$d.3 * d.1.3 * p3.3  ) +
    respvec$cy1.y2.y3 * ( -1/TIn$p011^2 * LgTRI$d.3 * LgTRI$p.3.01 * LgTRI$d.2 * LgTRI$p.2.01 + 1/TIn$p011 * LgTRI$d.3 * d.1.3 * p3.3.c ) -
    respvec$cy1.y2.cy3 * ( -1/TIn$p010^2 * LgTRI$d.3 * LgTRI$p.3.01 * LgTRI$d.2 * LgTRI$p.2.00 + 1/TIn$p010 * LgTRI$d.3 * d.1.3 * p3.3.c ) +
    respvec$cy1.cy2.cy3 * ( -1/TIn$p000^2 * LgTRI$d.3 * LgTRI$p.3.00 * LgTRI$d.2 * LgTRI$p.2.00 + 1/TIn$p000 * LgTRI$d.3 * d.1.3 * p3.3.c ) -
    respvec$cy1.cy2.y3 * ( -1/TIn$p001^2 * LgTRI$d.3 * LgTRI$p.3.00 * LgTRI$d.2 * LgTRI$p.2.01 + 1/TIn$p001 * LgTRI$d.3 * d.1.3 * p3.3.c ) +
    respvec$y1.cy2.cy3 * ( -1/TIn$p100^2 * LgTRI$d.3 * LgTRI$p.3.10 * LgTRI$d.2 * LgTRI$p.2.10 + 1/TIn$p100 * LgTRI$d.3 * d.1.3 * p3.3  ) -
    respvec$y1.cy2.y3 * ( -1/TIn$p101^2 * LgTRI$d.3 * LgTRI$p.3.10 * LgTRI$d.2 * LgTRI$p.2.11 + 1/TIn$p101 * LgTRI$d.3 * d.1.3 * p3.3   ) 
  
  #####
  # ! #
  ######################
  ## This bit is new  ##                 
  ######################
  
  d2l.de2.e3 <- d2l.dF2.F3 * LgTRI$dF2.de2 * LgTRI$dF3.de3  
  
  ###########################################################
  
  
  #####
  # ! #
  ########################################################################
  ## I replaced TIn$eta1 with TIn$mar1. Same for TIn$eta2 and TIn$eta3  ##                 
  ########################################################################
  
  d12 <- dnorm( (TIn$mar3 - LgTRI$mean.12)/LgTRI$sd.12 )
  d13 <- dnorm( (TIn$mar2 - LgTRI$mean.13)/LgTRI$sd.13 )
  d23 <- dnorm( (TIn$mar1 - LgTRI$mean.23)/LgTRI$sd.23 )
  
  ###########################################################
  
  d2l.dF1.theta12 <- respvec$y1.y2.y3 * ( -1/TIn$p111^2 * LgTRI$d.1 * LgTRI$p.1.11 * LgTRI$d11.12 * LgTRI$p12.g   + 1/TIn$p111 * ( LgTRI$d11.12 * (TIn$theta12*TIn$mar2  - TIn$mar1)/(1 - TIn$theta12^2) * LgTRI$p12.g   - LgTRI$d11.12 * d12/LgTRI$sd.12 * ((TIn$theta13 - TIn$theta12 * TIn$theta23)/(1 - TIn$theta12^2)) )) +
    respvec$y1.y2.cy3 * ( -1/TIn$p110^2 * LgTRI$d.1 * LgTRI$p.1.10 * LgTRI$d11.12 * LgTRI$p12.g.c + 1/TIn$p110 * ( LgTRI$d11.12 * (TIn$theta12*TIn$mar2  - TIn$mar1)/(1 - TIn$theta12^2) * LgTRI$p12.g.c + LgTRI$d11.12 * d12/LgTRI$sd.12 * ((TIn$theta13 - TIn$theta12 * TIn$theta23)/(1 - TIn$theta12^2)) )) -
    respvec$cy1.y2.y3 * (  1/TIn$p011^2 * LgTRI$d.1 * LgTRI$p.1.11 * LgTRI$d11.12 * LgTRI$p12.g   + 1/TIn$p011 * ( LgTRI$d11.12 * (TIn$theta12*TIn$mar2  - TIn$mar1)/(1 - TIn$theta12^2) * LgTRI$p12.g   - LgTRI$d11.12 * d12/LgTRI$sd.12 * ((TIn$theta13 - TIn$theta12 * TIn$theta23)/(1 - TIn$theta12^2)) )) -
    respvec$cy1.y2.cy3 * (  1/TIn$p010^2 * LgTRI$d.1 * LgTRI$p.1.10 * LgTRI$d11.12 * LgTRI$p12.g.c + 1/TIn$p010 * ( LgTRI$d11.12 * (TIn$theta12*TIn$mar2  - TIn$mar1)/(1 - TIn$theta12^2) * LgTRI$p12.g.c + LgTRI$d11.12 * d12/LgTRI$sd.12 * ((TIn$theta13 - TIn$theta12 * TIn$theta23)/(1 - TIn$theta12^2)) )) +
    respvec$cy1.cy2.cy3 * (  1/TIn$p000^2 * LgTRI$d.1 * LgTRI$p.1.00 * LgTRI$d11.12 * LgTRI$p12.g.c + 1/TIn$p000 * ( LgTRI$d11.12 * (TIn$theta12*TIn$mar2  - TIn$mar1)/(1 - TIn$theta12^2) * LgTRI$p12.g.c + LgTRI$d11.12 * d12/LgTRI$sd.12 * ((TIn$theta13 - TIn$theta12 * TIn$theta23)/(1 - TIn$theta12^2)) )) +
    respvec$cy1.cy2.y3 * (  1/TIn$p001^2 * LgTRI$d.1 * LgTRI$p.1.01 * LgTRI$d11.12 * LgTRI$p12.g   + 1/TIn$p001 * ( LgTRI$d11.12 * (TIn$theta12*TIn$mar2  - TIn$mar1)/(1 - TIn$theta12^2) * LgTRI$p12.g   - LgTRI$d11.12 * d12/LgTRI$sd.12 * ((TIn$theta13 - TIn$theta12 * TIn$theta23)/(1 - TIn$theta12^2)) )) -
    respvec$y1.cy2.cy3 * ( -1/TIn$p100^2 * LgTRI$d.1 * LgTRI$p.1.00 * LgTRI$d11.12 * LgTRI$p12.g.c + 1/TIn$p100 * ( LgTRI$d11.12 * (TIn$theta12*TIn$mar2  - TIn$mar1)/(1 - TIn$theta12^2) * LgTRI$p12.g.c + LgTRI$d11.12 * d12/LgTRI$sd.12 * ((TIn$theta13 - TIn$theta12 * TIn$theta23)/(1 - TIn$theta12^2)) )) -
    respvec$y1.cy2.y3 * ( -1/TIn$p101^2 * LgTRI$d.1 * LgTRI$p.1.01 * LgTRI$d11.12 * LgTRI$p12.g   + 1/TIn$p101 * ( LgTRI$d11.12 * (TIn$theta12*TIn$mar2  - TIn$mar1)/(1 - TIn$theta12^2) * LgTRI$p12.g   - LgTRI$d11.12 * d12/LgTRI$sd.12 * ((TIn$theta13 - TIn$theta12 * TIn$theta23)/(1 - TIn$theta12^2)) )) 
  
  #####
  # ! #
  ######################
  ## This bit is new  ##                 
  ######################
  
  d2l.de1.theta12 <- d2l.dF1.theta12 * LgTRI$dF1.de1
  
  ###########################################################
  
  d2l.dF1.theta13 <- respvec$y1.y2.y3 * ( -1/TIn$p111^2 * LgTRI$d.1 * LgTRI$p.1.11 * LgTRI$d11.13 * LgTRI$p13.g   + 1/TIn$p111 * ( LgTRI$d11.13 * (TIn$theta13 * TIn$mar3 - TIn$mar1)/(1 - TIn$theta13^2) * LgTRI$p13.g   - LgTRI$d11.13 * d13/LgTRI$sd.13 * ((TIn$theta12 - TIn$theta13 * TIn$theta23)/(1 - TIn$theta13^2)) )) -
    respvec$y1.y2.cy3   * ( -1/TIn$p110^2 * LgTRI$d.1 * LgTRI$p.1.10 * LgTRI$d11.13 * LgTRI$p13.g   + 1/TIn$p110 * ( LgTRI$d11.13 * (TIn$theta13 * TIn$mar3 - TIn$mar1)/(1 - TIn$theta13^2) * LgTRI$p13.g   - LgTRI$d11.13 * d13/LgTRI$sd.13 * ((TIn$theta12 - TIn$theta13 * TIn$theta23)/(1 - TIn$theta13^2)) )) -
    respvec$cy1.y2.y3   * (  1/TIn$p011^2 * LgTRI$d.1 * LgTRI$p.1.11 * LgTRI$d11.13 * LgTRI$p13.g   + 1/TIn$p011 * ( LgTRI$d11.13 * (TIn$theta13 * TIn$mar3 - TIn$mar1)/(1 - TIn$theta13^2) * LgTRI$p13.g   - LgTRI$d11.13 * d13/LgTRI$sd.13 * ((TIn$theta12 - TIn$theta13 * TIn$theta23)/(1 - TIn$theta13^2)) )) +
    respvec$cy1.y2.cy3  * (  1/TIn$p010^2 * LgTRI$d.1 * LgTRI$p.1.10 * LgTRI$d11.13 * LgTRI$p13.g   + 1/TIn$p010 * ( LgTRI$d11.13 * (TIn$theta13 * TIn$mar3 - TIn$mar1)/(1 - TIn$theta13^2) * LgTRI$p13.g   - LgTRI$d11.13 * d13/LgTRI$sd.13 * ((TIn$theta12 - TIn$theta13 * TIn$theta23)/(1 - TIn$theta13^2)) )) +
    respvec$cy1.cy2.cy3 * (  1/TIn$p000^2 * LgTRI$d.1 * LgTRI$p.1.00 * LgTRI$d11.13 * LgTRI$p13.g.c + 1/TIn$p000 * ( LgTRI$d11.13 * (TIn$theta13 * TIn$mar3 - TIn$mar1)/(1 - TIn$theta13^2) * LgTRI$p13.g.c + LgTRI$d11.13 * d13/LgTRI$sd.13 * ((TIn$theta12 - TIn$theta13 * TIn$theta23)/(1 - TIn$theta13^2)) )) -
    respvec$cy1.cy2.y3  * (  1/TIn$p001^2 * LgTRI$d.1 * LgTRI$p.1.01 * LgTRI$d11.13 * LgTRI$p13.g.c + 1/TIn$p001 * ( LgTRI$d11.13 * (TIn$theta13 * TIn$mar3 - TIn$mar1)/(1 - TIn$theta13^2) * LgTRI$p13.g.c + LgTRI$d11.13 * d13/LgTRI$sd.13 * ((TIn$theta12 - TIn$theta13 * TIn$theta23)/(1 - TIn$theta13^2)) )) -
    respvec$y1.cy2.cy3  * ( -1/TIn$p100^2 * LgTRI$d.1 * LgTRI$p.1.00 * LgTRI$d11.13 * LgTRI$p13.g.c + 1/TIn$p100 * ( LgTRI$d11.13 * (TIn$theta13 * TIn$mar3 - TIn$mar1)/(1 - TIn$theta13^2) * LgTRI$p13.g.c + LgTRI$d11.13 * d13/LgTRI$sd.13 * ((TIn$theta12 - TIn$theta13 * TIn$theta23)/(1 - TIn$theta13^2)) )) +
    respvec$y1.cy2.y3   * ( -1/TIn$p101^2 * LgTRI$d.1 * LgTRI$p.1.01 * LgTRI$d11.13 * LgTRI$p13.g.c + 1/TIn$p101 * ( LgTRI$d11.13 * (TIn$theta13 * TIn$mar3 - TIn$mar1)/(1 - TIn$theta13^2) * LgTRI$p13.g.c + LgTRI$d11.13 * d13/LgTRI$sd.13 * ((TIn$theta12 - TIn$theta13 * TIn$theta23)/(1 - TIn$theta13^2)) )) 
  
  #####
  # ! #
  ######################
  ## This bit is new  ##                 
  ######################
  
  d2l.de1.theta13 <- d2l.dF1.theta13 * LgTRI$dF1.de1
  
  ###########################################################
  
  d2l.dF1.theta23 <- respvec$y1.y2.y3 * ( -1/TIn$p111^2 * LgTRI$d.1 * LgTRI$p.1.11 * LgTRI$d11.23 * LgTRI$p23.g   + 1/TIn$p111 * LgTRI$d11.23 * d23/LgTRI$sd.23) -
    respvec$y1.y2.cy3  * ( -1/TIn$p110^2 * LgTRI$d.1 * LgTRI$p.1.10 * LgTRI$d11.23 * LgTRI$p23.g   + 1/TIn$p110 * LgTRI$d11.23 * d23/LgTRI$sd.23 ) -
    respvec$cy1.y2.y3 * ( -1/TIn$p011^2 * LgTRI$d.1 * LgTRI$p.1.11 * LgTRI$d11.23 * LgTRI$p23.g.c + 1/TIn$p011 * LgTRI$d11.23 * d23/LgTRI$sd.23 ) +
    respvec$cy1.y2.cy3 * ( -1/TIn$p010^2 * LgTRI$d.1 * LgTRI$p.1.10 * LgTRI$d11.23 * LgTRI$p23.g.c + 1/TIn$p010 * LgTRI$d11.23 * d23/LgTRI$sd.23 ) -
    respvec$cy1.cy2.cy3 * ( -1/TIn$p000^2 * LgTRI$d.1 * LgTRI$p.1.00 * LgTRI$d11.23 * LgTRI$p23.g.c + 1/TIn$p000 * LgTRI$d11.23 * d23/LgTRI$sd.23 ) +
    respvec$cy1.cy2.y3 * ( -1/TIn$p001^2 * LgTRI$d.1 * LgTRI$p.1.01 * LgTRI$d11.23 * LgTRI$p23.g.c + 1/TIn$p001 * LgTRI$d11.23 * d23/LgTRI$sd.23 ) +
    respvec$y1.cy2.cy3 * ( -1/TIn$p100^2 * LgTRI$d.1 * LgTRI$p.1.00 * LgTRI$d11.23 * LgTRI$p23.g   + 1/TIn$p100 * LgTRI$d11.23 * d23/LgTRI$sd.23 ) -
    respvec$y1.cy2.y3 * ( -1/TIn$p101^2 * LgTRI$d.1 * LgTRI$p.1.01 * LgTRI$d11.23 * LgTRI$p23.g   + 1/TIn$p101 * LgTRI$d11.23 * d23/LgTRI$sd.23 ) 
  
  #####
  # ! #
  ######################
  ## This bit is new  ##                 
  ######################
  
  d2l.de1.theta23 <- d2l.dF1.theta23 * LgTRI$dF1.de1
  
  ###########################################################
  
  d2l.dF2.theta12 <- respvec$y1.y2.y3 * ( -1/TIn$p111^2 * LgTRI$d.2 * LgTRI$p.2.11 * LgTRI$d11.12 * LgTRI$p12.g   + 1/TIn$p111 * (LgTRI$d11.12 * (TIn$theta12*TIn$mar1 - TIn$mar2)/(1 - TIn$theta12^2) * LgTRI$p12.g   - LgTRI$d11.12 * (d12/LgTRI$sd.12) * (TIn$theta23 - TIn$theta12 * TIn$theta13)/(1 - TIn$theta12^2) ) ) +
    respvec$y1.y2.cy3 * ( -1/TIn$p110^2 * LgTRI$d.2 * LgTRI$p.2.10 * LgTRI$d11.12 * LgTRI$p12.g.c + 1/TIn$p110 * (LgTRI$d11.12 * (TIn$theta12 * TIn$mar1 - TIn$mar2)/(1 - TIn$theta12^2) * LgTRI$p12.g.c + LgTRI$d11.12 * (d12/LgTRI$sd.12) * (TIn$theta23 - TIn$theta12 * TIn$theta13)/(1 - TIn$theta12^2) ) ) -
    respvec$cy1.y2.y3 * ( -1/TIn$p011^2 * LgTRI$d.2 * LgTRI$p.2.01 * LgTRI$d11.12 * LgTRI$p12.g   + 1/TIn$p011 * (LgTRI$d11.12 * (TIn$theta12 * TIn$mar1 - TIn$mar2)/(1 - TIn$theta12^2) * LgTRI$p12.g   - LgTRI$d11.12 * (d12/LgTRI$sd.12) * (TIn$theta23 - TIn$theta12 * TIn$theta13)/(1 - TIn$theta12^2) ) ) -
    respvec$cy1.y2.cy3 * ( -1/TIn$p010^2 * LgTRI$d.2 * LgTRI$p.2.00 * LgTRI$d11.12 * LgTRI$p12.g.c + 1/TIn$p010 * (LgTRI$d11.12 * (TIn$theta12 * TIn$mar1 - TIn$mar2)/(1 - TIn$theta12^2) * LgTRI$p12.g.c + LgTRI$d11.12 * (d12/LgTRI$sd.12) * (TIn$theta23 - TIn$theta12 * TIn$theta13)/(1 - TIn$theta12^2) ) ) +
    respvec$cy1.cy2.cy3 * (  1/TIn$p000^2 * LgTRI$d.2 * LgTRI$p.2.00 * LgTRI$d11.12 * LgTRI$p12.g.c + 1/TIn$p000 * (LgTRI$d11.12 * (TIn$theta12 * TIn$mar1 - TIn$mar2)/(1 - TIn$theta12^2) * LgTRI$p12.g.c + LgTRI$d11.12 * (d12/LgTRI$sd.12) * (TIn$theta23 - TIn$theta12 * TIn$theta13)/(1 - TIn$theta12^2) ) ) +
    respvec$cy1.cy2.y3 * (  1/TIn$p001^2 * LgTRI$d.2 * LgTRI$p.2.01 * LgTRI$d11.12 * LgTRI$p12.g   + 1/TIn$p001 * (LgTRI$d11.12 * (TIn$theta12 * TIn$mar1 - TIn$mar2)/(1 - TIn$theta12^2) * LgTRI$p12.g   - LgTRI$d11.12 * (d12/LgTRI$sd.12) * (TIn$theta23 - TIn$theta12 * TIn$theta13)/(1 - TIn$theta12^2) ) ) -
    respvec$y1.cy2.cy3 * (  1/TIn$p100^2 * LgTRI$d.2 * LgTRI$p.2.10 * LgTRI$d11.12 * LgTRI$p12.g.c + 1/TIn$p100 * (LgTRI$d11.12 * (TIn$theta12 * TIn$mar1 - TIn$mar2)/(1 - TIn$theta12^2) * LgTRI$p12.g.c + LgTRI$d11.12 * (d12/LgTRI$sd.12) * (TIn$theta23 - TIn$theta12 * TIn$theta13)/(1 - TIn$theta12^2) ) ) -
    respvec$y1.cy2.y3 * (  1/TIn$p101^2 * LgTRI$d.2 * LgTRI$p.2.11 * LgTRI$d11.12 * LgTRI$p12.g   + 1/TIn$p101 * (LgTRI$d11.12 * (TIn$theta12 * TIn$mar1 - TIn$mar2)/(1 - TIn$theta12^2) * LgTRI$p12.g   - LgTRI$d11.12 * (d12/LgTRI$sd.12) * (TIn$theta23 - TIn$theta12 * TIn$theta13)/(1 - TIn$theta12^2) ) )
  
  #####
  # ! #
  ######################
  ## This bit is new  ##                 
  ######################
  
  d2l.de2.theta12 <- d2l.dF2.theta12 * LgTRI$dF2.de2
  
  ###########################################################
  
  d2l.dF2.theta13 <- respvec$y1.y2.y3 * ( -1/TIn$p111^2 * LgTRI$d.2 * LgTRI$p.2.11 * LgTRI$d11.13 * LgTRI$p13.g   + 1/TIn$p111 * LgTRI$d11.13 * d13/LgTRI$sd.13 ) -
    respvec$y1.y2.cy3 * ( -1/TIn$p110^2 * LgTRI$d.2 * LgTRI$p.2.10 * LgTRI$d11.13 * LgTRI$p13.g   + 1/TIn$p110 * LgTRI$d11.13 * d13/LgTRI$sd.13 )-
    respvec$cy1.y2.y3 * ( -1/TIn$p011^2 * LgTRI$d.2 * LgTRI$p.2.01 * LgTRI$d11.13 * LgTRI$p13.g   + 1/TIn$p011 * LgTRI$d11.13 * d13/LgTRI$sd.13 ) +
    respvec$cy1.y2.cy3 * ( -1/TIn$p010^2 * LgTRI$d.2 * LgTRI$p.2.00 * LgTRI$d11.13 * LgTRI$p13.g   + 1/TIn$p010 * LgTRI$d11.13 * d13/LgTRI$sd.13 ) -
    respvec$cy1.cy2.cy3 * ( -1/TIn$p000^2 * LgTRI$d.2 * LgTRI$p.2.00 * LgTRI$d11.13 * LgTRI$p13.g.c + 1/TIn$p000 * LgTRI$d11.13 * d13/LgTRI$sd.13 ) +
    respvec$cy1.cy2.y3 * ( -1/TIn$p001^2 * LgTRI$d.2 * LgTRI$p.2.01 * LgTRI$d11.13 * LgTRI$p13.g.c + 1/TIn$p001 * LgTRI$d11.13 * d13/LgTRI$sd.13 ) +
    respvec$y1.cy2.cy3 * ( -1/TIn$p100^2 * LgTRI$d.2 * LgTRI$p.2.10 * LgTRI$d11.13 * LgTRI$p13.g.c + 1/TIn$p100 * LgTRI$d11.13 * d13/LgTRI$sd.13 ) -
    respvec$y1.cy2.y3 * ( -1/TIn$p101^2 * LgTRI$d.2 * LgTRI$p.2.11 * LgTRI$d11.13 * LgTRI$p13.g.c + 1/TIn$p101 * LgTRI$d11.13 * d13/LgTRI$sd.13 ) 
  
  #####
  # ! #
  ######################
  ## This bit is new  ##                 
  ######################
  
  d2l.de2.theta13 <- d2l.dF2.theta13 * LgTRI$dF2.de2
  
  ###########################################################
  
  d2l.dF2.theta23 <- respvec$y1.y2.y3 * ( -1/TIn$p111^2 * LgTRI$d.2 * LgTRI$p.2.11 * LgTRI$d11.23 * LgTRI$p23.g   + 1/TIn$p111 * (LgTRI$d11.23 * (TIn$theta23 * TIn$mar3 - TIn$mar2)/(1 - TIn$theta23^2) * LgTRI$p23.g   - LgTRI$d11.23 * d23/LgTRI$sd.23 * (TIn$theta12 - TIn$theta13 * TIn$theta23)/(1 - TIn$theta23^2) ) ) -
    respvec$y1.y2.cy3 * ( -1/TIn$p110^2 * LgTRI$d.2 * LgTRI$p.2.10 * LgTRI$d11.23 * LgTRI$p23.g   + 1/TIn$p110 * (LgTRI$d11.23 * (TIn$theta23 * TIn$mar3 - TIn$mar2)/(1 - TIn$theta23^2) * LgTRI$p23.g   - LgTRI$d11.23 * d23/LgTRI$sd.23 * (TIn$theta12 - TIn$theta13 * TIn$theta23)/(1 - TIn$theta23^2) ) ) +
    respvec$cy1.y2.y3 * ( -1/TIn$p011^2 * LgTRI$d.2 * LgTRI$p.2.01 * LgTRI$d11.23 * LgTRI$p23.g.c + 1/TIn$p011 * (LgTRI$d11.23 * (TIn$theta23 * TIn$mar3 - TIn$mar2)/(1 - TIn$theta23^2) * LgTRI$p23.g.c + LgTRI$d11.23 * d23/LgTRI$sd.23 * (TIn$theta12 - TIn$theta13 * TIn$theta23)/(1 - TIn$theta23^2) ) ) -
    respvec$cy1.y2.cy3 * ( -1/TIn$p010^2 * LgTRI$d.2 * LgTRI$p.2.00 * LgTRI$d11.23 * LgTRI$p23.g.c + 1/TIn$p010 * (LgTRI$d11.23 * (TIn$theta23 * TIn$mar3 - TIn$mar2)/(1 - TIn$theta23^2) * LgTRI$p23.g.c + LgTRI$d11.23 * d23/LgTRI$sd.23 * (TIn$theta12 - TIn$theta13 * TIn$theta23)/(1 - TIn$theta23^2) ) ) +
    respvec$cy1.cy2.cy3 * (  1/TIn$p000^2 * LgTRI$d.2 * LgTRI$p.2.00 * LgTRI$d11.23 * LgTRI$p23.g.c + 1/TIn$p000 * (LgTRI$d11.23 * (TIn$theta23 * TIn$mar3 - TIn$mar2)/(1 - TIn$theta23^2) * LgTRI$p23.g.c + LgTRI$d11.23 * d23/LgTRI$sd.23 * (TIn$theta12 - TIn$theta13 * TIn$theta23)/(1 - TIn$theta23^2) ) ) -
    respvec$cy1.cy2.y3 * (  1/TIn$p001^2 * LgTRI$d.2 * LgTRI$p.2.01 * LgTRI$d11.23 * LgTRI$p23.g.c + 1/TIn$p001 * (LgTRI$d11.23 * (TIn$theta23 * TIn$mar3 - TIn$mar2)/(1 - TIn$theta23^2) * LgTRI$p23.g.c + LgTRI$d11.23 * d23/LgTRI$sd.23 * (TIn$theta12 - TIn$theta13 * TIn$theta23)/(1 - TIn$theta23^2) ) ) +
    respvec$y1.cy2.cy3 * (  1/TIn$p100^2 * LgTRI$d.2 * LgTRI$p.2.10 * LgTRI$d11.23 * LgTRI$p23.g   + 1/TIn$p100 * (LgTRI$d11.23 * (TIn$theta23 * TIn$mar3 - TIn$mar2)/(1 - TIn$theta23^2) * LgTRI$p23.g   - LgTRI$d11.23 * d23/LgTRI$sd.23 * (TIn$theta12 - TIn$theta13 * TIn$theta23)/(1 - TIn$theta23^2) ) ) -
    respvec$y1.cy2.y3 * (  1/TIn$p101^2 * LgTRI$d.2 * LgTRI$p.2.11 * LgTRI$d11.23 * LgTRI$p23.g   + 1/TIn$p101 * (LgTRI$d11.23 * (TIn$theta23 * TIn$mar3 - TIn$mar2)/(1 - TIn$theta23^2) * LgTRI$p23.g   - LgTRI$d11.23 * d23/LgTRI$sd.23 * (TIn$theta12 - TIn$theta13 * TIn$theta23)/(1 - TIn$theta23^2) ) ) 
  
  #####
  # ! #
  ######################
  ## This bit is new  ##                 
  ######################
  
  d2l.de2.theta23 <- d2l.dF2.theta23 * LgTRI$dF2.de2
  
  ###########################################################
  
  
  d2l.dF3.theta12 <- respvec$y1.y2.y3 * ( -1/TIn$p111^2 * LgTRI$d.3 * LgTRI$p.3.11 * LgTRI$d11.12 * LgTRI$p12.g   + 1/TIn$p111 * LgTRI$d11.12 * d12/LgTRI$sd.12 ) -
    respvec$y1.y2.cy3 * ( -1/TIn$p110^2 * LgTRI$d.3 * LgTRI$p.3.11 * LgTRI$d11.12 * LgTRI$p12.g.c + 1/TIn$p110 * LgTRI$d11.12 * d12/LgTRI$sd.12 ) -
    respvec$cy1.y2.y3 * ( -1/TIn$p011^2 * LgTRI$d.3 * LgTRI$p.3.01 * LgTRI$d11.12 * LgTRI$p12.g   + 1/TIn$p011 * LgTRI$d11.12 * d12/LgTRI$sd.12 ) +
    respvec$cy1.y2.cy3 * ( -1/TIn$p010^2 * LgTRI$d.3 * LgTRI$p.3.01 * LgTRI$d11.12 * LgTRI$p12.g.c + 1/TIn$p010 * LgTRI$d11.12 * d12/LgTRI$sd.12 ) -
    respvec$cy1.cy2.cy3 * ( -1/TIn$p000^2 * LgTRI$d.3 * LgTRI$p.3.00 * LgTRI$d11.12 * LgTRI$p12.g.c + 1/TIn$p000 * LgTRI$d11.12 * d12/LgTRI$sd.12 ) +
    respvec$cy1.cy2.y3 * ( -1/TIn$p001^2 * LgTRI$d.3 * LgTRI$p.3.00 * LgTRI$d11.12 * LgTRI$p12.g   + 1/TIn$p001 * LgTRI$d11.12 * d12/LgTRI$sd.12 ) +
    respvec$y1.cy2.cy3 * ( -1/TIn$p100^2 * LgTRI$d.3 * LgTRI$p.3.10 * LgTRI$d11.12 * LgTRI$p12.g.c + 1/TIn$p100 * LgTRI$d11.12 * d12/LgTRI$sd.12 ) -
    respvec$y1.cy2.y3 * ( -1/TIn$p101^2 * LgTRI$d.3 * LgTRI$p.3.10 * LgTRI$d11.12 * LgTRI$p12.g   + 1/TIn$p101 * LgTRI$d11.12 * d12/LgTRI$sd.12 ) 
  
  #####
  # ! #
  ######################
  ## This bit is new  ##                 
  ######################
  
  d2l.de3.theta12 <- d2l.dF3.theta12 * LgTRI$dF3.de3
  
  ###########################################################
  
  d2l.dF3.theta13 <- respvec$y1.y2.y3 * ( -1/TIn$p111^2 * LgTRI$d.3 * LgTRI$p.3.11 * LgTRI$d11.13 * LgTRI$p13.g   + 1/TIn$p111 * (LgTRI$d11.13 * ((TIn$theta13*TIn$mar1 - TIn$mar3)/(1 - TIn$theta13^2)) * LgTRI$p13.g   - LgTRI$d11.13 * (d13/LgTRI$sd.13) * (TIn$theta23 - TIn$theta12 * TIn$theta13)/(1 - TIn$theta13^2) ) ) -
    respvec$y1.y2.cy3 * (  1/TIn$p110^2 * LgTRI$d.3 * LgTRI$p.3.11 * LgTRI$d11.13 * LgTRI$p13.g   + 1/TIn$p110 * (LgTRI$d11.13 * ((TIn$theta13 * TIn$mar1 - TIn$mar3)/(1 - TIn$theta13^2)) * LgTRI$p13.g   - LgTRI$d11.13 * (d13/LgTRI$sd.13) * (TIn$theta23 - TIn$theta12 * TIn$theta13)/(1 - TIn$theta13^2) ) ) -
    respvec$cy1.y2.y3 * ( -1/TIn$p011^2 * LgTRI$d.3 * LgTRI$p.3.01 * LgTRI$d11.13 * LgTRI$p13.g   + 1/TIn$p011 * (LgTRI$d11.13 * ((TIn$theta13 * TIn$mar1 - TIn$mar3)/(1 - TIn$theta13^2)) * LgTRI$p13.g   - LgTRI$d11.13 * (d13/LgTRI$sd.13) * (TIn$theta23 - TIn$theta12 * TIn$theta13)/(1 - TIn$theta13^2) ) ) +
    respvec$cy1.y2.cy3 * (  1/TIn$p010^2 * LgTRI$d.3 * LgTRI$p.3.01 * LgTRI$d11.13 * LgTRI$p13.g   + 1/TIn$p010 * (LgTRI$d11.13 * ((TIn$theta13 * TIn$mar1 - TIn$mar3)/(1 - TIn$theta13^2)) * LgTRI$p13.g   - LgTRI$d11.13 * (d13/LgTRI$sd.13) * (TIn$theta23 - TIn$theta12 * TIn$theta13)/(1 - TIn$theta13^2) ) ) +
    respvec$cy1.cy2.cy3 * (  1/TIn$p000^2 * LgTRI$d.3 * LgTRI$p.3.00 * LgTRI$d11.13 * LgTRI$p13.g.c + 1/TIn$p000 * (LgTRI$d11.13 * ((TIn$theta13 * TIn$mar1 - TIn$mar3)/(1 - TIn$theta13^2)) * LgTRI$p13.g.c + LgTRI$d11.13 * (d13/LgTRI$sd.13) * (TIn$theta23 - TIn$theta12 * TIn$theta13)/(1 - TIn$theta13^2) ) ) -
    respvec$cy1.cy2.y3 * ( -1/TIn$p001^2 * LgTRI$d.3 * LgTRI$p.3.00 * LgTRI$d11.13 * LgTRI$p13.g.c + 1/TIn$p001 * (LgTRI$d11.13 * ((TIn$theta13 * TIn$mar1 - TIn$mar3)/(1 - TIn$theta13^2)) * LgTRI$p13.g.c + LgTRI$d11.13 * (d13/LgTRI$sd.13) * (TIn$theta23 - TIn$theta12 * TIn$theta13)/(1 - TIn$theta13^2) ) ) -
    respvec$y1.cy2.cy3 * (  1/TIn$p100^2 * LgTRI$d.3 * LgTRI$p.3.10 * LgTRI$d11.13 * LgTRI$p13.g.c + 1/TIn$p100 * (LgTRI$d11.13 * ((TIn$theta13 * TIn$mar1 - TIn$mar3)/(1 - TIn$theta13^2)) * LgTRI$p13.g.c + LgTRI$d11.13 * (d13/LgTRI$sd.13) * (TIn$theta23 - TIn$theta12 * TIn$theta13)/(1 - TIn$theta13^2) ) ) +
    respvec$y1.cy2.y3 * ( -1/TIn$p101^2 * LgTRI$d.3 * LgTRI$p.3.10 * LgTRI$d11.13 * LgTRI$p13.g.c + 1/TIn$p101 * (LgTRI$d11.13 * ((TIn$theta13 * TIn$mar1 - TIn$mar3)/(1 - TIn$theta13^2)) * LgTRI$p13.g.c + LgTRI$d11.13 * (d13/LgTRI$sd.13) * (TIn$theta23 - TIn$theta12 * TIn$theta13)/(1 - TIn$theta13^2) ) ) 
  
  #####
  # ! #
  ######################
  ## This bit is new  ##                 
  ######################
  
  d2l.de3.theta13 <- d2l.dF3.theta13 * LgTRI$dF3.de3
  
  ###########################################################
  
  d2l.dF3.theta23 <- respvec$y1.y2.y3 * ( -1/TIn$p111^2 * LgTRI$d.3 * LgTRI$p.3.11 * LgTRI$d11.23 * LgTRI$p23.g   + 1/TIn$p111 * (LgTRI$d11.23 * ((TIn$theta23*TIn$mar2  - TIn$mar3)/(1 - TIn$theta23^2)) * LgTRI$p23.g   - LgTRI$d11.23 * d23/LgTRI$sd.23 * (TIn$theta13 - TIn$theta12 * TIn$theta23)/(1 - TIn$theta23^2) ) ) -
    respvec$y1.y2.cy3 * (  1/TIn$p110^2 * LgTRI$d.3 * LgTRI$p.3.11 * LgTRI$d11.23 * LgTRI$p23.g   + 1/TIn$p110 * (LgTRI$d11.23 * ((TIn$theta23*TIn$mar2  - TIn$mar3)/(1 - TIn$theta23^2)) * LgTRI$p23.g   - LgTRI$d11.23 * d23/LgTRI$sd.23 * (TIn$theta13 - TIn$theta12*TIn$theta23)/(1 - TIn$theta23^2) ) ) +
    respvec$cy1.y2.y3 * ( -1/TIn$p011^2 * LgTRI$d.3 * LgTRI$p.3.01 * LgTRI$d11.23 * LgTRI$p23.g.c + 1/TIn$p011 * (LgTRI$d11.23 * ((TIn$theta23*TIn$mar2  - TIn$mar3)/(1 - TIn$theta23^2)) * LgTRI$p23.g.c + LgTRI$d11.23 * d23/LgTRI$sd.23 * (TIn$theta13 - TIn$theta12*TIn$theta23)/(1 - TIn$theta23^2) ) ) -
    respvec$cy1.y2.cy3 * (  1/TIn$p010^2 * LgTRI$d.3 * LgTRI$p.3.01 * LgTRI$d11.23 * LgTRI$p23.g.c + 1/TIn$p010 * (LgTRI$d11.23 * ((TIn$theta23*TIn$mar2  - TIn$mar3)/(1 - TIn$theta23^2)) * LgTRI$p23.g.c + LgTRI$d11.23 * d23/LgTRI$sd.23 * (TIn$theta13 - TIn$theta12 * TIn$theta23)/(1 - TIn$theta23^2) ) ) +
    respvec$cy1.cy2.cy3 * (  1/TIn$p000^2 * LgTRI$d.3 * LgTRI$p.3.00 * LgTRI$d11.23 * LgTRI$p23.g.c + 1/TIn$p000 * (LgTRI$d11.23 * ((TIn$theta23*TIn$mar2  - TIn$mar3)/(1 - TIn$theta23^2)) * LgTRI$p23.g.c + LgTRI$d11.23 * d23/LgTRI$sd.23 * (TIn$theta13 - TIn$theta12 * TIn$theta23)/(1 - TIn$theta23^2) ) ) -
    respvec$cy1.cy2.y3 * ( -1/TIn$p001^2 * LgTRI$d.3 * LgTRI$p.3.00 * LgTRI$d11.23 * LgTRI$p23.g.c + 1/TIn$p001 * (LgTRI$d11.23 * ((TIn$theta23*TIn$mar2  - TIn$mar3)/(1 - TIn$theta23^2)) * LgTRI$p23.g.c + LgTRI$d11.23 * d23/LgTRI$sd.23 * (TIn$theta13 - TIn$theta12 * TIn$theta23)/(1 - TIn$theta23^2) ) ) +
    respvec$y1.cy2.cy3 * (  1/TIn$p100^2 * LgTRI$d.3 * LgTRI$p.3.10 * LgTRI$d11.23 * LgTRI$p23.g   + 1/TIn$p100 * (LgTRI$d11.23 * ((TIn$theta23*TIn$mar2  - TIn$mar3)/(1 - TIn$theta23^2)) * LgTRI$p23.g   - LgTRI$d11.23 * d23/LgTRI$sd.23 * (TIn$theta13 - TIn$theta12 * TIn$theta23)/(1 - TIn$theta23^2) ) ) -
    respvec$y1.cy2.y3 * ( -1/TIn$p101^2 * LgTRI$d.3 * LgTRI$p.3.10 * LgTRI$d11.23 * LgTRI$p23.g   + 1/TIn$p101 * (LgTRI$d11.23 * ((TIn$theta23*TIn$mar2  - TIn$mar3)/(1 - TIn$theta23^2)) * LgTRI$p23.g   - LgTRI$d11.23 * d23/LgTRI$sd.23 * (TIn$theta13 - TIn$theta12 * TIn$theta23)/(1 - TIn$theta23^2) ) ) 
  
  #####
  # ! #
  ######################
  ## This bit is new  ##                 
  ######################
  
  d2l.de3.theta23 <- d2l.dF3.theta23 * LgTRI$dF3.de3
  
  ###########################################################
  
  
  #####
  # ! #
  ########################################################################
  ## I replaced TIn$eta1 with TIn$mar1. Same for TIn$eta2 and TIn$eta3  ##                 
  ########################################################################
  
  dd11.12.dtheta12 <- (LgTRI$d11.12/(1 - TIn$theta12^2)) * ( ( (TIn$theta12 * TIn$mar2 - TIn$mar1) * (TIn$theta12 * TIn$mar1 - TIn$mar2)/(1 - TIn$theta12^2) ) + TIn$theta12)
  dd11.13.dtheta13 <- (LgTRI$d11.13/(1 - TIn$theta13^2)) * ( ( (TIn$theta13 * TIn$mar3 - TIn$mar1) * (TIn$theta13 * TIn$mar1 - TIn$mar3)/(1 - TIn$theta13^2) ) + TIn$theta13)
  dd11.23.dtheta23 <- (LgTRI$d11.23/(1 - TIn$theta23^2)) * ( ( (TIn$theta23 * TIn$mar3 - TIn$mar2) * (TIn$theta23 * TIn$mar2 - TIn$mar3)/(1 - TIn$theta23^2) ) + TIn$theta23)
  
  dmean12.dtheta12 <- ( (-TIn$mar1 * TIn$theta23 - TIn$mar2 * TIn$theta13) * (1 - TIn$theta12^2) + 2 * TIn$theta12 * (TIn$mar1 * (TIn$theta13 - TIn$theta12 * TIn$theta23) + TIn$mar2 * (TIn$theta23 - TIn$theta12 * TIn$theta13)))/(1 - TIn$theta12^2)^2
  dmean13.dtheta13 <- ( (-TIn$mar1 * TIn$theta23 - TIn$mar3 * TIn$theta12) * (1 - TIn$theta13^2) + 2 * TIn$theta13 * (TIn$mar1 * (TIn$theta12 - TIn$theta13 * TIn$theta23) + TIn$mar3 * (TIn$theta23 - TIn$theta12 * TIn$theta13)))/(1 - TIn$theta13^2)^2
  dmean23.dtheta23 <- ( (-TIn$mar2 * TIn$theta13 - TIn$mar3 * TIn$theta12) * (1 - TIn$theta23^2) + 2 * TIn$theta23 * (TIn$mar2 * (TIn$theta12 - TIn$theta13 * TIn$theta23) + TIn$mar3 * (TIn$theta13 - TIn$theta12 * TIn$theta23)))/(1 - TIn$theta23^2)^2
  
  dmean13.dtheta12 <- ( TIn$mar1 - TIn$mar3 * TIn$theta13 )/( 1 - TIn$theta13^2 )
  dmean23.dtheta12 <- ( TIn$mar2 - TIn$mar3 * TIn$theta23 )/( 1 - TIn$theta23^2 )
  dmean23.dtheta13 <- ( TIn$mar3 - TIn$mar2 * TIn$theta23 )/( 1 - TIn$theta23^2 )
  
  #############################################################################################################
  
  dvar12.dtheta12 <- ( 2 * (TIn$theta13 * TIn$theta23 - TIn$theta12) * (1 - TIn$theta12^2) + 2 * TIn$theta12 * deno )/(1 - TIn$theta12^2)^2
  dvar13.dtheta13 <- ( 2 * (TIn$theta12 * TIn$theta23 - TIn$theta13) * (1 - TIn$theta13^2) + 2 * TIn$theta13 * deno )/(1 - TIn$theta13^2)^2
  dvar23.dtheta23 <- ( 2 * (TIn$theta12 * TIn$theta13 - TIn$theta23) * (1 - TIn$theta23^2) + 2 * TIn$theta23 * deno )/(1 - TIn$theta23^2)^2
  
  dvar13.dtheta12 <- ( 2 * (TIn$theta13 * TIn$theta23 - TIn$theta12 ) )/(1 - TIn$theta13^2)
  dvar23.dtheta12 <- ( 2 * (TIn$theta13 * TIn$theta23 - TIn$theta12 ) )/(1 - TIn$theta23^2)
  dvar23.dtheta13 <- ( 2 * (TIn$theta12 * TIn$theta23 - TIn$theta13 ) )/(1 - TIn$theta23^2)
  
  
  d2l.dtheta12.theta12 <- respvec$y1.y2.y3 * ( -1/TIn$p111^2 * (LgTRI$d11.12 * LgTRI$p12.g   )^2 + 1/TIn$p111 * (( dd11.12.dtheta12 * LgTRI$p12.g  ) - ((d12/LgTRI$sd.12) * (dmean12.dtheta12  + (((TIn$mar3 - LgTRI$mean.12)/(2 * LgTRI$sd.12^2)) * dvar12.dtheta12)) * LgTRI$d11.12) )) +
    respvec$y1.y2.cy3 * ( -1/TIn$p110^2   * (LgTRI$d11.12 * LgTRI$p12.g.c )^2 + 1/TIn$p110 * (( dd11.12.dtheta12 * LgTRI$p12.g.c) + ((d12/LgTRI$sd.12) * (dmean12.dtheta12  + (((TIn$mar3 - LgTRI$mean.12)/(2 * LgTRI$sd.12^2)) * dvar12.dtheta12)) * LgTRI$d11.12) )) -
    respvec$cy1.y2.y3 * (  1/TIn$p011^2   * (LgTRI$d11.12 * LgTRI$p12.g   )^2 + 1/TIn$p011 * (( dd11.12.dtheta12 * LgTRI$p12.g  ) - ((d12/LgTRI$sd.12) * (dmean12.dtheta12  + (((TIn$mar3 - LgTRI$mean.12)/(2 * LgTRI$sd.12^2)) * dvar12.dtheta12)) * LgTRI$d11.12) )) -
    respvec$cy1.y2.cy3 * (  1/TIn$p010^2  * (LgTRI$d11.12 * LgTRI$p12.g.c )^2 + 1/TIn$p010 * (( dd11.12.dtheta12 * LgTRI$p12.g.c) + ((d12/LgTRI$sd.12) * (dmean12.dtheta12  + (((TIn$mar3 - LgTRI$mean.12)/(2 * LgTRI$sd.12^2)) * dvar12.dtheta12)) * LgTRI$d11.12) )) +
    respvec$cy1.cy2.cy3 * ( -1/TIn$p000^2 * (LgTRI$d11.12 * LgTRI$p12.g.c )^2 + 1/TIn$p000 * (( dd11.12.dtheta12 * LgTRI$p12.g.c) + ((d12/LgTRI$sd.12) * (dmean12.dtheta12  + (((TIn$mar3 - LgTRI$mean.12)/(2 * LgTRI$sd.12^2)) * dvar12.dtheta12)) * LgTRI$d11.12) )) +
    respvec$cy1.cy2.y3 * ( -1/TIn$p001^2  * (LgTRI$d11.12 * LgTRI$p12.g   )^2 + 1/TIn$p001 * (( dd11.12.dtheta12 * LgTRI$p12.g  ) - ((d12/LgTRI$sd.12) * (dmean12.dtheta12  + (((TIn$mar3 - LgTRI$mean.12)/(2 * LgTRI$sd.12^2)) * dvar12.dtheta12)) * LgTRI$d11.12) )) -
    respvec$y1.cy2.cy3 * (  1/TIn$p100^2  * (LgTRI$d11.12 * LgTRI$p12.g.c )^2 + 1/TIn$p100 * (( dd11.12.dtheta12 * LgTRI$p12.g.c) + ((d12/LgTRI$sd.12) * (dmean12.dtheta12  + (((TIn$mar3 - LgTRI$mean.12)/(2 * LgTRI$sd.12^2)) * dvar12.dtheta12)) * LgTRI$d11.12) )) -
    respvec$y1.cy2.y3 * (  1/TIn$p101^2   * (LgTRI$d11.12 * LgTRI$p12.g   )^2 + 1/TIn$p101 * (( dd11.12.dtheta12 * LgTRI$p12.g  ) - ((d12/LgTRI$sd.12) * (dmean12.dtheta12  + (((TIn$mar3 - LgTRI$mean.12)/(2 * LgTRI$sd.12^2)) * dvar12.dtheta12)) * LgTRI$d11.12) )) 
  
  d2l.dtheta13.theta13 <- respvec$y1.y2.y3 * ( -1/TIn$p111^2 * (LgTRI$d11.13 * LgTRI$p13.g)^2 + 1/TIn$p111 * (( dd11.13.dtheta13 * LgTRI$p13.g  ) - ((d13/LgTRI$sd.13) * (dmean13.dtheta13 + (((TIn$mar2 - LgTRI$mean.13)/(2 * LgTRI$sd.13^2)) * dvar13.dtheta13)) * LgTRI$d11.13) )) -
    respvec$y1.y2.cy3 * (  1/TIn$p110^2 * (LgTRI$d11.13 * LgTRI$p13.g)^2 + 1/TIn$p110 * ((dd11.13.dtheta13 * LgTRI$p13.g  ) - ((d13/LgTRI$sd.13) * (dmean13.dtheta13 + (((TIn$mar2 - LgTRI$mean.13)/(2 * LgTRI$sd.13^2)) * dvar13.dtheta13)) * LgTRI$d11.13) )) -
    respvec$cy1.y2.y3 * (  1/TIn$p011^2 * (LgTRI$d11.13 * LgTRI$p13.g)^2 + 1/TIn$p011 * ((dd11.13.dtheta13 * LgTRI$p13.g  ) - ((d13/LgTRI$sd.13) * (dmean13.dtheta13 + (((TIn$mar2 - LgTRI$mean.13)/(2 * LgTRI$sd.13^2)) * dvar13.dtheta13)) * LgTRI$d11.13) )) +
    respvec$cy1.y2.cy3 * ( -1/TIn$p010^2 * (LgTRI$d11.13 * LgTRI$p13.g)^2 + 1/TIn$p010 * ((dd11.13.dtheta13 * LgTRI$p13.g  ) - ((d13/LgTRI$sd.13) * (dmean13.dtheta13 + (((TIn$mar2 - LgTRI$mean.13)/(2 * LgTRI$sd.13^2)) * dvar13.dtheta13)) * LgTRI$d11.13) )) +
    respvec$cy1.cy2.cy3 * ( -1/TIn$p000^2 * (LgTRI$d11.13 * LgTRI$p13.g.c )^2 + 1/TIn$p000 * ((dd11.13.dtheta13 * LgTRI$p13.g.c) + ((d13/LgTRI$sd.13) * (dmean13.dtheta13 + (((TIn$mar2 - LgTRI$mean.13)/(2 * LgTRI$sd.13^2)) * dvar13.dtheta13)) * LgTRI$d11.13) )) -
    respvec$cy1.cy2.y3 * (  1/TIn$p001^2 * (LgTRI$d11.13 * LgTRI$p13.g.c )^2 + 1/TIn$p001 * ((dd11.13.dtheta13 * LgTRI$p13.g.c) + ((d13/LgTRI$sd.13) * (dmean13.dtheta13 + (((TIn$mar2 - LgTRI$mean.13)/(2 * LgTRI$sd.13^2)) * dvar13.dtheta13)) * LgTRI$d11.13) )) -
    respvec$y1.cy2.cy3 * (  1/TIn$p100^2 * (LgTRI$d11.13 * LgTRI$p13.g.c )^2 + 1/TIn$p100 * ((dd11.13.dtheta13 * LgTRI$p13.g.c) + ((d13/LgTRI$sd.13) * (dmean13.dtheta13 + (((TIn$mar2 - LgTRI$mean.13)/(2 * LgTRI$sd.13^2)) * dvar13.dtheta13)) * LgTRI$d11.13) )) +
    respvec$y1.cy2.y3 * ( -1/TIn$p101^2 * (LgTRI$d11.13 * LgTRI$p13.g.c )^2 + 1/TIn$p101 * ((dd11.13.dtheta13 * LgTRI$p13.g.c) + ((d13/LgTRI$sd.13) * (dmean13.dtheta13 + (((TIn$mar2 - LgTRI$mean.13)/(2 * LgTRI$sd.13^2)) * dvar13.dtheta13)) * LgTRI$d11.13) ))   
  
  d2l.dtheta23.theta23 <- respvec$y1.y2.y3 * ( -1/TIn$p111^2 * (LgTRI$d11.23 * LgTRI$p23.g   )^2 + 1/TIn$p111 * ((dd11.23.dtheta23 * LgTRI$p23.g  ) - ((d23/LgTRI$sd.23) * (dmean23.dtheta23 + ((TIn$mar1 - LgTRI$mean.23)/(2*LgTRI$sd.23^2) * dvar23.dtheta23)) * LgTRI$d11.23) )) -
    respvec$y1.y2.cy3 * (  1/TIn$p110^2 * (LgTRI$d11.23 * LgTRI$p23.g   )^2 + 1/TIn$p110 * ((dd11.23.dtheta23 * LgTRI$p23.g  ) - ((d23/LgTRI$sd.23) * (dmean23.dtheta23 + ((TIn$mar1 - LgTRI$mean.23)/(2*LgTRI$sd.23^2) * dvar23.dtheta23)) * LgTRI$d11.23) )) +
    respvec$cy1.y2.y3 * ( -1/TIn$p011^2 * (LgTRI$d11.23 * LgTRI$p23.g.c )^2 + 1/TIn$p011 * ((dd11.23.dtheta23 * LgTRI$p23.g.c) + ((d23/LgTRI$sd.23) * (dmean23.dtheta23 + ((TIn$mar1 - LgTRI$mean.23)/(2*LgTRI$sd.23^2) * dvar23.dtheta23)) * LgTRI$d11.23) )) -
    respvec$cy1.y2.cy3 * (  1/TIn$p010^2 * (LgTRI$d11.23 * LgTRI$p23.g.c )^2 + 1/TIn$p010 * ((dd11.23.dtheta23 * LgTRI$p23.g.c) + ((d23/LgTRI$sd.23) * (dmean23.dtheta23 + ((TIn$mar1 - LgTRI$mean.23)/(2*LgTRI$sd.23^2) * dvar23.dtheta23)) * LgTRI$d11.23) )) +
    respvec$cy1.cy2.cy3 * ( -1/TIn$p000^2 * (LgTRI$d11.23 * LgTRI$p23.g.c )^2 + 1/TIn$p000 * ((dd11.23.dtheta23 * LgTRI$p23.g.c) + ((d23/LgTRI$sd.23) * (dmean23.dtheta23 + ((TIn$mar1 - LgTRI$mean.23)/(2*LgTRI$sd.23^2) * dvar23.dtheta23)) * LgTRI$d11.23) )) -
    respvec$cy1.cy2.y3 * (  1/TIn$p001^2 * (LgTRI$d11.23 * LgTRI$p23.g.c )^2 + 1/TIn$p001 * ((dd11.23.dtheta23 * LgTRI$p23.g.c) + ((d23/LgTRI$sd.23) * (dmean23.dtheta23 + ((TIn$mar1 - LgTRI$mean.23)/(2*LgTRI$sd.23^2) * dvar23.dtheta23)) * LgTRI$d11.23) )) +
    respvec$y1.cy2.cy3 * ( -1/TIn$p100^2 * (LgTRI$d11.23 * LgTRI$p23.g   )^2 + 1/TIn$p100 * ((dd11.23.dtheta23 * LgTRI$p23.g  ) - ((d23/LgTRI$sd.23) * (dmean23.dtheta23 + ((TIn$mar1 - LgTRI$mean.23)/(2*LgTRI$sd.23^2) * dvar23.dtheta23)) * LgTRI$d11.23) )) -
    respvec$y1.cy2.y3 * (  1/TIn$p101^2 * (LgTRI$d11.23 * LgTRI$p23.g   )^2 + 1/TIn$p101 * ((dd11.23.dtheta23 * LgTRI$p23.g  ) - ((d23/LgTRI$sd.23) * (dmean23.dtheta23 + ((TIn$mar1 - LgTRI$mean.23)/(2*LgTRI$sd.23^2) * dvar23.dtheta23)) * LgTRI$d11.23) )) 
  
  d2l.dtheta12.theta13 <- respvec$y1.y2.y3 * ( -1/TIn$p111^2 * LgTRI$d11.12 * LgTRI$p12.g   * LgTRI$d11.13 * LgTRI$p13.g   - 1/TIn$p111 * LgTRI$d11.13 * (d13/LgTRI$sd.13) * (dmean13.dtheta12 + ((TIn$mar2-LgTRI$mean.13)/(2*LgTRI$sd.13^2) * dvar13.dtheta12) )) -
    respvec$y1.y2.cy3 * ( -1/TIn$p110^2 * LgTRI$d11.12 * LgTRI$p12.g.c * LgTRI$d11.13 * LgTRI$p13.g   - 1/TIn$p110 * LgTRI$d11.13 * (d13/LgTRI$sd.13) * (dmean13.dtheta12 + ((TIn$mar2-LgTRI$mean.13)/(2*LgTRI$sd.13^2) * dvar13.dtheta12) )) -
    respvec$cy1.y2.y3 * (  1/TIn$p011^2 * LgTRI$d11.12 * LgTRI$p12.g   * LgTRI$d11.13 * LgTRI$p13.g   - 1/TIn$p011 * LgTRI$d11.13 * (d13/LgTRI$sd.13) * (dmean13.dtheta12 + ((TIn$mar2-LgTRI$mean.13)/(2*LgTRI$sd.13^2) * dvar13.dtheta12) )) +
    respvec$cy1.y2.cy3 * (  1/TIn$p010^2 * LgTRI$d11.12 * LgTRI$p12.g.c * LgTRI$d11.13 * LgTRI$p13.g   - 1/TIn$p010 * LgTRI$d11.13 * (d13/LgTRI$sd.13) * (dmean13.dtheta12 + ((TIn$mar2-LgTRI$mean.13)/(2*LgTRI$sd.13^2) * dvar13.dtheta12) )) +
    respvec$cy1.cy2.cy3 * ( -1/TIn$p000^2 * LgTRI$d11.12 * LgTRI$p12.g.c * LgTRI$d11.13 * LgTRI$p13.g.c + 1/TIn$p000 * LgTRI$d11.13 * (d13/LgTRI$sd.13) * (dmean13.dtheta12 + ((TIn$mar2-LgTRI$mean.13)/(2*LgTRI$sd.13^2) * dvar13.dtheta12) )) -
    respvec$cy1.cy2.y3 * ( -1/TIn$p001^2 * LgTRI$d11.12 * LgTRI$p12.g   * LgTRI$d11.13 * LgTRI$p13.g.c + 1/TIn$p001 * LgTRI$d11.13 * (d13/LgTRI$sd.13) * (dmean13.dtheta12 + ((TIn$mar2-LgTRI$mean.13)/(2*LgTRI$sd.13^2) * dvar13.dtheta12) )) -
    respvec$y1.cy2.cy3 * (  1/TIn$p100^2 * LgTRI$d11.12 * LgTRI$p12.g.c * LgTRI$d11.13 * LgTRI$p13.g.c + 1/TIn$p100 * LgTRI$d11.13 * (d13/LgTRI$sd.13) * (dmean13.dtheta12 + ((TIn$mar2-LgTRI$mean.13)/(2*LgTRI$sd.13^2) * dvar13.dtheta12) )) +
    respvec$y1.cy2.y3 * (  1/TIn$p101^2 * LgTRI$d11.12 * LgTRI$p12.g   * LgTRI$d11.13 * LgTRI$p13.g.c + 1/TIn$p101 * LgTRI$d11.13 * (d13/LgTRI$sd.13) * (dmean13.dtheta12 + ((TIn$mar2-LgTRI$mean.13)/(2*LgTRI$sd.13^2) * dvar13.dtheta12) ))
  
  d2l.dtheta12.theta23 <- respvec$y1.y2.y3 * ( -1/TIn$p111^2 * LgTRI$d11.12 * LgTRI$p12.g   * LgTRI$d11.23 * LgTRI$p23.g   - 1/TIn$p111 * LgTRI$d11.23 * (d23/LgTRI$sd.23) * (dmean23.dtheta12 + ((TIn$mar1-LgTRI$mean.23)/(2*LgTRI$sd.23^2) * dvar23.dtheta12) )) -
    respvec$y1.y2.cy3 * ( -1/TIn$p110^2 * LgTRI$d11.12 * LgTRI$p12.g.c * LgTRI$d11.23 * LgTRI$p23.g   - 1/TIn$p110 * LgTRI$d11.23 * (d23/LgTRI$sd.23) * (dmean23.dtheta12 + ((TIn$mar1-LgTRI$mean.23)/(2*LgTRI$sd.23^2) * dvar23.dtheta12) )) +
    respvec$cy1.y2.y3 * (  1/TIn$p011^2 * LgTRI$d11.12 * LgTRI$p12.g   * LgTRI$d11.23 * LgTRI$p23.g.c + 1/TIn$p011 * LgTRI$d11.23 * (d23/LgTRI$sd.23) * (dmean23.dtheta12 + ((TIn$mar1-LgTRI$mean.23)/(2*LgTRI$sd.23^2) * dvar23.dtheta12) )) -
    respvec$cy1.y2.cy3 * (  1/TIn$p010^2 * LgTRI$d11.12 * LgTRI$p12.g.c * LgTRI$d11.23 * LgTRI$p23.g.c + 1/TIn$p010 * LgTRI$d11.23 * (d23/LgTRI$sd.23) * (dmean23.dtheta12 + ((TIn$mar1-LgTRI$mean.23)/(2*LgTRI$sd.23^2) * dvar23.dtheta12) )) +
    respvec$cy1.cy2.cy3 * ( -1/TIn$p000^2 * LgTRI$d11.12 * LgTRI$p12.g.c * LgTRI$d11.23 * LgTRI$p23.g.c + 1/TIn$p000 * LgTRI$d11.23 * (d23/LgTRI$sd.23) * (dmean23.dtheta12 + ((TIn$mar1-LgTRI$mean.23)/(2*LgTRI$sd.23^2) * dvar23.dtheta12) )) -
    respvec$cy1.cy2.y3 * ( -1/TIn$p001^2 * LgTRI$d11.12 * LgTRI$p12.g   * LgTRI$d11.23 * LgTRI$p23.g.c + 1/TIn$p001 * LgTRI$d11.23 * (d23/LgTRI$sd.23) * (dmean23.dtheta12 + ((TIn$mar1-LgTRI$mean.23)/(2*LgTRI$sd.23^2) * dvar23.dtheta12) )) +
    respvec$y1.cy2.cy3 * (  1/TIn$p100^2 * LgTRI$d11.12 * LgTRI$p12.g.c * LgTRI$d11.23 * LgTRI$p23.g   - 1/TIn$p100 * LgTRI$d11.23 * (d23/LgTRI$sd.23) * (dmean23.dtheta12 + ((TIn$mar1-LgTRI$mean.23)/(2*LgTRI$sd.23^2) * dvar23.dtheta12) )) -
    respvec$y1.cy2.y3 * (  1/TIn$p101^2 * LgTRI$d11.12 * LgTRI$p12.g   * LgTRI$d11.23 * LgTRI$p23.g   - 1/TIn$p101 * LgTRI$d11.23 * (d23/LgTRI$sd.23) * (dmean23.dtheta12 + ((TIn$mar1-LgTRI$mean.23)/(2*LgTRI$sd.23^2) * dvar23.dtheta12) ))
  
  d2l.dtheta13.theta23 <- respvec$y1.y2.y3 * ( -1/TIn$p111^2 * LgTRI$d11.13 * LgTRI$p13.g   * LgTRI$d11.23 * LgTRI$p23.g   - 1/TIn$p111 * LgTRI$d11.23 * (d23/LgTRI$sd.23) * (dmean23.dtheta13 + ((TIn$mar1-LgTRI$mean.23)/(2 * LgTRI$sd.23^2) * dvar23.dtheta13) )) -
    respvec$y1.y2.cy3 * (  1/TIn$p110^2 * LgTRI$d11.13 * LgTRI$p13.g   * LgTRI$d11.23 * LgTRI$p23.g   - 1/TIn$p110 * LgTRI$d11.23 * (d23/LgTRI$sd.23) * (dmean23.dtheta13 + ((TIn$mar1-LgTRI$mean.23)/(2 * LgTRI$sd.23^2) * dvar23.dtheta13) )) +
    respvec$cy1.y2.y3 * (  1/TIn$p011^2 * LgTRI$d11.13 * LgTRI$p13.g   * LgTRI$d11.23 * LgTRI$p23.g.c + 1/TIn$p011 * LgTRI$d11.23 * (d23/LgTRI$sd.23) * (dmean23.dtheta13 + ((TIn$mar1-LgTRI$mean.23)/(2 * LgTRI$sd.23^2) * dvar23.dtheta13) )) -
    respvec$cy1.y2.cy3 * ( -1/TIn$p010^2 * LgTRI$d11.13 * LgTRI$p13.g   * LgTRI$d11.23 * LgTRI$p23.g.c + 1/TIn$p010 * LgTRI$d11.23 * (d23/LgTRI$sd.23) * (dmean23.dtheta13 + ((TIn$mar1-LgTRI$mean.23)/(2 * LgTRI$sd.23^2) * dvar23.dtheta13) )) +
    respvec$cy1.cy2.cy3 * ( -1/TIn$p000^2 * LgTRI$d11.13 * LgTRI$p13.g.c * LgTRI$d11.23 * LgTRI$p23.g.c + 1/TIn$p000 * LgTRI$d11.23 * (d23/LgTRI$sd.23) * (dmean23.dtheta13 + ((TIn$mar1-LgTRI$mean.23)/(2 * LgTRI$sd.23^2) * dvar23.dtheta13) )) -
    respvec$cy1.cy2.y3 * (  1/TIn$p001^2 * LgTRI$d11.13 * LgTRI$p13.g.c * LgTRI$d11.23 * LgTRI$p23.g.c + 1/TIn$p001 * LgTRI$d11.23 * (d23/LgTRI$sd.23) * (dmean23.dtheta13 + ((TIn$mar1-LgTRI$mean.23)/(2 * LgTRI$sd.23^2) * dvar23.dtheta13) )) +
    respvec$y1.cy2.cy3 * (  1/TIn$p100^2 * LgTRI$d11.13 * LgTRI$p13.g.c * LgTRI$d11.23 * LgTRI$p23.g   - 1/TIn$p100 * LgTRI$d11.23 * (d23/LgTRI$sd.23) * (dmean23.dtheta13 + ((TIn$mar1-LgTRI$mean.23)/(2 * LgTRI$sd.23^2) * dvar23.dtheta13) )) -
    respvec$y1.cy2.y3 * ( -1/TIn$p101^2 * LgTRI$d11.13 * LgTRI$p13.g.c * LgTRI$d11.23 * LgTRI$p23.g   - 1/TIn$p101 * LgTRI$d11.23 * (d23/LgTRI$sd.23) * (dmean23.dtheta13 + ((TIn$mar1-LgTRI$mean.23)/(2 * LgTRI$sd.23^2) * dvar23.dtheta13) )) 
  
  
  d2l.dbe1.be1 <- crossprod( VC$X1 * c(VC$weights*d2l.de1.e1), VC$X1 )
  d2l.dbe2.be2 <- crossprod( VC$X2 * c(VC$weights*d2l.de2.e2), VC$X2 )
  d2l.dbe3.be3 <- crossprod( VC$X3 * c(VC$weights*d2l.de3.e3), VC$X3 )
  d2l.dbe1.be2 <- crossprod( VC$X1 * c(VC$weights*d2l.de1.e2), VC$X2 )
  d2l.dbe1.be3 <- crossprod( VC$X1 * c(VC$weights*d2l.de1.e3), VC$X3 )
  d2l.dbe2.be3 <- crossprod( VC$X2 * c(VC$weights*d2l.de2.e3), VC$X3 )
  
  if(VC$Chol == FALSE){
    dtheta12.dtheta12.st <- 4 * exp( 2 * TIn$theta12.st )/( exp(2 * TIn$theta12.st) + 1 )^2
    dtheta13.dtheta13.st <- 4 * exp( 2 * TIn$theta13.st )/( exp(2 * TIn$theta13.st) + 1 )^2
    dtheta23.dtheta23.st <- 4 * exp( 2 * TIn$theta23.st )/( exp(2 * TIn$theta23.st) + 1 )^2
    
    d2theta12.theta12.st <- (8 * exp( 2 * TIn$theta12.st ) - 8 * exp( 4 * TIn$theta12.st ))/( exp(2 * TIn$theta12.st) + 1 )^3
    d2theta13.theta13.st <- (8 * exp( 2 * TIn$theta13.st ) - 8 * exp( 4 * TIn$theta13.st ))/( exp(2 * TIn$theta13.st) + 1 )^3
    d2theta23.theta23.st <- (8 * exp( 2 * TIn$theta23.st ) - 8 * exp( 4 * TIn$theta23.st ))/( exp(2 * TIn$theta23.st) + 1 )^3
    
    
    d2l.dbe1.theta12.st <- colSums(c(VC$weights*d2l.de1.theta12* dtheta12.dtheta12.st) * VC$X1) 
    d2l.dbe1.theta13.st <- colSums(c(VC$weights*d2l.de1.theta13* dtheta13.dtheta13.st) * VC$X1) 
    d2l.dbe1.theta23.st <- colSums(c(VC$weights*d2l.de1.theta23* dtheta23.dtheta23.st) * VC$X1) 
    
    d2l.dbe2.theta12.st <- colSums(c(VC$weights*d2l.de2.theta12* dtheta12.dtheta12.st) * VC$X2) 
    d2l.dbe2.theta13.st <- colSums(c(VC$weights*d2l.de2.theta13* dtheta13.dtheta13.st) * VC$X2) 
    d2l.dbe2.theta23.st <- colSums(c(VC$weights*d2l.de2.theta23* dtheta23.dtheta23.st) * VC$X2) 
    
    d2l.dbe3.theta12.st <- colSums(c(VC$weights*d2l.de3.theta12* dtheta12.dtheta12.st) * VC$X3) 
    d2l.dbe3.theta13.st <- colSums(c(VC$weights*d2l.de3.theta13* dtheta13.dtheta13.st) * VC$X3) 
    d2l.dbe3.theta23.st <- colSums(c(VC$weights*d2l.de3.theta23* dtheta23.dtheta23.st) * VC$X3) 
    
    d2l.dtheta12.st <- sum( VC$weights*(d2l.dtheta12.theta12 * (dtheta12.dtheta12.st)^2  + LgTRI$dl.dtheta12 * d2theta12.theta12.st) )
    d2l.dtheta13.st <- sum( VC$weights*(d2l.dtheta13.theta13 * (dtheta13.dtheta13.st)^2  + LgTRI$dl.dtheta13 * d2theta13.theta13.st) )
    d2l.dtheta23.st <- sum( VC$weights*(d2l.dtheta23.theta23 * (dtheta23.dtheta23.st)^2  + LgTRI$dl.dtheta23 * d2theta23.theta23.st) )
    
    d2l.dtheta12.st.theta13.st <- sum(VC$weights*d2l.dtheta12.theta13 * dtheta12.dtheta12.st * dtheta13.dtheta13.st )
    d2l.dtheta12.st.theta23.st <- sum(VC$weights*d2l.dtheta12.theta23 * dtheta12.dtheta12.st * dtheta23.dtheta23.st )
    d2l.dtheta13.st.theta23.st <- sum(VC$weights*d2l.dtheta13.theta23 * dtheta13.dtheta13.st * dtheta23.dtheta23.st )
    
  }
  
  if(VC$Chol == TRUE){
    
    # #### When I compare corrs with regressors and corrs without regressors: ####
    # TIn$theta12.st <- c(rep(TIn$theta12.st, VC$n))
    # TIn$theta13.st <- c(rep(TIn$theta13.st, VC$n))
    # TIn$theta23.st <- c(rep(TIn$theta23.st, VC$n))
    
    # TIn$theta12.st <- theta12.st
    # TIn$theta13.st <- theta13.st
    # TIn$theta23.st <- theta23.st
    
    ########################################################
    
    d2l.de1.theta <- matrix(c ( VC$weights*d2l.de1.theta12,
                                VC$weights*d2l.de1.theta13,
                                VC$weights*d2l.de1.theta23), VC$n, 3)
    
    d2l.de2.theta <- matrix(c ( VC$weights*d2l.de2.theta12,
                                VC$weights*d2l.de2.theta13,
                                VC$weights*d2l.de2.theta23), VC$n, 3)
    
    d2l.de3.theta <- matrix(c ( VC$weights*d2l.de3.theta12,
                                VC$weights*d2l.de3.theta13,
                                VC$weights*d2l.de3.theta23), VC$n, 3)
    
    dth12.dth12.st <- 1/(1 + TIn$theta12.st^2)^(3/2)
    dth13.dth13.st <- (1 + TIn$theta23.st^2)/(1 + TIn$theta13.st^2 + 
                                                TIn$theta23.st^2)^(3/2)
    dth13.dth23.st <- -(TIn$theta13.st * TIn$theta23.st)/(1 + 
                                                            TIn$theta13.st^2 + TIn$theta23.st^2)^(3/2)
    dth23.dth12.st <- TIn$theta13.st/sqrt((1 + TIn$theta12.st^2) * 
                                            (1 + TIn$theta13.st^2 + TIn$theta23.st^2)) - (TIn$theta12.st * 
                                                                                            (TIn$theta12.st * TIn$theta13.st + TIn$theta23.st))/((1 + 
                                                                                                                                                    TIn$theta12.st^2)^(3/2) * sqrt(1 + TIn$theta13.st^2 + 
                                                                                                                                                                                     TIn$theta23.st^2))
    dth23.dth13.st <- TIn$theta12.st/sqrt((1 + TIn$theta12.st^2) * 
                                            (1 + TIn$theta13.st^2 + TIn$theta23.st^2)) - (TIn$theta13.st * 
                                                                                            (TIn$theta12.st * TIn$theta13.st + TIn$theta23.st))/(sqrt(1 + 
                                                                                                                                                        TIn$theta12.st^2) * (1 + TIn$theta13.st^2 + TIn$theta23.st^2)^(3/2))
    dth23.dth23.st <- 1/sqrt((1 + TIn$theta12.st^2) * (1 + TIn$theta13.st^2 + TIn$theta23.st^2)) - (TIn$theta23.st * (TIn$theta12.st * TIn$theta13.st + TIn$theta23.st))/(sqrt(1 + 
                                                                                                                                                                                 TIn$theta12.st^2) * (1 + TIn$theta13.st^2 + TIn$theta23.st^2)^(3/2))
    
    
    d2th12.dth12.st.th12.st <- -(3 * TIn$theta12.st)/(1 + TIn$theta12.st^2)^(5/2)
    d2th13.dth13.st.th13.st<- -(3 * TIn$theta13.st * (1 + TIn$theta23.st^2))/(1 + TIn$theta13.st^2 + TIn$theta23.st^2)^(5/2)
    d2th13.dth23.st.th13.st <- (2 * TIn$theta23.st)/(1 + TIn$theta13.st^2 + TIn$theta23.st^2)^(3/2) - (3 * TIn$theta23.st * (1 + TIn$theta23.st^2))/(1 + TIn$theta13.st^2 + TIn$theta23.st^2)^(5/2)
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
    
    
    if(is.null(VC$X4)){ 
      
      dth12.dth13.st <- 0
      dth12.dth23.st <- 0
      dth13.dth12.st <- 0
      
      dtheta.theta.st <- matrix( c(  dth12.dth12.st,  dth13.dth12.st, dth23.dth12.st,
                                     dth12.dth13.st,  dth13.dth13.st, dth23.dth13.st,
                                     dth12.dth23.st,  dth13.dth23.st, dth23.dth23.st ), 3 , 3)
      
      d2l.dbe1.theta.st <- t(VC$X1) %*%  d2l.de1.theta %*% dtheta.theta.st 
      d2l.dbe2.theta.st <- t(VC$X2) %*%  d2l.de2.theta %*% dtheta.theta.st 
      d2l.dbe3.theta.st <- t(VC$X3) %*%  d2l.de3.theta %*% dtheta.theta.st 
      
      d2l.dbe1.theta12.st <- d2l.dbe1.theta.st[, 1]
      d2l.dbe1.theta13.st <- d2l.dbe1.theta.st[, 2]
      d2l.dbe1.theta23.st <- d2l.dbe1.theta.st[, 3]
      
      d2l.dbe2.theta12.st <- d2l.dbe2.theta.st[, 1]
      d2l.dbe2.theta13.st <- d2l.dbe2.theta.st[, 2]
      d2l.dbe2.theta23.st <- d2l.dbe2.theta.st[, 3]
      
      d2l.dbe3.theta12.st <- d2l.dbe3.theta.st[, 1]
      d2l.dbe3.theta13.st <- d2l.dbe3.theta.st[, 2]
      d2l.dbe3.theta23.st <- d2l.dbe3.theta.st[, 3]
      
      ############################
      ## d2l.dtheta.st.theta.st ## 
      ############################
      d2th12.dth13.st.th12.st <- 0
      d2th12.dth23.st.th12.st <- 0
      d2th12.dth12.st.th13.st <- 0
      d2th12.dth13.st.th13.st <- 0
      d2th12.dth23.st.th13.st <- 0
      d2th12.dth12.st.th23.st <- 0
      d2th12.dth13.st.th23.st <- 0
      d2th12.dth23.st.th23.st <- 0
      d2th13.dth12.st.th12.st <- 0
      d2th13.dth12.st.th23.st <- 0
      d2th13.dth13.st.th12.st <- 0
      d2th13.dth23.st.th12.st <- 0
      d2th13.dth12.st.th13.st <- 0    

      d2theta12.theta12.st <- matrix( c(          d2th12.dth12.st.th12.st,  d2th13.dth12.st.th12.st, d2th23.dth12.st.th12.st,
                                                  d2th12.dth12.st.th13.st,  d2th13.dth12.st.th13.st, d2th23.dth12.st.th13.st,
                                                  d2th12.dth12.st.th23.st,  d2th13.dth12.st.th23.st, d2th23.dth12.st.th23.st ), 3 , 3)
      
      d2theta13.theta13.st <- matrix( c(          d2th12.dth13.st.th12.st,  d2th13.dth13.st.th12.st, d2th23.dth13.st.th12.st,
                                                  d2th12.dth13.st.th13.st,  d2th13.dth13.st.th13.st, d2th23.dth13.st.th13.st,
                                                  d2th12.dth13.st.th23.st,  d2th13.dth13.st.th23.st, d2th23.dth13.st.th23.st ), 3 , 3)
      
      d2theta23.theta23.st <- matrix( c(          d2th12.dth23.st.th12.st,  d2th13.dth23.st.th12.st, d2th23.dth23.st.th12.st,
                                                  d2th12.dth23.st.th13.st,  d2th13.dth23.st.th13.st, d2th23.dth23.st.th13.st,
                                                  d2th12.dth23.st.th23.st,  d2th13.dth23.st.th23.st, d2th23.dth23.st.th23.st ), 3 , 3)
      
      d2l.theta <- matrix( c(          sum(VC$weights*d2l.dtheta12.theta12),  sum(VC$weights*d2l.dtheta12.theta13), sum(VC$weights*d2l.dtheta12.theta23),
                                       sum(VC$weights*d2l.dtheta12.theta13),  sum(VC$weights*d2l.dtheta13.theta13), sum(VC$weights*d2l.dtheta13.theta23),
                                       sum(VC$weights*d2l.dtheta12.theta23),  sum(VC$weights*d2l.dtheta13.theta23), sum(VC$weights*d2l.dtheta23.theta23) ), 3 , 3)
      
      dl.dtheta <- matrix(c ( LgTRI$dl.dtheta12,
                              LgTRI$dl.dtheta13,
                              LgTRI$dl.dtheta23), VC$n, 3)
      
      d2l.dtheta.st <- t(dtheta.theta.st) %*% d2l.theta %*% dtheta.theta.st + matrix(c(colSums(dl.dtheta) %*% d2theta12.theta12.st,
                                                                                       colSums(dl.dtheta) %*% d2theta13.theta13.st,
                                                                                       colSums(dl.dtheta) %*% d2theta23.theta23.st), 3, 3)
      
      d2l.dtheta12.st <- d2l.dtheta.st[1, 1]
      d2l.dtheta13.st <- d2l.dtheta.st[2, 2]
      d2l.dtheta23.st <- d2l.dtheta.st[3, 3]
      
      d2l.dtheta12.st.theta13.st <- d2l.dtheta.st[1, 2]
      d2l.dtheta12.st.theta23.st <- d2l.dtheta.st[1, 3]
      d2l.dtheta13.st.theta23.st <- d2l.dtheta.st[2, 3]
      
      
    }
    
    
    if(!is.null(VC$X4)){  
      
      dth12.dth13.st <- rep(0, VC$n)
      dth12.dth23.st <- rep(0, VC$n)
      dth13.dth12.st <- rep(0, VC$n)
      
      mattheta12 <- matrix(c(dth12.dth12.st,
                             dth12.dth13.st,
                             dth12.dth23.st),
                           VC$n, 3)
      
      mattheta13 <- matrix(c(dth13.dth12.st,
                             dth13.dth13.st,
                             dth13.dth23.st),
                           VC$n, 3)
      
      mattheta23 <- matrix(c(dth23.dth12.st,
                             dth23.dth13.st,
                             dth23.dth23.st),
                           VC$n, 3)
      
      d2l.dbe1.theta12.st <- crossprod(VC$X1 * (d2l.de1.theta[, 1] * mattheta12)[,1], VC$X4) +
        crossprod(VC$X1 * (d2l.de1.theta[, 2] * mattheta13)[,1], VC$X4) +
        crossprod(VC$X1 * (d2l.de1.theta[, 3] * mattheta23)[,1], VC$X4)
      
      d2l.dbe1.theta13.st <- crossprod(VC$X1 * (d2l.de1.theta[, 1] * mattheta12)[,2], VC$X5) +
        crossprod(VC$X1 * (d2l.de1.theta[, 2] * mattheta13)[,2], VC$X5) +
        crossprod(VC$X1 * (d2l.de1.theta[, 3] * mattheta23)[,2], VC$X5)
      
      d2l.dbe1.theta23.st <- crossprod(VC$X1 * (d2l.de1.theta[, 1] * mattheta12)[,3], VC$X6) +
        crossprod(VC$X1 * (d2l.de1.theta[, 2] * mattheta13)[,3], VC$X6) +
        crossprod(VC$X1 * (d2l.de1.theta[, 3] * mattheta23)[,3], VC$X6)       
      
      d2l.dbe2.theta12.st <- crossprod(VC$X2 * (d2l.de2.theta[, 1] * mattheta12)[,1], VC$X4) +
        crossprod(VC$X2 * (d2l.de2.theta[, 2] * mattheta13)[,1], VC$X4) +
        crossprod(VC$X2 * (d2l.de2.theta[, 3] * mattheta23)[,1], VC$X4)
      
      d2l.dbe2.theta13.st <- crossprod(VC$X2 * (d2l.de2.theta[, 1] * mattheta12)[,2], VC$X5) +
        crossprod(VC$X2 * (d2l.de2.theta[, 2] * mattheta13)[,2], VC$X5) +
        crossprod(VC$X2 * (d2l.de2.theta[, 3] * mattheta23)[,2], VC$X5)
      
      d2l.dbe2.theta23.st <- crossprod(VC$X2 * (d2l.de2.theta[, 1] * mattheta12)[,3], VC$X6) +
        crossprod(VC$X2 * (d2l.de2.theta[, 2] * mattheta13)[,3], VC$X6) +
        crossprod(VC$X2 * (d2l.de2.theta[, 3] * mattheta23)[,3], VC$X6) 
      
      d2l.dbe3.theta12.st <- crossprod(VC$X3 * (d2l.de3.theta[, 1] * mattheta12)[,1], VC$X4) +
        crossprod(VC$X3 * (d2l.de3.theta[, 2] * mattheta13)[,1], VC$X4) +
        crossprod(VC$X3 * (d2l.de3.theta[, 3] * mattheta23)[,1], VC$X4)
      
      d2l.dbe3.theta13.st <- crossprod(VC$X3 * (d2l.de3.theta[, 1] * mattheta12)[,2], VC$X5) +
        crossprod(VC$X3 * (d2l.de3.theta[, 2] * mattheta13)[,2], VC$X5) +
        crossprod(VC$X3 * (d2l.de3.theta[, 3] * mattheta23)[,2], VC$X5)
      
      d2l.dbe3.theta23.st <- crossprod(VC$X3 * (d2l.de3.theta[, 1] * mattheta12)[,3], VC$X6) +
        crossprod(VC$X3 * (d2l.de3.theta[, 2] * mattheta13)[,3], VC$X6) +
        crossprod(VC$X3 * (d2l.de3.theta[, 3] * mattheta23)[,3], VC$X6) 
      
      
      ############################
      ## d2l.dtheta.st.theta.st ## 
      ############################
      
      d2th12.dth13.st.th12.st <- rep(0, VC$n)
      d2th12.dth23.st.th12.st <- rep(0, VC$n)
      d2th12.dth12.st.th13.st <- rep(0, VC$n)
      d2th12.dth13.st.th13.st <- rep(0, VC$n)
      d2th12.dth23.st.th13.st <- rep(0, VC$n)
      d2th12.dth12.st.th23.st <- rep(0, VC$n)
      d2th12.dth13.st.th23.st <- rep(0, VC$n)
      d2th12.dth23.st.th23.st <- rep(0, VC$n)
      
      d2th13.dth12.st.th12.st <- rep(0, VC$n)
      d2th13.dth12.st.th23.st <- rep(0, VC$n)
      d2th13.dth13.st.th12.st <- rep(0, VC$n)
      d2th13.dth23.st.th12.st <- rep(0, VC$n)
      d2th13.dth12.st.th13.st <- rep(0, VC$n)

      
      # d2theta12.theta12.st <- matrix( c(          d2th12.dth12.st.th12.st,  d2th13.dth12.st.th12.st, d2th23.dth12.st.th12.st,
      #                                             d2th12.dth12.st.th13.st,  d2th13.dth12.st.th13.st, d2th23.dth12.st.th13.st,
      #                                             d2th12.dth12.st.th23.st,  d2th13.dth12.st.th23.st, d2th23.dth12.st.th23.st ), 3 * VC$n , 3 * VC$n)
      # 
      # d2theta13.theta13.st <- matrix( c(          d2th12.dth13.st.th12.st,  d2th13.dth13.st.th12.st, d2th23.dth13.st.th12.st,
      #                                             d2th12.dth13.st.th13.st,  d2th13.dth13.st.th13.st, d2th23.dth13.st.th13.st,
      #                                             d2th12.dth13.st.th23.st,  d2th13.dth13.st.th23.st, d2th23.dth13.st.th23.st ), 3 * VC$n , 3 * VC$n)
      # 
      # d2theta23.theta23.st <- matrix( c(          d2th12.dth23.st.th12.st,  d2th13.dth23.st.th12.st, d2th23.dth23.st.th12.st,
      #                                             d2th12.dth23.st.th13.st,  d2th13.dth23.st.th13.st, d2th23.dth23.st.th13.st,
      #                                             d2th12.dth23.st.th23.st,  d2th13.dth23.st.th23.st, d2th23.dth23.st.th23.st ), 3 * VC$n , 3 * VC$n)
      # 
      d2l.theta <- matrix( c(          VC$weights*d2l.dtheta12.theta12,  VC$weights*d2l.dtheta12.theta13, VC$weights*d2l.dtheta12.theta23,
                                       VC$weights*d2l.dtheta12.theta13,  VC$weights*d2l.dtheta13.theta13, VC$weights*d2l.dtheta13.theta23,
                                       VC$weights*d2l.dtheta12.theta23,  VC$weights*d2l.dtheta13.theta23, VC$weights*d2l.dtheta23.theta23 ), 3 * VC$n , 3 )
      
      
      ### Part A ### 
      
      dtheta.theta.st <- matrix(c(dth12.dth12.st, dth13.dth12.st, dth23.dth12.st, 
                                  dth12.dth13.st, dth13.dth13.st, dth23.dth13.st, 
                                  dth12.dth23.st, dth13.dth23.st, dth23.dth23.st), 
                                3 * VC$n, 3 )
      
      tmattheta12rep <- matrix(c(dth12.dth12.st, dth12.dth12.st, dth12.dth12.st,
                                 dth13.dth12.st, dth13.dth12.st, dth13.dth12.st,
                                 dth23.dth12.st, dth23.dth12.st, dth23.dth12.st), 3 * VC$n, 3)
      
      
      tmattheta13rep <- matrix(c(dth12.dth13.st, dth12.dth13.st, dth12.dth13.st,
                                 dth13.dth13.st, dth13.dth13.st, dth13.dth13.st,
                                 dth23.dth13.st, dth23.dth13.st, dth23.dth13.st), 3 * VC$n, 3)
      
      tmattheta23rep <- matrix(c(dth12.dth23.st, dth12.dth23.st, dth12.dth23.st,
                                 dth13.dth23.st, dth13.dth23.st, dth13.dth23.st,
                                 dth23.dth23.st, dth23.dth23.st, dth23.dth23.st), 3 * VC$n, 3)
      
      
      d2l.dtheta.st.A1 <- cbind(crossprod(rbind(VC$X4, VC$X4, VC$X4) * t(tmattheta12rep * d2l.theta)[1,], dtheta.theta.st[,1] * rbind(VC$X4, VC$X4, VC$X4)),
                                crossprod(rbind(VC$X4, VC$X4, VC$X4) * t(tmattheta12rep * d2l.theta)[1,], dtheta.theta.st[,2] * rbind(VC$X5, VC$X5, VC$X5)),
                                crossprod(rbind(VC$X4, VC$X4, VC$X4) * t(tmattheta12rep * d2l.theta)[1,], dtheta.theta.st[,3] * rbind(VC$X6, VC$X6, VC$X6))) +
        cbind(crossprod(rbind(VC$X4, VC$X4, VC$X4) * t(tmattheta12rep * d2l.theta)[2,], dtheta.theta.st[,1] * rbind(VC$X4, VC$X4, VC$X4)),
              crossprod(rbind(VC$X4, VC$X4, VC$X4) * t(tmattheta12rep * d2l.theta)[2,], dtheta.theta.st[,2] * rbind(VC$X5, VC$X5, VC$X5)),
              crossprod(rbind(VC$X4, VC$X4, VC$X4) * t(tmattheta12rep * d2l.theta)[2,], dtheta.theta.st[,3] * rbind(VC$X6, VC$X6, VC$X6))) +
        cbind(crossprod(rbind(VC$X4, VC$X4, VC$X4) * t(tmattheta12rep * d2l.theta)[3,], dtheta.theta.st[,1] * rbind(VC$X4, VC$X4, VC$X4)),
              crossprod(rbind(VC$X4, VC$X4, VC$X4) * t(tmattheta12rep * d2l.theta)[3,], dtheta.theta.st[,2] * rbind(VC$X5, VC$X5, VC$X5)),
              crossprod(rbind(VC$X4, VC$X4, VC$X4) * t(tmattheta12rep * d2l.theta)[3,], dtheta.theta.st[,3] * rbind(VC$X6, VC$X6, VC$X6)))
      
      
      d2l.dtheta.st.A2 <- cbind(crossprod(rbind(VC$X5, VC$X5, VC$X5) * t(tmattheta13rep * d2l.theta)[1,], dtheta.theta.st[,1] * rbind(VC$X4, VC$X4, VC$X4)),
                                crossprod(rbind(VC$X5, VC$X5, VC$X5) * t(tmattheta13rep * d2l.theta)[1,], dtheta.theta.st[,2] * rbind(VC$X5, VC$X5, VC$X5)),
                                crossprod(rbind(VC$X5, VC$X5, VC$X5) * t(tmattheta13rep * d2l.theta)[1,], dtheta.theta.st[,3] * rbind(VC$X6, VC$X6, VC$X6))) +
        cbind(crossprod(rbind(VC$X5, VC$X5, VC$X5) * t(tmattheta13rep * d2l.theta)[2,], dtheta.theta.st[,1] * rbind(VC$X4, VC$X4, VC$X4)),
              crossprod(rbind(VC$X5, VC$X5, VC$X5) * t(tmattheta13rep * d2l.theta)[2,], dtheta.theta.st[,2] * rbind(VC$X5, VC$X5, VC$X5)),
              crossprod(rbind(VC$X5, VC$X5, VC$X5) * t(tmattheta13rep * d2l.theta)[2,], dtheta.theta.st[,3] * rbind(VC$X6, VC$X6, VC$X6))) +
        cbind(crossprod(rbind(VC$X5, VC$X5, VC$X5) * t(tmattheta13rep * d2l.theta)[3,], dtheta.theta.st[,1] * rbind(VC$X4, VC$X4, VC$X4)),
              crossprod(rbind(VC$X5, VC$X5, VC$X5) * t(tmattheta13rep * d2l.theta)[3,], dtheta.theta.st[,2] * rbind(VC$X5, VC$X5, VC$X5)),
              crossprod(rbind(VC$X5, VC$X5, VC$X5) * t(tmattheta13rep * d2l.theta)[3,], dtheta.theta.st[,3] * rbind(VC$X6, VC$X6, VC$X6)))
      
      d2l.dtheta.st.A3 <- cbind(crossprod(rbind(VC$X6, VC$X6, VC$X6) * t(tmattheta23rep * d2l.theta)[1,], dtheta.theta.st[,1] * rbind(VC$X4, VC$X4, VC$X4)),
                                crossprod(rbind(VC$X6, VC$X6, VC$X6) * t(tmattheta23rep * d2l.theta)[1,], dtheta.theta.st[,2] * rbind(VC$X5, VC$X5, VC$X5)),
                                crossprod(rbind(VC$X6, VC$X6, VC$X6) * t(tmattheta23rep * d2l.theta)[1,], dtheta.theta.st[,3] * rbind(VC$X6, VC$X6, VC$X6))) +
        cbind(crossprod(rbind(VC$X6, VC$X6, VC$X6) * t(tmattheta23rep * d2l.theta)[2,], dtheta.theta.st[,1] * rbind(VC$X4, VC$X4, VC$X4)),
              crossprod(rbind(VC$X6, VC$X6, VC$X6) * t(tmattheta23rep * d2l.theta)[2,], dtheta.theta.st[,2] * rbind(VC$X5, VC$X5, VC$X5)),
              crossprod(rbind(VC$X6, VC$X6, VC$X6) * t(tmattheta23rep * d2l.theta)[2,], dtheta.theta.st[,3] * rbind(VC$X6, VC$X6, VC$X6))) +
        cbind(crossprod(rbind(VC$X6, VC$X6, VC$X6) * t(tmattheta23rep * d2l.theta)[3,], dtheta.theta.st[,1] * rbind(VC$X4, VC$X4, VC$X4)),
              crossprod(rbind(VC$X6, VC$X6, VC$X6) * t(tmattheta23rep * d2l.theta)[3,], dtheta.theta.st[,2] * rbind(VC$X5, VC$X5, VC$X5)),
              crossprod(rbind(VC$X6, VC$X6, VC$X6) * t(tmattheta23rep * d2l.theta)[3,], dtheta.theta.st[,3] * rbind(VC$X6, VC$X6, VC$X6)))
      
      
      
      d2l.dtheta.st.A <- rbind(d2l.dtheta.st.A1, d2l.dtheta.st.A2, d2l.dtheta.st.A3)
      
      ###################################################################################################
      
      ### Part B ### 
      
      d2theta12.theta12.st1 <- matrix( c( d2th12.dth12.st.th12.st,  d2th13.dth12.st.th12.st, d2th23.dth12.st.th12.st), 3 * VC$n , 1)
      d2theta12.theta12.st2 <- matrix( c( d2th12.dth12.st.th13.st,  d2th13.dth12.st.th13.st, d2th23.dth12.st.th13.st), 3 * VC$n , 1)
      d2theta12.theta12.st3 <- matrix( c( d2th12.dth12.st.th23.st,  d2th13.dth12.st.th23.st, d2th23.dth12.st.th23.st), 3 * VC$n , 1)
      
      d2theta13.theta13.st1 <- matrix( c( d2th12.dth13.st.th12.st,  d2th13.dth13.st.th12.st, d2th23.dth13.st.th12.st), 3 * VC$n , 1)
      d2theta13.theta13.st2 <- matrix( c( d2th12.dth13.st.th13.st,  d2th13.dth13.st.th13.st, d2th23.dth13.st.th13.st), 3 * VC$n , 1)
      d2theta13.theta13.st3 <- matrix( c( d2th12.dth13.st.th23.st,  d2th13.dth13.st.th23.st, d2th23.dth13.st.th23.st), 3 * VC$n , 1)
      
      d2theta23.theta23.st1 <- matrix( c( d2th12.dth23.st.th12.st,  d2th13.dth23.st.th12.st, d2th23.dth23.st.th12.st), 3 * VC$n , 1)
      d2theta23.theta23.st2 <- matrix( c( d2th12.dth23.st.th13.st,  d2th13.dth23.st.th13.st, d2th23.dth23.st.th13.st), 3 * VC$n , 1)
      d2theta23.theta23.st3 <- matrix( c( d2th12.dth23.st.th23.st,  d2th13.dth23.st.th23.st, d2th23.dth23.st.th23.st), 3 * VC$n , 1)
      
      
      dl.dtheta <- matrix(c ( LgTRI$dl.dtheta12,
                              LgTRI$dl.dtheta13,
                              LgTRI$dl.dtheta23), VC$n, 3)
      

      d2l.dtheta12.st <- d2l.dtheta.st.A[1:VC$X4.d2,1:VC$X4.d2] + crossprod(rbind(VC$X4, VC$X4, VC$X4) * c(d2theta12.theta12.st1), rbind(VC$X4, VC$X4, VC$X4) * c(dl.dtheta))
      d2l.dtheta13.st <- d2l.dtheta.st.A[(VC$X4.d2+1):(VC$X4.d2+VC$X5.d2),(VC$X4.d2+1):(VC$X4.d2+VC$X5.d2)] + crossprod(rbind(VC$X5, VC$X5, VC$X5) * c(d2theta13.theta13.st2), rbind(VC$X5, VC$X5, VC$X5) * c(dl.dtheta))
      d2l.dtheta23.st <- d2l.dtheta.st.A[(VC$X4.d2+VC$X5.d2+1):(VC$X4.d2+VC$X5.d2+VC$X6.d2), (VC$X4.d2+VC$X5.d2+1):(VC$X4.d2+VC$X5.d2+VC$X6.d2)] + crossprod(rbind(VC$X6, VC$X6, VC$X6) * c(d2theta23.theta23.st3), rbind(VC$X6, VC$X6, VC$X6) * c(dl.dtheta))
        
      d2l.dtheta12.st.theta13.st <- d2l.dtheta.st.A[1:VC$X4.d2, (VC$X4.d2+1):(VC$X4.d2+VC$X5.d2)] + crossprod(rbind(VC$X4, VC$X4, VC$X4) * c(d2theta13.theta13.st1), rbind(VC$X5, VC$X5, VC$X5) * c(dl.dtheta))
      d2l.dtheta12.st.theta23.st <- d2l.dtheta.st.A[1:VC$X4.d2, (VC$X4.d2+VC$X5.d2+1):(VC$X4.d2+VC$X5.d2+VC$X6.d2)] + crossprod(rbind(VC$X4, VC$X4, VC$X4) * c(d2theta23.theta23.st1), rbind(VC$X6, VC$X6, VC$X6) * c(dl.dtheta))
      d2l.dtheta13.st.theta23.st <- d2l.dtheta.st.A[(VC$X4.d2+1):(VC$X4.d2+VC$X5.d2), (VC$X4.d2+VC$X5.d2+1):(VC$X4.d2+VC$X5.d2+VC$X6.d2)] + crossprod(rbind(VC$X5, VC$X5, VC$X5) * c(d2theta23.theta23.st2), rbind(VC$X6, VC$X6, VC$X6) * c(dl.dtheta))
      
      
      
    }
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
