bcontSurvGBINROY <- function(params, respvec, VC, ps, AT = FALSE){

p1 <- p2 <- p3 <- pdf1 <- pdf2 <- pdf3 <- NA
etad <- etas1 <- etas2 <- etas3 <- l.ln <- l.par <- NULL 

rotConst1 <- rotConst2 <- 1

##################

    params1 <- params[1:VC$X1.d2]
    params2 <- params[(VC$X1.d2 + 1):(VC$X1.d2 + VC$X2.d2)] 
    params3 <- params[(VC$X1.d2 + VC$X2.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2)] 
    
    params2[VC$mono.sm.pos2] <- exp( params2[VC$mono.sm.pos2] ) 
    params3[VC$mono.sm.pos3] <- exp( params3[VC$mono.sm.pos3] ) 
    
    eta1 <- VC$X1%*%params1
    eta2 <- VC$X2%*%params2
    eta3 <- VC$X3%*%params3
    
    Xd2P <- VC$Xd2%*%params2
    indNeq2 <- as.numeric(Xd2P < 0)
    Xd2P <- ifelse(Xd2P < VC$min.dn, VC$min.dn, Xd2P) # not really needed 
    
    Xd3P <- VC$Xd3%*%params3
    indNeq3 <- as.numeric(Xd3P < 0)
    Xd3P <- ifelse(Xd3P < VC$min.dn, VC$min.dn, Xd3P)     
  
    teta.st1 <- VC$X4%*%params[(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2)]
    teta.st2 <- VC$X5%*%params[(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2 + VC$X5.d2)]
 
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
    
##################

if(VC$BivD1 %in% c("GAL180","C180","J180","G180","GAL90","C90","J90","G90","GAL270","C270","J270","G270") ) rotConst1 <- -1
if(VC$BivD2 %in% c("GAL180","C180","J180","G180","GAL90","C90","J90","G90","GAL270","C270","J270","G270") ) rotConst2 <- -1

##################

  pd1 <-  probm(eta1, VC$margins[1], bc = TRUE, only.pr = FALSE, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr) 
  pd2 <- probmS(eta2, VC$margins[2], min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr) 
  pd3 <- probmS(eta3, VC$margins[3], min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr) 

  
  p1   <- 1 - pd1$pr # Pr(Y_1 = 0) as this is the set up we follow (like in bin-cont case)
  
  derp1.dereta1      <- pd1$derp1.dereta1      # includes - sign already
  der2p1.dereta1eta1 <- pd1$der2p1.dereta1eta1
  

  p2       <- pd2$pr 
  dS2eta2  <- pd2$dS
  d2S2eta2 <- pd2$d2S 
  d3S2eta2 <- pd2$d3S
  
  p3       <- pd3$pr 
  dS2eta3  <- pd3$dS
  d2S2eta3 <- pd3$d2S 
  d3S2eta3 <- pd3$d3S  
  
  
  
  der.par2 <- der2.par2 <- params2

  der.par2[ -c( VC$mono.sm.pos2 )] <- 1
  der2.par2[-c( VC$mono.sm.pos2 )] <- 0

  der2eta2dery2b2   <- t(t(VC$Xd2)*der.par2) 
  der3eta2dery2b2b2 <- t(t(VC$Xd2)*der2.par2) 

  dereta2derb2    <- t(t(VC$X2)*der.par2)
  der2eta2derb2   <- t(t(VC$X2)*der2.par2)  
  
  
  
  der.par3 <- der2.par3 <- params3

  der.par3[ -c( VC$mono.sm.pos3 )] <- 1
  der2.par3[-c( VC$mono.sm.pos3 )] <- 0

  der2eta3dery3b3   <- t(t(VC$Xd3)*der.par3) 
  der3eta3dery3b3b3 <- t(t(VC$Xd3)*der2.par3) 

  dereta3derb3    <- t(t(VC$X3)*der.par3)
  der2eta3derb3   <- t(t(VC$X3)*der2.par3)  
  

##################

dH1.0   <- copgHs(p1[VC$inde0], p2, eta1=NULL, eta2=NULL, teta1, teta.st1, Cop1, VC$dof, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
BITS1.0 <- copgHsCont(p1[VC$inde0], p2, teta1, teta.st1, Cop1, Cont = TRUE, par2 = VC$dof, nu.st = log(VC$dof - 2))
  
  c.copula2.be1be2.0 <- dH1.0$c.copula2.be1be2
  c.copula.be1.0     <- dH1.0$c.copula.be1
  c.copula.be2.0     <- dH1.0$c.copula.be2
  c.copula.theta.0   <- dH1.0$c.copula.theta
  c.copula.thet.0    <- dH1.0$c.copula.thet 
  c.copula2.be1.0    <- dH1.0$c.copula2.be1   
  c.copula2.be2.0    <- dH1.0$c.copula2.be2 
  c.copula2.be2th.0  <- dH1.0$c.copula2.be2th
  c.copula2.be1t.0   <- dH1.0$c.copula2.be1t 
  c.copula2.be2t.0   <- dH1.0$c.copula2.be2t                           
  bit1.th2ATE.0      <- dH1.0$bit1.th2ATE
  
  der2h.derp1p1.0              <- BITS1.0$der2h.derp1p1
  der2h.derp1p2.0              <- BITS1.0$der2h.derp1p2                            
  der2h.derp1teta.0            <- BITS1.0$der2h.derp1teta   
  derteta.derteta.st.0         <- BITS1.0$derteta.derteta.st   
  der2h.derteta.teta.st.0      <- BITS1.0$der2h.derteta.teta.st   
  der2teta.derteta.stteta.st.0 <- BITS1.0$der2teta.derteta.stteta.st   
  der2h.derp2p2.0              <- BITS1.0$der2h.derp2p2   
  der2h.derp2teta.0            <- BITS1.0$der2h.derp2teta   
                                                                  
  der2h.derp1teta.st.0  <- der2h.derp1teta.0*derteta.derteta.st.0  


dH1.1   <- copgHs(p1[VC$inde1], p3, eta1=NULL, eta2=NULL, teta2, teta.st2, Cop2, VC$dof, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
BITS1.1 <- copgHsCont(p1[VC$inde1], p3, teta2, teta.st2, Cop2, Cont = TRUE, par2 = VC$dof, nu.st = log(VC$dof - 2))
  
  c.copula2.be1be2.1 <- dH1.1$c.copula2.be1be2
  c.copula.be1.1     <- dH1.1$c.copula.be1
  c.copula.be2.1     <- dH1.1$c.copula.be2
  c.copula.theta.1   <- dH1.1$c.copula.theta
  c.copula.thet.1    <- dH1.1$c.copula.thet 
  c.copula2.be1.1    <- dH1.1$c.copula2.be1   
  c.copula2.be2.1    <- dH1.1$c.copula2.be2 
  c.copula2.be2th.1  <- dH1.1$c.copula2.be2th
  c.copula2.be1t.1   <- dH1.1$c.copula2.be1t 
  c.copula2.be2t.1   <- dH1.1$c.copula2.be2t                           
  bit1.th2ATE.1      <- dH1.1$bit1.th2ATE
  
  der2h.derp1p1.1              <- BITS1.1$der2h.derp1p1
  der2h.derp1p2.1              <- BITS1.1$der2h.derp1p2                            
  der2h.derp1teta.1            <- BITS1.1$der2h.derp1teta   
  derteta.derteta.st.1         <- BITS1.1$derteta.derteta.st   
  der2h.derteta.teta.st.1      <- BITS1.1$der2h.derteta.teta.st   
  der2teta.derteta.stteta.st.1 <- BITS1.1$der2teta.derteta.stteta.st   
  der2h.derp2p2.1              <- BITS1.1$der2h.derp2p2   
  der2h.derp2teta.1            <- BITS1.1$der2h.derp2teta   
                                                                  
  der2h.derp1teta.st.1  <- der2h.derp1teta.1*derteta.derteta.st.1   
  
  
  
  
  p00.0 <- mm(BiCDF(p1[VC$inde0], p2, nC1, teta1, VC$dof), min.pr = VC$min.pr, max.pr = VC$max.pr  ) 
  p01.0 <- c.copula.be2.0*-dS2eta2*Xd2P  
  
  p00.1 <- mm(BiCDF(p1[VC$inde1], p3, nC2, teta2, VC$dof), min.pr = VC$min.pr, max.pr = VC$max.pr  ) 
  p10.1 <- p3 - p00.1 
                                              

  p01.1 <- c.copula.be2.1*-dS2eta3*Xd3P  
  p11.1 <- -dS2eta3*Xd3P - p01.1
   
#################
#################
   
# check whether I need these


p00.0 <- ifelse(p00.0 < VC$min.pr, VC$min.pr, p00.0)  
p01.0 <- ifelse(p01.0 < VC$min.pr, VC$min.pr, p01.0)
p10.1 <- ifelse(p10.1 < VC$min.pr, VC$min.pr, p10.1)
p11.1 <- ifelse(p11.1 < VC$min.pr, VC$min.pr, p11.1)


#################
#################

l.par[VC$y10.y20] <- log(p00.0[VC$y10.y20R]) 
l.par[VC$y10.y21] <- log(p01.0[VC$y10.y21R]) 
l.par[VC$y11.y30] <- log(p10.1[VC$y11.y30R]) 
l.par[VC$y11.y31] <- log(p11.1[VC$y11.y31R]) 
  
l.par <- VC$weights*l.par

res <- -sum(l.par)


#########################################################################################################################################
# SCORE
#########################################################################################################################################

p00.cop.be1.derp1eta1.0                 <-  c.copula.be1.0*derp1.dereta1[VC$inde0]
p00.cop.be1.derp1eta1.1                 <-  c.copula.be1.1*derp1.dereta1[VC$inde1]   

p01.cop2.be1be2.dSeta.etatt.derp1eta1.0 <- -c.copula2.be1be2.0*dS2eta2*Xd2P*derp1.dereta1[VC$inde0]   
p01.cop2.be1be2.dSeta.etatt.derp1eta1.1 <- -c.copula2.be1be2.1*dS2eta3*Xd3P*derp1.dereta1[VC$inde1]   
   
   
p00.cop.be2.dSeta.0 <- c(c.copula.be2.0*dS2eta2) 
p00.cop.be2.dSeta.1 <- c(c.copula.be2.1*dS2eta3)   

p01.Gbe2.bit.0      <- -c(c.copula2.be2.0*dS2eta2^2*Xd2P + c.copula.be2.0*d2S2eta2*Xd2P)*dereta2derb2 - c(c.copula.be2.0*dS2eta2)*der2eta2dery2b2
p01.Gbe2.bit.1      <- -c(c.copula2.be2.1*dS2eta3^2*Xd3P + c.copula.be2.1*d2S2eta3*Xd3P)*dereta3derb3 - c(c.copula.be2.1*dS2eta3)*der2eta3dery3b3

p01.cop2.b2th.0 <- -dS2eta2*Xd2P*c.copula2.be2th.0
p01.cop2.b2th.1 <- -dS2eta3*Xd3P*c.copula2.be2th.1


ind.dl.dbe.11.1 <- c(VC$weights[VC$inde1]*p11.1^-1) 
ind.dl.dbe.10.1 <- c(VC$weights[VC$inde1]*p10.1^-1)
ind.dl.dbe.01.0 <- c(VC$weights[VC$inde0]*p01.0^-1)
ind.dl.dbe.00.0 <- c(VC$weights[VC$inde0]*p00.0^-1)

#############

dl.dbe1 <- NA
  
dl.dbe1[VC$y11.y31] <- (ind.dl.dbe.11.1*-p01.cop2.be1be2.dSeta.etatt.derp1eta1.1)[VC$y11.y31R]
dl.dbe1[VC$y11.y30] <- (ind.dl.dbe.10.1*-p00.cop.be1.derp1eta1.1		)[VC$y11.y30R]
dl.dbe1[VC$y10.y21] <- (ind.dl.dbe.01.0*p01.cop2.be1be2.dSeta.etatt.derp1eta1.0	)[VC$y10.y21R]
dl.dbe1[VC$y10.y20] <- (ind.dl.dbe.00.0*p00.cop.be1.derp1eta1.0			)[VC$y10.y20R]

dl.dbe1 <- c(dl.dbe1)*VC$X1 
dl.dbe1 <- -colSums(dl.dbe1)

#############

dl.dbe2 <- rbind(  (ind.dl.dbe.01.0*p01.Gbe2.bit.0)[VC$y10.y21R, ],
                   (ind.dl.dbe.00.0*p00.cop.be2.dSeta.0*dereta2derb2)[VC$y10.y20R, ]   
                )  
                
dl.dbe2 <- -colSums(dl.dbe2)

 
#############

dl.dbe3 <- rbind(  (ind.dl.dbe.11.1*( -c(d2S2eta3*Xd3P)*dereta3derb3 - c(dS2eta3)*der2eta3dery3b3 - p01.Gbe2.bit.1) )[VC$y11.y31R, ],
                   (ind.dl.dbe.10.1*c(dS2eta3 - p00.cop.be2.dSeta.1)*dereta3derb3      )[VC$y11.y30R, ]   
                )   
                
dl.dbe3 <- -colSums(dl.dbe3)

   
#############
  
dl.dteta.st1 <- dl.dteta.st2 <- NA  
  
dl.dteta.st1[VC$y10.y21R] <- (ind.dl.dbe.01.0* p01.cop2.b2th.0)[VC$y10.y21R]   
dl.dteta.st1[VC$y10.y20R] <- (ind.dl.dbe.00.0* c.copula.theta.0)[VC$y10.y20R] 
 
dl.dteta.st1 <- dl.dteta.st1*VC$X4
   
dl.dteta.st2[VC$y11.y31R] <- (ind.dl.dbe.11.1*-p01.cop2.b2th.1)[VC$y11.y31R] 
dl.dteta.st2[VC$y11.y30R] <- (ind.dl.dbe.10.1*-c.copula.theta.1)[VC$y11.y30R]
   
dl.dteta.st2 <- dl.dteta.st2*VC$X5


dl.dteta.st1 <- -colSums( dl.dteta.st1)  
dl.dteta.st2 <- -colSums( dl.dteta.st2)     


##########################


G <- c( dl.dbe1, dl.dbe2, dl.dbe3, dl.dteta.st1, dl.dteta.st2 )

        
#################################################################################################
# HESSIAN
#################################################################################################

  
p00.der2.eta1.0 <-  c.copula2.be1.0*derp1.dereta1[VC$inde0]^2 + c.copula.be1.0*der2p1.dereta1eta1[VC$inde0]
p01.der2.eta1.0 <- -der2h.derp1p1.0*derp1.dereta1[VC$inde0]^2*dS2eta2*Xd2P - c.copula2.be1be2.0*der2p1.dereta1eta1[VC$inde0]*dS2eta2*Xd2P

p00.der2.eta1.1 <-  c.copula2.be1.1*derp1.dereta1[VC$inde1]^2 + c.copula.be1.1*der2p1.dereta1eta1[VC$inde1]
p01.der2.eta1.1 <- -der2h.derp1p1.1*derp1.dereta1[VC$inde1]^2*dS2eta3*Xd3P - c.copula2.be1be2.1*der2p1.dereta1eta1[VC$inde1]*dS2eta3*Xd3P

be1.be1 <- NA

be1.be1[VC$y10.y20] <- ( p00.der2.eta1.0/p00.0 - p00.cop.be1.derp1eta1.0^2/p00.0^2                   )[VC$y10.y20R]  
be1.be1[VC$y10.y21] <- ( p01.der2.eta1.0/p01.0 - p01.cop2.be1be2.dSeta.etatt.derp1eta1.0^2/p01.0^2   )[VC$y10.y21R]
be1.be1[VC$y11.y30] <- (-p00.der2.eta1.1/p10.1 - p00.cop.be1.derp1eta1.1^2/p10.1^2                   )[VC$y11.y30R]
be1.be1[VC$y11.y31] <- (-p01.der2.eta1.1/p11.1 - p01.cop2.be1be2.dSeta.etatt.derp1eta1.1^2/p11.1^2   )[VC$y11.y31R]

be1.be1 <- c(VC$weights*be1.be1) 

be1.be1 <- -crossprod(VC$X1*be1.be1, VC$X1)  
  
  
  
#######################################################################################################################################

# 
 
p00.ind1.0 <- VC$weights[VC$inde0]*p00.0^-1*as.numeric(VC$y10.y20R)
p00.ind2.0 <- VC$weights[VC$inde0]*p00.0^-2*as.numeric(VC$y10.y20R)
 
p00.der2be2be2.0   <- crossprod(c(p00.ind1.0*c.copula2.be2.0*dS2eta2^2)*dereta2derb2, dereta2derb2) + crossprod(c(p00.ind1.0*c.copula.be2.0*d2S2eta2)*dereta2derb2, dereta2derb2) + diag( colSums( t( t(c(p00.ind1.0*c.copula.be2.0*dS2eta2)*VC$X2)*der2.par2 ) ) )  
p00.der2be2be2.2.0 <- crossprod(c(p00.ind2.0*p00.cop.be2.dSeta.0^2)*dereta2derb2, dereta2derb2)
 

p01.ind1.0  <- VC$weights[VC$inde0]*p01.0^-1*as.numeric(VC$y10.y21R)

p01.der2be2be2.0 <- - crossprod(c(p01.ind1.0*der2h.derp2p2.0*dS2eta2^3*Xd2P)*dereta2derb2, dereta2derb2)          -
                    crossprod(c(p01.ind1.0*2*c.copula2.be2.0*dS2eta2*d2S2eta2*Xd2P)*dereta2derb2,dereta2derb2)  -
                    diag( colSums( t( t(c(p01.ind1.0*c.copula2.be2.0*dS2eta2^2*Xd2P)*VC$X2)*der2.par2 ) ) )     - 
                    crossprod(c(p01.ind1.0*c.copula2.be2.0*dS2eta2^2)*dereta2derb2, der2eta2dery2b2)            - #         
                    crossprod(c(p01.ind1.0*c.copula2.be2.0*dS2eta2*d2S2eta2*Xd2P)*dereta2derb2, dereta2derb2)   - 
                    crossprod(c(p01.ind1.0*c.copula.be2.0*d3S2eta2*Xd2P)*dereta2derb2, dereta2derb2)            -
                    diag( colSums( t( t(c(p01.ind1.0*c.copula.be2.0*d2S2eta2*Xd2P)*VC$X2)*der2.par2 ) ) )       -
                    crossprod(c(p01.ind1.0*c.copula.be2.0*d2S2eta2)*dereta2derb2, der2eta2dery2b2)              - #
                    crossprod(c(p01.ind1.0*c.copula2.be2.0*dS2eta2^2)*der2eta2dery2b2, dereta2derb2)            - # this appears twice but crossprod swapped
                    crossprod(c(p01.ind1.0*c.copula.be2.0*d2S2eta2)*der2eta2dery2b2, dereta2derb2)              - # this appears twice but crossprod swapped
                    diag( colSums( t( t(c(p01.ind1.0*c.copula.be2.0*dS2eta2)*VC$Xd2)*der2.par2 ) ) )  
                    
p01.Gbe2.bit2.0 <- c(ind.dl.dbe.01.0/VC$weights[VC$inde0])*as.numeric(VC$y10.y21R)*p01.Gbe2.bit.0
p01.Gbe2.bit2.0 <- crossprod(VC$weights[VC$inde0]*p01.Gbe2.bit2.0, p01.Gbe2.bit2.0)

                           
be2.be2 <- p01.der2be2be2.0 - p01.Gbe2.bit2.0 +  
           p00.der2be2be2.0 - p00.der2be2be2.2.0   
 
 
be2.be2 <- -be2.be2 


#######################################################################################################################################



p10.ind1.1 <- VC$weights[VC$inde1]*p10.1^-1*as.numeric(VC$y11.y30R)
p10.ind2.1 <- VC$weights[VC$inde1]*p10.1^-2*as.numeric(VC$y11.y30R)

p00.der2be2be2.10.1 <- crossprod(c(p10.ind1.1*c.copula2.be2.1*dS2eta3^2)*dereta3derb3, dereta3derb3) + crossprod(c(p10.ind1.1*c.copula.be2.1*d2S2eta3)*dereta3derb3, dereta3derb3) + diag( colSums( t( t(c(p10.ind1.1*c.copula.be2.1*dS2eta3)*VC$X3)*der2.par3 ) ) )  
p10.der2be2be2.10.1 <- crossprod(c(p10.ind2.1*c(dS2eta3 - p00.cop.be2.dSeta.1)^2)*dereta3derb3, dereta3derb3)

#

p11.ind1.1 <- VC$weights[VC$inde1]*p11.1^-1*as.numeric(VC$y11.y31R)
 
p01.der2be2be2.11.1 <- - crossprod(c(p11.ind1.1*der2h.derp2p2.1*dS2eta3^3*Xd3P)*dereta3derb3, dereta3derb3)          -
                       crossprod(c(p11.ind1.1*2*c.copula2.be2.1*dS2eta3*d2S2eta3*Xd3P)*dereta3derb3,dereta3derb3)  -
                       diag( colSums( t( t(c(p11.ind1.1*c.copula2.be2.1*dS2eta3^2*Xd3P)*VC$X3)*der2.par3 ) ) )     - 
                       crossprod(c(p11.ind1.1*c.copula2.be2.1*dS2eta3^2)*dereta3derb3, der2eta3dery3b3)            - #         
                       crossprod(c(p11.ind1.1*c.copula2.be2.1*dS2eta3*d2S2eta3*Xd3P)*dereta3derb3, dereta3derb3)   - 
                       crossprod(c(p11.ind1.1*c.copula.be2.1*d3S2eta3*Xd3P)*dereta3derb3, dereta3derb3)            -
                       diag( colSums( t( t(c(p11.ind1.1*c.copula.be2.1*d2S2eta3*Xd3P)*VC$X3)*der2.par3 ) ) )       -
                       crossprod(c(p11.ind1.1*c.copula.be2.1*d2S2eta3)*dereta3derb3, der2eta3dery3b3)              - #
                       crossprod(c(p11.ind1.1*c.copula2.be2.1*dS2eta3^2)*der2eta3dery3b3, dereta3derb3)            - # this appears twice but crossprod swapped
                       crossprod(c(p11.ind1.1*c.copula.be2.1*d2S2eta3)*der2eta3dery3b3, dereta3derb3)              - # this appears twice but crossprod swapped
                       diag( colSums( t( t(c(p11.ind1.1*c.copula.be2.1*dS2eta3)*VC$Xd3)*der2.par3 ) ) )  
 
 
P11.BIT1.1 <- -crossprod(c(p11.ind1.1*d3S2eta3*Xd3P)*dereta3derb3, dereta3derb3) - 
             diag( colSums( t( t(c(p11.ind1.1*d2S2eta3*Xd3P)*VC$X3)*der2.par3 ) ) ) - 
             crossprod(c(p11.ind1.1*d2S2eta3)*dereta3derb3, der2eta3dery3b3) - 
             crossprod(c(p11.ind1.1*d2S2eta3)*der2eta3dery3b3, dereta3derb3) - # swap here as well
             diag( colSums( t( t(c(p11.ind1.1*dS2eta3)*VC$Xd3)*der2.par3 ) ) ) - 
             p01.der2be2be2.11.1 
 
 
p01.Gbe2.bit3.1 <- c(ind.dl.dbe.11.1/VC$weights[VC$inde1])*as.numeric(VC$y11.y31R)*(-c(d2S2eta3*Xd3P)*dereta3derb3 - c(dS2eta3)*der2eta3dery3b3 - p01.Gbe2.bit.1) 
 
p01.Gbe2.bit3.1 <- crossprod(VC$weights[VC$inde1]*p01.Gbe2.bit3.1, p01.Gbe2.bit3.1)  
 
                                 
be3.be3 <- P11.BIT1.1 - p01.Gbe2.bit3.1 +  
           crossprod(c(p10.ind1.1*d2S2eta3)*dereta3derb3, dereta3derb3) + 
           diag( colSums( t( t(c(p10.ind1.1*dS2eta3)*VC$X3)*der2.par3 ) ) ) - 
           p00.der2be2be2.10.1 - p10.der2be2be2.10.1 
   
be3.be3 <- -be3.be3 

#######################################################################################################################################


p01.der2.derthstar.0 <- -der2h.derteta.teta.st.0*derteta.derteta.st.0^2*dS2eta2*Xd2P - rotConst1*c.copula2.be2t.0*der2teta.derteta.stteta.st.0*dS2eta2*Xd2P
p00.der2.derthstar.0 <- bit1.th2ATE.0*derteta.derteta.st.0^2 + rotConst1*c.copula.thet.0*der2teta.derteta.stteta.st.0

th1.th1 <- NA
th1.th1[VC$y10.y20R] <- (  p00.der2.derthstar.0*p00.0^-1 - c.copula.theta.0^2*p00.0^-2  )[VC$y10.y20R] 
th1.th1[VC$y10.y21R] <- (  p01.der2.derthstar.0*p01.0^-1 -  p01.cop2.b2th.0^2*p01.0^-2  )[VC$y10.y21R]

th1.th1 <- -crossprod(VC$X4*c(VC$weights[VC$inde0]*th1.th1), VC$X4)  


####################


p00.der2.derthstar.1 <- bit1.th2ATE.1*derteta.derteta.st.1^2 + rotConst2*c.copula.thet.1*der2teta.derteta.stteta.st.1
p01.der2.derthstar.1 <- -der2h.derteta.teta.st.1*derteta.derteta.st.1^2*dS2eta3*Xd3P - rotConst2*c.copula2.be2t.1*der2teta.derteta.stteta.st.1*dS2eta3*Xd3P

th2.th2 <- NA
th2.th2[VC$y11.y30R] <- ( -p00.der2.derthstar.1*p10.1^-1 - c.copula.theta.1^2*p10.1^-2 )[VC$y11.y30R]
th2.th2[VC$y11.y31R] <- ( -p01.der2.derthstar.1*p11.1^-1 -  p01.cop2.b2th.1^2*p11.1^-2 )[VC$y11.y31R]

th2.th2 <- -crossprod(VC$X5*c(VC$weights[VC$inde1]*th2.th2), VC$X5) 



#######################################################################################################################################
 
p00.derbe1b2.B1.0 <- c(c.copula2.be1be2.0*dS2eta2*derp1.dereta1[VC$inde0])*dereta2derb2 
p00.derbe1b2.B2.0 <- c(p00.cop.be1.derp1eta1.0)*dereta2derb2 
 
p01.derbe1b2.B1.0 <- c(-der2h.derp1p2.0*dS2eta2^2*Xd2P*derp1.dereta1[VC$inde0])*dereta2derb2 
p01.derbe1b2.B2.0 <- c(c.copula2.be1be2.0*derp1.dereta1[VC$inde0]*d2S2eta2*Xd2P)*dereta2derb2
p01.derbe1b2.B3.0 <- c(c.copula2.be1be2.0*derp1.dereta1[VC$inde0]*dS2eta2)*der2eta2dery2b2
 
be1.be2 <- as.numeric(VC$y10.y21R)*( c( p01.0^-1)*p01.derbe1b2.B1.0 - c(p01.0^-1)*p01.derbe1b2.B2.0 - c(p01.0^-1)*p01.derbe1b2.B3.0 - p01.Gbe2.bit.0*c(p01.0^-2*p01.cop2.be1be2.dSeta.etatt.derp1eta1.0)  ) + 
           as.numeric(VC$y10.y20R)*( c( p00.0^-1)*p00.derbe1b2.B1.0 - c(p00.0^-2*p00.cop.be2.dSeta.0)*p00.derbe1b2.B2.0 )  

be1.be2 <- -crossprod(VC$X1[VC$inde0, ], c(VC$weights[VC$inde0])*be1.be2) 
 
 

##########################




p00.derbe1b2.B1.1 <- c(c.copula2.be1be2.1*dS2eta3*derp1.dereta1[VC$inde1])*dereta3derb3 
p00.derbe1b2.B2.1 <- c(p00.cop.be1.derp1eta1.1)*dereta3derb3 
 
p01.derbe1b2.B1.1 <- c(-der2h.derp1p2.1*dS2eta3^2*Xd3P*derp1.dereta1[VC$inde1])*dereta3derb3 
p01.derbe1b2.B2.1 <- c(c.copula2.be1be2.1*derp1.dereta1[VC$inde1]*d2S2eta3*Xd3P)*dereta3derb3
p01.derbe1b2.B3.1 <- c(c.copula2.be1be2.1*derp1.dereta1[VC$inde1]*dS2eta3)*der2eta3dery3b3
 
be1.be3 <- as.numeric(VC$y11.y31R)*( c( p11.1^-2)*(-(c(p11.1)*p01.derbe1b2.B1.1 - c(p11.1)*p01.derbe1b2.B2.1 - c(p11.1)*p01.derbe1b2.B3.1) - c(-p01.cop2.be1be2.dSeta.etatt.derp1eta1.1)*(-c(d2S2eta3*Xd3P)*dereta3derb3 - c(dS2eta3)*der2eta3dery3b3 - p01.Gbe2.bit.1) )  ) +
           as.numeric(VC$y11.y30R)*( c(-p10.1^-1)*p00.derbe1b2.B1.1 - c(p10.1^-2*-(dS2eta3 - p00.cop.be2.dSeta.1))*p00.derbe1b2.B2.1 )

be1.be3 <- -crossprod(VC$X1[VC$inde1, ], c(VC$weights[VC$inde1])*be1.be3) 
 


############################

p00.der2.eta1.th.0 <- rotConst1*c.copula2.be1t.0*derteta.derteta.st.0*derp1.dereta1[VC$inde0] 
p01.der2.eta1.th.0 <- -der2h.derp1teta.0*derteta.derteta.st.0*derp1.dereta1[VC$inde0]*dS2eta2*Xd2P


be1.th1 <- as.numeric(VC$y10.y21R)*(  p01.der2.eta1.th.0*p01.0^-1 - p01.0^-2*p01.cop2.be1be2.dSeta.etatt.derp1eta1.0*p01.cop2.b2th.0   ) + 
           as.numeric(VC$y10.y20R)*(  p00.der2.eta1.th.0*p00.0^-1 - p00.0^-2*p00.cop.be1.derp1eta1.0*c.copula.theta.0   ) 
 
be1.th1 <- -crossprod(VC$X1[VC$inde0, ]*c(VC$weights[VC$inde0]*be1.th1), VC$X4) 



############################



p00.der2.eta1.th.1 <- rotConst2*c.copula2.be1t.1*derteta.derteta.st.1*derp1.dereta1[VC$inde1]
p01.der2.eta1.th.1 <- -der2h.derp1teta.1*derteta.derteta.st.1*derp1.dereta1[VC$inde1]*dS2eta3*Xd3P


be1.th2 <- as.numeric(VC$y11.y31R)*( -p01.der2.eta1.th.1*p11.1^-1 - p11.1^-2*p01.cop2.be1be2.dSeta.etatt.derp1eta1.1*p01.cop2.b2th.1   ) +
           as.numeric(VC$y11.y30R)*( -p00.der2.eta1.th.1*p10.1^-1 - p10.1^-2*p00.cop.be1.derp1eta1.1*c.copula.theta.1   )                
 
be1.th2 <- -crossprod(VC$X1[VC$inde1, ]*c(VC$weights[VC$inde1]*be1.th2), VC$X5)




####################################################################################################################################### 


ind.r.01.0 <- c(VC$weights[VC$inde0]*as.numeric(VC$y10.y21R))
ind.r.00.0 <- c(VC$weights[VC$inde0]*as.numeric(VC$y10.y20R))

p00.be2.rho.0 <- c(rotConst1*c.copula2.be2t.0*derteta.derteta.st.0*dS2eta2)*dereta2derb2
p01.be2.rho.0 <- c(-der2h.derp2teta.0*derteta.derteta.st.0*dS2eta2^2*Xd2P)*dereta2derb2 - c(rotConst1*c.copula2.be2t.0*derteta.derteta.st.0*d2S2eta2*Xd2P)*dereta2derb2 - c(rotConst1*c.copula2.be2t.0*derteta.derteta.st.0*dS2eta2)*der2eta2dery2b2   

be2.th1 <- crossprod( ind.r.00.0*( c(p00.0^-1)*p00.be2.rho.0 - c(p00.0^-2*c.copula.theta.0*p00.cop.be2.dSeta.0)*dereta2derb2 ), VC$X4) + 
	   crossprod( ind.r.01.0*( c(p01.0^-1)*p01.be2.rho.0 - c(p01.0^-2*p01.cop2.b2th.0)*p01.Gbe2.bit.0 ), VC$X4)

be2.th1 <- -be2.th1

##########################################################



ind.r.11.1 <- c(VC$weights[VC$inde1]*as.numeric(VC$y11.y31R))
ind.r.10.1 <- c(VC$weights[VC$inde1]*as.numeric(VC$y11.y30R))

p00.be2.rho.1 <- c(rotConst2*c.copula2.be2t.1*derteta.derteta.st.1*dS2eta3)*dereta3derb3
p01.be2.rho.1 <- c(-der2h.derp2teta.1*derteta.derteta.st.1*dS2eta3^2*Xd3P)*dereta3derb3 - c(rotConst2*c.copula2.be2t.1*derteta.derteta.st.1*d2S2eta3*Xd3P)*dereta3derb3 - c(rotConst2*c.copula2.be2t.1*derteta.derteta.st.1*dS2eta3)*der2eta3dery3b3   

be3.th2 <- crossprod( ind.r.10.1*( -c(p10.1^-1)*p00.be2.rho.1), VC$X5) - crossprod( c(ind.r.10.1*p10.1^-2*(dS2eta3 - p00.cop.be2.dSeta.1)*(-c.copula.theta.1))*dereta3derb3, VC$X5) +
	   crossprod( ind.r.11.1*( c(-p11.1^-1)*p01.be2.rho.1 - c(-p01.cop2.b2th.1*p11.1^-2)*(-c(d2S2eta3*Xd3P)*dereta3derb3 - c(dS2eta3)*der2eta3dery3b3 - p01.Gbe2.bit.1) ), VC$X5)


be3.th2 <- -be3.th2

##########################################################

  
  
  be2.be3 <- matrix(0, dim(VC$X2)[2],       dim(VC$X3)[2])
  be2.th2 <- matrix(0, dim(VC$X2)[2],       dim(VC$X5)[2])  
  be3.th1 <- matrix(0, dim(VC$X3)[2],       dim(VC$X4)[2])
  th1.th2 <- matrix(0, dim(VC$X4)[2],       dim(VC$X5)[2])
  

  

  H <- rbind( cbind(   be1.be1 ,   be1.be2 ,   be1.be3 ,   be1.th1,  be1.th2 ), 
              cbind( t(be1.be2),   be2.be2 ,   be2.be3 ,   be2.th1,  be2.th2 ), 
              cbind( t(be1.be3), t(be2.be3),   be3.be3 ,   be3.th1,  be3.th2 ),
              cbind( t(be1.th1), t(be2.th1), t(be3.th1),   th1.th1,  th1.th2 ),
              cbind( t(be1.th2), t(be2.th2), t(be3.th2), t(th1.th2), th2.th2 )               
            )



# rotations done previously??? bother with that??


#########################################################################


if(VC$extra.regI == "pC") H <- regH(H, type = 1)
   
S.h  <- ps$S.h                                 # hess
S.h1 <- 0.5*crossprod(params, ps$S.h)%*%params # lik
S.h2 <- S.h%*%params                           # grad   
  
  S.res <- res
  res   <- S.res + S.h1
  G     <- G + S.h2
  H     <- H + S.h  
      
          
if(VC$extra.regI == "sED") H <- regH(H, type = 2)   
  

         list(value=res, gradient=G, hessian=H, S.h=S.h, S.h1=S.h1, S.h2=S.h2, 
              l=S.res, l.ln = l.ln, l.par=l.par, ps = ps, 
              eta1=eta1, eta2=eta2, etad=etad, etas1 = 1, etas2 = 1, 
              BivD=VC$BivD, p1 = p1, p2 = p2, p2 = p3, pdf1 = p1, pdf2 = -dS2eta2,  pdf2 = -dS2eta3,          
              dl.dbe1          = NULL,       
              dl.dbe2          = NULL,       
              dl.dteta.st      = NULL, 
              Cop1 = Cop1, Cop2 = Cop2, teta1 = teta1, teta2 = teta2,
              indNeq2 = indNeq2, indNeq3 = indNeq3) 
              
  }
  

  
