bcontSurvGBIN <- function(params, respvec, VC, ps, AT = FALSE){

p1 <- p2 <- pdf1 <- pdf2 <- c.copula.be2 <- c.copula.be1 <- c.copula2.be1be2 <- NA
etad <- etas1 <- etas2 <- l.ln <- NULL 

    rotConst <- 1

    params1 <- params[1:VC$X1.d2]
    params2 <- params[(VC$X1.d2 + 1):(VC$X1.d2 + VC$X2.d2)] 

    params2[VC$mono.sm.pos2] <- exp( params2[VC$mono.sm.pos2] ) 

    eta1 <- VC$X1%*%params1
    eta2 <- VC$X2%*%params2
    
    etad <- etas1 <- etas2 <- l.ln <- NULL 
     
    Xd2P <- VC$Xd2%*%params2
    indNeq2 <- as.numeric(Xd2P < 0)
    Xd2P <- ifelse(Xd2P < VC$min.dn, VC$min.dn, Xd2P) 
  
  
if( is.null(VC$X3) ){
        X3 <- matrix(1, VC$n, 1)
        teta.st <- etad <- params[(VC$X1.d2 + VC$X2.d2 + 1)]
                    }

if( !is.null(VC$X3) ){
        X3 <- VC$X3
        teta.st <- etad <- X3%*%params[(VC$X1.d2 + VC$X2.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2)]
                      }

         
##################
## Transformations
##################
     
resT  <- teta.tr(VC, teta.st)

teta.st1 <- teta.st2 <- teta.st <- resT$teta.st
teta1 <- teta2 <- teta <- resT$teta 

##################

Cop1 <- Cop2 <- VC$BivD 
nC1 <- nC2 <- VC$nC

teta.ind1 <- as.logical(c(1,0,round(runif(VC$n-2))) ) 
teta.ind2 <- teta.ind1 == FALSE 


if(!(VC$BivD %in% VC$BivD2) && length(teta.st) > 1){

teta.st1 <- teta.st[teta.ind1]
teta.st2 <- teta.st[teta.ind2]

teta1 <- teta[teta.ind1]
teta2 <- teta[teta.ind2]

}

###
 
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

##################

  pd1 <-  probm(eta1, VC$margins[1], bc = TRUE, only.pr = FALSE, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr) 
  pd2 <- probmS(eta2, VC$margins[2], min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr) 
  
  p1   <- 1 - pd1$pr # Pr(Y_1 = 0) as this is the set up we follow (like in bin-cont case)
  
  derp1.dereta1      <- pd1$derp1.dereta1      # includes - sign already
  der2p1.dereta1eta1 <- pd1$der2p1.dereta1eta1
  

  p2 <- pd2$pr 
  dS2eta2  <- pd2$dS
  d2S2eta2 <- pd2$d2S 
  d3S2eta2 <- pd2$d3S

##################

  if( length(teta1) != 0) dH1 <- copgHs(p1[teta.ind1], p2[teta.ind1], eta1=NULL, eta2=NULL, teta1, teta.st1, Cop1, VC$dof, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
  if( length(teta2) != 0) dH2 <- copgHs(p1[teta.ind2], p2[teta.ind2], eta1=NULL, eta2=NULL, teta2, teta.st2, Cop2, VC$dof, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
  
  if( length(teta1) != 0) BITS1 <- copgHsCont(p1[teta.ind1], p2[teta.ind1], teta1, teta.st1, Cop1, Cont = TRUE, par2 = VC$dof, nu.st = log(VC$dof - 2))
  if( length(teta2) != 0) BITS2 <- copgHsCont(p1[teta.ind2], p2[teta.ind2], teta2, teta.st2, Cop2, Cont = TRUE, par2 = VC$dof, nu.st = log(VC$dof - 2)) 
  
  c.copula2.be1be2 <- c.copula.be1 <- c.copula.be2 <- p00 <- p10 <- c.copula.theta <- c.copula.thet <- bit1.th2ATE <- NA
  c.copula2.be1 <- c.copula2.be2 <- c.copula2.be2th <- c.copula2.be1t <- c.copula2.be2t <- NA

  
  if( length(teta1) != 0){ c.copula2.be1be2[teta.ind1] <- dH1$c.copula2.be1be2
                           c.copula.be1[teta.ind1]     <- dH1$c.copula.be1
                           c.copula.be2[teta.ind1]     <- dH1$c.copula.be2
                           c.copula.theta[teta.ind1]   <- dH1$c.copula.theta
                           c.copula.thet[teta.ind1]    <- dH1$c.copula.thet 
			   c.copula2.be1[teta.ind1]    <- dH1$c.copula2.be1   
			   c.copula2.be2[teta.ind1]    <- dH1$c.copula2.be2 
			   c.copula2.be2th[teta.ind1]  <- dH1$c.copula2.be2th
			   c.copula2.be1t[teta.ind1]   <- dH1$c.copula2.be1t 
			   c.copula2.be2t[teta.ind1]   <- dH1$c.copula2.be2t                           
                           bit1.th2ATE[teta.ind1]      <- dH1$bit1.th2ATE
                           p00[teta.ind1] <- mm(BiCDF(p1[teta.ind1], p2[teta.ind1], nC1, teta1, VC$dof), min.pr = VC$min.pr, max.pr = VC$max.pr  ) 
                           p10[teta.ind1] <- p2[teta.ind1] - p00[teta.ind1] 
                           
                         }
  
  if( length(teta2) != 0){ c.copula2.be1be2[teta.ind2] <- dH2$c.copula2.be1be2  
                           c.copula.be1[teta.ind2]     <- dH2$c.copula.be1
                           c.copula.be2[teta.ind2]     <- dH2$c.copula.be2
                           c.copula.theta[teta.ind2]   <- dH2$c.copula.theta
                           c.copula.thet[teta.ind2]    <- dH2$c.copula.thet
			   c.copula2.be1[teta.ind2]    <- dH2$c.copula2.be1   
			   c.copula2.be2[teta.ind2]    <- dH2$c.copula2.be2 
			   c.copula2.be2th[teta.ind2]  <- dH2$c.copula2.be2th
			   c.copula2.be1t[teta.ind2]   <- dH2$c.copula2.be1t 
			   c.copula2.be2t[teta.ind2]   <- dH2$c.copula2.be2t                           
                           bit1.th2ATE[teta.ind2]      <- dH2$bit1.th2ATE
                           p00[teta.ind2] <- mm(BiCDF(p1[teta.ind2], p2[teta.ind2], nC2, teta2, VC$dof), min.pr = VC$min.pr, max.pr = VC$max.pr  ) 
                           p10[teta.ind2] <- p2[teta.ind2] - p00[teta.ind2]
                         }
 

   der2h.derp1p1 <- der2h.derp1p2 <- der2h.derp1teta <- derteta.derteta.st <- der2h.derteta.teta.st <- der2teta.derteta.stteta.st <- NA
   der2h.derp2p2 <- der2h.derp2teta <- NA
   
   if( length(teta1) != 0){ 
   
                           der2h.derp1p1[teta.ind1]              <- BITS1$der2h.derp1p1
			   der2h.derp1p2[teta.ind1]              <- BITS1$der2h.derp1p2                            
			   der2h.derp1teta[teta.ind1]            <- BITS1$der2h.derp1teta   
			   derteta.derteta.st[teta.ind1]         <- BITS1$derteta.derteta.st   
			   der2h.derteta.teta.st[teta.ind1]      <- BITS1$der2h.derteta.teta.st   
			   der2teta.derteta.stteta.st[teta.ind1] <- BITS1$der2teta.derteta.stteta.st   
			   der2h.derp2p2[teta.ind1]              <- BITS1$der2h.derp2p2   
			   der2h.derp2teta[teta.ind1]            <- BITS1$der2h.derp2teta   
 
                          }
                          
                          
   if( length(teta2) != 0){ 
   
                           der2h.derp1p1[teta.ind2]              <- BITS2$der2h.derp1p1
			   der2h.derp1p2[teta.ind2]              <- BITS2$der2h.derp1p2   
			   der2h.derp1teta[teta.ind2]            <- BITS2$der2h.derp1teta    
			   derteta.derteta.st[teta.ind2]         <- BITS2$derteta.derteta.st   
			   der2h.derteta.teta.st[teta.ind2]      <- BITS2$der2h.derteta.teta.st     
			   der2teta.derteta.stteta.st[teta.ind2] <- BITS2$der2teta.derteta.stteta.st  
			   der2h.derp2p2[teta.ind2]              <- BITS2$der2h.derp2p2  
			   der2h.derp2teta[teta.ind2]            <- BITS2$der2h.derp2teta  
                           
                          }                          
   
   
der2h.derp1teta.st  <- der2h.derp1teta*derteta.derteta.st
   
#################
#################
   
p01 <- c.copula.be2*-dS2eta2*Xd2P
p11 <- -dS2eta2*Xd2P - p01

p00 <- ifelse(p00 < 0, VC$min.pr, p00)  
p01 <- ifelse(p01 < 0, VC$min.pr, p01)
p11 <- ifelse(p11 < 0, VC$min.pr, p11)
p10 <- ifelse(p10 < 0, VC$min.pr, p10)


der.par2 <- der2.par2 <- params2

der.par2[ -c( VC$mono.sm.pos2 )] <- 1
der2.par2[-c( VC$mono.sm.pos2 )] <- 0

der2eta2dery2b2   <- t(t(VC$Xd2)*der.par2) 
der3eta2dery2b2b2 <- t(t(VC$Xd2)*der2.par2) 

dereta2derb2    <- t(t(VC$X2)*der.par2)
der2eta2derb2   <- t(t(VC$X2)*der2.par2)

#################
#################

l.par <- VC$weights*( VC$c00*log(p00) +
                      VC$c01*log(p01) +    
                      VC$c11*log(p11) +   
                      VC$c10*log(p10)
                     )
 
res <- -sum(l.par)

#########################################################################################################################################
# SCORE
#########################################################################################################################################

p00.cop.be1.derp1eta1                 <-  c.copula.be1*derp1.dereta1   
p01.cop2.be1be2.dSeta.etatt.derp1eta1 <- -c.copula2.be1be2*dS2eta2*Xd2P*derp1.dereta1   
   
   
p00.cop.be2.dSeta <- c.copula.be2*dS2eta2   
p01.Gbe2.bit      <- -c(c.copula2.be2*dS2eta2^2*Xd2P + c.copula.be2*d2S2eta2*Xd2P)*dereta2derb2 - c(c.copula.be2*dS2eta2)*der2eta2dery2b2

p01.cop2.b2th <- -dS2eta2*Xd2P*c.copula2.be2th   


###

ind.dl.dbe.11 <- VC$weights*VC$c11*p11^-1
ind.dl.dbe.10 <- VC$weights*VC$c10*p10^-1
ind.dl.dbe.01 <- VC$weights*VC$c01*p01^-1
ind.dl.dbe.00 <- VC$weights*VC$c00*p00^-1

###

dl.dbe1 <- c(ind.dl.dbe.11*-p01.cop2.be1be2.dSeta.etatt.derp1eta1 +     
             ind.dl.dbe.10*-p00.cop.be1.derp1eta1                 + 
             ind.dl.dbe.01*p01.cop2.be1be2.dSeta.etatt.derp1eta1  +
             ind.dl.dbe.00*p00.cop.be1.derp1eta1 )*VC$X1 
   
dl.dbe1 <- -colSums(dl.dbe1)
   
###
 
dl.dbe2 <- c(ind.dl.dbe.11)*(-c(d2S2eta2*Xd2P)*dereta2derb2 - c(dS2eta2)*der2eta2dery2b2 - p01.Gbe2.bit) +
           c(ind.dl.dbe.10)*c(dS2eta2 - p00.cop.be2.dSeta)*dereta2derb2                                  + 
           c(ind.dl.dbe.01)*p01.Gbe2.bit                                                                 +   
           c(ind.dl.dbe.00)*c(p00.cop.be2.dSeta)*dereta2derb2  
     
dl.dbe2 <- -colSums(dl.dbe2)
   
###
    
dl.dteta.st <- c(ind.dl.dbe.11*-p01.cop2.b2th  +
    		 ind.dl.dbe.10*-c.copula.theta +    
    		 ind.dl.dbe.01* p01.cop2.b2th  +   
    		 ind.dl.dbe.00* c.copula.theta  )*X3
   
dl.dteta.st <- -colSums( dl.dteta.st)  
   


G <- c( dl.dbe1, dl.dbe2, dl.dteta.st )

        
#################################################################################################
# HESSIAN
#################################################################################################

  
p00.der2.eta1 <-  c.copula2.be1*derp1.dereta1^2 + c.copula.be1*der2p1.dereta1eta1
p01.der2.eta1 <- -der2h.derp1p1*derp1.dereta1^2*dS2eta2*Xd2P - c.copula2.be1be2*der2p1.dereta1eta1*dS2eta2*Xd2P
   
be1.be1 <- VC$c00*( p00.der2.eta1/p00 - p00.cop.be1.derp1eta1^2/p00^2)                 + 
           VC$c10*(-p00.der2.eta1/p10 - p00.cop.be1.derp1eta1^2/p10^2)                 + 
           VC$c01*( p01.der2.eta1/p01 - p01.cop2.be1be2.dSeta.etatt.derp1eta1^2/p01^2) +
           VC$c11*(-p01.der2.eta1/p11 - p01.cop2.be1be2.dSeta.etatt.derp1eta1^2/p11^2) 
  
be1.be1 <- -crossprod(VC$X1*c(VC$weights*be1.be1), VC$X1)  
  
  
#######################################################################################################################################

# 
 
p00.ind1 <- VC$weights*VC$c00*p00^-1
p00.ind2 <- VC$weights*VC$c00*p00^-2
 
p00.der2be2be2   <- crossprod(c(p00.ind1*c.copula2.be2*dS2eta2^2)*dereta2derb2, dereta2derb2) + crossprod(c(p00.ind1*c.copula.be2*d2S2eta2)*dereta2derb2, dereta2derb2) + diag( colSums( t( t(c(p00.ind1*c.copula.be2*dS2eta2)*VC$X2)*der2.par2 ) ) )  
p00.der2be2be2.2 <- crossprod(c(p00.ind2*p00.cop.be2.dSeta^2)*dereta2derb2, dereta2derb2)
 
#

p10.ind1 <- VC$weights*VC$c10*p10^-1
p10.ind2 <- VC$weights*VC$c10*p10^-2

p00.der2be2be2.10 <- crossprod(c(p10.ind1*c.copula2.be2*dS2eta2^2)*dereta2derb2, dereta2derb2) + crossprod(c(p10.ind1*c.copula.be2*d2S2eta2)*dereta2derb2, dereta2derb2) + diag( colSums( t( t(c(p10.ind1*c.copula.be2*dS2eta2)*VC$X2)*der2.par2 ) ) )  
p10.der2be2be2.10 <- crossprod(c(p10.ind2*c(dS2eta2 - p00.cop.be2.dSeta)^2)*dereta2derb2, dereta2derb2)

#

p01.ind1  <- VC$weights*VC$c01*p01^-1

p01.der2be2be2 <- - crossprod(c(p01.ind1*der2h.derp2p2*dS2eta2^3*Xd2P)*dereta2derb2, dereta2derb2)          -
                    crossprod(c(p01.ind1*2*c.copula2.be2*dS2eta2*d2S2eta2*Xd2P)*dereta2derb2,dereta2derb2)  -
                    diag( colSums( t( t(c(p01.ind1*c.copula2.be2*dS2eta2^2*Xd2P)*VC$X2)*der2.par2 ) ) )     - 
                    crossprod(c(p01.ind1*c.copula2.be2*dS2eta2^2)*dereta2derb2, der2eta2dery2b2)            - #         
                    crossprod(c(p01.ind1*c.copula2.be2*dS2eta2*d2S2eta2*Xd2P)*dereta2derb2, dereta2derb2)   - 
                    crossprod(c(p01.ind1*c.copula.be2*d3S2eta2*Xd2P)*dereta2derb2, dereta2derb2)            -
                    diag( colSums( t( t(c(p01.ind1*c.copula.be2*d2S2eta2*Xd2P)*VC$X2)*der2.par2 ) ) )       -
                    crossprod(c(p01.ind1*c.copula.be2*d2S2eta2)*dereta2derb2, der2eta2dery2b2)              - #
                    crossprod(c(p01.ind1*c.copula2.be2*dS2eta2^2)*der2eta2dery2b2, dereta2derb2)            - # this appears twice but crossprod swapped
                    crossprod(c(p01.ind1*c.copula.be2*d2S2eta2)*der2eta2dery2b2, dereta2derb2)              - # this appears twice but crossprod swapped
                    diag( colSums( t( t(c(p01.ind1*c.copula.be2*dS2eta2)*VC$Xd2)*der2.par2 ) ) )  
                    
p01.Gbe2.bit2 <- c(ind.dl.dbe.01/VC$weights)*p01.Gbe2.bit
p01.Gbe2.bit2 <- crossprod(VC$weights*p01.Gbe2.bit2, p01.Gbe2.bit2)


p11.ind1 <- VC$weights*VC$c11*p11^-1
 
p01.der2be2be2.11 <- - crossprod(c(p11.ind1*der2h.derp2p2*dS2eta2^3*Xd2P)*dereta2derb2, dereta2derb2)          -
                       crossprod(c(p11.ind1*2*c.copula2.be2*dS2eta2*d2S2eta2*Xd2P)*dereta2derb2,dereta2derb2)  -
                       diag( colSums( t( t(c(p11.ind1*c.copula2.be2*dS2eta2^2*Xd2P)*VC$X2)*der2.par2 ) ) )     - 
                       crossprod(c(p11.ind1*c.copula2.be2*dS2eta2^2)*dereta2derb2, der2eta2dery2b2)            - #         
                       crossprod(c(p11.ind1*c.copula2.be2*dS2eta2*d2S2eta2*Xd2P)*dereta2derb2, dereta2derb2)   - 
                       crossprod(c(p11.ind1*c.copula.be2*d3S2eta2*Xd2P)*dereta2derb2, dereta2derb2)            -
                       diag( colSums( t( t(c(p11.ind1*c.copula.be2*d2S2eta2*Xd2P)*VC$X2)*der2.par2 ) ) )       -
                       crossprod(c(p11.ind1*c.copula.be2*d2S2eta2)*dereta2derb2, der2eta2dery2b2)              - #
                       crossprod(c(p11.ind1*c.copula2.be2*dS2eta2^2)*der2eta2dery2b2, dereta2derb2)            - # this appears twice but crossprod swapped
                       crossprod(c(p11.ind1*c.copula.be2*d2S2eta2)*der2eta2dery2b2, dereta2derb2)              - # this appears twice but crossprod swapped
                       diag( colSums( t( t(c(p11.ind1*c.copula.be2*dS2eta2)*VC$Xd2)*der2.par2 ) ) )  
 
 
P11.BIT1 <- -crossprod(c(p11.ind1*d3S2eta2*Xd2P)*dereta2derb2, dereta2derb2) - 
             diag( colSums( t( t(c(p11.ind1*d2S2eta2*Xd2P)*VC$X2)*der2.par2 ) ) ) - 
             crossprod(c(p11.ind1*d2S2eta2)*dereta2derb2, der2eta2dery2b2) - 
             crossprod(c(p11.ind1*d2S2eta2)*der2eta2dery2b2, dereta2derb2) - # swap here as well
             diag( colSums( t( t(c(p11.ind1*dS2eta2)*VC$Xd2)*der2.par2 ) ) ) - 
             p01.der2be2be2.11 
 
 
p01.Gbe2.bit3 <- c(ind.dl.dbe.11/VC$weights)*(-c(d2S2eta2*Xd2P)*dereta2derb2 - c(dS2eta2)*der2eta2dery2b2 - p01.Gbe2.bit) 
 
p01.Gbe2.bit3 <- crossprod(VC$weights*p01.Gbe2.bit3, p01.Gbe2.bit3)  
 
                                 
be2.be2 <- P11.BIT1 - p01.Gbe2.bit3 +
           p01.der2be2be2 - p01.Gbe2.bit2 +  
           crossprod(c(p10.ind1*d2S2eta2)*dereta2derb2, dereta2derb2) + diag( colSums( t( t(c(p10.ind1*dS2eta2)*VC$X2)*der2.par2 ) ) ) - p00.der2be2be2.10 - p10.der2be2be2.10  + 
           p00.der2be2be2 - p00.der2be2be2.2   
 
 
be2.be2 <- -be2.be2 

############################################

#######################################################################################################################################
 
p00.derbe1b2.B1 <- c(c.copula2.be1be2*dS2eta2*derp1.dereta1)*dereta2derb2 
p00.derbe1b2.B2 <- c(p00.cop.be1.derp1eta1)*dereta2derb2 
 
p01.derbe1b2.B1 <- c(-der2h.derp1p2*dS2eta2^2*Xd2P*derp1.dereta1)*dereta2derb2 
p01.derbe1b2.B2 <- c(c.copula2.be1be2*derp1.dereta1*d2S2eta2*Xd2P)*dereta2derb2
p01.derbe1b2.B3 <- c(c.copula2.be1be2*derp1.dereta1*dS2eta2)*der2eta2dery2b2
 
be1.be2 <- VC$c11*( c( p11^-2)*(-(c(p11)*p01.derbe1b2.B1 - c(p11)*p01.derbe1b2.B2 - c(p11)*p01.derbe1b2.B3) - c(-p01.cop2.be1be2.dSeta.etatt.derp1eta1)*(-c(d2S2eta2*Xd2P)*dereta2derb2 - c(dS2eta2)*der2eta2dery2b2 - p01.Gbe2.bit) )  ) +
           VC$c01*( c( p01^-1)*p01.derbe1b2.B1 - c(p01^-1)*p01.derbe1b2.B2 - c(p01^-1)*p01.derbe1b2.B3 - p01.Gbe2.bit*c(p01^-2*p01.cop2.be1be2.dSeta.etatt.derp1eta1)  ) + 
           VC$c10*( c(-p10^-1)*p00.derbe1b2.B1 - c(p10^-2*-(dS2eta2 - p00.cop.be2.dSeta))*p00.derbe1b2.B2 ) +
           VC$c00*( c( p00^-1)*p00.derbe1b2.B1 - c(p00^-2*p00.cop.be2.dSeta)*p00.derbe1b2.B2 )  

be1.be2 <- -crossprod(VC$X1, c(VC$weights)*be1.be2) 
 

  
#######################################################################################################################################

if(VC$BivD %in% c("GAL180","C180","J180","G180","GAL90","C90","J90","G90","GAL270","C270","J270","G270") ) rotConst <- -1
if(VC$BivD %in% VC$BivD2) rotConst <- VC$my.env$signind

#######################################################################################################################################

p00.der2.derthstar <- bit1.th2ATE*derteta.derteta.st^2 + rotConst*c.copula.thet*der2teta.derteta.stteta.st
p01.der2.derthstar <- -der2h.derteta.teta.st*derteta.derteta.st^2*dS2eta2*Xd2P - rotConst*c.copula2.be2t*der2teta.derteta.stteta.st*dS2eta2*Xd2P

rho.rho <- VC$c11*( -p01.der2.derthstar*p11^-1 -  p01.cop2.b2th^2*p11^-2  ) +
           VC$c01*(  p01.der2.derthstar*p01^-1 -  p01.cop2.b2th^2*p01^-2  ) + 
           VC$c10*( -p00.der2.derthstar*p10^-1 - c.copula.theta^2*p10^-2  ) +
           VC$c00*(  p00.der2.derthstar*p00^-1 - c.copula.theta^2*p00^-2  ) 

rho.rho <- -crossprod(X3*c(VC$weights*rho.rho), X3) 
  
#######################################################################################################################################
 
p00.der2.eta1.th <- rotConst*c.copula2.be1t*derteta.derteta.st*derp1.dereta1 
p01.der2.eta1.th <- -der2h.derp1teta*derteta.derteta.st*derp1.dereta1*dS2eta2*Xd2P


be1.rho <- VC$c11*( -p01.der2.eta1.th*p11^-1 - p11^-2*p01.cop2.be1be2.dSeta.etatt.derp1eta1*p01.cop2.b2th   ) +
           VC$c01*(  p01.der2.eta1.th*p01^-1 - p01^-2*p01.cop2.be1be2.dSeta.etatt.derp1eta1*p01.cop2.b2th   ) + 
           VC$c10*( -p00.der2.eta1.th*p10^-1 - p10^-2*p00.cop.be1.derp1eta1*c.copula.theta   )                +
           VC$c00*(  p00.der2.eta1.th*p00^-1 - p00^-2*p00.cop.be1.derp1eta1*c.copula.theta   ) 
 
be1.rho <- -crossprod(VC$X1*c(VC$weights*be1.rho), X3) 

####################################################################################################################################### 


ind.r.11 <- c(VC$weights*VC$c11)
ind.r.10 <- c(VC$weights*VC$c10)
ind.r.01 <- c(VC$weights*VC$c01)
ind.r.00 <- c(VC$weights*VC$c00)

p00.be2.rho <- c(rotConst*c.copula2.be2t*derteta.derteta.st*dS2eta2)*dereta2derb2
p01.be2.rho <- c(-der2h.derp2teta*derteta.derteta.st*dS2eta2^2*Xd2P)*dereta2derb2 - c(rotConst*c.copula2.be2t*derteta.derteta.st*d2S2eta2*Xd2P)*dereta2derb2 - c(rotConst*c.copula2.be2t*derteta.derteta.st*dS2eta2)*der2eta2dery2b2   

be2.rho <- crossprod( ind.r.00*( c(p00^-1)*p00.be2.rho - c(p00^-2*c.copula.theta*p00.cop.be2.dSeta)*dereta2derb2 ), X3) + 
	   crossprod( ind.r.10*( -c(p10^-1)*p00.be2.rho), X3) - crossprod( c(ind.r.10*p10^-2*(dS2eta2 - p00.cop.be2.dSeta)*(-c.copula.theta))*dereta2derb2, X3) +
	   crossprod( ind.r.01*( c(p01^-1)*p01.be2.rho - c(p01^-2*p01.cop2.b2th)*p01.Gbe2.bit ), X3) +
	   crossprod( ind.r.11*( c(-p11^-1)*p01.be2.rho - c(-p01.cop2.b2th*p11^-2)*(-c(d2S2eta2*Xd2P)*dereta2derb2 - c(dS2eta2)*der2eta2dery2b2 - p01.Gbe2.bit) ), X3)

be2.rho <- -be2.rho

##########################################################

    H <- rbind( cbind( be1.be1        ,   be1.be2      ,      be1.rho    ), 
                cbind( t(be1.be2)     ,   be2.be2      ,      be2.rho    ), 
                cbind( t(be1.rho)     ,   t(be2.rho)   ,      rho.rho    ) ) 
                
########################################################################

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
              BivD=VC$BivD, p1 = p1, p2 = p2, pdf1 = p1, pdf2 = -dS2eta2,          
              c.copula.be2 = c.copula.be2,
              c.copula.be1 = c.copula.be1,
              c.copula2.be1be2 = c.copula2.be1be2, 
              dl.dbe1          = NULL,       
              dl.dbe2          = NULL,       
              dl.dteta.st      = NULL, 
              teta.ind2 = teta.ind2, teta.ind1 = teta.ind1,
              Cop1 = Cop1, Cop2 = Cop2, teta1 = teta1, teta2 = teta2,
              indNeq2 = indNeq2,
              k1 = VC$my.env$k1, k2 = VC$my.env$k2) 
              
  }
  

  
