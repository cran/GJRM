bcontSurvGcont2Surv <- function(params, respvec, VC, ps, AT = FALSE){
p1 <- p2 <- pdf1 <- pdf2 <- c.copula.be2 <- c.copula.be1 <- c.copula2.be1be2 <- NA

    rotConst <- 1

    etad <- etas1 <- etas2 <- l.ln <- NULL 

    params1 <- params[1:VC$X1.d2]
    params2 <- params[(VC$X1.d2 + 1):(VC$X1.d2 + VC$X2.d2)] 

    params2[VC$mono.sm.pos2] <- exp( params2[VC$mono.sm.pos2] ) 

    eta1 <- eta.tr(VC$X1%*%params1, VC$margins[1])
    eta2 <- VC$X2%*%params2
    
    etad <- etas1 <- etas2 <- l.ln <- NULL 
    
    Xd2P <- VC$Xd2%*%params2
    Xd2P <- ifelse(Xd2P < VC$min.dn, VC$min.dn, Xd2P ) # safety thing, never used as per model definition
    
if(is.null(VC$X3)){  
    X3 <- X4 <- matrix(1, VC$n, 1)
    sigma21.st <- etas1 <- params[(VC$X1.d2 + VC$X2.d2 + 1)]
    teta.st    <- etad  <- params[(VC$X1.d2 + VC$X2.d2 + 2)]
} 
  
  
if(!is.null(VC$X3)){  
    X3 <- VC$X3; X4 <- VC$X4
    sigma21.st <- etas1 <- X3%*%params[(VC$X1.d2 + VC$X2.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2)]
    teta.st    <- etad  <- X4%*%params[(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2)]
}   
 
 

##################
## Transformations
##################
   
sstr1      <- esp.tr(sigma21.st, VC$margins[1])  
sigma21.st <- sstr1$vrb.st 
sigma21    <- sstr1$vrb 
    
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

  dHs1 <- distrHs(respvec$y1, eta1, sigma21, sigma21.st, nu = 1, nu.st = 1, margin2 = VC$margins[1], naive = FALSE, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr, left.trunc = VC$left.trunc1)
  pdf1 <- dHs1$pdf2
  p1   <- dHs1$p2
  derp1.dersigma21.st   <- dHs1$derp2.dersigma.st
  derpdf1.dereta1       <- dHs1$derpdf2.dereta2 
  derp1.dereta1         <- dHs1$derp2.dereta2
  derpdf1.dersigma21.st <- dHs1$derpdf2.dersigma2.st 

  der2pdf1.dereta1              <- dHs1$der2pdf2.dereta2
  der2pdf1.dersigma21.st2       <- dHs1$der2pdf2.dersigma2.st2
  der2p1.dereta1eta1            <- dHs1$der2p2.dereta2eta2
  der2p1.dersigma21.st2         <- dHs1$der2p2.dersigma2.st2
  der2pdf1.dereta1dersigma21.st <- dHs1$der2pdf2.dereta2dersigma2.st
  der2p1.dereta1dersigma21.st   <- dHs1$der2p2.dereta2dersigma2.st

  pd2      <- probmS(eta2, VC$margins[2], min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr) 
  p2       <- pd2$pr 
  dS2eta2  <- pd2$dS
  d2S2eta2 <- pd2$d2S
  d3S2eta2 <- pd2$d3S

##################

  if( length(teta1) != 0) dH1 <- copgHs(p1[teta.ind1], p2[teta.ind1], eta1=NULL, eta2=NULL, teta1, teta.st1, Cop1, VC$dof, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
  if( length(teta2) != 0) dH2 <- copgHs(p1[teta.ind2], p2[teta.ind2], eta1=NULL, eta2=NULL, teta2, teta.st2, Cop2, VC$dof, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
  
  c.copula2.be1be2 <- c.copula.be1 <- c.copula.be2 <- p00 <- c.copula.theta <- c.copula.thet <- bit1.th2ATE <- NA
  
  if( length(teta1) != 0){ c.copula2.be1be2[teta.ind1] <- dH1$c.copula2.be1be2
                           c.copula.be1[teta.ind1]     <- dH1$c.copula.be1
                           c.copula.theta[teta.ind1]   <- dH1$c.copula.theta
                           c.copula.thet[teta.ind1]    <- dH1$c.copula.thet 
                           bit1.th2ATE[teta.ind1]      <- dH1$bit1.th2ATE
                         }
                         
# CHECK THAT WE NEED ALL BITS ABOVE                         
  
  if( length(teta2) != 0){ c.copula2.be1be2[teta.ind2] <- dH2$c.copula2.be1be2  
                           c.copula.be1[teta.ind2]     <- dH2$c.copula.be1
                           c.copula.theta[teta.ind2]   <- dH2$c.copula.theta
                           c.copula.thet[teta.ind2]    <- dH2$c.copula.thet
                           bit1.th2ATE[teta.ind2]      <- dH2$bit1.th2ATE
                         }
 

l.par <- VC$weights*( VC$c11*( log(c.copula2.be1be2) + log(pdf1) + log(-dS2eta2) + log(Xd2P) ) + 
                      VC$c10*( log(c.copula.be1) + log(pdf1) ) 
                     )
 
res <- -sum(l.par)

##################

der.par1 <- der2.par1 <- params1; der.par2 <- der2.par2 <- params2

der.par2[-c( VC$mono.sm.pos2 )]  <- 1
der2.par2[-c( VC$mono.sm.pos2 )] <- 0

der2eta2dery2b2 <- t(t(VC$Xd2)*der.par2) 
dereta2derb2    <- t(t(VC$X2)*der.par2)





##################
##################


 if( length(teta1) != 0) BITS1 <- copgHsCont(p1[teta.ind1], p2[teta.ind1], teta1, teta.st1, Cop1, Cont = TRUE, par2 = VC$dof, nu.st = log(VC$dof - 2))
 if( length(teta2) != 0) BITS2 <- copgHsCont(p1[teta.ind2], p2[teta.ind2], teta2, teta.st2, Cop2, Cont = TRUE, par2 = VC$dof, nu.st = log(VC$dof - 2)) 
 

   der2h.derp1p1 <- NA
   if( length(teta1) != 0) der2h.derp1p1[teta.ind1]   <- BITS1$der2h.derp1p1
   if( length(teta2) != 0) der2h.derp1p1[teta.ind2]   <- BITS2$der2h.derp1p1 
   
   
   der2h.derp1p2 <- NA
   if( length(teta1) != 0) der2h.derp1p2[teta.ind1] <- BITS1$der2h.derp1p2 
   if( length(teta2) != 0) der2h.derp1p2[teta.ind2] <- BITS2$der2h.derp1p2 
   
   der2h.derp1teta            <- NA
   derteta.derteta.st         <- NA
   if( length(teta1) != 0) der2h.derp1teta[teta.ind1]    <- BITS1$der2h.derp1teta
   if( length(teta2) != 0) der2h.derp1teta[teta.ind2]    <- BITS2$der2h.derp1teta   
   if( length(teta1) != 0) derteta.derteta.st[teta.ind1] <- BITS1$derteta.derteta.st
   if( length(teta2) != 0) derteta.derteta.st[teta.ind2] <- BITS2$derteta.derteta.st 

   der2h.derp1teta.st  <- der2h.derp1teta * derteta.derteta.st # new bit
   
   derc.dersigma21.st  <- der2h.derp1p1 * derp1.dersigma21.st
   derc.dereta1        <- der2h.derp1p1 * derp1.dereta1 

   
#################
#################
   
c.copula2.be1 <- c.copula2.be2 <- c.copula2.be1th <- c.copula2.be2th <- bit1.th2 <- c.copula2.be1t <- c.copula2.be2t <- NA

# CHECK THAT WE NEED ALL this     

if( length(teta1) != 0){
 
c.copula2.be1[teta.ind1]    <- dH1$c.copula2.be1   
c.copula2.be2[teta.ind1]    <- dH1$c.copula2.be2 
c.copula2.be1th[teta.ind1]  <- dH1$c.copula2.be1th 
c.copula2.be2th[teta.ind1]  <- dH1$c.copula2.be2th
c.copula2.be1t[teta.ind1]   <- dH1$c.copula2.be1t 
c.copula2.be2t[teta.ind1]   <- dH1$c.copula2.be2t
bit1.th2[teta.ind1]         <- dH1$bit1.th2

}


if( length(teta2) != 0){

c.copula2.be1[teta.ind2]    <- dH2$c.copula2.be1   
c.copula2.be2[teta.ind2]    <- dH2$c.copula2.be2 
c.copula2.be1th[teta.ind2]  <- dH2$c.copula2.be1th 
c.copula2.be2th[teta.ind2]  <- dH2$c.copula2.be2th
c.copula2.be1t[teta.ind2]   <- dH2$c.copula2.be1t 
c.copula2.be2t[teta.ind2]   <- dH2$c.copula2.be2t
bit1.th2[teta.ind2]         <- dH2$bit1.th2

}

#################
#################
   

dl.dbe1 <- c(-VC$weights*(   
   
   VC$c11*(  derpdf1.dereta1/pdf1 + derc.dereta1/c.copula2.be1be2  ) +   
   
   VC$c10*( c.copula.be1^-1 * c.copula2.be1 * derp1.dereta1 + derpdf1.dereta1 * pdf1^-1) 
    
  ))*VC$X1
   
   
dl.dbe1 <- colSums(dl.dbe1)
   

 
dl.dbe2 <- -VC$weights*(   
   
      VC$c11*( c(c.copula2.be1be2^-1*der2h.derp1p2*dS2eta2 + dS2eta2^-1*d2S2eta2)*dereta2derb2  +  
               c(Xd2P^-1)*der2eta2dery2b2 ) 
       +
   
      VC$c10*( c(c.copula.be1^-1*c.copula2.be1be2*dS2eta2)*dereta2derb2 ) 
       
  )   
   
   
dl.dbe2 <- colSums(dl.dbe2)
   
   
 

   
dl.dsigma21.st <- c(-VC$weights*(   
   
   
   VC$c11*(  derpdf1.dersigma21.st/pdf1 + derc.dersigma21.st/c.copula2.be1be2  ) +
   
   VC$c10*( c.copula.be1^-1 * c.copula2.be1 * derp1.dersigma21.st + derpdf1.dersigma21.st * pdf1^-1   ) 
   
   
   ))*X3
    

dl.dsigma21.st <- colSums(dl.dsigma21.st) 
 
 
 
dl.dteta.st    <- c(-VC$weights*(   
   
   VC$c11*(  der2h.derp1teta.st/c.copula2.be1be2  ) +
   
   VC$c10*( c.copula.be1^-1 * c.copula2.be1th   )  
   
   ))*X4
   
   
dl.dteta.st <- colSums( dl.dteta.st) 
   



G <- c( dl.dbe1, dl.dbe2, dl.dsigma21.st, dl.dteta.st ) 

#round(as.numeric(G + VC$S.h%*%params),4)
#round(as.numeric(gr),4)
        
#################################################################################################

der2c.derrho.derrho    <- NA
der2c.derp1.derp1      <- NA
der2c.derp2.derp2      <- NA
der2c.derp1.derp2      <- NA 
der2c.derp1.derrho     <- NA
der2c.derp2.derrho     <- NA
der2teta.derteta.stteta.st <- NA 

if( length(teta1) != 0){ der2c.derrho.derrho[teta.ind1] <- BITS1$der2c.derrho.derrho
 der2c.derp1.derp1[teta.ind1]      <- BITS1$der2c.derp1.derp1  
 der2c.derp2.derp2[teta.ind1]      <- BITS1$der2c.derp2.derp2  
 der2c.derp1.derp2[teta.ind1]      <- BITS1$der2c.derp1.derp2  
 der2c.derp1.derrho[teta.ind1]     <- BITS1$der2c.derp1.derrho 
 der2c.derp2.derrho[teta.ind1]     <- BITS1$der2c.derp2.derrho }


if( length(teta2) != 0){ der2c.derrho.derrho[teta.ind2] <- BITS2$der2c.derrho.derrho
 der2c.derp1.derp1[teta.ind2]      <- BITS2$der2c.derp1.derp1  
 der2c.derp2.derp2[teta.ind2]      <- BITS2$der2c.derp2.derp2  
 der2c.derp1.derp2[teta.ind2]      <- BITS2$der2c.derp1.derp2  
 der2c.derp1.derrho[teta.ind2]     <- BITS2$der2c.derp1.derrho 
 der2c.derp2.derrho[teta.ind2]     <- BITS2$der2c.derp2.derrho }

if( length(teta1) != 0) der2teta.derteta.stteta.st[teta.ind1] <- BITS1$der2teta.derteta.stteta.st
if( length(teta2) != 0) der2teta.derteta.stteta.st[teta.ind2] <- BITS2$der2teta.derteta.stteta.st 


########################

der3C.derp1p1p1 <- der3C.derp1tetateta <- der2h.derteta.teta.st <- der3C.p1p1teta <- der2h.derp2teta <- der2h.derp2p2 <- NA

if( length(teta1) != 0){der3C.derp1p1p1[teta.ind1]       <- BITS1$der3C.derp1p1p1 
                        der2h.derteta.teta.st[teta.ind1] <- BITS1$der2h.derteta.teta.st
                        der3C.derp1tetateta[teta.ind1]   <- BITS1$der3C.derp1tetateta
                        der3C.p1p1teta[teta.ind1]        <- BITS1$der3C.p1p1teta  
                        der2h.derp2teta[teta.ind1]       <- BITS1$der2h.derp2teta
                        der2h.derp2p2[teta.ind1]         <- BITS1$der2h.derp2p2
                        der2h.derp1teta[teta.ind1]       <- BITS1$der2h.derp1teta
                       }


if( length(teta2) != 0){der3C.derp1p1p1[teta.ind2]       <- BITS2$der3C.derp1p1p1
                        der2h.derteta.teta.st[teta.ind2] <- BITS2$der2h.derteta.teta.st
                        der3C.derp1tetateta[teta.ind2]   <- BITS2$der3C.derp1tetateta
                        der3C.p1p1teta[teta.ind2]        <- BITS2$der3C.p1p1teta
                        der2h.derp2teta[teta.ind2]       <- BITS2$der2h.derp2teta
                        der2h.derp2p2[teta.ind2]         <- BITS2$der2h.derp2p2
                        der2h.derp1teta[teta.ind2]       <- BITS2$der2h.derp1teta
                       }


########################
########################                                                                                                 





           
  d2l.be1.be1      <-  -VC$weights*( 
  
   VC$c11*(  (der2pdf1.dereta1 * pdf1 - derpdf1.dereta1^2) / pdf1^2  + 
                       ((der2c.derp1.derp1 * derp1.dereta1^2 + der2h.derp1p1 * der2p1.dereta1eta1) * c.copula2.be1be2 + -derc.dereta1^2) /c.copula2.be1be2^2   ) +
   
   VC$c10*( -(c.copula.be1)^-2 * c.copula2.be1^2 * derp1.dereta1^2 + c.copula.be1^-1 * der3C.derp1p1p1 * derp1.dereta1^2 + c.copula.be1^-1 * c.copula2.be1 * der2p1.dereta1eta1 + #-der2p1.dereta1eta1 
            der2pdf1.dereta1 * pdf1^-1 +  derpdf1.dereta1^2*-pdf1^-2 ) 
  
  )             
            
be1.be1 <- crossprod(VC$X1*c(d2l.be1.be1), VC$X1)            
            
            
            
    
be2.be2 <-  -( 
  
            crossprod(VC$weights*VC$c11*c(-c.copula2.be1be2^-2*der2h.derp1p2^2*dS2eta2^2 + c.copula2.be1be2^-1*der2c.derp2.derp2*dS2eta2^2 + c.copula2.be1be2^-1*der2h.derp1p2*d2S2eta2 -dS2eta2^-2*d2S2eta2^2 + dS2eta2^-1*d3S2eta2)*dereta2derb2, dereta2derb2)  +  

            diag( colSums( t( t(VC$weights*VC$c11*c(c.copula2.be1be2^-1*der2h.derp1p2*dS2eta2 + dS2eta2^-1*d2S2eta2)*VC$X2)*der2.par2 ) ) ) +
            
            crossprod(VC$weights*VC$c11*c(-Xd2P^-2)*der2eta2dery2b2, der2eta2dery2b2)  +  

            diag( colSums( t( t(VC$weights*VC$c11*c(Xd2P^-1)*VC$Xd2)*der2.par2 ) ) )        +
         
            
                        
            crossprod(VC$weights*VC$c10*c(-c.copula.be1^-2*c.copula2.be1be2^2*dS2eta2^2 + c.copula.be1^-1*der2h.derp1p2*dS2eta2^2 + c.copula.be1^-1*c.copula2.be1be2*d2S2eta2)*dereta2derb2, dereta2derb2)  +  

            diag( colSums( t( t(VC$weights*VC$c10*c( c.copula.be1^-1*c.copula2.be1be2*dS2eta2  )*VC$X2)*der2.par2 ) ) ) 
               
               
            )


be1.be2 <- -( 

crossprod(VC$weights*VC$c11*c((-c.copula2.be1be2^-2*der2h.derp1p2*der2h.derp1p1 + c.copula2.be1be2^-1*der2c.derp1.derp2)*derp1.dereta1*dS2eta2)*VC$X1, dereta2derb2)  +  

crossprod(VC$weights*VC$c10*c((-c.copula.be1^-2*c.copula2.be1be2*c.copula2.be1 + c.copula.be1^-1*der2h.derp1p1)*derp1.dereta1*dS2eta2)*VC$X1, dereta2derb2)    

)

 
if(VC$BivD %in% c("GAL180","C180","J180","G180","GAL90","C90","J90","G90","GAL270","C270","J270","G270") ) rotConst <- -1
if(VC$BivD %in% VC$BivD2) rotConst <- VC$my.env$signind


d2l.rho.rho <- -( 

VC$weights*VC$c11*( -c.copula2.be1be2^-2*der2h.derp1teta^2*derteta.derteta.st^2 + c.copula2.be1be2^-1*der2c.derrho.derrho*derteta.derteta.st^2 + c.copula2.be1be2^-1*der2h.derp1teta*der2teta.derteta.stteta.st) +

VC$weights*VC$c10*( -c.copula.be1^-2*c.copula2.be1t^2*derteta.derteta.st^2 + c.copula.be1^-1*der3C.derp1tetateta*derteta.derteta.st^2 + rotConst*c.copula.be1^-1*c.copula2.be1t*der2teta.derteta.stteta.st ) 


) 

rho.rho <- crossprod(X4*c(d2l.rho.rho), X4)



be1.rho <- -( 

crossprod(VC$weights*VC$c11*c((-c.copula2.be1be2^-2*der2h.derp1p1*der2h.derp1teta + c.copula2.be1be2^-1*der2c.derp1.derrho)*derp1.dereta1*derteta.derteta.st)*VC$X1, X4)  +  

crossprod(VC$weights*VC$c10*c((rotConst*-c.copula.be1^-2*c.copula2.be1*c.copula2.be1t + c.copula.be1^-1*der3C.p1p1teta)*derp1.dereta1*derteta.derteta.st)*VC$X1, X4)    


)



be2.rho <- -( 

crossprod(VC$weights*VC$c11*c((-c.copula2.be1be2^-2*der2h.derp1p2*der2h.derp1teta + c.copula2.be1be2^-1*der2c.derp2.derrho)*dS2eta2*derteta.derteta.st)*dereta2derb2, X4)  +  

crossprod(VC$weights*VC$c10*c((rotConst*-c.copula.be1^-2*c.copula2.be1be2*c.copula2.be1t + c.copula.be1^-1*der2h.derp1teta)*dS2eta2*derteta.derteta.st)*dereta2derb2, X4)    

)


                            
d2l.sigma21.sigma21  <- -VC$weights*( 
  
   VC$c11*(  (der2pdf1.dersigma21.st2 * pdf1 - derpdf1.dersigma21.st^2) / pdf1^2   + 
                          ((der2c.derp1.derp1*derp1.dersigma21.st^2+der2h.derp1p1*der2p1.dersigma21.st2)*c.copula2.be1be2 - derc.dersigma21.st^2) /c.copula2.be1be2^2  ) +
   
   VC$c10*( -(c.copula.be1)^-2 * c.copula2.be1^2 * derp1.dersigma21.st^2 + c.copula.be1^-1 * der3C.derp1p1p1 * derp1.dersigma21.st^2 + c.copula.be1^-1 * c.copula2.be1 * der2p1.dersigma21.st2 +
            der2pdf1.dersigma21.st2 * pdf1^-1 +  derpdf1.dersigma21.st^2*-pdf1^-2    )  
  ) 
  
  
sigma21.sigma21 <- crossprod(X3*c(d2l.sigma21.sigma21),X3)     
  


d2l.be1.sigma21  <- -VC$weights*( 
 
   VC$c11*(  (der2pdf1.dereta1dersigma21.st * pdf1 - derpdf1.dereta1*derpdf1.dersigma21.st) / pdf1^2   + 
                       ((der2c.derp1.derp1 * derp1.dereta1* derp1.dersigma21.st+der2h.derp1p1*der2p1.dereta1dersigma21.st) * c.copula2.be1be2 - derc.dereta1*derc.dersigma21.st) /c.copula2.be1be2^2  ) +
   
   VC$c10*( -(c.copula.be1)^-2 * c.copula2.be1^2 * derp1.dereta1*derp1.dersigma21.st + c.copula.be1^-1 * der3C.derp1p1p1 * derp1.dereta1*derp1.dersigma21.st + c.copula.be1^-1 * c.copula2.be1 * der2p1.dereta1dersigma21.st +
            der2pdf1.dereta1dersigma21.st * pdf1^-1 +  derpdf1.dereta1*derpdf1.dersigma21.st*-pdf1^-2    ) 
  )
  
be1.sigma21 <- crossprod(VC$X1*c(d2l.be1.sigma21), X3)      



dl.dbe2 <- -VC$weights*(   
   
      VC$c11*( c(c.copula2.be1be2^-1*der2h.derp1p2*dS2eta2 + dS2eta2^-1*d2S2eta2)*dereta2derb2 ) 
       +
   
      VC$c10*( c(c.copula.be1^-1*c.copula2.be1be2*dS2eta2)*dereta2derb2 ) 
       
  )   


be2.sigma21 <- -(   
   
      crossprod( VC$weights*VC$c11*( c(-(c.copula2.be1be2)^-2*derc.dersigma21.st*der2h.derp1p2*dS2eta2 + c.copula2.be1be2^-1*der2c.derp1.derp2*derp1.dersigma21.st*dS2eta2 )*dereta2derb2 ), X3) 
      +
      crossprod( VC$weights*VC$c10*( c(-(c.copula.be1)^-2*c.copula2.be1*derp1.dersigma21.st*c.copula2.be1be2*dS2eta2 + c.copula.be1^-1*der2h.derp1p1*derp1.dersigma21.st*dS2eta2)*dereta2derb2 ) , X3)
       
  )  


d2l.rho.sigma21 <-  -VC$weights*( 

   VC$c11*(  (der2c.derp1.derrho*derp1.dersigma21.st*derteta.derteta.st*c.copula2.be1be2 - derc.dersigma21.st*der2h.derp1teta *derteta.derteta.st)/c.copula2.be1be2^2  ) +
   
   VC$c10*( -(c.copula.be1)^-2 * c.copula2.be1th *  c.copula2.be1 * derp1.dersigma21.st + c.copula.be1^-1 * der3C.p1p1teta *derteta.derteta.st*derp1.dersigma21.st   )  
  
  )
  
  
 rho.sigma21 <- crossprod(X3*c(d2l.rho.sigma21), X4)
  



 
  
 
    H <- rbind( cbind( be1.be1        ,   be1.be2      ,   be1.sigma21,     be1.rho    ), 
                cbind( t(be1.be2)     ,   be2.be2      ,   be2.sigma21,     be2.rho    ), 
                cbind( t(be1.sigma21) , t(be2.sigma21) ,   sigma21.sigma21, rho.sigma21),
                cbind( t(be1.rho)     ,   t(be2.rho)   ,   t(rho.sigma21),  rho.rho    ) ) 
                
########################################################################

if(VC$extra.regI == "pC") H <- regH(H, type = 1)
   
  S.h  <- ps$S.h  
  
  if( length(S.h) != 1){
  
  S.h1 <- 0.5*crossprod(params,S.h)%*%params # as(S.h, "sparseMatrix") not sure about advantage...
  S.h2 <- S.h%*%params
  
  } else S.h <- S.h1 <- S.h2 <- 0   
  
  S.res <- res
  res   <- S.res + S.h1
  G     <- G + S.h2
  H     <- H + S.h  
  
     
if(VC$extra.regI == "sED") H <- regH(H, type = 2)   
  

         list(value=res, gradient=G, hessian=H, S.h=S.h, S.h1=S.h1, S.h2=S.h2, 
              l=S.res, l.ln = l.ln, l.par=l.par, ps = ps, 
              eta1=eta1, eta2=eta2, etad=etad, etas1 = sigma21.st, etas2 = 1, 
              BivD=VC$BivD, p1 = p1, p2 = p2, pdf1 = pdf1, pdf2 = -dS2eta2,          
              c.copula.be2 = c.copula.be2,
              c.copula.be1 = c.copula.be1,
              c.copula2.be1be2 = c.copula2.be1be2, 
              dl.dbe1          = NULL,       
              dl.dbe2          = NULL,       
              dl.dteta.st      = NULL, 
              teta.ind2 = teta.ind2, teta.ind1 = teta.ind1,
              Cop1 = Cop1, Cop2 = Cop2, teta1 = teta1, teta2 = teta2) 
              
              

  }
  

  
