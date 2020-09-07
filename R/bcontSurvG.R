bcontSurvG <- function(params, respvec, VC, ps, AT = FALSE){
p1 <- p2 <- pdf1 <- pdf2 <- c.copula.be2 <- c.copula.be1 <- c.copula2.be1be2 <- NA

    monP <- monP1 <- k1 <- k2 <- 0; Veq1 <- Veq2 <- list()  
    monP2 <- matrix(0, length(params),length(params))
   
    rotConst <- 1

    params1 <- params[1:VC$X1.d2]
    params2 <- params[(VC$X1.d2 + 1):(VC$X1.d2 + VC$X2.d2)] 

    params1[VC$mono.sm.pos1] <- exp( params1[VC$mono.sm.pos1] )
    params2[VC$mono.sm.pos2] <- exp( params2[VC$mono.sm.pos2] ) 

    eta1 <- VC$X1%*%params1
    eta2 <- VC$X2%*%params2
    
    etad <- etas1 <- etas2 <- l.ln <- NULL 
    
    Xd1P <- VC$Xd1%*%params1  
    Xd2P <- VC$Xd2%*%params2
  

  
    etad <- etas1 <- etas2 <- l.ln <- NULL 
 
if( is.null(VC$X3) ){
        X3 <- matrix(1, VC$n, 1)
        teta.st <- etad <- params[(VC$X1.d2 + VC$X2.d2 + 1)]
                    }

if( !is.null(VC$X3) ){
        X3 <- VC$X3
        teta.st <- etad <- X3%*%params[(VC$X1.d2 + VC$X2.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2)]
                      }

##################

    indNeq1 <- as.numeric(Xd1P < 0)
    indNeq2 <- as.numeric(Xd2P < 0)
  
    Xd1P <- ifelse(Xd1P < VC$min.dn, VC$min.dn, Xd1P ) 
    Xd2P <- ifelse(Xd2P < VC$min.dn, VC$min.dn, Xd2P ) 

##################

    if( any(indNeq1 == TRUE) && !is.null(VC$indexTeq1) ){
        
         k1 <- VC$my.env$k1

         for(i in 1:length(VC$pos.pbeq1)){

           Veq1[[i]] <- as.numeric(diff(params1[ VC$pos.pbeq1[[i]] ]) < 0)
           monP2[ VC$pos.pbeq1[[i]], VC$pos.pbeq1[[i]] ] <- k1*t(VC$Deq1[[i]]*Veq1[[i]])%*%VC$Deq1[[i]]
 
                                       }      
         VC$my.env$k1 <- k1*2
 
   }
   
    if( any(indNeq2 == TRUE) && !is.null(VC$indexTeq2) ){
        
         k2 <- VC$my.env$k2

         for(i in 1:length(VC$pos.pbeq2)){

           Veq2[[i]] <- as.numeric(diff(params2[ VC$pos.pbeq2[[i]] ]) < 0)
           monP2[ (VC$X1.d2 + VC$pos.pbeq2[[i]]), (VC$X1.d2 + VC$pos.pbeq2[[i]]) ] <- k2*t(VC$Deq2[[i]]*Veq2[[i]])%*%VC$Deq2[[i]]
 
                                       }      
         VC$my.env$k2 <- k2*2
 
   }
   

     monP  <- 0.5*crossprod(params, monP2)%*%params 
     monP1 <- monP2%*%params 
     
         
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

if(VC$BivD %in% VC$BivD2[1:4])  teta.ind1 <- ifelse(VC$my.env$signind*teta > exp(VC$zerov), TRUE, FALSE)
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

  pd1 <- probmS(eta1, VC$margins[1], min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr) 
  pd2 <- probmS(eta2, VC$margins[2], min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr) 
  
  p1 <- pd1$pr
  p2 <- pd2$pr 
  
  dS1eta1 <- pd1$dS
  dS2eta2 <- pd2$dS
  
  d2S1eta1 <- pd1$d2S
  d2S2eta2 <- pd2$d2S
  
  d3S1eta1 <- pd1$d3S
  d3S2eta2 <- pd2$d3S

##################

  if( length(teta1) != 0) dH1 <- copgHs(p1[teta.ind1], p2[teta.ind1], eta1=NULL, eta2=NULL, teta1, teta.st1, Cop1, VC$dof, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
  if( length(teta2) != 0) dH2 <- copgHs(p1[teta.ind2], p2[teta.ind2], eta1=NULL, eta2=NULL, teta2, teta.st2, Cop2, VC$dof, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
  
  c.copula2.be1be2 <- c.copula.be1 <- c.copula.be2 <- p00 <- c.copula.theta <- c.copula.thet <- bit1.th2ATE <- NA
  
  if( length(teta1) != 0){ c.copula2.be1be2[teta.ind1] <- dH1$c.copula2.be1be2
                           c.copula.be1[teta.ind1]     <- dH1$c.copula.be1
                           c.copula.be2[teta.ind1]     <- dH1$c.copula.be2
                           c.copula.theta[teta.ind1]   <- dH1$c.copula.theta
                           c.copula.thet[teta.ind1]    <- dH1$c.copula.thet 
                           bit1.th2ATE[teta.ind1]      <- dH1$bit1.th2ATE
                           p00[teta.ind1] <- mm(BiCDF(p1[teta.ind1], p2[teta.ind1], nC1, teta1, VC$dof), min.pr = VC$min.pr, max.pr = VC$max.pr  ) 
                         }
  
  if( length(teta2) != 0){ c.copula2.be1be2[teta.ind2] <- dH2$c.copula2.be1be2  
                           c.copula.be1[teta.ind2]     <- dH2$c.copula.be1
                           c.copula.be2[teta.ind2]     <- dH2$c.copula.be2
                           c.copula.theta[teta.ind2]   <- dH2$c.copula.theta
                           c.copula.thet[teta.ind2]    <- dH2$c.copula.thet
                           bit1.th2ATE[teta.ind2]      <- dH2$bit1.th2ATE
                           p00[teta.ind2] <- mm(BiCDF(p1[teta.ind2], p2[teta.ind2], nC2, teta2, VC$dof), min.pr = VC$min.pr, max.pr = VC$max.pr  ) 
                         }
 

l.par <- VC$weights*( VC$c11*( log(c.copula2.be1be2) + log(-dS1eta1) + log(-dS2eta2) + log(Xd1P) + log(Xd2P) ) + 
                      VC$c10*( log(c.copula.be1) + log(-dS1eta1) + log(Xd1P) ) + 
                      VC$c01*( log(c.copula.be2) + log(-dS2eta2) + log(Xd2P) ) + 
                      VC$c00*log(p00) 
                     )
 
res <- -sum(l.par)

##################

der.par1 <- der2.par1 <- params1; der.par2 <- der2.par2 <- params2

der.par1[-c( VC$mono.sm.pos1 )] <- 1
der.par2[-c( VC$mono.sm.pos2 )] <- 1

der2.par1[-c( VC$mono.sm.pos1 )] <- 0
der2.par2[-c( VC$mono.sm.pos2 )] <- 0

der2eta1dery1b1 <- t(t(VC$Xd1)*der.par1)
der2eta2dery2b2 <- t(t(VC$Xd2)*der.par2) 

dereta1derb1    <- t(t(VC$X1)*der.par1)
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
   
#################
#################
   
c.copula2.be1 <- c.copula2.be2 <- c.copula2.be1th <- c.copula2.be2th <- bit1.th2 <- c.copula2.be1t <- c.copula2.be2t <- NA

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
   

dl.dbe1 <- -VC$weights*(   
   
      VC$c11*( c(c.copula2.be1be2^-1*der2h.derp1p1*dS1eta1 + dS1eta1^-1*d2S1eta1)*dereta1derb1 +  
               c(Xd1P^-1)*der2eta1dery1b1 
             ) 
       +
   
      VC$c10*( c(c.copula.be1^-1*c.copula2.be1*dS1eta1 + dS1eta1^-1*d2S1eta1)*dereta1derb1 +
               c(Xd1P^-1)*der2eta1dery1b1  
                
             ) 
       + 
   
      VC$c01*( c(c.copula.be2^-1*c.copula2.be1be2*dS1eta1)*dereta1derb1 
      
             ) 
       +
   
      VC$c00*( c( p00^-1*c.copula.be1*dS1eta1)*dereta1derb1   
      
             )
   
  )
   
   
dl.dbe1 <- colSums(dl.dbe1)
   

 
dl.dbe2 <- -VC$weights*(   
   
      VC$c11*( c(c.copula2.be1be2^-1*der2h.derp1p2*dS2eta2 + dS2eta2^-1*d2S2eta2)*dereta2derb2  +  
               c(Xd2P^-1)*der2eta2dery2b2 
             ) 
       +
   
      VC$c10*( c(c.copula.be1^-1*c.copula2.be1be2*dS2eta2)*dereta2derb2
       
             ) 
       + 
   
      VC$c01*( c(c.copula.be2^-1*c.copula2.be2*dS2eta2 + dS2eta2^-1*d2S2eta2)*dereta2derb2 +
               c(Xd2P^-1)*der2eta2dery2b2
      
             ) 
       +
   
      VC$c00*( c( p00^-1*c.copula.be2*dS2eta2)*dereta2derb2   
      
             )
   
  )   
   
   
dl.dbe2 <- colSums(dl.dbe2)
   
   
 
dl.dteta.st    <- -VC$weights*(   
   
   VC$c11*(  der2h.derp1teta.st/c.copula2.be1be2  ) +
   
   VC$c10*( c.copula.be1^-1 * c.copula2.be1th   ) + 
   
   VC$c01*( c.copula.be2^-1 * c.copula2.be2th   ) +
   
   VC$c00*(  c.copula.theta/p00  )
   
   )*X3
   
   
dl.dteta.st <- colSums( dl.dteta.st) 
   

G <- c( dl.dbe1, dl.dbe2, dl.dteta.st ) 

        
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
   
# change position of weights   
   
be1.be1 <-  -( 
  
            crossprod(VC$weights*VC$c11*c(-c.copula2.be1be2^-2*der2h.derp1p1^2*dS1eta1^2 + c.copula2.be1be2^-1*der2c.derp1.derp1*dS1eta1^2 + c.copula2.be1be2^-1*der2h.derp1p1*d2S1eta1 -dS1eta1^-2*d2S1eta1^2 + dS1eta1^-1*d3S1eta1)*dereta1derb1, dereta1derb1)  +  

            diag( colSums( t( t(VC$weights*VC$c11*c(c.copula2.be1be2^-1*der2h.derp1p1*dS1eta1 + dS1eta1^-1*d2S1eta1)*VC$X1)*der2.par1 ) ) ) +
            
            crossprod(VC$weights*VC$c11*c(-Xd1P^-2)*der2eta1dery1b1, der2eta1dery1b1)  +  

            diag( colSums( t( t(VC$weights*VC$c11*c(Xd1P^-1)*VC$Xd1)*der2.par1 ) ) )        +
         
         
         
         
            crossprod(VC$weights*VC$c10*c(-c.copula.be1^-2*c.copula2.be1^2*dS1eta1^2 + c.copula.be1^-1*der3C.derp1p1p1*dS1eta1^2 + c.copula.be1^-1*c.copula2.be1*d2S1eta1 -dS1eta1^-2*d2S1eta1^2 + dS1eta1^-1*d3S1eta1)*dereta1derb1, dereta1derb1)  +  

            diag( colSums( t( t(VC$weights*VC$c10*c(c.copula.be1^-1*c.copula2.be1*dS1eta1 + dS1eta1^-1*d2S1eta1)*VC$X1)*der2.par1 ) ) ) +
            
            crossprod(VC$weights*VC$c10*c(-Xd1P^-2)*der2eta1dery1b1, der2eta1dery1b1)  +  

            diag( colSums( t( t(VC$weights*VC$c10*c(Xd1P^-1)*VC$Xd1)*der2.par1 ) ) )        + 
       
       


            
            crossprod(VC$weights*VC$c01*c(-c.copula.be2^-2*c.copula2.be1be2^2*dS1eta1^2 + c.copula.be2^-1*der2h.derp1p1*dS1eta1^2 + c.copula.be2^-1*c.copula2.be1be2*d2S1eta1)*dereta1derb1, dereta1derb1)  +  

            diag( colSums( t( t(VC$weights*VC$c01*c( c.copula.be2^-1*c.copula2.be1be2*dS1eta1  )*VC$X1)*der2.par1 ) ) ) +
                 
       
       

            crossprod(VC$weights*VC$c00*c(-p00^-2*c.copula.be1^2*dS1eta1^2 + p00^-1*c.copula2.be1*dS1eta1^2 + p00^-1*c.copula.be1*d2S1eta1)*dereta1derb1, dereta1derb1)  +  

            diag( colSums( t( t(VC$weights*VC$c00*c( p00^-1*c.copula.be1*dS1eta1  )*VC$X1)*der2.par1 ) ) ) 
            
            
            
            )
                                


be2.be2 <-  -( 
  
            crossprod(VC$weights*VC$c11*c(-c.copula2.be1be2^-2*der2h.derp1p2^2*dS2eta2^2 + c.copula2.be1be2^-1*der2c.derp2.derp2*dS2eta2^2 + c.copula2.be1be2^-1*der2h.derp1p2*d2S2eta2 -dS2eta2^-2*d2S2eta2^2 + dS2eta2^-1*d3S2eta2)*dereta2derb2, dereta2derb2)  +  

            diag( colSums( t( t(VC$weights*VC$c11*c(c.copula2.be1be2^-1*der2h.derp1p2*dS2eta2 + dS2eta2^-1*d2S2eta2)*VC$X2)*der2.par2 ) ) ) +
            
            crossprod(VC$weights*VC$c11*c(-Xd2P^-2)*der2eta2dery2b2, der2eta2dery2b2)  +  

            diag( colSums( t( t(VC$weights*VC$c11*c(Xd2P^-1)*VC$Xd2)*der2.par2 ) ) )        +
         
            
                        
            crossprod(VC$weights*VC$c10*c(-c.copula.be1^-2*c.copula2.be1be2^2*dS2eta2^2 + c.copula.be1^-1*der2h.derp1p2*dS2eta2^2 + c.copula.be1^-1*c.copula2.be1be2*d2S2eta2)*dereta2derb2, dereta2derb2)  +  

            diag( colSums( t( t(VC$weights*VC$c10*c( c.copula.be1^-1*c.copula2.be1be2*dS2eta2  )*VC$X2)*der2.par2 ) ) ) +
               
               
                              
            crossprod(VC$weights*VC$c01*c(-c.copula.be2^-2*c.copula2.be2^2*dS2eta2^2 + c.copula.be2^-1*der2h.derp2p2*dS2eta2^2 + c.copula.be2^-1*c.copula2.be2*d2S2eta2 -dS2eta2^-2*d2S2eta2^2 + dS2eta2^-1*d3S2eta2)*dereta2derb2, dereta2derb2)  +  

            diag( colSums( t( t(VC$weights*VC$c01*c(c.copula.be2^-1*c.copula2.be2*dS2eta2 + dS2eta2^-1*d2S2eta2)*VC$X2)*der2.par2 ) ) ) +
            
            crossprod(VC$weights*VC$c01*c(-Xd2P^-2)*der2eta2dery2b2, der2eta2dery2b2)  +  

            diag( colSums( t( t(VC$weights*VC$c01*c(Xd2P^-1)*VC$Xd2)*der2.par2 ) ) )        +        
       


            crossprod(VC$weights*VC$c00*c(-p00^-2*c.copula.be2^2*dS2eta2^2 + p00^-1*c.copula2.be2*dS2eta2^2 + p00^-1*c.copula.be2*d2S2eta2)*dereta2derb2, dereta2derb2)  +  

            diag( colSums( t( t(VC$weights*VC$c00*c( p00^-1*c.copula.be2*dS2eta2  )*VC$X2)*der2.par2 ) ) ) 
            
            
            
            )



be1.be2 <- -( 

crossprod(VC$weights*VC$c11*c((-c.copula2.be1be2^-2*der2h.derp1p2*der2h.derp1p1 + c.copula2.be1be2^-1*der2c.derp1.derp2)*dS1eta1*dS2eta2)*dereta1derb1, dereta2derb2)  +  

crossprod(VC$weights*VC$c10*c((-c.copula.be1^-2*c.copula2.be1be2*c.copula2.be1 + c.copula.be1^-1*der2h.derp1p1)*dS1eta1*dS2eta2)*dereta1derb1, dereta2derb2)  +  

crossprod(VC$weights*VC$c01*c((-c.copula.be2^-2*c.copula2.be1be2*c.copula2.be2 + c.copula.be2^-1*der2h.derp1p2)*dS1eta1*dS2eta2)*dereta1derb1, dereta2derb2)  +  

crossprod(VC$weights*VC$c00*c((-p00^-2*c.copula.be2*c.copula.be1 + p00^-1*c.copula2.be1be2)*dS1eta1*dS2eta2)*dereta1derb1, dereta2derb2)   

)

 
if(VC$BivD %in% c("C180","J180","G180","C90","J90","G90","C270","J270","G270") ) rotConst <- -1
if(VC$BivD %in% VC$BivD2) rotConst <- VC$my.env$signind


d2l.rho.rho <- -( 

VC$weights*VC$c11*( -c.copula2.be1be2^-2*der2h.derp1teta^2*derteta.derteta.st^2 + c.copula2.be1be2^-1*der2c.derrho.derrho*derteta.derteta.st^2 + c.copula2.be1be2^-1*der2h.derp1teta*der2teta.derteta.stteta.st) +

VC$weights*VC$c10*( -c.copula.be1^-2*c.copula2.be1t^2*derteta.derteta.st^2 + c.copula.be1^-1*der3C.derp1tetateta*derteta.derteta.st^2 + rotConst*c.copula.be1^-1*c.copula2.be1t*der2teta.derteta.stteta.st ) +

VC$weights*VC$c01*( -c.copula.be2^-2*c.copula2.be2t^2*derteta.derteta.st^2 + c.copula.be2^-1*der2h.derteta.teta.st*derteta.derteta.st^2 + rotConst*c.copula.be2^-1*c.copula2.be2t*der2teta.derteta.stteta.st ) + 

VC$weights*VC$c00*( -p00^-2*c.copula.thet^2*derteta.derteta.st^2 + p00^-1*bit1.th2ATE*derteta.derteta.st^2 + rotConst*p00^-1*c.copula.thet*der2teta.derteta.stteta.st )

) 

rho.rho <- crossprod(X3*c(d2l.rho.rho), X3)



be1.rho <- -( 

crossprod(VC$weights*VC$c11*c((-c.copula2.be1be2^-2*der2h.derp1p1*der2h.derp1teta + c.copula2.be1be2^-1*der2c.derp1.derrho)*dS1eta1*derteta.derteta.st)*dereta1derb1, X3)  +  

crossprod(VC$weights*VC$c10*c((rotConst*-c.copula.be1^-2*c.copula2.be1*c.copula2.be1t + c.copula.be1^-1*der3C.p1p1teta)*dS1eta1*derteta.derteta.st)*dereta1derb1, X3)  +  

crossprod(VC$weights*VC$c01*c((rotConst*-c.copula.be2^-2*c.copula2.be1be2*c.copula2.be2t + c.copula.be2^-1*der2h.derp1teta)*dS1eta1*derteta.derteta.st)*dereta1derb1, X3)  +  

crossprod(VC$weights*VC$c00*c(rotConst*(-p00^-2*c.copula.be1*c.copula.thet + p00^-1*c.copula2.be1t)*dS1eta1*derteta.derteta.st)*dereta1derb1, X3)   

)



be2.rho <- -( 

crossprod(VC$weights*VC$c11*c((-c.copula2.be1be2^-2*der2h.derp1p2*der2h.derp1teta + c.copula2.be1be2^-1*der2c.derp2.derrho)*dS2eta2*derteta.derteta.st)*dereta2derb2, X3)  +  

crossprod(VC$weights*VC$c10*c((rotConst*-c.copula.be1^-2*c.copula2.be1be2*c.copula2.be1t + c.copula.be1^-1*der2h.derp1teta)*dS2eta2*derteta.derteta.st)*dereta2derb2, X3)  +  

crossprod(VC$weights*VC$c01*c((rotConst*-c.copula.be2^-2*c.copula2.be2*c.copula2.be2t + c.copula.be2^-1*der2h.derp2teta)*dS2eta2*derteta.derteta.st)*dereta2derb2, X3)  +  

crossprod(VC$weights*VC$c00*c(rotConst*(-p00^-2*c.copula.be2*c.copula.thet + p00^-1*c.copula2.be2t)*dS2eta2*derteta.derteta.st)*dereta2derb2, X3)   

)


 
  
 
    H <- rbind( cbind( be1.be1        ,   be1.be2      ,      be1.rho    ), 
                cbind( t(be1.be2)     ,   be2.be2      ,      be2.rho    ), 
                cbind( t(be1.rho)     ,   t(be2.rho)   ,      rho.rho    ) ) 
                
########################################################################

if(VC$extra.regI == "pC") H <- regH(H, type = 1)
   
S.h  <- ps$S.h + monP2                                # hess
S.h1 <- 0.5*crossprod(params, ps$S.h)%*%params + monP # lik
S.h2 <- S.h%*%params + monP1                          # grad   
  
  S.res <- res
  res   <- S.res + S.h1
  G     <- G + S.h2
  H     <- H + S.h  
      
      
      
      
if(VC$extra.regI == "sED") H <- regH(H, type = 2)   
  

         list(value=res, gradient=G, hessian=H, S.h=S.h, S.h1=S.h1, S.h2=S.h2, 
              l=S.res, l.ln = l.ln, l.par=l.par, ps = ps, 
              eta1=eta1, eta2=eta2, etad=etad, etas1 = 1, etas2 = 1, 
              BivD=VC$BivD,               p1 = p1, p2 = p2, pdf1 = -dS1eta1, pdf2 = -dS2eta2,          
              c.copula.be2 = c.copula.be2,
              c.copula.be1 = c.copula.be1,
              c.copula2.be1be2 = c.copula2.be1be2, 
              dl.dbe1          = NULL,       
              dl.dbe2          = NULL,       
              dl.dteta.st      = NULL, 
              teta.ind2 = teta.ind2, teta.ind1 = teta.ind1,
              Cop1 = Cop1, Cop2 = Cop2, teta1 = teta1, teta2 = teta2,
              indNeq1 = indNeq1, indNeq2 = indNeq2,
              Veq1 = Veq1, Veq2 = Veq2, 
              k1 = VC$my.env$k1, k2 = VC$my.env$k2, monP2 = monP2) 
              
  }
  

  
