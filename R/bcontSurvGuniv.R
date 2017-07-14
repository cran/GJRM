bcontSurvGuniv <- function(params, respvec, VC, ps, AT = FALSE){

    params1 <- params[1:VC$X1.d2]
    params1[VC$mono.sm.pos] <- exp( params1[VC$mono.sm.pos] )

    eta1 <- VC$X1%*%params1
    Xd1P <- VC$Xd1%*%params1  
    
    Xd1P <- ifelse(Xd1P < 1e-08, 1e-08, Xd1P ) # safety thing, never used as per model definition

    etad <- etas1 <- etas2 <- l.ln <- NULL 
  
##################
    
pd1  <- probmS(eta1, VC$margins[1]) 
  
p1       <- pd1$pr
dS1eta1  <- pd1$dS
d2S1eta1 <- pd1$d2S
d3S1eta1 <- pd1$d3S

##################


l.par <- VC$weights*( VC$cens*( log(-dS1eta1) + log(Xd1P) ) + (1 - VC$cens)*log(p1) )    
res   <- -sum(l.par)

##################

der.par1 <- der2.par1 <- params1

 der.par1[-c( VC$mono.sm.pos )] <- 1
der2.par1[-c( VC$mono.sm.pos )] <- 0

der2eta1dery1b1 <- t(t(VC$Xd1)*der.par1)
dereta1derb1    <- t(t(VC$X1)*der.par1)

##################

dl.dbe1 <- -VC$weights*(   
   
   VC$cens*( c((dS1eta1*Xd1P)^-1)*(c(d2S1eta1*Xd1P)*dereta1derb1 + c(dS1eta1)*der2eta1dery1b1) ) + (1 - VC$cens)*c(p1^-1*dS1eta1)*dereta1derb1
   
  )
   
G <- colSums(dl.dbe1)
   
########################
########################                                                                                                 
   
   
H <-  -(   
   
   crossprod(c(VC$weights*VC$cens*(-dS1eta1^-2*d2S1eta1^2 + dS1eta1^-1*d3S1eta1))*dereta1derb1, dereta1derb1) +
    
   diag( colSums( t( t(  c(VC$weights*VC$cens*dS1eta1^-1*d2S1eta1)*VC$X1)*der2.par1 ) ) ) + 
   
   diag( colSums( t( t(VC$weights*VC$cens*c(Xd1P^-1)*VC$Xd1)*der2.par1 ) ) )  +
    
   crossprod(VC$weights*VC$cens*c(-Xd1P^-2)*der2eta1dery1b1, der2eta1dery1b1) +  
    
   crossprod(c(VC$weights*(1 - VC$cens)*(-p1^-2*dS1eta1^2+p1^-1*d2S1eta1))*dereta1derb1, dereta1derb1) +
   
   diag( colSums( t( t(c(VC$weights*(1 - VC$cens)*p1^-1*dS1eta1)*VC$X1)*der2.par1 ) ) ) 


  )
  
  
                
########################################################################

if(VC$extra.regI == "pC") H <- regH(H, type = 1)
   
  S.h  <- ps$S.h  
  
  if( length(S.h) != 1){
  
  S.h1 <- 0.5*crossprod(params,S.h)%*%params 
  S.h2 <- S.h%*%params
  
  } else S.h <- S.h1 <- S.h2 <- 0   
  
  S.res <- res
  res   <- S.res + S.h1
  G     <- G + S.h2
  H     <- H + S.h  
            
if(VC$extra.regI == "sED") H <- regH(H, type = 2)   
  

         list(value=res, gradient=G, hessian=H, S.h=S.h, S.h1=S.h1, S.h2=S.h2, 
              l=S.res, l.ln = l.ln, l.par=l.par, ps = ps, 
              eta1=eta1, 
              p1 = p1,
              dl.dbe1          = NULL,       
              dl.dbe2          = NULL,       
              dl.dteta.st      = NULL) 
              
  }
  

  
