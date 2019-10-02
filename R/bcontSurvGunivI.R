bcontSurvGunivI <- function(params, respvec, VC, ps, AT = FALSE){

monP <- monP1 <- monP2 <- k <- 0; V <- list()

etad <- etas1 <- etas2 <- l.ln <- NULL 

    params1 <- params[1:VC$X1.d2]
    params1[VC$mono.sm.pos] <- exp( params1[VC$mono.sm.pos] )

    eta1 <- VC$X1%*%params1
    eta2 <- VC$X2%*%params1
     
    Xd1P <- VC$Xd1%*%params1  
   
    indN <- as.numeric(Xd1P < 0) 

    Xd1P <- ifelse(Xd1P < sqrt(.Machine$double.eps), sqrt(.Machine$double.eps), Xd1P ) # safety check 

##################
    
pd1  <- probmS(eta1, VC$margins[1]) 
pd2  <- probmS(eta2, VC$margins[1]) 
  
p1       <- pd1$pr
dS1eta1  <- pd1$dS
d2S1eta1 <- pd1$d2S
d3S1eta1 <- pd1$d3S

p2       <- pd2$pr
dS2eta2  <- pd2$dS
d2S2eta2 <- pd2$d2S
d3S2eta2 <- pd2$d3S


##################


l.par <- VC$weights*( VC$cens*( log(-dS1eta1) + log(Xd1P) ) + (1 - VC$cens)*log(mm(p1-p2)) )   
res   <- -sum(l.par)

##################

der.par1 <- der2.par1 <- params1

 der.par1[-c( VC$mono.sm.pos )] <- 1
der2.par1[-c( VC$mono.sm.pos )] <- 0

der2eta1dery1b1 <- t(t(VC$Xd1)*der.par1)

dereta1derb1    <- t(t(VC$X1)*der.par1)
dereta2derb1    <- t(t(VC$X2)*der.par1)


##################

dl.dbe1 <- -VC$weights*(   

   VC$cens*( c((dS1eta1*Xd1P)^-1)*(c(d2S1eta1*Xd1P)*dereta1derb1 + c(dS1eta1)*der2eta1dery1b1) ) + 
   
   (1 - VC$cens)*c(mm(p1-p2))^-1*(c(dS1eta1)*dereta1derb1-c(dS2eta2)*dereta2derb1)    
     
                       )
   
G <- colSums(dl.dbe1)
   
########################
########################                                                                                                 
   
   
 H <-  -(   

   crossprod(c(VC$weights*VC$cens*(-dS1eta1^-2*d2S1eta1^2 + dS1eta1^-1*d3S1eta1))*dereta1derb1, dereta1derb1) +
    
   diag( colSums( t( t(  c(VC$weights*VC$cens*dS1eta1^-1*d2S1eta1)*VC$X1)*der2.par1 ) ) ) + 
   
   diag( colSums( t( t(VC$weights*VC$cens*c(Xd1P^-1)*VC$Xd1)*der2.par1 ) ) )  +
    
   crossprod(VC$weights*VC$cens*c(-Xd1P^-2)*der2eta1dery1b1, der2eta1dery1b1) +  


  crossprod(c(VC$weights*(1 - VC$cens)*(mm(p1-p2)^-1*d2S1eta1))*dereta1derb1, dereta1derb1) +      ### 

  crossprod(c(VC$weights*(1 - VC$cens)*(-mm(p1-p2)^-1*d2S2eta2))*dereta2derb1, dereta2derb1) +  
 
   
     diag( colSums( t( t(c(VC$weights*(1 - VC$cens)*mm(p1-p2)^-1*dS1eta1)*VC$X1)*der2.par1 ) ) ) + 
  
     diag( colSums( t( t(c(VC$weights*(1 - VC$cens)*-mm(p1-p2)^-1*dS2eta2)*VC$X2)*der2.par1 ) ) ) +
     

   crossprod(c(VC$weights*(1 - VC$cens)*(-mm(p1-p2)^-2*dS1eta1^2))*dereta1derb1, dereta1derb1) +
     
   crossprod(c(VC$weights*(1 - VC$cens)*(-mm(p1-p2)^-2*dS2eta2^2))*dereta2derb1, dereta2derb1) +
     
     
   crossprod(c(VC$weights*(1 - VC$cens)*(mm(p1-p2)^-2*dS1eta1*dS2eta2))*dereta1derb1, dereta2derb1) +
   
   crossprod(c(VC$weights*(1 - VC$cens)*(mm(p1-p2)^-2*dS1eta1*dS2eta2))*dereta2derb1, dereta1derb1) 


  )
  
                
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
  

         list(value=res, gradient=G, hessian=H, S.h=S.h, S.h1=S.h1, S.h2=S.h2, indN = indN, V = V, 
              l=S.res, l.ln = l.ln, l.par=l.par, ps = ps, k = VC$my.env$k, monP2 = monP2, params1 = params1,
              eta1=eta1, 
              p1 = p1,
              dl.dbe1          = NULL,       
              dl.dbe2          = NULL,       
              dl.dteta.st      = NULL) 
              
  }
  