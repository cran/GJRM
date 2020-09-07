bcontSurvGunivInform <- function(params, respvec, VC, ps, AT = FALSE){
p1 <- p2 <- pdf1 <- pdf2 <- c.copula.be2 <- c.copula.be1 <- c.copula2.be1be2 <- NA

etad <- etas1 <- etas2 <- l.ln <- NULL

Xi <- Xi <- VC$Xi  
X1ni <- VC$X1ni
X2ni <- VC$X2ni
Gmat12 <- VC$Gmat12
H      <- VC$Hmat12

inde.inf1 <- VC$inde.inf1
inde.inf2 <- VC$inde.inf2

    params1 <- params[1:VC$X1.d2]
    params1[VC$mono.sm.pos] <- exp( params1[VC$mono.sm.pos] )

    eta1 <- VC$X1%*%params1
    Xd1P <- VC$Xd1%*%params1  
    Xd1P <- ifelse(Xd1P < VC$min.dn, VC$min.dn, Xd1P )
        
    params2 <- params[-c(1:VC$X1.d2)]
    
    
    inde.b2 <- VC$X1.d2 + (1:length(params2))
    inde.b1 <- 1:VC$X1.d2
    inde.b1 <- inde.b1[-c(inde.inf1)]
    inde.b0 <- inde.inf1
    
    if(!is.null(VC$infsetupR$pcv)){
    
          for(i in 1:VC$infsetupR$lp1) params2 <- append(params2, params1[c(VC$infsetupR$Lposs2PAR[[i]])], (sort(VC$infsetupR$poss2PAR)[i]-1))   
    
                                  }
    
    
    
    if(!is.null(VC$infsetupR$smo.pos1)){
    
          for(i in 1:length(VC$infsetupR$poss2)) params2 <- append(params2, params1[c(VC$infsetupR$Lposs2[[i]])], (sort(VC$infsetupR$poss2)[i]-1)) 
    
                                       }
    
    
    params2[VC$mono.sm.pos2] <- exp( params2[VC$mono.sm.pos2]) # mono is same as mono1 anyway

    eta2 <- VC$X2%*%params2
    Xd2P <- VC$Xd2%*%params2

    Xd2P <- ifelse(Xd2P < VC$min.dn, VC$min.dn, Xd2P )

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

l.par <- VC$weights*(

              log(p1) + log(p2) + VC$cens*log( Xd1P*-dS1eta1/p1) + (1-VC$cens)*log( Xd2P*-dS2eta2/p2)  

                    )  

res   <- -sum(l.par)

##################

######
######

der.par1 <- der2.par1 <- params1#[-c(inde.inf1)]
der.par2 <- der2.par2 <- params2#[-c(inde.inf2)]

der.par1[-c( VC$mono.sm.pos  )] <- 1
der.par2[-c( VC$mono.sm.pos2 )] <- 1


der2.par1[-c( VC$mono.sm.pos )]  <- 0
der2.par2[-c( VC$mono.sm.pos2 )] <- 0


der.par1  <- der.par1[-c(inde.inf1)] 
der.par2  <- der.par2[-c(inde.inf2)] 
der2.par1 <- der2.par1[-c(inde.inf1)] 
der2.par2 <- der2.par2[-c(inde.inf2)] # new part that was creating a mismatch! fixed now


der2eta1dery1b1  <- t(t(VC$Xd1[,-c(inde.inf1)])*der.par1)
der2eta2dery2b2  <- t(t(VC$Xd2[,-c(inde.inf2)])*der.par2)  


der3eta1dery1b12 <- t(t(VC$Xd1[,-c(inde.inf1)])*der2.par1)  # NEW, should be like above #
der3eta2dery2b22 <- t(t(VC$Xd2[,-c(inde.inf2)])*der2.par2)  # NEW, should be like above #


der2eta1dery1b0  <- der2eta2dery2b0 <- 0


der3eta1dery1b02 <- der3eta2dery2b02  <- 0 # NEW # ok fine
der3eta1dery1b0b1<- der3eta2dery2b0b2 <- 0 # NEW # ok fine


dereta1derb1 <- t(t(X1ni)*der.par1)
dereta2derb2 <- t(t(X2ni)*der.par2)  

der2eta1derb12 <- t(t(X1ni)*der2.par1) # NEW #
der2eta2derb22 <- t(t(X2ni)*der2.par2) # NEW #


dereta1derb0 <- dereta2derb0 <- Xi


 der2eta1derb02 <-  der2eta2derb02 <- 0 # NEW # ok fine
der2eta1derb0b1 <- der2eta2derb0b2 <- 0 # NEW # ok fine


##################


dl.dbe0 <- as.matrix( -VC$weights*(  c(dS1eta1/p1)*dereta1derb0 + c(dS2eta2/p2)*dereta2derb0  +    
                       VC$cens*( dereta1derb0*c(d2S1eta1/dS1eta1 -  dS1eta1/p1)) +
                       (1-VC$cens)*( dereta2derb0*c(d2S2eta2/dS2eta2 -dS2eta2/p2) ) ) )

dl.dbe1 <- as.matrix( -VC$weights*(  c(dS1eta1/p1)*dereta1derb1 + VC$cens*( dereta1derb1*c(d2S1eta1/dS1eta1 -  dS1eta1/p1) + der2eta1dery1b1*c(Xd1P^-1)) ) )

dl.dbe2 <- as.matrix( -VC$weights*(  c(dS2eta2/p2)*dereta2derb2 + (1-VC$cens)*( dereta2derb2*c(d2S2eta2/dS2eta2 -  dS2eta2/p2) + der2eta2dery2b2*c(Xd2P^-1)) ) )


Gmat12[, inde.inf1]     <- dl.dbe0
Gmat12[, -c(inde.inf1)] <- dl.dbe1

G <- c( colSums(Gmat12), colSums(dl.dbe2) )

########################
########################    


b0b0 <- -(

  crossprod(c(VC$weights* ( (d2S1eta1/p1 - dS1eta1^2/p1^2 ) + VC$cens*(d3S1eta1/dS1eta1 - d2S1eta1^2/dS1eta1^2 - d2S1eta1/p1 + dS1eta1^2/p1^2  ) ) ) *dereta1derb0, dereta1derb0)

+ crossprod(c(VC$weights*((d2S2eta2/p2 - dS2eta2^2/p2^2 )+ (1-VC$cens)*(d3S2eta2/dS2eta2 - d2S2eta2^2/dS2eta2^2 - d2S2eta2/p2 + dS2eta2^2/p2^2  ) ) )*dereta2derb0, dereta2derb0)  
 
)


b0b1 <- - crossprod(c(VC$weights*((d2S1eta1/p1 - dS1eta1^2/p1^2 )+ VC$cens*(d3S1eta1/dS1eta1 - d2S1eta1^2/dS1eta1^2 - d2S1eta1/p1 + dS1eta1^2/p1^2  ) ) )*dereta1derb0, dereta1derb1) 

b0b2 <- - crossprod(c(VC$weights*( (d2S2eta2/p2 - dS2eta2^2/p2^2 )+ (1-VC$cens)*(d3S2eta2/dS2eta2 - d2S2eta2^2/dS2eta2^2 - d2S2eta2/p2 + dS2eta2^2/p2^2  ) ) )*dereta2derb0, dereta2derb2)  


b1b1 <-  -(   
   
   crossprod(c(VC$weights*( VC$cens*(-dS1eta1^-2*d2S1eta1^2 + dS1eta1^-1*d3S1eta1) + (1 - VC$cens)*(-p1^-2*dS1eta1^2+p1^-1*d2S1eta1)))*dereta1derb1, dereta1derb1) +
    
   diag( colSums( t( t(  c(VC$weights*VC$cens*dS1eta1^-1*d2S1eta1)*X1ni)*der2.par1 ) ) ) +
   
   diag( colSums( t( t(VC$weights*VC$cens*c(Xd1P^-1)*VC$Xd1[,-c(inde.inf1)])*der2.par1 ) ) )  +
    
   crossprod(VC$weights*VC$cens*c(-Xd1P^-2)*der2eta1dery1b1, der2eta1dery1b1) +  
    
   diag( colSums( t( t(c(VC$weights*(1 - VC$cens)*p1^-1*dS1eta1)*X1ni)*der2.par1 ) ) )

  )


b2b2 <-  -(   
   
   crossprod(c(VC$weights*( (1-VC$cens)*(-dS2eta2^-2*d2S2eta2^2 + dS2eta2^-1*d3S2eta2) + VC$cens*(-p2^-2*dS2eta2^2+p2^-1*d2S2eta2 ) ) )*dereta2derb2, dereta2derb2) +
    
   diag( colSums( t( t(  c(VC$weights*(1-VC$cens)*dS2eta2^-1*d2S2eta2)*X2ni)*der2.par2 ) ) ) +
   
   diag( colSums( t( t(VC$weights*(1-VC$cens)*c(Xd2P^-1)*VC$Xd2[,-c(inde.inf2)])*der2.par2 ) ) )  +
    
   crossprod(VC$weights*(1-VC$cens)*c(-Xd2P^-2)*der2eta2dery2b2, der2eta2dery2b2) +
   
   diag( colSums( t( t(c(VC$weights*(VC$cens)*p2^-1*dS2eta2)*X2ni)*der2.par2 ) ) )

)


b1b2 <- matrix(0, length(params1[-c(inde.inf1)]), length(params[-c(1:VC$X1.d2)]))  
 
 
H[inde.b0, inde.b0] <- b0b0  
H[inde.b0, inde.b1] <- b0b1
H[inde.b0, inde.b2] <- b0b2
H[inde.b1, inde.b0] <- t(b0b1)
H[inde.b1, inde.b1] <- b1b1
H[inde.b1, inde.b2] <- b1b2
H[inde.b2, inde.b0] <- t(b0b2)
H[inde.b2, inde.b1] <- t(b1b2)
H[inde.b2, inde.b2] <- b2b2

    
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
              l=S.res, l.ln = l.ln, l.par=l.par, ps = ps, params1 = params1, params2 = params2,
              eta1=eta1, eta2 = eta2, 
                                   p1 = p1, p2 = p2, pdf1 = -dS1eta1, pdf2 = -dS2eta2,          
	      	                           c.copula.be2 = c.copula.be2,
	      	                           c.copula.be1 = c.copula.be1,
              c.copula2.be1be2 = c.copula2.be1be2,
              dl.dbe1          = NULL,       
              dl.dbe2          = NULL,       
              dl.dteta.st      = NULL)
              
  }
  
