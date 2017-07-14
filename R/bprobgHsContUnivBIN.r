bprobgHsContUnivBIN <- function(params, respvec, VC, ps, AT = FALSE){

weights <- VC$weights

l.lnun <- NULL

eta2 <- VC$X1%*%params
pd   <- probm(eta2, VC$margins[1], only.pr = FALSE, tau = VC$gev.par)
y    <- respvec$y1

#########################################################################
tauetaIND    <- pd$tauetaIND==FALSE

pr           <- pd$pr[tauetaIND]
d.n          <- pd$d.n[tauetaIND]
der2p.dereta <- pd$der2p.dereta[tauetaIND]

X1 <- VC$X1[tauetaIND,]
y  <- y[tauetaIND]
weights <- weights[tauetaIND]
########################################################################################################

l.par <- weights*( y*log(pr) + (1-y)*log(1-pr) )
res   <- -sum(l.par)

########################################################################################################
 
dl.dbe    <- weights*(  ( y/pr - (1-y)/(1-pr) )*d.n )                 
d2l.be.be <- weights*(  ( y/pr - (1-y)/(1-pr) )*der2p.dereta + d.n^2*( -y/pr^2 - (1-y)/(1-pr)^2  )  )

########################################################################################################
     
  G <- -c( colSums( c(dl.dbe)*X1 ) )     
  H <- -crossprod(X1*c(d2l.be.be),X1) 

if(VC$extra.regI == "pC") H <- regH(H, type = 1)
  
  S.h <- ps$S.h  

  if( length(S.h) != 1 ){
  
  S.h1 <- 0.5*crossprod(params,S.h)%*%params
  S.h2 <- S.h%*%params
  
  } else S.h <- S.h1 <- S.h2 <- 0   

  
  S.res <- res
  res   <- S.res + S.h1
  G     <- G + S.h2
  H     <- H + S.h 
        
if(VC$extra.regI == "sED") H <- regH(H, type = 2) 
  
list(value=res, gradient=G, hessian=H, S.h=S.h, S.h1=S.h1, S.h2=S.h2, l=S.res, l.lnun = l.lnun, 
     l.par=l.par, ps = ps, sigma2.st = NULL,
     etas1 = NULL, eta1 = eta2, 
     BivD=VC$BivD, eta2 = eta2, sigma2 = NULL, nu = NULL, tauetaIND = tauetaIND)      


}


