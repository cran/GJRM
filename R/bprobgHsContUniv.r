bprobgHsContUniv <- function(params, respvec, VC, ps, AT = FALSE){

epsilon <- 0.0000001 # 0.9999999 0.0001 # sqrt(.Machine$double.eps)

weights <- VC$weights

l.lnun <- bcorR <- NULL

eta2 <- VC$X1%*%params[1:VC$X1.d2] # this is eta1 but changed to eta2 for convenience
eta2 <- eta.tr(eta2, VC$margins[1]) # eta2 is xi for GP and there is no restriction


if( !(VC$margins[1] %in% c(VC$m1d)) ){ 

if(is.null(VC$X2))  sigma2.st <- params[(VC$X1.d2 + 1)] 
if(!is.null(VC$X2)) sigma2.st <- VC$X2%*%params[(VC$X1.d2+1):(VC$X1.d2+VC$X2.d2)]

}

if(VC$margins[1] %in% VC$m1d) sigma2.st <- 0


    sstr1 <- esp.tr(sigma2.st, VC$margins[1])  # is sigma for GP
    sigma2.st <- sstr1$vrb.st 
    sigma2    <- sstr1$vrb 


if(VC$surv == TRUE) naiveind <- FALSE else naiveind <- TRUE  

if(VC$margins[1] %in% VC$m2)            dHs  <-      distrHs(respvec$y1, eta2, sigma2, sigma2.st, nu = 1, nu.st = 1, margin2=VC$margins[1], naive = naiveind)
if(VC$margins[1] %in% c(VC$m1d,VC$m2d)) dHs  <- distrHsDiscr(respvec$y1, eta2, sigma2, sigma2.st, nu = 1, nu.st = 1, margin2=VC$margins[1], naive = TRUE, y2m = VC$y1m)

#########################################################################
#########################################################################


pdf2                         <- dHs$pdf2
derpdf2.dereta2              <- dHs$derpdf2.dereta2 
derpdf2.dersigma2.st         <- dHs$derpdf2.dersigma2.st  
der2pdf2.dereta2             <- dHs$der2pdf2.dereta2
der2pdf2.dersigma2.st2       <- dHs$der2pdf2.dersigma2.st2         
der2pdf2.dereta2dersigma2.st <- dHs$der2pdf2.dereta2dersigma2.st  

if(VC$surv == TRUE){

p2                           <- dHs$p2
derp2.dereta2                <- dHs$derp2.dereta2 
derp2.dersigma2.st           <- dHs$derp2.dersigma.st  
der2p2.dereta2               <- dHs$der2p2.dereta2eta2
der2p2.dersigma2.st2         <- dHs$der2p2.dersigma2.st2         
der2p2.dereta2dersigma2.st   <- dHs$der2p2.dereta2dersigma2.st 

}


########################################################################################################

if(VC$robust == FALSE){

if(VC$surv == FALSE) l.par <- weights*log(pdf2)
if(VC$surv == TRUE)  l.par <- weights*(VC$cens*log(pdf2) + (1-VC$cens)*log(1-p2)  )


res   <- -sum(l.par)
d.psi <- 1

bcorR <- list(b = 0, bp = 0, bs = 0)

}else{


VC$y <- respvec$y1  # do I need this?

if(VC$margins[1] %in% c("NBI", "NBII","NBIa", "NBIIa","PIG","PO","ZTP") ) bcorR <- bcorrecDiscr(VC, params) else{




     #if(is.null(VC$my.env$lB)){ bb <- bounds(VC, params, lo = 1000, tol = VC$tol.rc)
     #                           VC$my.env$lB <- bb$lB
     #                           VC$my.env$uB <- bb$uB
     #                      }
 
     #bcorR <- bcorrec(VC, params, VC$my.env$lB, VC$my.env$uB)
     
     
     #if(VC$r.type == "a") 
     bcorR <- bcorrec(VC, params)
     #if(VC$r.type == "n") bcorR <- bcorrec2(VC, params)
 

}




l.par1    <- log(pdf2)
Robj.lpar <- llpsi(l.par1, VC$rc) # check what happens here with 0 log lik
psi       <- Robj.lpar$psi
d.psi     <- Robj.lpar$d.psi
d2.psi    <- Robj.lpar$d2.psi 
l.par     <- psi  
res       <- -( sum(weights*l.par) - bcorR$b )


}



  
########################################################################################################
 
if(VC$surv == FALSE){ 
 
dl.dbe0       <- derpdf2.dereta2/pdf2
dl.dsigma.st0 <- derpdf2.dersigma2.st/pdf2  

}


if(VC$surv == TRUE){ 
 
dl.dbe0       <- VC$cens*derpdf2.dereta2/pdf2       + (1-VC$cens)*-derp2.dereta2/(1-p2) 
dl.dsigma.st0 <- VC$cens*derpdf2.dersigma2.st/pdf2  + (1-VC$cens)*-derp2.dersigma2.st/(1-p2)  

}


dl.dbe       <- weights*d.psi*dl.dbe0            
dl.dsigma.st <- weights*d.psi*dl.dsigma.st0

                    
########################################################################################################
      
if(VC$robust == TRUE){     
      
d2l.be.be0        <- (der2pdf2.dereta2*pdf2 - (derpdf2.dereta2)^2)/pdf2^2 
d2l.sigma.sigma0  <- (der2pdf2.dersigma2.st2*pdf2-(derpdf2.dersigma2.st)^2)/pdf2^2 
d2l.be.sigma0     <- (der2pdf2.dereta2dersigma2.st*pdf2 - derpdf2.dereta2*derpdf2.dersigma2.st)/pdf2^2 
 
d2l.be.be         <- weights*( d2.psi*dl.dbe0^2             + d.psi*d2l.be.be0 )  
d2l.sigma.sigma   <- weights*( d2.psi*dl.dsigma.st0^2       + d.psi*d2l.sigma.sigma0 )
d2l.be.sigma      <- weights*( d2.psi*dl.dsigma.st0*dl.dbe0 + d.psi*d2l.be.sigma0 )
 
} 

if(VC$robust == FALSE){   



if(VC$surv == FALSE){ 
      
d2l.be.be        <- weights*( (der2pdf2.dereta2*pdf2 - (derpdf2.dereta2)^2)/pdf2^2 )
d2l.sigma.sigma  <- weights*( (der2pdf2.dersigma2.st2*pdf2-(derpdf2.dersigma2.st)^2)/pdf2^2 ) 
d2l.be.sigma     <- weights*( (der2pdf2.dereta2dersigma2.st*pdf2 - derpdf2.dereta2*derpdf2.dersigma2.st)/pdf2^2 )

}


if(VC$surv == TRUE){ 
      
d2l.be.be        <- weights*(     VC$cens*( (der2pdf2.dereta2*pdf2 - (derpdf2.dereta2)^2)/pdf2^2 ) +
                               
                              (1-VC$cens)*( -(der2p2.dereta2*(1-p2) + derp2.dereta2^2)/(1-p2)^2 ) 
                            )

d2l.sigma.sigma  <- weights*(     VC$cens*( (der2pdf2.dersigma2.st2*pdf2-(derpdf2.dersigma2.st)^2)/pdf2^2 ) + 

                              (1-VC$cens)*( -(der2p2.dersigma2.st2*(1-p2) + derp2.dersigma2.st^2)/(1-p2)^2 )
                            )

d2l.be.sigma     <- weights*(    VC$cens*( (der2pdf2.dereta2dersigma2.st*pdf2 - derpdf2.dereta2*derpdf2.dersigma2.st)/pdf2^2 ) +

                             (1-VC$cens)*( -(der2p2.dereta2dersigma2.st*(1-p2) + derp2.dereta2*derp2.dersigma2.st)/(1-p2)^2 )
                            )

}


 
} 
 
########################################################################################################
     





if( !(VC$margins[1] %in% c(VC$m1d)) ){ ##

if( is.null(VC$X2) ){

  G   <- -c( colSums( c(dl.dbe)*VC$X1 ) ,
            sum( dl.dsigma.st )
          ) + bcorR$bp 
          
                
  be.be    <- crossprod(VC$X1*c(d2l.be.be),VC$X1)
  be.sigma <- t(t(rowSums(t(VC$X1*c(d2l.be.sigma))))) 

  H <- -rbind( cbind( be.be      , be.sigma             ), 
              cbind( t(be.sigma), sum(d2l.sigma.sigma) )  
            )    +  bcorR$bs  
            
                    }


if( !is.null(VC$X2) ){

  G   <- -c( colSums( c(dl.dbe)*VC$X1        ) ,
            colSums( c(dl.dsigma.st)*VC$X2  )
          ) + bcorR$bp 
                
  be.be    <- crossprod(VC$X1*c(d2l.be.be),VC$X1)
  be.sigma <- crossprod(VC$X1*c(d2l.be.sigma),VC$X2)
  si.sigma <- crossprod(VC$X2*c(d2l.sigma.sigma),VC$X2)  

  H <- -rbind( cbind( be.be      , be.sigma  ), 
              cbind( t(be.sigma), si.sigma  )  
            ) +  bcorR$bs
                     }
                                         
                     

} ##


if( VC$margins[1] %in% c(VC$m1d) ){

  G <- -c( colSums( c(dl.dbe)*VC$X1 ) )     +  bcorR$bp
  H <- -crossprod(VC$X1*c(d2l.be.be),VC$X1) +  bcorR$bs

}




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
  
  
  
if( VC$margins[1] == "LN"){
  
  dHs1 <- distrHsAT(exp(respvec$y1), eta2, sigma2, 1, margin2=VC$margins[1])
  if(VC$surv == FALSE) l.lnun <- -( sum(weights*log(dHs1$pdf2)) - bcorR$b )
  if(VC$surv == TRUE)  l.lnun <- -sum(weights*(VC$cens*log(dHs1$pdf2)+(1-VC$cens)*log(1-dHs1$p2) )      )
    
  }  
  



if(VC$margins[1] %in% c("GP","DGP")) indx <- dHs$indx else indx <- 1 



list(value=res, gradient=G, hessian=H, S.h=S.h, S.h1=S.h1, S.h2=S.h2, l=S.res, l.lnun = l.lnun, 
     l.par=l.par, ps = ps, sigma2.st = sigma2.st,
     etas1 = sigma2.st, eta1 = eta2, bcorR = bcorR,  
     BivD=VC$BivD, eta2 = eta2, sigma2 = sigma2, nu = NULL, d.psi = d.psi,
     dl.dbe = dl.dbe, dl.dsigma.st = dl.dsigma.st, indx = indx)      


}


