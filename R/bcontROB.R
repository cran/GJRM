bcontROB <- function(params, respvec, VC, ps, AT = FALSE){
p1 <- p2 <- pdf1 <- pdf2 <- c.copula.be2 <- c.copula.be1 <- c.copula2.be1be2 <- NA

    eta1 <- VC$X1%*%params[1:VC$X1.d2]
    eta2 <- VC$X2%*%params[(VC$X1.d2 + 1):(VC$X1.d2 + VC$X2.d2)]
    etad <- etas1 <- etas2 <- l.ln <- NULL 
  
      
  if(is.null(VC$X3)){  
    sigma21.st <- etas1 <- params[(VC$X1.d2 + VC$X2.d2 + 1)]
    sigma22.st <- etas2 <- params[(VC$X1.d2 + VC$X2.d2 + 2)]
    teta.st    <- etad  <- params[(VC$X1.d2 + VC$X2.d2 + 3)]
  } 
  
  
  if(!is.null(VC$X3)){  
    sigma21.st <- etas1 <- VC$X3%*%params[(VC$X1.d2 + VC$X2.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2)]
    sigma22.st <- etas2 <- VC$X4%*%params[(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2)]
    teta.st    <- etad  <- VC$X5%*%params[(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2 + VC$X5.d2)]
  }  
  
  
##################
## Transformations
##################
  
    sstr1 <- esp.tr(sigma21.st, VC$margins[1])  
    sstr2 <- esp.tr(sigma22.st, VC$margins[2])  
  
    sigma21.st <- sstr1$vrb.st 
    sigma22.st <- sstr2$vrb.st 
    
    sigma21    <- sstr1$vrb 
    sigma22    <- sstr2$vrb 
    
    eta1 <- eta.tr(eta1, VC$margins[1])
    eta2 <- eta.tr(eta2, VC$margins[2])
    
resT    <- teta.tr(VC, teta.st)
teta.st <- resT$teta.st
teta    <- resT$teta
    
##################
##################

  dHs1 <- distrHs(respvec$y1, eta1, sigma21, sigma21.st, nu = 1, nu.st = 1, margin2=VC$margins[1], naive = FALSE, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
  dHs2 <- distrHs(respvec$y2, eta2, sigma22, sigma22.st, nu = 1, nu.st = 1, margin2=VC$margins[2], naive = FALSE, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)

  pdf1 <- dHs1$pdf2
  pdf2 <- dHs2$pdf2

  p1 <- dHs1$p2
  p2 <- dHs2$p2
  
  dH <- copgHsAT(p1, p2, teta, VC$BivD, Ln = TRUE, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr, par2 = VC$dof)

  c.copula2.be1be2 <- dH$c.copula2.be1be2
  
  l.par1 <- log(pdf1) + log(pdf2) + log(c.copula2.be1be2) 
  
  Robj.lpar <- llpsi(l.par1, VC$tc)

  psi    <- Robj.lpar$psi
  d.psi  <- Robj.lpar$d.psi
  d2.psi <- Robj.lpar$d2.psi 
 
  l.par <- VC$weights*( psi )
 
##################
##################

 derpdf1.dereta1              <- dHs1$derpdf2.dereta2 
 derpdf1.dersigma21.st        <- dHs1$derpdf2.dersigma2.st 
 
 derpdf2.dereta2              <- dHs2$derpdf2.dereta2 
 derpdf2.dersigma22.st        <- dHs2$derpdf2.dersigma2.st 
 
 derp1.dereta1                <- dHs1$derp2.dereta2
 derp1.dersigma21.st          <- dHs1$derp2.dersigma.st
 
 derp2.dereta2                <- dHs2$derp2.dereta2
 derp2.dersigma22.st          <- dHs2$derp2.dersigma.st 
 
 
 BITS <- copgHsCont(p1, p2, teta, teta.st, VC$BivD, Cont = TRUE, par2 = VC$dof, nu.st = log(VC$dof - 2))
 
 
   der2h.derp1p1              <- BITS$der2h.derp1p1
   derc.dereta1               <- der2h.derp1p1 * derp1.dereta1 
   derc.dersigma21.st         <- der2h.derp1p1 * derp1.dersigma21.st
   der2h.derp1p2              <- BITS$der2h.derp1p2  
   derc.dereta2               <- der2h.derp1p2 * derp2.dereta2    
   derc.dersigma22.st         <- der2h.derp1p2 * derp2.dersigma22.st
   der2h.derp1teta            <- BITS$der2h.derp1teta                                     
   derteta.derteta.st         <- BITS$derteta.derteta.st 
   der2h.derp1teta.st         <- der2h.derp1teta * derteta.derteta.st # new bit
   
   
 
   dl.dbe1.0        <- derpdf1.dereta1/pdf1 + derc.dereta1/c.copula2.be1be2                 
   dl.dbe2.0        <- derpdf2.dereta2/pdf2 + derc.dereta2/c.copula2.be1be2              
   dl.dsigma21.st.0 <- derpdf1.dersigma21.st/pdf1 + derc.dersigma21.st/c.copula2.be1be2  
   dl.dsigma22.st.0 <- derpdf2.dersigma22.st/pdf2 + derc.dersigma22.st/c.copula2.be1be2                
   dl.dteta.st.0    <- der2h.derp1teta.st/c.copula2.be1be2                               
   

   dl.dbe1        <- VC$weights*( d.psi*dl.dbe1.0          )   
   dl.dbe2        <- VC$weights*( d.psi*dl.dbe2.0          )
   dl.dsigma21.st <- VC$weights*( d.psi*dl.dsigma21.st.0   )
   dl.dsigma22.st <- VC$weights*( d.psi*dl.dsigma22.st.0   )              
   dl.dteta.st    <- VC$weights*( d.psi*dl.dteta.st.0      )   
  

               
#################################################################################################

der2c.derrho.derrho    <- BITS$der2c.derrho.derrho
der2c.derp1.derp1      <- BITS$der2c.derp1.derp1  
der2c.derp2.derp2      <- BITS$der2c.derp2.derp2  
der2c.derp1.derp2      <- BITS$der2c.derp1.derp2  
der2c.derp1.derrho     <- BITS$der2c.derp1.derrho 
der2c.derp2.derrho     <- BITS$der2c.derp2.derrho 

der2pdf1.dereta1 <- dHs1$der2pdf2.dereta2
der2pdf2.dereta2 <- dHs2$der2pdf2.dereta2

der2pdf1.dersigma21.st2  <- dHs1$der2pdf2.dersigma2.st2
der2pdf2.dersigma22.st2  <- dHs2$der2pdf2.dersigma2.st2

der2p1.dereta1eta1 <- dHs1$der2p2.dereta2eta2
der2p2.dereta2eta2 <- dHs2$der2p2.dereta2eta2
 
der2p1.dersigma21.st2 <-  dHs1$der2p2.dersigma2.st2
der2p2.dersigma22.st2 <-  dHs2$der2p2.dersigma2.st2
 
der2pdf1.dereta1dersigma21.st <-  dHs1$der2pdf2.dereta2dersigma2.st
der2pdf2.dereta2dersigma22.st <-  dHs2$der2pdf2.dereta2dersigma2.st
 
der2p1.dereta1dersigma21.st <-  dHs1$der2p2.dereta2dersigma2.st
der2p2.dereta2dersigma22.st <-  dHs2$der2p2.dereta2dersigma2.st


der2teta.derteta.stteta.st <- BITS$der2teta.derteta.stteta.st 

                  
  d2l.be1.be1.0         <- (der2pdf1.dereta1 * pdf1 - derpdf1.dereta1^2) / pdf1^2 +((der2c.derp1.derp1 * derp1.dereta1^2 + der2h.derp1p1 * der2p1.dereta1eta1) * c.copula2.be1be2 - derc.dereta1^2) /c.copula2.be1be2^2                    
  d2l.be2.be2.0         <- (der2pdf2.dereta2 * pdf2 - derpdf2.dereta2^2) / pdf2^2   + ((der2c.derp2.derp2 * derp2.dereta2^2 + der2h.derp1p2 * der2p2.dereta2eta2) * c.copula2.be1be2 - derc.dereta2^2) /c.copula2.be1be2^2                 
  d2l.rho.rho.0         <- ((der2c.derrho.derrho*derteta.derteta.st^2 + der2h.derp1teta *der2teta.derteta.stteta.st)*c.copula2.be1be2 -der2h.derp1teta.st^2) /c.copula2.be1be2^2                                         
  d2l.sigma21.sigma21.0 <- (der2pdf1.dersigma21.st2 * pdf1 - derpdf1.dersigma21.st^2) / pdf1^2   + ((der2c.derp1.derp1*derp1.dersigma21.st^2+der2h.derp1p1*der2p1.dersigma21.st2)*c.copula2.be1be2 - derc.dersigma21.st^2) /c.copula2.be1be2^2                  
  d2l.sigma22.sigma22.0 <- (der2pdf2.dersigma22.st2 * pdf2 - derpdf2.dersigma22.st^2) / pdf2^2   + ((der2c.derp2.derp2*derp2.dersigma22.st^2+der2h.derp1p2*der2p2.dersigma22.st2)*c.copula2.be1be2 - derc.dersigma22.st^2) /c.copula2.be1be2^2 
  d2l.be1.be2.0         <- (der2c.derp1.derp2*derp1.dereta1*derp2.dereta2* c.copula2.be1be2 - derc.dereta1*derc.dereta2) /c.copula2.be1be2^2 
  d2l.be1.rho.0         <- (der2c.derp1.derrho *derp1.dereta1*derteta.derteta.st* c.copula2.be1be2 - derc.dereta1*der2h.derp1teta*derteta.derteta.st)  /c.copula2.be1be2^2 
  d2l.be2.rho.0         <- (der2c.derp2.derrho *derp2.dereta2*derteta.derteta.st* c.copula2.be1be2 - derc.dereta2*der2h.derp1teta*derteta.derteta.st) /c.copula2.be1be2^2 
  d2l.be1.sigma21.0     <- (der2pdf1.dereta1dersigma21.st * pdf1 - derpdf1.dereta1*derpdf1.dersigma21.st) / pdf1^2   +  ((der2c.derp1.derp1 * derp1.dereta1* derp1.dersigma21.st+der2h.derp1p1*der2p1.dereta1dersigma21.st) * c.copula2.be1be2 - derc.dereta1*derc.dersigma21.st) /c.copula2.be1be2^2                   
  d2l.be2.sigma22.0     <- (der2pdf2.dereta2dersigma22.st * pdf2 - derpdf2.dereta2*derpdf2.dersigma22.st) / pdf2^2   + ((der2c.derp2.derp2 * derp2.dereta2* derp2.dersigma22.st+der2h.derp1p2*der2p2.dereta2dersigma22.st) * c.copula2.be1be2 - derc.dereta2*derc.dersigma22.st) /c.copula2.be1be2^2                                                   
  d2l.be2.sigma21.0     <- (der2c.derp1.derp2*derp2.dereta2*derp1.dersigma21.st* c.copula2.be1be2 - derc.dereta2*derc.dersigma21.st)/c.copula2.be1be2^2
  d2l.be1.sigma22.0     <- (der2c.derp1.derp2*derp1.dereta1*derp2.dersigma22.st* c.copula2.be1be2 - derc.dereta1*derc.dersigma22.st)/c.copula2.be1be2^2
  d2l.rho.sigma21.0     <- (der2c.derp1.derrho*derp1.dersigma21.st*derteta.derteta.st*c.copula2.be1be2 - derc.dersigma21.st*der2h.derp1teta *derteta.derteta.st)/c.copula2.be1be2^2 
  d2l.rho.sigma22.0     <- (der2c.derp2.derrho*derp2.dersigma22.st*derteta.derteta.st*c.copula2.be1be2 - derc.dersigma22.st*der2h.derp1teta *derteta.derteta.st)/c.copula2.be1be2^2 
  d2l.sigma21.sigma22.0 <- (der2c.derp1.derp2*derp1.dersigma21.st*derp2.dersigma22.st* c.copula2.be1be2 - derc.dersigma21.st*derc.dersigma22.st ) /c.copula2.be1be2^2 
  
  
  d2l.be1.be1         <- -VC$weights*( d2.psi*dl.dbe1.0^2		        +d.psi*d2l.be1.be1.0        )
  d2l.be2.be2         <- -VC$weights*( d2.psi*dl.dbe2.0^2		        +d.psi*d2l.be2.be2.0        )
  d2l.rho.rho         <- -VC$weights*( d2.psi*dl.dteta.st.0^2		        +d.psi*d2l.rho.rho.0        )
  d2l.sigma21.sigma21 <- -VC$weights*( d2.psi*dl.dsigma21.st.0^2		+d.psi*d2l.sigma21.sigma21.0)
  d2l.sigma22.sigma22 <- -VC$weights*( d2.psi*dl.dsigma22.st.0^2		+d.psi*d2l.sigma22.sigma22.0)
  d2l.be1.be2         <- -VC$weights*( d2.psi*dl.dbe1.0*dl.dbe2.0		+d.psi*d2l.be1.be2.0        )
  d2l.be1.rho         <- -VC$weights*( d2.psi*dl.dbe1.0*dl.dteta.st.0		+d.psi*d2l.be1.rho.0        )
  d2l.be2.rho         <- -VC$weights*( d2.psi*dl.dbe2.0*dl.dteta.st.0		+d.psi*d2l.be2.rho.0        )
  d2l.be1.sigma21     <- -VC$weights*( d2.psi*dl.dbe1.0*dl.dsigma21.st.0	+d.psi*d2l.be1.sigma21.0    )
  d2l.be2.sigma22     <- -VC$weights*( d2.psi*dl.dbe2.0*dl.dsigma22.st.0	+d.psi*d2l.be2.sigma22.0    )
  d2l.be2.sigma21     <- -VC$weights*( d2.psi*dl.dbe2.0*dl.dsigma21.st.0	+d.psi*d2l.be2.sigma21.0    )
  d2l.be1.sigma22     <- -VC$weights*( d2.psi*dl.dbe1.0*dl.dsigma22.st.0	+d.psi*d2l.be1.sigma22.0    )
  d2l.rho.sigma21     <- -VC$weights*( d2.psi*dl.dteta.st.0*dl.dsigma21.st.0	+d.psi*d2l.rho.sigma21.0    )
  d2l.rho.sigma22     <- -VC$weights*( d2.psi*dl.dteta.st.0*dl.dsigma22.st.0	+d.psi*d2l.rho.sigma22.0    )
  d2l.sigma21.sigma22 <- -VC$weights*( d2.psi*dl.dsigma21.st.0*dl.dsigma22.st.0	+d.psi*d2l.sigma21.sigma22.0)
  

                 
 

if( is.null(VC$X3) ){



  G   <- -c(colSums( c(dl.dbe1)*VC$X1 ),
            colSums( c(dl.dbe2)*VC$X2 ),
            sum( dl.dsigma21.st ),
            sum( dl.dsigma22.st ),
            sum( dl.dteta.st )            )



  be1.be1 <- crossprod(VC$X1*c(d2l.be1.be1),VC$X1)
  be2.be2 <- crossprod(VC$X2*c(d2l.be2.be2),VC$X2)
  be1.be2 <- crossprod(VC$X1*c(d2l.be1.be2),VC$X2)
  be1.rho <-   t(t(rowSums(t(VC$X1*c(d2l.be1.rho)))))
  be1.sigma21 <- t(t(rowSums(t(VC$X1*c(d2l.be1.sigma21))))) 
  be1.sigma22 <- t(t(rowSums(t(VC$X1*c(d2l.be1.sigma22))))) 
  be2.rho <-   t(t(rowSums(t(VC$X2*c(d2l.be2.rho)))))
  be2.sigma21 <- t(t(rowSums(t(VC$X2*c(d2l.be2.sigma21))))) 
  be2.sigma22 <- t(t(rowSums(t(VC$X2*c(d2l.be2.sigma22)))))

  H <- rbind( cbind( be1.be1    ,   be1.be2    ,   be1.sigma21,    be1.sigma22,       be1.rho  ), 
              cbind( t(be1.be2) ,   be2.be2    ,   be2.sigma21,    be2.sigma22,        be2.rho  ), 
              cbind( t(be1.sigma21) , t(be2.sigma21) , sum(d2l.sigma21.sigma21), sum(d2l.sigma21.sigma22), sum(d2l.rho.sigma21)),
              cbind( t(be1.sigma22) , t(be2.sigma22) , sum(d2l.sigma21.sigma22), sum(d2l.sigma22.sigma22), sum(d2l.rho.sigma22)),
              cbind( t(be1.rho) ,   t(be2.rho) ,   sum(d2l.rho.sigma21), sum(d2l.rho.sigma22), sum(d2l.rho.rho) ) 
              
              ) 

}



if( !is.null(VC$X3) ){



G   <- -c( colSums(       c(dl.dbe1)*VC$X1 ) ,
           colSums(       c(dl.dbe2)*VC$X2 ) ,
           colSums(c(dl.dsigma21.st)*VC$X3 ) ,
           colSums(c(dl.dsigma22.st)*VC$X4 ) ,
           colSums(   c(dl.dteta.st)*VC$X5 )  )  

                 
    be1.be1         <- crossprod(VC$X1*c(d2l.be1.be1),VC$X1)
    be2.be2         <- crossprod(VC$X2*c(d2l.be2.be2),VC$X2)
    be1.be2         <- crossprod(VC$X1*c(d2l.be1.be2),VC$X2)
    be1.rho         <- crossprod(VC$X1*c(d2l.be1.rho),VC$X5)    
    be2.rho         <- crossprod(VC$X2*c(d2l.be2.rho),VC$X5)   
    be1.sigma21     <- crossprod(VC$X1*c(d2l.be1.sigma21),VC$X3)      
    be1.sigma22     <- crossprod(VC$X1*c(d2l.be1.sigma22),VC$X4)      
    be2.sigma21     <- crossprod(VC$X2*c(d2l.be2.sigma21),VC$X3)   
    be2.sigma22     <- crossprod(VC$X2*c(d2l.be2.sigma22),VC$X4)   
    sigma21.sigma21 <- crossprod(VC$X3*c(d2l.sigma21.sigma21),VC$X3)     
    sigma21.sigma22 <- crossprod(VC$X3*c(d2l.sigma21.sigma22),VC$X4)
    rho.sigma21     <- crossprod(VC$X3*c(d2l.rho.sigma21),VC$X5)
    sigma22.sigma22 <- crossprod(VC$X4*c(d2l.sigma22.sigma22),VC$X4)
    rho.sigma22     <- crossprod(VC$X4*c(d2l.rho.sigma22),VC$X5)   
    rho.rho         <- crossprod(VC$X5*c(d2l.rho.rho),VC$X5)    
    
    
    H <- rbind( cbind( be1.be1        ,   be1.be2      ,   be1.sigma21,      be1.sigma22,       be1.rho    ), 
                cbind( t(be1.be2)     ,   be2.be2      ,   be2.sigma21,      be2.sigma22,       be2.rho    ), 
                cbind( t(be1.sigma21) , t(be2.sigma21) ,   sigma21.sigma21,  sigma21.sigma22,   rho.sigma21),
                cbind( t(be1.sigma22) , t(be2.sigma22) , t(sigma21.sigma22), sigma22.sigma22,   rho.sigma22),
                cbind( t(be1.rho)     ,   t(be2.rho)   ,   t(rho.sigma21),   t(rho.sigma22),    rho.rho    ) ) 
                
                 

}



res <- -sum(l.par)

 
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
  
  
  
if( VC$margins[1] == "LN" || VC$margins[2] == "LN"){
  
  
  if(VC$margins[1] == "LN") dHs1 <- distrHsAT(exp(respvec$y1), eta1, sigma21, 1, margin2=VC$margins[1], min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
  if(VC$margins[2] == "LN") dHs2 <- distrHsAT(exp(respvec$y2), eta2, sigma22, 1, margin2=VC$margins[2], min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)

  pdf1 <- dHs1$pdf2
  pdf2 <- dHs2$pdf2

  p1 <- dHs1$p2
  p2 <- dHs2$p2
  
  dH <- copgHsAT(p1, p2, teta, VC$BivD, Ln = TRUE, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr, par2 = VC$dof)

  c.copula2.be1be2 <- dH$c.copula2.be1be2
  
  l.ln <- -sum( VC$weights*( log(pdf1) + log(pdf2) + log(c.copula2.be1be2) ) )

  }
  
  
  
  

         list(value=res, gradient=G, hessian=H, S.h=S.h, S.h1=S.h1, S.h2=S.h2, l=S.res, l.ln = l.ln, l.par=l.par, ps = ps, 
              eta1=eta1, eta2=eta2, etad=etad, etas1 = etas1, etas2 = etas2, 
              BivD=VC$BivD,               p1 = p1, p2 = p2, pdf1 = pdf1, pdf2 = pdf2,          
              c.copula.be2 = c.copula.be2,
              c.copula.be1 = c.copula.be1,
              c.copula2.be1be2 = c.copula2.be1be2, 
              dl.dbe1          =dl.dbe1,       
              dl.dbe2          =dl.dbe2,       
              dl.dsigma21.st   =dl.dsigma21.st,
              dl.dsigma22.st   =dl.dsigma22.st,
              dl.dteta.st      =dl.dteta.st) 
              
              

  }
  

  
