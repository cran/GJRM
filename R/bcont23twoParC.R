bcont23twoParC <- function(params, respvec, VC, ps, AT = FALSE){

    eta1 <- VC$X1%*%params[1:VC$X1.d2]
    eta2 <- VC$X2%*%params[(VC$X1.d2 + 1):(VC$X1.d2 + VC$X2.d2)]
    nu <- etad <- etas1 <- etas2 <- etan <- etan1 <- etan2 <- NULL 
  
    epsilon <- 0.0000001 
    max.p   <- 0.9999999
    
    
  if(is.null(VC$X3)){  
    sigma21.st <- etas1 <- params[(VC$X1.d2 + VC$X2.d2 + 1)]
    sigma22.st <- etas2 <- params[(VC$X1.d2 + VC$X2.d2 + 2)]
    nu2.st     <- etan2 <- params[(VC$X1.d2 + VC$X2.d2 + 3)]   
    nu.st      <- etan  <- params[(VC$X1.d2 + VC$X2.d2 + 4)] 
    teta.st    <- etad  <- params[(VC$X1.d2 + VC$X2.d2 + 5)]
  } 
  
  
  if(!is.null(VC$X3)){  
    sigma21.st <- etas1 <- VC$X3%*%params[(VC$X1.d2 + VC$X2.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2)]
    sigma22.st <- etas2 <- VC$X4%*%params[(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2)]
    nu2.st     <- etan2 <- VC$X5%*%params[(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2 + VC$X5.d2)]
    nu.st      <- etan  <- VC$X6%*%params[(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2 + VC$X5.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2 + VC$X5.d2 + VC$X6.d2 )]
    teta.st    <- etad  <- VC$X7%*%params[(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2 + VC$X5.d2 + VC$X6.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2 + VC$X5.d2 + VC$X6.d2 + VC$X7.d2)]
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
    

sstr2 <- esp.tr(nu2.st, VC$margins[2])  
nu2.st <- sstr2$vrb.st 
nu2    <- sstr2$vrb 

nu.stt <- dof.tr(nu.st) 
nu.st  <- nu.stt$var.st    
nu     <- nu.stt$vao 

eta1 <- eta.tr(eta1, VC$margins[1])
eta2 <- eta.tr(eta2, VC$margins[2])
    
  
resT    <- teta.tr(VC, teta.st)
teta.st <- resT$teta.st
teta    <- resT$teta  
  
##################
##################

  dHs1 <- distrHs(respvec$y1, eta1, sigma21, sigma21.st,   1,      1, margin2=VC$margins[1], naive = FALSE)
  dHs2 <- distrHs(respvec$y2, eta2, sigma22, sigma22.st, nu2, nu2.st, margin2=VC$margins[2], naive = FALSE)

  pdf1 <- dHs1$pdf2
  pdf2 <- dHs2$pdf2

  p1 <- dHs1$p2
  p2 <- dHs2$p2
  
  dH <- copgHsAT(p1, p2, teta, VC$BivD, Ln = TRUE, par2 = nu)

  c.copula2.be1be2 <- dH$c.copula2.be1be2
  
  l.par <- VC$weights*( log(pdf1) + log(pdf2) + log(c.copula2.be1be2) )
 
 
##################
##################

 derpdf1.dereta1              <- dHs1$derpdf2.dereta2 
 derpdf1.dersigma21.st        <- dHs1$derpdf2.dersigma2.st 
 derpdf1.dernu1.st            <- dHs1$derpdf2.dernu.st                  
 
 derpdf2.dereta2              <- dHs2$derpdf2.dereta2 
 derpdf2.dersigma22.st        <- dHs2$derpdf2.dersigma2.st 
 derpdf2.dernu2.st            <- dHs2$derpdf2.dernu.st
 
 derp1.dereta1                <- dHs1$derp2.dereta2
 derp1.dersigma21.st          <- dHs1$derp2.dersigma.st
 derp1.dernu1.st              <- dHs1$derp2.nu.st
 
 derp2.dereta2                <- dHs2$derp2.dereta2
 derp2.dersigma22.st          <- dHs2$derp2.dersigma.st
 derp2.dernu2.st              <- dHs2$derp2.nu.st
 
 
 BITS <- copgHsCont(p1, p2, teta, teta.st, VC$BivD, Cont = TRUE, par2 = nu, nu.st = nu.st)
 
 
   der2h.derp1p1              <- BITS$der2h.derp1p1
   derc.dereta1               <- der2h.derp1p1 * derp1.dereta1 
   derc.dersigma21.st         <- der2h.derp1p1 * derp1.dersigma21.st
   der2h.derp1p2              <- BITS$der2h.derp1p2  
   derc.dereta2               <- der2h.derp1p2 * derp2.dereta2    
   derc.dersigma22.st         <- der2h.derp1p2 * derp2.dersigma22.st
   der2h.derp1teta            <- BITS$der2h.derp1teta                                     
   derteta.derteta.st         <- BITS$derteta.derteta.st 
   der2h.derp1teta.st         <- der2h.derp1teta * derteta.derteta.st # new bit
   
   dernu.dernu.st             <- BITS$dernu.dernu.st
   der2h.derp1nu              <- BITS$der2h.derp1nu
   der2h.derp1nu.st           <- der2h.derp1nu * dernu.dernu.st
   
   derc.dernu1.st             <- der2h.derp1p1 * derp1.dernu1.st
   derc.dernu2.st             <- der2h.derp1p2 * derp2.dernu2.st
   
   
 
   dl.dbe1        <- VC$weights*( derpdf1.dereta1/pdf1 + derc.dereta1/c.copula2.be1be2              )   
   dl.dbe2        <- VC$weights*( derpdf2.dereta2/pdf2 + derc.dereta2/c.copula2.be1be2              )
   dl.dsigma21.st <- VC$weights*( derpdf1.dersigma21.st/pdf1 + derc.dersigma21.st/c.copula2.be1be2  )
   dl.dsigma22.st <- VC$weights*( derpdf2.dersigma22.st/pdf2 + derc.dersigma22.st/c.copula2.be1be2  )
   dl.dnu2.st     <- VC$weights*( derpdf2.dernu2.st/pdf2 + derc.dernu2.st/c.copula2.be1be2  )
   dl.dnu.st      <- VC$weights*( der2h.derp1nu.st/c.copula2.be1be2)   
   dl.dteta.st    <- VC$weights*( der2h.derp1teta.st/c.copula2.be1be2 )
  


               
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

der2pdf1.dernu1.st2  <- dHs1$der2pdf2.dernu.st2
der2pdf2.dernu2.st2  <- dHs2$der2pdf2.dernu.st2

der2p1.dereta1eta1 <- dHs1$der2p2.dereta2eta2
der2p2.dereta2eta2 <- dHs2$der2p2.dereta2eta2
 
der2p1.dersigma21.st2 <-  dHs1$der2p2.dersigma2.st2
der2p2.dersigma22.st2 <-  dHs2$der2p2.dersigma2.st2

der2p1.dernu1.st2 <-  dHs1$der2p2.dernu.st2 
der2p2.dernu2.st2 <-  dHs2$der2p2.dernu.st2 
 
der2pdf1.dereta1dersigma21.st <-  dHs1$der2pdf2.dereta2dersigma2.st
der2pdf2.dereta2dersigma22.st <-  dHs2$der2pdf2.dereta2dersigma2.st

der2pdf1.dereta1dernu1.st <-  dHs1$der2pdf2.dereta2dernu.st
der2pdf2.dereta2dernu2.st <-  dHs2$der2pdf2.dereta2dernu.st
 
der2p1.dereta1dersigma21.st <-  dHs1$der2p2.dereta2dersigma2.st
der2p2.dereta2dersigma22.st <-  dHs2$der2p2.dereta2dersigma2.st

der2p1.dereta1dernu1.st <-  dHs1$der2p2.dereta2dernu.st
der2p2.dereta2dernu2.st <-  dHs2$der2p2.dereta2dernu.st


der2pdf1.dersigma21.stdernu1.st <-  dHs1$der2pdf2.sigma2.st2dernu.st
der2pdf2.dersigma22.stdernu2.st <-  dHs2$der2pdf2.sigma2.st2dernu.st


der2p1.dersigma21.stdernu1.st <-  dHs1$der2p2.dersigma2.stdernu.st
der2p2.dersigma22.stdernu2.st <-  dHs2$der2p2.dersigma2.stdernu.st

der2teta.derteta.stteta.st <- BITS$der2teta.derteta.stteta.st 
der2nu.dernu.stnu.st       <- BITS$der2nu.dernu.stnu.st

der2c.derrho.dernu <- BITS$der2c.derrho.dernu  
der2c.dernu.dernu  <- BITS$der2c.dernu.dernu


der2c.derp1.dernu  <- BITS$der2c.derp1.dernu
der2c.derp2.dernu  <- BITS$der2c.derp2.dernu

                               
                               
  d2l.be1.be1      <-  -VC$weights*( (der2pdf1.dereta1 * pdf1 - derpdf1.dereta1^2) / pdf1^2   + 
                       ((der2c.derp1.derp1 * derp1.dereta1^2 + der2h.derp1p1 * der2p1.dereta1eta1) * c.copula2.be1be2 - derc.dereta1^2) /c.copula2.be1be2^2 )                   
  d2l.be2.be2      <-  -VC$weights*( (der2pdf2.dereta2 * pdf2 - derpdf2.dereta2^2) / pdf2^2   + 
                       ((der2c.derp2.derp2 * derp2.dereta2^2 + der2h.derp1p2 * der2p2.dereta2eta2) * c.copula2.be1be2 - derc.dereta2^2) /c.copula2.be1be2^2 )
                  
  d2l.rho.rho      <-   -VC$weights*( ((der2c.derrho.derrho*derteta.derteta.st^2 + der2h.derp1teta *der2teta.derteta.stteta.st)*c.copula2.be1be2 -
                          der2h.derp1teta.st^2) /c.copula2.be1be2^2  ) 
                                      
  d2l.sigma21.sigma21  <- -VC$weights*( (der2pdf1.dersigma21.st2 * pdf1 - derpdf1.dersigma21.st^2) / pdf1^2   + 
                          ((der2c.derp1.derp1*derp1.dersigma21.st^2+der2h.derp1p1*der2p1.dersigma21.st2)*c.copula2.be1be2 - derc.dersigma21.st^2) /c.copula2.be1be2^2 )                  
  d2l.sigma22.sigma22  <- -VC$weights*( (der2pdf2.dersigma22.st2 * pdf2 - derpdf2.dersigma22.st^2) / pdf2^2   + 
                          ((der2c.derp2.derp2*derp2.dersigma22.st^2+der2h.derp1p2*der2p2.dersigma22.st2)*c.copula2.be1be2 - derc.dersigma22.st^2) /c.copula2.be1be2^2 )
  
        
  d2l.nu2.nu2  <- -VC$weights*( (der2pdf2.dernu2.st2 * pdf2 - derpdf2.dernu2.st^2) / pdf2^2   + 
                            ((der2c.derp2.derp2*derp2.dernu2.st^2+der2h.derp1p2*der2p2.dernu2.st2)*c.copula2.be1be2 - derc.dernu2.st^2) /c.copula2.be1be2^2 )

  
  d2l.be1.be2    <- -VC$weights*( (der2c.derp1.derp2*derp1.dereta1*derp2.dereta2* c.copula2.be1be2 - derc.dereta1*derc.dereta2) /c.copula2.be1be2^2 )
  
  d2l.be1.rho   <- -VC$weights*( (der2c.derp1.derrho *derp1.dereta1*derteta.derteta.st* c.copula2.be1be2 - derc.dereta1*der2h.derp1teta*derteta.derteta.st)  /c.copula2.be1be2^2 )
  
d2l.be2.rho   <- -VC$weights*( (der2c.derp2.derrho *derp2.dereta2*derteta.derteta.st* c.copula2.be1be2 - derc.dereta2*der2h.derp1teta*derteta.derteta.st) /c.copula2.be1be2^2 )

  d2l.be1.sigma21  <- -VC$weights*( (der2pdf1.dereta1dersigma21.st * pdf1 - derpdf1.dereta1*derpdf1.dersigma21.st) / pdf1^2   + 
                       ((der2c.derp1.derp1 * derp1.dereta1* derp1.dersigma21.st+der2h.derp1p1*der2p1.dereta1dersigma21.st) * c.copula2.be1be2 - derc.dereta1*derc.dersigma21.st) /c.copula2.be1be2^2 )
                   
  d2l.be2.sigma22 <- -VC$weights*( (der2pdf2.dereta2dersigma22.st * pdf2 - derpdf2.dereta2*derpdf2.dersigma22.st) / pdf2^2   + 
                      ((der2c.derp2.derp2 * derp2.dereta2* derp2.dersigma22.st+der2h.derp1p2*der2p2.dereta2dersigma22.st) * c.copula2.be1be2 - derc.dereta2*derc.dersigma22.st) /c.copula2.be1be2^2 )
  
  
                     
  d2l.be2.nu2 <-  -VC$weights*( (der2pdf2.dereta2dernu2.st * pdf2 - derpdf2.dereta2*derpdf2.dernu2.st) / pdf2^2   + 
                         ((der2c.derp2.derp2 * derp2.dereta2* derp2.dernu2.st+der2h.derp1p2*der2p2.dereta2dernu2.st) * c.copula2.be1be2 - derc.dereta2*derc.dernu2.st) /c.copula2.be1be2^2 )
   
   
   
  d2l.be2.sigma21 <- -VC$weights*( (der2c.derp1.derp2*derp2.dereta2*derp1.dersigma21.st* c.copula2.be1be2 - derc.dereta2*derc.dersigma21.st)/c.copula2.be1be2^2)
  
  d2l.be1.sigma22 <- -VC$weights*( (der2c.derp1.derp2*derp1.dereta1*derp2.dersigma22.st* c.copula2.be1be2 - derc.dereta1*derc.dersigma22.st)/c.copula2.be1be2^2)
  
 d2l.be1.nu2 <- -VC$weights*( (der2c.derp1.derp2*derp1.dereta1*derp2.dernu2.st* c.copula2.be1be2 - derc.dereta1*derc.dernu2.st)/c.copula2.be1be2^2)
  
  
  d2l.rho.sigma21 <-  -VC$weights*( (der2c.derp1.derrho*derp1.dersigma21.st*derteta.derteta.st*c.copula2.be1be2 - derc.dersigma21.st*der2h.derp1teta *derteta.derteta.st)/c.copula2.be1be2^2 )
  
d2l.rho.sigma22 <-  -VC$weights*( (der2c.derp2.derrho*derp2.dersigma22.st*derteta.derteta.st*c.copula2.be1be2 - derc.dersigma22.st*der2h.derp1teta *derteta.derteta.st)/c.copula2.be1be2^2 )
  


d2l.rho.nu2 <-  -VC$weights*( (der2c.derp2.derrho*derp2.dernu2.st*derteta.derteta.st*c.copula2.be1be2 - derc.dernu2.st*der2h.derp1teta *derteta.derteta.st)/c.copula2.be1be2^2 )
  
  d2l.sigma21.sigma22 <- -VC$weights*( (der2c.derp1.derp2*derp1.dersigma21.st*derp2.dersigma22.st* c.copula2.be1be2 - derc.dersigma21.st*derc.dersigma22.st ) /c.copula2.be1be2^2 )
  

d2l.sigma21.nu2 <- -VC$weights*( (der2c.derp1.derp2*derp1.dersigma21.st*derp2.dernu2.st* c.copula2.be1be2 - derc.dersigma21.st*derc.dernu2.st)/c.copula2.be1be2^2)
                 
 
                  
  d2l.sigma22.nu2 <-  -VC$weights*( (der2pdf2.dersigma22.stdernu2.st * pdf2 - derpdf2.dersigma22.st*derpdf2.dernu2.st) / pdf2^2   + 
                          ((der2c.derp2.derp2 * derp2.dersigma22.st* derp2.dernu2.st+der2h.derp1p2*der2p2.dersigma22.stdernu2.st) * c.copula2.be1be2 - derc.dersigma22.st*derc.dernu2.st) /c.copula2.be1be2^2 )
   





d2l.nu.nu   <-  -VC$weights*( ((der2c.dernu.dernu*dernu.dernu.st^2 + der2h.derp1nu *der2nu.dernu.stnu.st)*c.copula2.be1be2 - der2h.derp1nu.st^2) /c.copula2.be1be2^2 )
d2l.rho.nu  <-  -VC$weights*( (der2c.derrho.dernu*derteta.derteta.st* dernu.dernu.st*c.copula2.be1be2 - der2h.derp1teta.st*der2h.derp1nu.st) /c.copula2.be1be2^2 )


d2l.nu.sigma21 <-  -VC$weights*( (der2c.derp1.dernu*derp1.dersigma21.st*dernu.dernu.st*c.copula2.be1be2 - derc.dersigma21.st*der2h.derp1nu *dernu.dernu.st)/c.copula2.be1be2^2 )

d2l.nu.sigma22 <-  -VC$weights*( (der2c.derp2.dernu*derp2.dersigma22.st*dernu.dernu.st*c.copula2.be1be2 - derc.dersigma22.st*der2h.derp1nu *dernu.dernu.st)/c.copula2.be1be2^2 )


d2l.nu.nu2 <-  -VC$weights*( (der2c.derp2.dernu*derp2.dernu2.st*dernu.dernu.st*c.copula2.be1be2 - derc.dernu2.st*der2h.derp1nu *dernu.dernu.st)/c.copula2.be1be2^2 )

d2l.be1.nu   <- -VC$weights*( (der2c.derp1.dernu *derp1.dereta1*dernu.dernu.st* c.copula2.be1be2 - derc.dereta1*der2h.derp1nu*dernu.dernu.st)  /c.copula2.be1be2^2 )

d2l.be2.nu   <- -VC$weights*( (der2c.derp2.dernu *derp2.dereta2*dernu.dernu.st* c.copula2.be1be2 - derc.dereta2*der2h.derp1nu*dernu.dernu.st) /c.copula2.be1be2^2 )




if( is.null(VC$X3) ){



  G   <- -c(colSums( c(dl.dbe1)*VC$X1 ),
            colSums( c(dl.dbe2)*VC$X2 ),
            sum( dl.dsigma21.st ),
            sum( dl.dsigma22.st ),
            sum( dl.dnu2.st ),
            sum(dl.dnu.st),
            sum( dl.dteta.st )            )



  be1.be1     <- crossprod(VC$X1*c(d2l.be1.be1),VC$X1)
  be2.be2     <- crossprod(VC$X2*c(d2l.be2.be2),VC$X2)
  be1.be2     <- crossprod(VC$X1*c(d2l.be1.be2),VC$X2)
  be1.rho     <- t(t(rowSums(t(VC$X1*c(d2l.be1.rho)))))
  be1.sigma21 <- t(t(rowSums(t(VC$X1*c(d2l.be1.sigma21))))) 
  be1.sigma22 <- t(t(rowSums(t(VC$X1*c(d2l.be1.sigma22))))) 
  be1.nu2     <- t(t(rowSums(t(VC$X1*c(d2l.be1.nu2))))) 
  be2.rho     <- t(t(rowSums(t(VC$X2*c(d2l.be2.rho)))))
  be2.sigma21 <- t(t(rowSums(t(VC$X2*c(d2l.be2.sigma21))))) 
  be2.sigma22 <- t(t(rowSums(t(VC$X2*c(d2l.be2.sigma22)))))
  be2.nu2     <- t(t(rowSums(t(VC$X2*c(d2l.be2.nu2)))))
  
  be1.nu     <- t(t(rowSums(t(VC$X1*c(d2l.be1.nu))))) 
  be2.nu     <- t(t(rowSums(t(VC$X2*c(d2l.be2.nu)))))


  
  
  

  H <- rbind( cbind(     be1.be1    , be1.be2        ,   be1.sigma21,            be1.sigma22,                           be1.nu2,              be1.nu, be1.rho  ), 
              cbind( t(be1.be2)     , be2.be2        ,   be2.sigma21,            be2.sigma22,                            be2.nu2,              be2.nu,  be2.rho  ), 
              cbind( t(be1.sigma21) , t(be2.sigma21) , sum(d2l.sigma21.sigma21), sum(d2l.sigma21.sigma22), sum(d2l.sigma21.nu2), sum(d2l.nu.sigma21), sum(d2l.rho.sigma21)),
              cbind( t(be1.sigma22) , t(be2.sigma22) , sum(d2l.sigma21.sigma22), sum(d2l.sigma22.sigma22), sum(d2l.sigma22.nu2), sum(d2l.nu.sigma22), sum(d2l.rho.sigma22)),
              cbind( t(be1.nu2)     , t(be2.nu2)     , sum(d2l.sigma21.nu2),     sum(d2l.sigma22.nu2),          sum(d2l.nu2.nu2),     sum(d2l.nu.nu2), sum(d2l.rho.nu2)),
              cbind( t(be1.nu)      ,    t(be2.nu)   ,    sum(d2l.nu.sigma21),   sum(d2l.nu.sigma22),           sum(d2l.nu.nu2),      sum(d2l.nu.nu),  sum(d2l.rho.nu) ),
              cbind( t(be1.rho)     ,   t(be2.rho)   ,   sum(d2l.rho.sigma21),   sum(d2l.rho.sigma22),          sum(d2l.rho.nu2),     sum(d2l.rho.nu), sum(d2l.rho.rho) ) 
              
              ) 

}



if( !is.null(VC$X3) ){



G   <- -c( colSums(       c(dl.dbe1)*VC$X1 ) ,
           colSums(       c(dl.dbe2)*VC$X2 ) ,
           colSums(c(dl.dsigma21.st)*VC$X3 ) ,
           colSums(c(dl.dsigma22.st)*VC$X4 ) ,
           colSums(    c(dl.dnu2.st)*VC$X5 ) ,
           colSums(     c(dl.dnu.st)*VC$X6 ) ,
           colSums(   c(dl.dteta.st)*VC$X7 )  )  

                 
    be1.be1         <- crossprod(VC$X1*c(d2l.be1.be1),VC$X1)
    be2.be2         <- crossprod(VC$X2*c(d2l.be2.be2),VC$X2)
    be1.be2         <- crossprod(VC$X1*c(d2l.be1.be2),VC$X2)
    be1.rho         <- crossprod(VC$X1*c(d2l.be1.rho),VC$X7)    
    be2.rho         <- crossprod(VC$X2*c(d2l.be2.rho),VC$X7)   
    be1.sigma21     <- crossprod(VC$X1*c(d2l.be1.sigma21),VC$X3)      
    be1.sigma22     <- crossprod(VC$X1*c(d2l.be1.sigma22),VC$X4) 
    
   
    be1.nu2     <- crossprod(VC$X1*c(d2l.be1.nu2),VC$X5)  
    

    be2.nu2   <- crossprod(VC$X2*c(d2l.be2.nu2),VC$X5)  
    
    be2.sigma21     <- crossprod(VC$X2*c(d2l.be2.sigma21),VC$X3)   
    be2.sigma22     <- crossprod(VC$X2*c(d2l.be2.sigma22),VC$X4)   
    
    sigma21.sigma21 <- crossprod(VC$X3*c(d2l.sigma21.sigma21),VC$X3)   
    

    nu2.nu2 <- crossprod(VC$X5*c(d2l.nu2.nu2),VC$X5)
    
    sigma21.sigma22 <- crossprod(VC$X3*c(d2l.sigma21.sigma22),VC$X4)
    



    rho.sigma21     <- crossprod(VC$X3*c(d2l.rho.sigma21),VC$X7)
    

    rho.nu2     <- crossprod(VC$X5*c(d2l.rho.nu2),VC$X7)

    sigma22.sigma22 <- crossprod(VC$X4*c(d2l.sigma22.sigma22),VC$X4)
    rho.sigma22     <- crossprod(VC$X4*c(d2l.rho.sigma22),VC$X7)   
    rho.rho         <- crossprod(VC$X7*c(d2l.rho.rho),VC$X7) 
    

    sigma21.nu2 <- crossprod(VC$X3*c(d2l.sigma21.nu2),VC$X5)

    sigma22.nu2 <- crossprod(VC$X4*c(d2l.sigma22.nu2),VC$X5)

    
   nu.sigma21  <- crossprod(VC$X3*c(d2l.nu.sigma21),VC$X6) 
   nu.sigma22  <- crossprod(VC$X4*c(d2l.nu.sigma22),VC$X6)
   nu.nu2      <- crossprod(VC$X5*c(d2l.nu.nu2),VC$X6)
   nu.nu       <- crossprod(VC$X6*c(d2l.nu.nu),VC$X6)
   rho.nu      <- crossprod(VC$X6*c(d2l.rho.nu),VC$X7)
  
   
      be1.nu     <- crossprod(VC$X1*c(d2l.be1.nu),VC$X6)
   be2.nu      <- crossprod(VC$X2*c(d2l.be2.nu),VC$X6)
  

    H <- rbind( cbind( be1.be1        ,   be1.be2      ,   be1.sigma21     , be1.sigma22        , be1.nu2    , be1.nu, be1.rho    ), 
                cbind( t(be1.be2)     ,   be2.be2      ,   be2.sigma21     , be2.sigma22       , be2.nu2    , be2.nu, be2.rho    ), 
                cbind( t(be1.sigma21) , t(be2.sigma21) ,   sigma21.sigma21 , sigma21.sigma22, sigma21.nu2, nu.sigma21, rho.sigma21),
                cbind( t(be1.sigma22) , t(be2.sigma22) , t(sigma21.sigma22), sigma22.sigma22,  sigma22.nu2, nu.sigma22, rho.sigma22),
                cbind( t(be1.nu2)     , t(be2.nu2)     , t(sigma21.nu2)    , t(sigma22.nu2) ,  nu2.nu2    , nu.nu2, rho.nu2    ),
                cbind( t(be1.nu)      ,    t(be2.nu)   ,    t(nu.sigma21),   t(nu.sigma22),     t(nu.nu2),     nu.nu,      rho.nu ),
                cbind( t(be1.rho)     ,   t(be2.rho)   , t(rho.sigma21)    , t(rho.sigma22) , t(rho.nu2) , t(rho.nu), rho.rho    )  ) 
                
                 

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
  
  
         list(value=res, gradient=G, hessian=H, S.h=S.h, S.h1=S.h1, S.h2=S.h2, l=S.res, l.par=l.par, ps = ps, 
              eta1=eta1, eta2=eta2, etad=etad, etas1 = etas1, etas2 = etas2, etan = etan, etan2 = etan2, 
              BivD=VC$BivD, p1 = p1, p2 = p2,
              dl.dbe1          =dl.dbe1,       
              dl.dbe2          =dl.dbe2,       
              dl.dsigma21.st   =dl.dsigma21.st,
              dl.dsigma22.st   =dl.dsigma22.st,
              dl.dnu2.st       =dl.dnu2.st,
              dl.dnu.st        =dl.dnu.st,
              dl.dteta.st      =dl.dteta.st) 
              

  }
  

  
