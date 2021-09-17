bprobgHsDiscr1SS <- function(params, respvec, VC, ps, AT = FALSE){
p1 <- p2 <- pdf1 <- pdf2 <- c.copula.be2 <- c.copula.be1 <- c.copula2.be1be2 <- NA

 
   eta1 <- VC$X1%*%params[1:VC$X1.d2]
   eta2 <- VC$X2%*%params[(VC$X1.d2+1):(VC$X1.d2+VC$X2.d2)]
   etad <- etas <- l.ln <- NULL 
 

   
 if(is.null(VC$X3))  teta.st <- etad <- params[(VC$X1.d2 + VC$X2.d2 + 1)]
 if(!is.null(VC$X3)) teta.st <- etad <- VC$X3%*%params[(VC$X1.d2+VC$X2.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2)]
 
   eta2 <- eta.tr(eta2, VC$margins[2])
     
   
  dHs <- distrHsDiscr(respvec$y2, eta2, 1, 1, nu = 1, nu.st = 1, margin2=VC$margins[2], naive = FALSE, y2m = VC$y2m, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
   

    pdf2                         <- dHs$pdf2
    p2                           <- dHs$p2 
    derpdf2.dereta2              <- dHs$derpdf2.dereta2 
    derp2.dereta2                <- dHs$derp2.dereta2
    der2p2.dereta2eta2           <- dHs$der2p2.dereta2eta2 
    der2pdf2.dereta2             <- dHs$der2pdf2.dereta2
     
     
     pd1 <- probm(eta1, VC$margins[1], bc = TRUE, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr) 
     p1  <- 1 - pd1$pr                           #   pnorm(-eta1), p(y1=0)
   
   
   ########################################################################################################  
     
   resT    <- teta.tr(VC, teta.st)
   teta.st <- resT$teta.st
   teta    <- resT$teta
       
   
   ########################################################################################################
   
   C1 <- mm(BiCDF(p1[VC$inde], p2,          VC$nC, teta, VC$dof), min.pr = VC$min.pr, max.pr = VC$max.pr  )
   C2 <- mm(BiCDF(p1[VC$inde], mm(p2-pdf2, min.pr = VC$min.pr, max.pr = VC$max.pr), VC$nC, teta, VC$dof), min.pr = VC$min.pr, max.pr = VC$max.pr  )
   
   A <- mm(C1 - C2, min.pr = VC$min.pr, max.pr = VC$max.pr)
   B <- mm( pdf2 - A, min.pr = VC$min.pr, max.pr = VC$max.pr)
   
   l.par1          <- log(p1)
   l.par1[VC$inde] <- log(B) 
   
   l.par <- VC$weights*l.par1    
     
   
   ########################################################################################################
    
   dH1 <- copgHs(p1[VC$inde], p2,          eta1=NULL, eta2=NULL, teta, teta.st, VC$BivD, VC$dof, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
   dH2 <- copgHs(p1[VC$inde], mm(p2-pdf2, min.pr = VC$min.pr, max.pr = VC$max.pr), eta1=NULL, eta2=NULL, teta, teta.st, VC$BivD, VC$dof, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr) 
    
   c.copula.be1.C1 <- dH1$c.copula.be1 
   c.copula.be1.C2 <- dH2$c.copula.be1 
   
   c.copula.be2.C1 <- dH1$c.copula.be2 
   c.copula.be2.C2 <- dH2$c.copula.be2 
   
   derp2m1.dereta2     <- derp2.dereta2 - derpdf2.dereta2
   
   c.copula.theta.C1 <- dH1$c.copula.theta 
   c.copula.theta.C2 <- dH2$c.copula.theta
   
   derp1.dereta1   <- pd1$derp1.dereta1    
   
   Cc <- c.copula.be1.C1 - c.copula.be1.C2 # mm(c.copula.be1.C1 - c.copula.be1.C2, min.pr = VC$min.pr, max.pr = VC$max.pr) 
   C  <- Cc*derp1.dereta1[VC$inde] 
   
   Cs    <- c.copula.theta.C1  - c.copula.theta.C2
   Cssb2 <- c.copula.be2.C1*derp2.dereta2 - c.copula.be2.C2*derp2m1.dereta2  
   
     dl.dbe11          <- 1/p1*derp1.dereta1 
     dl.dbe11[VC$inde] <- -C/B  
     dl.dbe1           <- VC$weights*dl.dbe11 
     dl.dteta.st  <- VC$weights[VC$inde]*(-Cs/B)                     
     dl.dbe2      <- VC$weights[VC$inde]*( (derpdf2.dereta2 - Cssb2)/B    )
    
    
   ######################################################################################################## 
    
     
     c.copula2.be1.C1 <- dH1$c.copula2.be1
     c.copula2.be1.C2 <- dH2$c.copula2.be1
     
     c.copula2.be2.C1 <- dH1$c.copula2.be2
     c.copula2.be2.C2 <- dH2$c.copula2.be2
     
     c.copula2.be1be2.C1 <- dH1$c.copula2.be1be2
     c.copula2.be1be2.C2 <- dH2$c.copula2.be1be2
     
     c.copula2.be2th.C1 <- dH1$c.copula2.be2th
     c.copula2.be2th.C2 <- dH2$c.copula2.be2th  
     
     
    
     der2p1.dereta1eta1 <- pd1$der2p1.dereta1eta1
     
     derC.dereta1 <- (c.copula2.be1.C1 - c.copula2.be1.C2)*derp1.dereta1[VC$inde]^2 + Cc*der2p1.dereta1eta1[VC$inde]
     
   
     c.copula2.theta.C1 <- dH1$bit1.th2ATE 
     c.copula2.theta.C2 <- dH2$bit1.th2ATE 
     
     c.copula.thet.C1 <- dH1$c.copula.thet # NO star
     c.copula.thet.C2 <- dH2$c.copula.thet
     
     derteta.derteta.st         <- dH1$derteta.derteta.st          
     der2teta.derteta.stteta.st <- dH1$der2teta.derteta.stteta.st
     
     derCs.dertheta.st <- (c.copula2.theta.C1 - c.copula2.theta.C2)*derteta.derteta.st^2 + (c.copula.thet.C1 - c.copula.thet.C2)*der2teta.derteta.stteta.st
     
     derA.dereta2 <- c.copula.be2.C1*derp2.dereta2   - c.copula.be2.C2*derp2m1.dereta2 
     derB.dereta2 <- derpdf2.dereta2 - derA.dereta2 
     
     der2p2m1.dereta2eta2 <- der2p2.dereta2eta2 - der2pdf2.dereta2
     derCssb2.dereta2 <- c.copula2.be2.C1*derp2.dereta2^2 + c.copula.be2.C1*der2p2.dereta2eta2 - (c.copula2.be2.C2*derp2m1.dereta2^2 + c.copula.be2.C2*der2p2m1.dereta2eta2)                                                                     
                  
     derC.dereta2 <- (c.copula2.be1be2.C1*derp2.dereta2 - c.copula2.be1be2.C2*derp2m1.dereta2)*derp1.dereta1[VC$inde]    
    
     c.copula2.be1th.C1 <- dH1$c.copula2.be1th 
     c.copula2.be1th.C2 <- dH2$c.copula2.be1th 
     
     derC.dertheta.st <- (c.copula2.be1th.C1 - c.copula2.be1th.C2)*derp1.dereta1[VC$inde] 
     
     derCs.dereta2      <- c.copula2.be2th.C1*derp2.dereta2 - c.copula2.be2th.C2*derp2m1.dereta2
     
     d2l.be1.be11          <- -1/p1^2*derp1.dereta1*derp1.dereta1 + 1/p1*der2p1.dereta1eta1  
     d2l.be1.be11[VC$inde] <- -C^2/B^2 - derC.dereta1/B    
     d2l.be1.be1           <- -VC$weights*d2l.be1.be11 
   
     d2l.rho.rho      <- -VC$weights[VC$inde]*(  -Cs^2/B^2  - derCs.dertheta.st/B )
     d2l.be2.be2      <- -VC$weights[VC$inde]*( -derB.dereta2/B^2*(derpdf2.dereta2-Cssb2) + ( der2pdf2.dereta2 - derCssb2.dereta2)/B  )
     d2l.be1.be2      <- -VC$weights[VC$inde]*( ( derB.dereta2/B^2)*C - derC.dereta2/B )
     d2l.be1.rho      <- -VC$weights[VC$inde]*(  -C*Cs/B^2 - derC.dertheta.st/B )
     d2l.be2.rho      <- -VC$weights[VC$inde]*(  derB.dereta2/B^2*Cs - derCs.dereta2/B )
     
  
  
  
  


if( is.null(VC$X3) ){

  be1.be1 <- crossprod(VC$X1*c(d2l.be1.be1),VC$X1)
  be2.be2 <- crossprod(VC$X2*c(d2l.be2.be2),VC$X2)
  be1.be2 <- crossprod(VC$X1[VC$inde,]*c(d2l.be1.be2),VC$X2)
  be1.rho <- t(t(rowSums(t(VC$X1[VC$inde,]*c(d2l.be1.rho)))))
  be2.rho <- t(t(rowSums(t(VC$X2*c(d2l.be2.rho)))))
  
  H <- rbind( cbind( be1.be1    , be1.be2    , be1.rho ), 
              cbind( t(be1.be2) , be2.be2    , be2.rho ), 
              cbind( t(be1.rho) , t(be2.rho) , sum(d2l.rho.rho) ) 
            ) 
            
           
         
  G   <- -c( colSums( c(dl.dbe1)*VC$X1 ) ,
             colSums( c(dl.dbe2)*VC$X2 ) ,
             sum( dl.dteta.st ) )
    
}

if( !is.null(VC$X3) ){





  be1.be1 <- crossprod(VC$X1*c(d2l.be1.be1),VC$X1)
  be2.be2 <- crossprod(VC$X2*c(d2l.be2.be2),VC$X2)
  be1.be2 <- crossprod(VC$X1[VC$inde,]*c(d2l.be1.be2),VC$X2)
  be1.rho <- crossprod(VC$X1[VC$inde,]*c(d2l.be1.rho),VC$X3)
  be2.rho <- crossprod(VC$X2*c(d2l.be2.rho),VC$X3)
  rho.rho <- crossprod(VC$X3*c(d2l.rho.rho),VC$X3)
  
  H <- rbind( cbind( be1.be1    , be1.be2    , be1.rho ), 
              cbind( t(be1.be2) , be2.be2    , be2.rho ), 
              cbind( t(be1.rho) , t(be2.rho) , rho.rho ) 
            ) 
            
           
  G   <- -c( colSums(      c(dl.dbe1)*VC$X1 ) ,
             colSums(      c(dl.dbe2)*VC$X2 ) ,
             colSums(  c(dl.dteta.st)*VC$X3 ) ) 
    
}

         res <- -sum(l.par)

if(VC$extra.regI == "pC" && VC$hess==FALSE) H <- regH(H, type = 1)
  
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
  
  



         list(value=res, gradient=G, hessian=H, S.h=S.h, S.h1=S.h1, S.h2=S.h2, l=S.res, l.par=l.par, ps = ps, etas = etas,
              eta1=eta1, eta2=eta2, etad=etad,
              dl.dbe1=dl.dbe1, dl.dbe2=dl.dbe2, dl.dteta.st = dl.dteta.st,
              BivD=VC$BivD,                             p1 = 1-p1, p2 = p2, pdf1 = pdf1, pdf2 = pdf2,          
	      	                    c.copula.be2 = c.copula.be2,
	      	                    c.copula.be1 = c.copula.be1,
              c.copula2.be1be2 = c.copula2.be1be2, theta.star = teta.st)      

}




     























