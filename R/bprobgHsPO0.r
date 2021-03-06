bprobgHsPO0 <- function(params, respvec, VC, ps){
p1 <- p2 <- pdf1 <- pdf2 <- c.copula.be2 <- c.copula.be1 <- c.copula2.be1be2 <- NA


  eta1 <- VC$X1%*%params[1:VC$X1.d2]
  eta2 <- VC$X2%*%params[(VC$X1.d2+1):(VC$X1.d2+VC$X2.d2)]

  pd1 <- probm(eta1, VC$margins[1], only.pr = FALSE, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
  pd2 <- probm(eta2, VC$margins[2], only.pr = FALSE, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
  
  p1 <- pd1$pr; d.n1 <- pd1$d.n 
  p2 <- pd2$pr; d.n2 <- pd2$d.n  




########################################################################################################  
  
teta.st <- teta <- 0  
  
p11 <-  mm(BiCDF(p1, p2, VC$nC, teta), min.pr = VC$min.pr, max.pr = VC$max.pr  )

########################################################################################################

  cp11 <- mm(1 - p11, min.pr = VC$min.pr, max.pr = VC$max.pr)
  
  l.par <- VC$weights*( respvec$y1*log(p11) + respvec$cy*log(cp11) )


c.copula.be1 <- p2                          
c.copula.be2 <- p1

derf1.dereta1 <- pd1$der2p.dereta
derf2.dereta2 <- pd2$der2p.dereta

bit1.b1b1 <- c.copula.be1*derf1.dereta1 
bit1.b2b2 <- c.copula.be2*derf2.dereta2
bit1.b1b2 <- d.n1*d.n2

  dl.dbe1 <-  VC$weights*d.n1*( respvec$y1*c.copula.be1/p11   - respvec$cy*c.copula.be1/cp11   )                                     
  dl.dbe2 <-  VC$weights*d.n2*( respvec$y1*c.copula.be2/p11   - respvec$cy*c.copula.be2/cp11   )
                                                  

if(VC$hess==TRUE){

  d2l.be1.be1  <- -VC$weights*(respvec$y1*(bit1.b1b1*p11-(c.copula.be1*d.n1)^2)/p11^2 -
                           respvec$cy*(bit1.b1b1*cp11+(c.copula.be1*d.n1)^2)/cp11^2 )

  d2l.be2.be2  <- -VC$weights*(respvec$y1*(bit1.b2b2*p11 - (c.copula.be2*d.n2)^2 )/p11^2 -
                           respvec$cy*(bit1.b2b2*cp11 + (c.copula.be2*d.n2)^2 )/cp11^2 )

  d2l.be1.be2  <- -VC$weights*(respvec$y1*(bit1.b1b2*p11-(c.copula.be1*d.n1*c.copula.be2*d.n2))/p11^2 -
                           respvec$cy*(bit1.b1b2*cp11+(c.copula.be1*d.n1*c.copula.be2*d.n2))/cp11^2)


}  
  
if(VC$hess==FALSE){
                          
  d2l.be1.be1  <- -VC$weights*((bit1.b1b1*p11-(c.copula.be1*d.n1)^2)/p11 -
                            (bit1.b1b1*cp11+(c.copula.be1*d.n1)^2)/cp11 )

  d2l.be2.be2  <- -VC$weights*((bit1.b2b2*p11 - (c.copula.be2*d.n2)^2 )/p11 -
                            (bit1.b2b2*cp11 + (c.copula.be2*d.n2)^2 )/cp11 )

  d2l.be1.be2  <- -VC$weights*((bit1.b1b2*p11-(c.copula.be1*d.n1*c.copula.be2*d.n2))/p11 -
                            (bit1.b1b2*cp11+(c.copula.be1*d.n1*c.copula.be2*d.n2))/cp11)

                         
}

       
if( is.null(VC$X3) ){

  be1.be1 <- crossprod(VC$X1*c(d2l.be1.be1),VC$X1)
  be2.be2 <- crossprod(VC$X2*c(d2l.be2.be2),VC$X2)
  be1.be2 <- crossprod(VC$X1*c(d2l.be1.be2),VC$X2)

  
  H <- rbind( cbind( be1.be1    , be1.be2 ), 
              cbind( t(be1.be2) , be2.be2 )    
            ) 
            
            
    
         
         G   <- -c( colSums( c(dl.dbe1)*VC$X1 ),
                    colSums( c(dl.dbe2)*VC$X2 )
                  )
    
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
   
 
     
    

         list(value=res, gradient=G, hessian=H, S.h=S.h, S.h1=S.h1, S.h2=S.h2, l=S.res, l.par=l.par, ps = ps,
              p11=p11, cp11=cp11, eta1=eta1, eta2=eta2,  
              dl.dbe1=dl.dbe1, dl.dbe2=dl.dbe2,
              #d2l.be1.be1=d2l.be1.be1, d2l.be2.be2=d2l.be2.be2, 
              #d2l.be1.be2=d2l.be1.be2,
              BivD=VC$BivD, p1 = p1, p2 = p2, pdf1 = d.n1, pdf2 = d.n2,          
	      	                    c.copula.be2 = c.copula.be2,
	      	                    c.copula.be1 = c.copula.be1,
              c.copula2.be1be2 = c.copula2.be1be2)      

}




     























