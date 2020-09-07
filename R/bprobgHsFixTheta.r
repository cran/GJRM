bprobgHsFixTheta <- function(params, respvec, VC, ps){

p1 <- p2 <- pdf1 <- pdf2 <- c.copula.be2 <- c.copula.be1 <- c.copula2.be1be2 <- NA

  eta1 <- VC$X1%*%params[1:VC$X1.d2]
  eta2 <- VC$X2%*%params[(VC$X1.d2+1):(VC$X1.d2+VC$X2.d2)]
  etad <- teta.st <- etad <- NULL 
  
  pd1 <- probm(eta1, VC$margins[1], only.pr = FALSE, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
  pd2 <- probm(eta2, VC$margins[2], only.pr = FALSE, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
  
  p1 <- pd1$pr; d.n1 <- pd1$d.n 
  p2 <- pd2$pr; d.n2 <- pd2$d.n  
  
  der2p1.dereta12 <- pd1$der2p.dereta
  der2p2.dereta22 <- pd2$der2p.dereta
                         
########################################################################################################  
  
p11 <- mm(BiCDF(p1, p2, VC$nC, VC$theta), min.pr = VC$min.pr, max.pr = VC$max.pr  )

########################################################################################################

  p10 <- mm(p1 - p11, min.pr = VC$min.pr, max.pr = VC$max.pr)
  p01 <- mm(p2 - p11, min.pr = VC$min.pr, max.pr = VC$max.pr)
  p00 <- mm(1 - p11 - p10 - p01, min.pr = VC$min.pr, max.pr = VC$max.pr)

  sall.p <- rowSums(cbind(p11, p10, p01, p00))
  
  p11 <- p11/sall.p 
  p10 <- p10/sall.p
  p01 <- p01/sall.p
  p00 <- p00/sall.p



  l.par <- VC$weights*( respvec$y1.y2*log(p11) + respvec$y1.cy2*log(p10) + respvec$cy1.y2*log(p01) + respvec$cy1.cy2*log(p00) )

dH <- copgHs(p1,p2,eta1=NULL,eta2=NULL,VC$theta,VC$theta,VC$BivD, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)

c.copula.be1   <- dH$c.copula.be1
c.copula.be2   <- dH$c.copula.be2
  
c.copula2.be1 <- dH$c.copula2.be1   
c.copula2.be2 <- dH$c.copula2.be2 


bit1.b1b1 <- c.copula2.be1*d.n1^2 + c.copula.be1*der2p1.dereta12                                                                                                                                
bit2.b1b1 <- -c.copula2.be1*d.n1^2  + (1-c.copula.be1)*der2p1.dereta12
bit3.b1b1 <- -bit1.b1b1
bit4.b1b1 <- -bit2.b1b1


bit1.b2b2 <- c.copula2.be2*d.n2^2 + c.copula.be2*der2p2.dereta22                                                                                                                             
bit2.b2b2 <- -bit1.b2b2
bit3.b2b2 <- -c.copula2.be2*d.n2^2 + (1-c.copula.be2)*der2p2.dereta22
bit4.b2b2 <- -bit3.b2b2

c.copula2.be1be2 <- dH$c.copula2.be1be2
bit1.b1b2 <- c.copula2.be1be2*d.n1*d.n2
bit2.b1b2 <- -bit1.b1b2
bit3.b1b2 <- -bit1.b1b2
bit4.b1b2 <- bit1.b1b2


  dl.dbe1 <-  VC$weights*d.n1*( (respvec$y1.y2*c.copula.be1/p11)  +
                      (respvec$y1.cy2*(1-c.copula.be1)/p10) +
                      (respvec$cy1.y2*c.copula.be1/(-p01)) +
                      (respvec$cy1.cy2*(c.copula.be1-1)/p00) )
                                
  dl.dbe2 <-  VC$weights*d.n2*( (respvec$y1.y2*c.copula.be2/p11)  +
                            (respvec$y1.cy2*c.copula.be2/(-p10)) +
                                (respvec$cy1.y2*(1-c.copula.be2)/(p01)) +
                                (respvec$cy1.cy2*(c.copula.be2-1)/p00) )



if(VC$hess==TRUE){

  d2l.be1.be1  <- -VC$weights*(respvec$y1.y2*(bit1.b1b1*p11-(c.copula.be1*d.n1)^2)/p11^2+
                              respvec$y1.cy2*(bit2.b1b1*p10-((1-c.copula.be1)*d.n1)^2)/p10^2+
                              respvec$cy1.y2*(bit3.b1b1*p01-(-c.copula.be1*d.n1)^2)/p01^2+
                              respvec$cy1.cy2*(bit4.b1b1*p00-((c.copula.be1-1)*d.n1)^2)/p00^2 )

  d2l.be2.be2  <- -VC$weights*(respvec$y1.y2*(bit1.b2b2*p11 - (c.copula.be2*d.n2)^2 )/p11^2+
                              respvec$y1.cy2*(bit2.b2b2*p10-(-c.copula.be2*d.n2)^2)/p10^2+
                              respvec$cy1.y2*(bit3.b2b2*p01- ((1-c.copula.be2)*d.n2)^2 )/p01^2+
                              respvec$cy1.cy2*(bit4.b2b2*p00-((c.copula.be2-1)*d.n2)^2)/p00^2 )

  d2l.be1.be2  <- -VC$weights*(respvec$y1.y2*(bit1.b1b2*p11-(c.copula.be1*d.n1*c.copula.be2*d.n2))/p11^2+
                              respvec$y1.cy2*(bit2.b1b2*p10-((1-c.copula.be1)*d.n1)*(-c.copula.be2*d.n2))/p10^2+
                              respvec$cy1.y2*(bit3.b1b2*p01-(-c.copula.be1*d.n1*(1-c.copula.be2)*d.n2))/p01^2+
                              respvec$cy1.cy2*(bit4.b1b2*p00-((c.copula.be1-1)*d.n1*((c.copula.be2-1)*d.n2)))/p00^2 )

                          
}     


if(VC$hess==FALSE && VC$end==0){

  d2l.be1.be1  <- -VC$weights*(  (bit1.b1b1*p11-(c.copula.be1*d.n1)^2)/p11+
                              (bit2.b1b1*p10-((1-c.copula.be1)*d.n1)^2)/p10+
                              (bit3.b1b1*p01-(-c.copula.be1*d.n1)^2)/p01+
                              (bit4.b1b1*p00-((c.copula.be1-1)*d.n1)^2)/p00 )

  d2l.be2.be2  <- -VC$weights*((bit1.b2b2*p11 - (c.copula.be2*d.n2)^2 )/p11+
                              (bit2.b2b2*p10-(-c.copula.be2*d.n2)^2)/p10+
                              (bit3.b2b2*p01- ((1-c.copula.be2)*d.n2)^2 )/p01+
                              (bit4.b2b2*p00-((c.copula.be2-1)*d.n2)^2)/p00 )

  d2l.be1.be2  <- -VC$weights*((bit1.b1b2*p11-(c.copula.be1*d.n1*c.copula.be2*d.n2))/p11+
                              (bit2.b1b2*p10-((1-c.copula.be1)*d.n1)*(-c.copula.be2*d.n2))/p10+
                              (bit3.b1b2*p01-(-c.copula.be1*d.n1*(1-c.copula.be2)*d.n2))/p01+
                              (bit4.b1b2*p00-((c.copula.be1-1)*d.n1*((c.copula.be2-1)*d.n2)))/p00 )


}  

if(VC$hess==FALSE && (VC$end==1 || VC$end==2)){

if(VC$end==1){

fi <- p11/p1
se <- p10/p1
th <- p01/(1-p1)
fo <- p00/(1-p1)

resp1 <- respvec$y1
resp2 <- respvec$y1
resp3 <- 1 - respvec$y1 
resp4 <- 1 - respvec$y1 


}

if(VC$end==2){

fi <- p11/p2
se <- p10/(1-p2)
th <- p01/p2
fo <- p00/(1-p2)

resp1 <- respvec$y2
resp2 <- 1 - respvec$y2
resp3 <- respvec$y2 
resp4 <- 1 - respvec$y2

}

  d2l.be1.be1  <- -VC$weights*(  fi*(bit1.b1b1*p11-(c.copula.be1*d.n1)^2)/p11^2*resp1+
                              se*(bit2.b1b1*p10-((1-c.copula.be1)*d.n1)^2)/p10^2*resp2+
                              th*(bit3.b1b1*p01-(-c.copula.be1*d.n1)^2)/p01^2*resp3+
                              fo*(bit4.b1b1*p00-((c.copula.be1-1)*d.n1)^2)/p00^2*resp4 )

  d2l.be2.be2  <- -VC$weights*(  fi*(bit1.b2b2*p11 - (c.copula.be2*d.n2)^2 )/p11^2*resp1+
                              se*(bit2.b2b2*p10-(-c.copula.be2*d.n2)^2)/p10^2*resp2+
                              th*(bit3.b2b2*p01- ((1-c.copula.be2)*d.n2)^2 )/p01^2*resp3+
                              fo*(bit4.b2b2*p00-((c.copula.be2-1)*d.n2)^2)/p00^2*resp4 )

  d2l.be1.be2  <- -VC$weights*(  fi*(bit1.b1b2*p11-(c.copula.be1*d.n1*c.copula.be2*d.n2))/p11^2*resp1+
                              se*(bit2.b1b2*p10-((1-c.copula.be1)*d.n1)*(-c.copula.be2*d.n2))/p10^2*resp2+
                              th*(bit3.b1b2*p01-(-c.copula.be1*d.n1*(1-c.copula.be2)*d.n2))/p01^2*resp3+
                              fo*(bit4.b1b2*p00-((c.copula.be1-1)*d.n1*((c.copula.be2-1)*d.n2)))/p00^2*resp4 )


}






  be1.be1 <- crossprod(VC$X1*c(d2l.be1.be1),VC$X1)
  be2.be2 <- crossprod(VC$X2*c(d2l.be2.be2),VC$X2)
  be1.be2 <- crossprod(VC$X1*c(d2l.be1.be2),VC$X2)

  
  H <- rbind( cbind( be1.be1    , be1.be2     ), 
              cbind( t(be1.be2) , be2.be2     ) 
            ) 
            
           
         
         G   <- -c( colSums( c(dl.dbe1)*VC$X1 ),
                    colSums( c(dl.dbe2)*VC$X2 )  )
    




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
              p11=p11, p10=p10, p01=p01, p00=p00, eta1=eta1, eta2=eta2, etad = etad, 
              dl.dbe1=dl.dbe1, dl.dbe2=dl.dbe2, dl.drho = NULL,
              BivD=VC$BivD,                             p1 = p1, p2 = p2, pdf1 = d.n1, pdf2 = d.n2,          
	      	                    c.copula.be2 = c.copula.be2,
	      	                    c.copula.be1 = c.copula.be1,
              c.copula2.be1be2 = c.copula2.be1be2, theta.star = 0)      

}




     























