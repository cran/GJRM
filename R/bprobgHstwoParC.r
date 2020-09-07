bprobgHstwoParC <- function(params, respvec, VC, ps, AT = FALSE){
p1 <- p2 <- pdf1 <- pdf2 <- c.copula.be2 <- c.copula.be1 <- c.copula2.be1be2 <- NA


  eta1 <- VC$X1%*%params[1:VC$X1.d2]
  eta2 <- VC$X2%*%params[(VC$X1.d2+1):(VC$X1.d2+VC$X2.d2)]
  nu <- etad <- etas1 <- etas2 <- etan <- etan1 <- etan2 <- NULL  
  
  pd1 <- probm(eta1, VC$margins[1], only.pr = FALSE, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
  pd2 <- probm(eta2, VC$margins[2], only.pr = FALSE, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
  
  p1 <- pd1$pr; d.n1 <- pd1$d.n 
  p2 <- pd2$pr; d.n2 <- pd2$d.n  
  
  
  der2p1.dereta12 <- pd1$der2p.dereta
  der2p2.dereta22 <- pd2$der2p.dereta
  
                         
  if(is.null(VC$X3)){
                 nu.st <- etan   <- params[(VC$X1.d2+VC$X2.d2+1)] 
                 teta.st <- etad <- params[(VC$X1.d2+VC$X2.d2+2)]
                     }
  
  
  if(!is.null(VC$X3)){
                 nu.st <- etan   <- VC$X3%*%params[(VC$X1.d2+VC$X2.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2)]
                 teta.st <- etad <- VC$X4%*%params[(VC$X1.d2+VC$X2.d2+VC$X3.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2+VC$X4.d2)]
                      }
 

########################################################################################################  
  
nu.stt <- dof.tr(nu.st) 
nu.st  <- nu.stt$var.st    
nu     <- nu.stt$vao  
  
  
resT    <- teta.tr(VC, teta.st)
teta.st <- resT$teta.st
teta    <- resT$teta 
    
p11 <- mm(BiCDF(p1, p2, VC$nC, teta, nu), min.pr = VC$min.pr, max.pr = VC$max.pr  )

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

dH <- copgHs(p1, p2, eta1=NULL, eta2=NULL, teta, teta.st, VC$BivD, nu, nu.st, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)

c.copula.be1   <- dH$c.copula.be1
c.copula.be2   <- dH$c.copula.be2
c.copula.theta <- dH$c.copula.theta 
c.copula.dof   <- dH$c.copula.dof.st

c.copula2.be1    <- dH$c.copula2.be1   
c.copula2.be2    <- dH$c.copula2.be2 
c.copula2.be1be2 <- dH$c.copula2.be1be2

c.copula2.be1dof <- dH$c.copula2.be1dof.st 
c.copula2.be2dof <- dH$c.copula2.be2dof.st
c.copula2.be1th  <- dH$c.copula2.be1th 

c.copula2.dof2     <- dH$c.copula2.dof2.st
c.copula2.tetadof  <- dH$c.copula2.thdof.st
c.copula2.be2th    <- dH$c.copula2.be2th



bit1.b1b1 <- c.copula2.be1*d.n1^2 + c.copula.be1*der2p1.dereta12                                                                                                                                
bit2.b1b1 <- -c.copula2.be1*d.n1^2  + (1-c.copula.be1)*der2p1.dereta12
bit3.b1b1 <- -bit1.b1b1
bit4.b1b1 <- -bit2.b1b1


bit1.b2b2 <- c.copula2.be2*d.n2^2 + c.copula.be2*der2p2.dereta22                                                                                                                             
bit2.b2b2 <- -bit1.b2b2
bit3.b2b2 <- -c.copula2.be2*d.n2^2 + (1-c.copula.be2)*der2p2.dereta22
bit4.b2b2 <- -bit3.b2b2

bit1.b1b2 <- c.copula2.be1be2*d.n1*d.n2
bit2.b1b2 <- -bit1.b1b2
bit3.b1b2 <- -bit1.b1b2
bit4.b1b2 <- bit1.b1b2

bit1.b1th <- c.copula2.be1th*d.n1
bit2.b1th <- -bit1.b1th 
bit3.b1th <- -bit1.b1th 
bit4.b1th <- bit1.b1th 

bit1.b2th <- c.copula2.be2th*d.n2
bit2.b2th <- -bit1.b2th 
bit3.b2th <- -bit1.b2th 
bit4.b2th <- bit1.b2th 

if(AT==TRUE) bit1.th2 <- dH$bit1.th2ATE else bit1.th2 <- dH$bit1.th2

bit2.th2 <- -bit1.th2 
bit3.th2 <- -bit1.th2 
bit4.th2 <- bit1.th2 




bit1.b1dof <- c.copula2.be1dof*d.n1
bit2.b1dof <- -bit1.b1dof 
bit3.b1dof <- -bit1.b1dof 
bit4.b1dof <- bit1.b1dof 

bit1.b2dof <- c.copula2.be2dof*d.n2
bit2.b2dof <- -bit1.b2dof 
bit3.b2dof <- -bit1.b2dof 
bit4.b2dof <- bit1.b2dof 


#Not sure about these below for the signs

bit1.dof2 <- c.copula2.dof2
bit2.dof2 <- -bit1.dof2 
bit3.dof2 <- -bit1.dof2 
bit4.dof2 <-  bit1.dof2


bit1.thdof <- c.copula2.tetadof
bit2.thdof <- -bit1.thdof 
bit3.thdof <- -bit1.thdof 
bit4.thdof <-  bit1.thdof 









##########################################################################################

  dl.dbe1 <-  VC$weights*d.n1*( (respvec$y1.y2*c.copula.be1/p11)  +
                      (respvec$y1.cy2*(1-c.copula.be1)/p10) +
                      (respvec$cy1.y2*c.copula.be1/(-p01)) +
                      (respvec$cy1.cy2*(c.copula.be1-1)/p00) )
                                
  dl.dbe2 <-  VC$weights*d.n2*( (respvec$y1.y2*c.copula.be2/p11)  +
                            (respvec$y1.cy2*c.copula.be2/(-p10)) +
                                (respvec$cy1.y2*(1-c.copula.be2)/(p01)) +
                                (respvec$cy1.cy2*(c.copula.be2-1)/p00) )

  dl.drho <- VC$weights*( respvec$y1.y2*c.copula.theta/p11+respvec$y1.cy2*(-c.copula.theta)/p10 + 
                       respvec$cy1.y2*(-c.copula.theta)/p01+respvec$cy1.cy2*c.copula.theta/p00 ) 

  dl.ddof <- VC$weights*( respvec$y1.y2*c.copula.dof/p11+respvec$y1.cy2*(-c.copula.dof)/p10 + 
                       respvec$cy1.y2*(-c.copula.dof)/p01+respvec$cy1.cy2*c.copula.dof/p00 )

##########################################################################################
add.b  <- 1
if(AT==TRUE){
    if(VC$BivD %in% c("N") )                                    add.b <- 1/cosh(teta.st)^2
    if(VC$BivD %in% c("C0", "C180","J0", "J180","G0", "G180") ) add.b <-  exp(teta.st)     
    if(VC$BivD %in% c("C90","C270","J90","J270","G90","G270") ) add.b <- -exp(teta.st)        
}
##########################################################################################


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

  d2l.be1.rho  <- -VC$weights*(respvec$y1.y2*(bit1.b1th*p11-(c.copula.be1*d.n1*c.copula.theta))/p11^2+
                              respvec$y1.cy2*(bit2.b1th*p10-((1-c.copula.be1)*d.n1)*(-c.copula.theta))/p10^2+
                             respvec$cy1.y2*(bit3.b1th*p01-(-c.copula.be1*d.n1)*(-c.copula.theta))/p01^2+
                              respvec$cy1.cy2*(bit4.b1th*p00-((c.copula.be1-1)*d.n1)*c.copula.theta)/p00^2 )/add.b

  d2l.be2.rho  <- -VC$weights*(respvec$y1.y2*(bit1.b2th*p11-(c.copula.be2*d.n2*c.copula.theta))/p11^2+
                              respvec$y1.cy2*(bit2.b2th*p10-(-c.copula.be2*d.n2)*(-c.copula.theta))/p10^2+
                             respvec$cy1.y2*(bit3.b2th*p01-((1-c.copula.be2)*d.n2)*(-c.copula.theta))/p01^2+
                              respvec$cy1.cy2*(bit4.b2th*p00-((c.copula.be2-1)*d.n2)*c.copula.theta)/p00^2 )/add.b
                              
  d2l.rho.rho  <- (-VC$weights*(   respvec$y1.y2*(bit1.th2*p11-( c.copula.theta/add.b)^2)/p11^2+
                               respvec$y1.cy2*(bit2.th2*p10-(-c.copula.theta/add.b)^2)/p10^2+
                               respvec$cy1.y2*(bit3.th2*p01-(-c.copula.theta/add.b)^2)/p01^2+
                              respvec$cy1.cy2*(bit4.th2*p00-( c.copula.theta/add.b)^2)/p00^2 ) )                            


  d2l.be1.dof  <- -VC$weights*(respvec$y1.y2*(bit1.b1dof*p11-(c.copula.be1*d.n1*c.copula.dof))/p11^2+
                              respvec$y1.cy2*(bit2.b1dof*p10-((1-c.copula.be1)*d.n1)*(-c.copula.dof))/p10^2+
                              respvec$cy1.y2*(bit3.b1dof*p01-(-c.copula.be1*d.n1)*(-c.copula.dof))/p01^2+
                              respvec$cy1.cy2*(bit4.b1dof*p00-((c.copula.be1-1)*d.n1)*c.copula.dof)/p00^2 )

  d2l.be2.dof  <- -VC$weights*(respvec$y1.y2*(bit1.b2dof*p11-(c.copula.be2*d.n2*c.copula.dof))/p11^2+
                               respvec$y1.cy2*(bit2.b2dof*p10-(-c.copula.be2*d.n2)*(-c.copula.dof))/p10^2+
                               respvec$cy1.y2*(bit3.b2dof*p01-((1-c.copula.be2)*d.n2)*(-c.copula.dof))/p01^2+
                               respvec$cy1.cy2*(bit4.b2dof*p00-((c.copula.be2-1)*d.n2)*c.copula.dof)/p00^2 )
                              
  d2l.dof.dof  <- (-VC$weights*(respvec$y1.y2*(bit1.dof2*p11-( c.copula.dof)^2)/p11^2+
                               respvec$y1.cy2*(bit2.dof2*p10-(-c.copula.dof)^2)/p10^2+
                               respvec$cy1.y2*(bit3.dof2*p01-(-c.copula.dof)^2)/p01^2+
                              respvec$cy1.cy2*(bit4.dof2*p00-( c.copula.dof)^2)/p00^2 ) )      

  d2l.rho.dof  <- (-VC$weights*(  respvec$y1.y2*(bit1.thdof*p11-( c.copula.theta/add.b *c.copula.dof ))/p11^2+
                                  respvec$y1.cy2*(bit2.thdof*p10-(-c.copula.theta/add.b * (-c.copula.dof)))/p10^2+
                                  respvec$cy1.y2*(bit3.thdof*p01-(-c.copula.theta/add.b * (-c.copula.dof)))/p01^2+
                                  respvec$cy1.cy2*(bit4.thdof*p00-( c.copula.theta/add.b * c.copula.dof )^2)/p00^2 ) )     







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

  d2l.be1.rho  <- -VC$weights*((bit1.b1th*p11-(c.copula.be1*d.n1*c.copula.theta))/p11+
                              (bit2.b1th*p10-((1-c.copula.be1)*d.n1)*(-c.copula.theta))/p10+
                              (bit3.b1th*p01-(-c.copula.be1*d.n1)*(-c.copula.theta))/p01+
                              (bit4.b1th*p00-((c.copula.be1-1)*d.n1)*c.copula.theta)/p00 )/add.b

  d2l.be2.rho  <- -VC$weights*(  (bit1.b2th*p11-(c.copula.be2*d.n2*c.copula.theta))/p11+
                              (bit2.b2th*p10-(-c.copula.be2*d.n2)*(-c.copula.theta))/p10+
                              (bit3.b2th*p01-((1-c.copula.be2)*d.n2)*(-c.copula.theta))/p01+
                              (bit4.b2th*p00-((c.copula.be2-1)*d.n2)*c.copula.theta)/p00 )/add.b

  d2l.rho.rho  <- -VC$weights*(  (bit1.th2*p11-( c.copula.theta/add.b)^2)/p11+
                              (bit2.th2*p10-(-c.copula.theta/add.b)^2)/p10+
                              (bit3.th2*p01-(-c.copula.theta/add.b)^2)/p01+
                              (bit4.th2*p00-( c.copula.theta/add.b)^2)/p00 )

  d2l.be1.dof  <- -VC$weights*(  (bit1.b1dof*p11-(c.copula.be1*d.n1*c.copula.dof))/p11+
                              (bit2.b1dof*p10-((1-c.copula.be1)*d.n1)*(-c.copula.dof))/p10+
                              (bit3.b1dof*p01-(-c.copula.be1*d.n1)*(-c.copula.dof))/p01+
                              (bit4.b1dof*p00-((c.copula.be1-1)*d.n1)*c.copula.dof)/p00 )

  d2l.be2.dof  <- -VC$weights*( (bit1.b2dof*p11-(c.copula.be2*d.n2*c.copula.dof))/p11+
                              (bit2.b2dof*p10-(-c.copula.be2*d.n2)*(-c.copula.dof))/p10+
                             (bit3.b2dof*p01-((1-c.copula.be2)*d.n2)*(-c.copula.dof))/p01+
                              (bit4.b2dof*p00-((c.copula.be2-1)*d.n2)*c.copula.dof)/p00 )
                              
  d2l.dof.dof  <- (-VC$weights*((bit1.dof2*p11-( c.copula.dof)^2)/p11+
                               (bit2.dof2*p10-(-c.copula.dof)^2)/p10+
                               (bit3.dof2*p01-(-c.copula.dof)^2)/p01+
                              (bit4.dof2*p00-( c.copula.dof)^2)/p00 ) )      

  d2l.rho.dof  <- (-VC$weights*(  (bit1.thdof*p11-( c.copula.theta/add.b *c.copula.dof ))/p11+
                                  (bit2.thdof*p10-(-c.copula.theta/add.b * (-c.copula.dof)))/p10+
                                  (bit3.thdof*p01-(-c.copula.theta/add.b * (-c.copula.dof)))/p01+
                                  (bit4.thdof*p00-( c.copula.theta/add.b * c.copula.dof )^2)/p00 ) )


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

  d2l.be1.rho  <- -VC$weights*(  fi*(bit1.b1th*p11-(c.copula.be1*d.n1*c.copula.theta))/p11^2*resp1+
                              se*(bit2.b1th*p10-((1-c.copula.be1)*d.n1)*(-c.copula.theta))/p10^2*resp2+
                              th*(bit3.b1th*p01-(-c.copula.be1*d.n1)*(-c.copula.theta))/p01^2*resp3+
                              fo*(bit4.b1th*p00-((c.copula.be1-1)*d.n1)*c.copula.theta)/p00^2*resp4 )/add.b

  d2l.be2.rho  <- -VC$weights*(  fi*(bit1.b2th*p11-(c.copula.be2*d.n2*c.copula.theta))/p11^2*resp1+
                              se*(bit2.b2th*p10-(-c.copula.be2*d.n2)*(-c.copula.theta))/p10^2*resp2+
                              th*(bit3.b2th*p01-((1-c.copula.be2)*d.n2)*(-c.copula.theta))/p01^2*resp3+
                              fo*(bit4.b2th*p00-((c.copula.be2-1)*d.n2)*c.copula.theta)/p00^2*resp4 )/add.b

  d2l.rho.rho  <- -VC$weights*(  fi*(bit1.th2*p11-( c.copula.theta/add.b)^2)/p11^2*resp1+
                              se*(bit2.th2*p10-(-c.copula.theta/add.b)^2)/p10^2*resp2+
                              th*(bit3.th2*p01-(-c.copula.theta/add.b)^2)/p01^2*resp3+
                              fo*(bit4.th2*p00-( c.copula.theta/add.b)^2)/p00^2*resp4 )



  d2l.be1.dof  <- -VC$weights*(  fi*(bit1.b1dof*p11-(c.copula.be1*d.n1*c.copula.dof))/p11^2*resp1+
                              se*(bit2.b1dof*p10-((1-c.copula.be1)*d.n1)*(-c.copula.dof))/p10^2*resp2+
                              th*(bit3.b1dof*p01-(-c.copula.be1*d.n1)*(-c.copula.dof))/p01^2*resp3+
                              fo*(bit4.b1dof*p00-((c.copula.be1-1)*d.n1)*c.copula.dof)/p00^2*resp4 )

  d2l.be2.dof  <- -VC$weights*(  fi*(bit1.b2dof*p11-(c.copula.be2*d.n2*c.copula.dof))/p11^2*resp1+
                              se*(bit2.b2dof*p10-(-c.copula.be2*d.n2)*(-c.copula.dof))/p10^2*resp2+
                              th*(bit3.b2dof*p01-((1-c.copula.be2)*d.n2)*(-c.copula.dof))/p01^2*resp3+
                              fo*(bit4.b2dof*p00-((c.copula.be2-1)*d.n2)*c.copula.dof)/p00^2*resp4 )

  d2l.dof.dof  <- -VC$weights*(  fi*(bit1.dof2*p11-( c.copula.dof)^2)/p11^2*resp1+
                              se*(bit2.dof2*p10-(-c.copula.dof)^2)/p10^2*resp2+
                              th*(bit3.dof2*p01-(-c.copula.dof)^2)/p01^2*resp3+
                              fo*(bit4.dof2*p00-( c.copula.dof)^2)/p00^2*resp4 )


  d2l.rho.dof  <- VC$weights*(  fi*(bit1.thdof*p11-( c.copula.theta/add.b * c.copula.dof))/p11^2*resp1+
                                se*(bit2.thdof*p10-(-c.copula.theta/add.b *(-c.copula.dof)))/p10^2*resp2+
                                th*(bit3.thdof*p01-(-c.copula.theta/add.b *(-c.copula.dof)))/p01^2*resp3+
                                fo*(bit4.thdof*p00-( c.copula.theta/add.b * c.copula.dof))/p00^2*resp4 )




}




if( is.null(VC$X3) ){

  be1.be1 <- crossprod(VC$X1*c(d2l.be1.be1),VC$X1)
  be2.be2 <- crossprod(VC$X2*c(d2l.be2.be2),VC$X2)
  be1.be2 <- crossprod(VC$X1*c(d2l.be1.be2),VC$X2)
  be1.rho <- t(t(rowSums(t(VC$X1*c(d2l.be1.rho)))))
  be2.rho <- t(t(rowSums(t(VC$X2*c(d2l.be2.rho)))))
  be1.dof <- t(t(rowSums(t(VC$X1*c(d2l.be1.dof)))))
  be2.dof <- t(t(rowSums(t(VC$X2*c(d2l.be2.dof)))))  
  
  
  H <- rbind( cbind( be1.be1    , be1.be2    , be1.dof, be1.rho ), 
              cbind( t(be1.be2) , be2.be2    , be2.dof, be2.rho ), 
              cbind( t(be1.dof) , t(be2.dof) , sum(d2l.dof.dof), sum(d2l.rho.dof) ),
              cbind( t(be1.rho) , t(be2.rho) , sum(d2l.rho.dof), sum(d2l.rho.rho) )
            ) 
            
           
         
         G   <- -c( colSums( c(dl.dbe1)*VC$X1 ),
                    colSums( c(dl.dbe2)*VC$X2 ),
                    sum( dl.ddof ), sum( dl.drho )  )
                                        
}

if( !is.null(VC$X3) ){

  be1.be1 <- crossprod(VC$X1*c(d2l.be1.be1),VC$X1)
  be2.be2 <- crossprod(VC$X2*c(d2l.be2.be2),VC$X2)
  be1.be2 <- crossprod(VC$X1*c(d2l.be1.be2),VC$X2)
  be1.rho <- crossprod(VC$X1*c(d2l.be1.rho),VC$X4)
  be2.rho <- crossprod(VC$X2*c(d2l.be2.rho),VC$X4)
  rho.rho <- crossprod(VC$X4*c(d2l.rho.rho),VC$X4)
  be1.dof <- crossprod(VC$X1*c(d2l.be1.dof),VC$X3)
  be2.dof <- crossprod(VC$X2*c(d2l.be2.dof),VC$X3)
  dof.dof <- crossprod(VC$X3*c(d2l.dof.dof),VC$X3)
  rho.dof <- crossprod(VC$X3*c(d2l.rho.dof),VC$X4)
  
  H <- rbind( cbind(   be1.be1    ,   be1.be2    ,   be1.dof   ,  be1.rho), 
              cbind( t(be1.be2)   ,   be2.be2    ,   be2.dof   ,  be2.rho), 
              cbind( t(be1.dof)   , t(be2.dof)   ,   dof.dof  ,  rho.dof),
              cbind( t(be1.rho)   , t(be2.rho)   , t(rho.dof) ,  rho.rho)
            )
            
           
         G   <- -c( colSums( c(dl.dbe1)*VC$X1 ),
                    colSums( c(dl.dbe2)*VC$X2 ),
                    colSums( c(dl.ddof)*VC$X3 ),
                    colSums( c(dl.drho)*VC$X4 )   )
    
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
              p11=p11, p10=p10, p01=p01, p00=p00, eta1=eta1, eta2=eta2, etad=etad, 
              dl.dbe1=dl.dbe1, dl.dbe2=dl.dbe2, dl.drho=dl.drho,
              BivD=VC$BivD, p1 = p1, p2 = p2, pdf1 = d.n1, pdf2 = d.n2,          
	      	                    c.copula.be2 = c.copula.be2,
	      	                    c.copula.be1 = c.copula.be1,
              c.copula2.be1be2 = c.copula2.be1be2, theta.star = teta.st,
              etan = etan, dl.dnu.st=dl.ddof)      

}




     























