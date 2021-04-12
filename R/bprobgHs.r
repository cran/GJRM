bprobgHs <- function(params, respvec, VC, ps, AT = FALSE){

p1 <- p2 <- pdf1 <- pdf2 <- c.copula.be2 <- c.copula.be1 <- c.copula2.be1be2 <- NA

  eta1 <- VC$X1%*%params[1:VC$X1.d2]
  eta2 <- VC$X2%*%params[(VC$X1.d2+1):(VC$X1.d2+VC$X2.d2)]
  etad <- NULL 
  
  pd1 <- probm(eta1, VC$margins[1], only.pr = FALSE, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
  pd2 <- probm(eta2, VC$margins[2], only.pr = FALSE, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
  
  p1 <- pd1$pr; d.n1 <- pd1$d.n 
  p2 <- pd2$pr; d.n2 <- pd2$d.n  
  
  
  der2p1.dereta12 <- pd1$der2p.dereta
  der2p2.dereta22 <- pd2$der2p.dereta
  
                         
  if(is.null(VC$X3))  teta.st <- etad <- params[(VC$X1.d2+VC$X2.d2+1)]
  if(!is.null(VC$X3)) teta.st <- etad <- VC$X3%*%params[(VC$X1.d2+VC$X2.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2)]
 

########################################################################################################  
  
resT    <- teta.tr(VC, teta.st)

teta.st1 <- teta.st2 <- teta.st <- resT$teta.st
teta1 <- teta2 <- teta <- resT$teta 
    
##################

Cop1 <- Cop2 <- VC$BivD 
nC1 <- nC2 <- VC$nC 


teta.ind1 <- as.logical(c(1,0,round(runif(VC$n-2))) ) 
teta.ind2 <- teta.ind1 == FALSE  


if(!(VC$BivD %in% VC$BivD2) && length(teta.st) > 1){

teta.st1 <- teta.st[teta.ind1]
teta.st2 <- teta.st[teta.ind2]

teta1 <- teta[teta.ind1]
teta2 <- teta[teta.ind2]

}

 
 
if(VC$BivD %in% VC$BivD2){

if(VC$BivD %in% VC$BivD2[c(1:4,13:16)])  teta.ind1 <- ifelse(VC$my.env$signind*teta > exp(VC$zerov), TRUE, FALSE)
if(VC$BivD %in% VC$BivD2[5:12]) teta.ind1 <- ifelse(VC$my.env$signind*teta > exp(VC$zerov) + 1, TRUE, FALSE) 
teta.ind2 <- teta.ind1 == FALSE 

VC$my.env$signind <- ifelse(teta.ind1 == TRUE,  1, -1) 

teta1 <-  teta[teta.ind1]
teta2 <- -teta[teta.ind2]

teta.st1 <- teta.st[teta.ind1]
teta.st2 <- teta.st[teta.ind2]

if(length(teta) == 1) teta.ind2 <- teta.ind1 <- rep(TRUE, VC$n)  

Cop1Cop2R <- Cop1Cop2(VC$BivD)
Cop1 <- Cop1Cop2R$Cop1
Cop2 <- Cop1Cop2R$Cop2

nC1 <- VC$ct[which(VC$ct[,1] == Cop1),2] 
nC2 <- VC$ct[which(VC$ct[,1] == Cop2),2]

} 



    
    
########################################################################################################

p11 <- NA

if( length(teta1) != 0) p11[teta.ind1] <- mm(BiCDF(p1[teta.ind1], p2[teta.ind1], nC1, teta1, VC$dof), min.pr = VC$min.pr, max.pr = VC$max.pr  )
if( length(teta2) != 0) p11[teta.ind2] <- mm(BiCDF(p1[teta.ind2], p2[teta.ind2], nC2, teta2, VC$dof), min.pr = VC$min.pr, max.pr = VC$max.pr  )

########################################################################################################

  # stuff below is pretty much a consistency check which will only be useful for those
  # odd model fits

  p10 <- mm(p1 - p11, min.pr = VC$min.pr, max.pr = VC$max.pr)
  p01 <- mm(p2 - p11, min.pr = VC$min.pr, max.pr = VC$max.pr)
  p00 <- mm(1 - p11 - p10 - p01, min.pr = VC$min.pr, max.pr = VC$max.pr)
  
  sall.p <- rowSums(cbind(p11, p10, p01, p00))
  
  p11 <- p11/sall.p 
  p10 <- p10/sall.p
  p01 <- p01/sall.p
  p00 <- p00/sall.p
  
  l.par <- VC$weights*( respvec$y1.y2*log(p11) + respvec$y1.cy2*log(p10) + respvec$cy1.y2*log(p01) + respvec$cy1.cy2*log(p00) )

########################################################################################################





c.copula.be1 <- c.copula.be2 <- c.copula.theta <- c.copula2.be1 <- c.copula2.be2 <- c.copula2.be1be2 <- c.copula2.be1th <- c.copula2.be2th <- bit1.th2 <- NA

if( length(teta1) != 0){

dH <- copgHs(p1[teta.ind1],p2[teta.ind1],eta1=NULL,eta2=NULL,teta1,teta.st1,Cop1, VC$dof, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)

c.copula.be1[teta.ind1]     <- dH$c.copula.be1
c.copula.be2[teta.ind1]     <- dH$c.copula.be2
c.copula.theta[teta.ind1]   <- dH$c.copula.theta 
c.copula2.be1[teta.ind1]    <- dH$c.copula2.be1   
c.copula2.be2[teta.ind1]    <- dH$c.copula2.be2 
c.copula2.be1be2[teta.ind1] <- dH$c.copula2.be1be2
c.copula2.be1th[teta.ind1]  <- dH$c.copula2.be1th 
c.copula2.be2th[teta.ind1]  <- dH$c.copula2.be2th

if(AT==TRUE) bit1.th2[teta.ind1] <- dH$bit1.th2ATE else bit1.th2[teta.ind1] <- dH$bit1.th2

}



if( length(teta2) != 0){

dH <- copgHs(p1[teta.ind2],p2[teta.ind2],eta1=NULL,eta2=NULL,teta2,teta.st2,Cop2, VC$dof, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)

c.copula.be1[teta.ind2]     <- dH$c.copula.be1
c.copula.be2[teta.ind2]     <- dH$c.copula.be2
c.copula.theta[teta.ind2]   <- dH$c.copula.theta 
c.copula2.be1[teta.ind2]    <- dH$c.copula2.be1   
c.copula2.be2[teta.ind2]    <- dH$c.copula2.be2 
c.copula2.be1be2[teta.ind2] <- dH$c.copula2.be1be2
c.copula2.be1th[teta.ind2]  <- dH$c.copula2.be1th 
c.copula2.be2th[teta.ind2]  <- dH$c.copula2.be2th

if(AT==TRUE) bit1.th2[teta.ind2] <- dH$bit1.th2ATE else bit1.th2[teta.ind2] <- dH$bit1.th2

}



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

bit2.th2 <- -bit1.th2 
bit3.th2 <- -bit1.th2 
bit4.th2 <- bit1.th2 

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


add.b  <- 1
if(AT==TRUE){
    if(VC$BivD %in% c("N") )                                    add.b <- 1/cosh(teta.st)^2
    if(VC$BivD %in% c("C0", "C180","GAL0", "GAL180","J0", "J180","G0", "G180") ) add.b <-  exp(teta.st)     
    if(VC$BivD %in% c("C90","C270","GAL90","GAL270","J90","J270","G90","G270") ) add.b <- -exp(teta.st)        
}



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
}




if( is.null(VC$X3) ){

  be1.be1 <- crossprod(VC$X1*c(d2l.be1.be1),VC$X1)
  be2.be2 <- crossprod(VC$X2*c(d2l.be2.be2),VC$X2)
  be1.be2 <- crossprod(VC$X1*c(d2l.be1.be2),VC$X2)
  be1.rho <- t(t(rowSums(t(VC$X1*c(d2l.be1.rho)))))
  be2.rho <- t(t(rowSums(t(VC$X2*c(d2l.be2.rho)))))
  
  H <- rbind( cbind( be1.be1    , be1.be2    , be1.rho ), 
              cbind( t(be1.be2) , be2.be2    , be2.rho ), 
              cbind( t(be1.rho) , t(be2.rho) , sum(d2l.rho.rho) ) 
            ) 
            
           
         
         G   <- -c( colSums( c(dl.dbe1)*VC$X1 ),
                    colSums( c(dl.dbe2)*VC$X2 ),
                    sum( dl.drho )  )
    
}

if( !is.null(VC$X3) ){

  be1.be1 <- crossprod(VC$X1*c(d2l.be1.be1),VC$X1)
  be2.be2 <- crossprod(VC$X2*c(d2l.be2.be2),VC$X2)
  be1.be2 <- crossprod(VC$X1*c(d2l.be1.be2),VC$X2)
  be1.rho <- crossprod(VC$X1*c(d2l.be1.rho),VC$X3)
  be2.rho <- crossprod(VC$X2*c(d2l.be2.rho),VC$X3)
  rho.rho <- crossprod(VC$X3*c(d2l.rho.rho),VC$X3)
  
  H <- rbind( cbind( be1.be1    , be1.be2    , be1.rho ), 
              cbind( t(be1.be2) , be2.be2    , be2.rho ), 
              cbind( t(be1.rho) , t(be2.rho) , rho.rho ) 
            ) 
            
           
         G   <- -c( colSums( c(dl.dbe1)*VC$X1 ),
                    colSums( c(dl.dbe2)*VC$X2 ),
                    colSums( c(dl.drho)*VC$X3 )  )
    
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
              BivD=VC$BivD, p1=p1, p2=p2, theta.star = teta.st,
              teta.ind2 = teta.ind2, teta.ind1 = teta.ind1,
              p1 = p1, p2 = p2, pdf1 = d.n1, pdf2 = d.n2,          
	                    c.copula.be2 = c.copula.be2,
	                    c.copula.be1 = c.copula.be1,
              c.copula2.be1be2 = c.copula2.be1be2, 
              Cop1 = Cop1, Cop2 = Cop2, teta1 = teta1, teta2 = teta2)      

}




     























