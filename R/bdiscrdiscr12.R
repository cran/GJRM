bdiscrdiscr12 <- function(params, respvec, VC, ps, AT = FALSE){
p1 <- p2 <- pdf1 <- pdf2 <- c.copula.be2 <- c.copula.be1 <- c.copula2.be1be2 <- NA

    eta1 <- VC$X1%*%params[1:VC$X1.d2]
    eta2 <- VC$X2%*%params[(VC$X1.d2 + 1):(VC$X1.d2 + VC$X2.d2)]
    etad <- etas1 <- etas2 <- l.ln <- NULL 
  
  if(is.null(VC$X3)){  
    sigma21.st <- etas1 <- 0
    sigma22.st <- etas2 <- params[(VC$X1.d2 + VC$X2.d2 + 1)]
    teta.st    <- etad  <- params[(VC$X1.d2 + VC$X2.d2 + 2)]
  } 
  
  
  if(!is.null(VC$X3)){  
    sigma21.st <- etas1 <- 0
    sigma22.st <- etas2 <- VC$X3%*%params[(VC$X1.d2 + VC$X2.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2)]
    teta.st    <- etad  <- VC$X4%*%params[(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2)]
  }  
  
  
##################
## Transformations
##################
  
    sstr1 <- 0
    sstr2 <- esp.tr(sigma22.st, VC$margins[2])  
  
    sigma21.st <- 0
    sigma22.st <- sstr2$vrb.st 
    
    sigma21    <- 0
    sigma22    <- sstr2$vrb 
    
    eta1 <- eta.tr(eta1, VC$margins[1])
    eta2 <- eta.tr(eta2, VC$margins[2])
    
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

##################

  dHs1 <- distrHsDiscr(respvec$y1, eta1, sigma21, sigma21.st, nu = 1, nu.st = 1, margin2=VC$margins[1], naive = FALSE, y2m = VC$y1m, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
  dHs2 <- distrHsDiscr(respvec$y2, eta2, sigma22, sigma22.st, nu = 1, nu.st = 1, margin2=VC$margins[2], naive = FALSE, y2m = VC$y2m, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)

  pdf1 <- dHs1$pdf2
  pdf2 <- dHs2$pdf2

  p1 <- dHs1$p2
  p2 <- dHs2$p2
  
C11 <- C01 <- C10 <- C00 <- NA

if( length(teta1) != 0){  
  C11[teta.ind1] <- mm(BiCDF(p1[teta.ind1], p2[teta.ind1], nC1, teta1, VC$dof), min.pr = VC$min.pr, max.pr = VC$max.pr  )
  C01[teta.ind1] <- mm(BiCDF(mm(p1[teta.ind1]-pdf1[teta.ind1], min.pr = VC$min.pr, max.pr = VC$max.pr), p2[teta.ind1], nC1, teta1, VC$dof), min.pr = VC$min.pr, max.pr = VC$max.pr  )
  C10[teta.ind1] <- mm(BiCDF(p1[teta.ind1], mm(p2[teta.ind1]-pdf2[teta.ind1], min.pr = VC$min.pr, max.pr = VC$max.pr), nC1, teta1, VC$dof), min.pr = VC$min.pr, max.pr = VC$max.pr  )
  C00[teta.ind1] <- mm(BiCDF(mm(p1[teta.ind1]-pdf1[teta.ind1], min.pr = VC$min.pr, max.pr = VC$max.pr), mm(p2[teta.ind1]-pdf2[teta.ind1], min.pr = VC$min.pr, max.pr = VC$max.pr), nC1, teta1, VC$dof), min.pr = VC$min.pr, max.pr = VC$max.pr  )
}  

if( length(teta2) != 0){
  C11[teta.ind2] <- mm(BiCDF(p1[teta.ind2], p2[teta.ind2], nC2, teta2, VC$dof), min.pr = VC$min.pr, max.pr = VC$max.pr  )
  C01[teta.ind2] <- mm(BiCDF(mm(p1[teta.ind2]-pdf1[teta.ind2], min.pr = VC$min.pr, max.pr = VC$max.pr), p2[teta.ind2], nC2, teta2, VC$dof), min.pr = VC$min.pr, max.pr = VC$max.pr  )
  C10[teta.ind2] <- mm(BiCDF(p1[teta.ind2], mm(p2[teta.ind2]-pdf2[teta.ind2], min.pr = VC$min.pr, max.pr = VC$max.pr), nC2, teta2, VC$dof), min.pr = VC$min.pr, max.pr = VC$max.pr  )
  C00[teta.ind2] <- mm(BiCDF(mm(p1[teta.ind2]-pdf1[teta.ind2], min.pr = VC$min.pr, max.pr = VC$max.pr), mm(p2[teta.ind2]-pdf2[teta.ind2], min.pr = VC$min.pr, max.pr = VC$max.pr), nC2, teta2, VC$dof), min.pr = VC$min.pr, max.pr = VC$max.pr  )
}

  E <- mm(C11 - C01 - C10 + C00, min.pr = VC$min.pr, max.pr = VC$max.pr) 


  l.par <- VC$weights*log(E)
  
  
##################



if( length(teta1) != 0){  

  dHC11F <- copgHs(p1[teta.ind1], p2[teta.ind1],                   eta1=NULL, eta2=NULL, teta1, teta.st1, Cop1, VC$dof, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
  dHC01F <- copgHs(mm(p1[teta.ind1]-pdf1[teta.ind1], min.pr = VC$min.pr, max.pr = VC$max.pr), p2[teta.ind1],          eta1=NULL, eta2=NULL, teta1, teta.st1, Cop1, VC$dof, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
  dHC10F <- copgHs(p1[teta.ind1], mm(p2[teta.ind1]-pdf2[teta.ind1], min.pr = VC$min.pr, max.pr = VC$max.pr),          eta1=NULL, eta2=NULL, teta1, teta.st1, Cop1, VC$dof, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
  dHC00F <- copgHs(mm(p1[teta.ind1]-pdf1[teta.ind1], min.pr = VC$min.pr, max.pr = VC$max.pr), mm(p2[teta.ind1]-pdf2[teta.ind1], min.pr = VC$min.pr, max.pr = VC$max.pr), eta1=NULL, eta2=NULL, teta1, teta.st1, Cop1, VC$dof, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
 
} 
 
 
 
if( length(teta2) != 0){  

  dHC11S <- copgHs(p1[teta.ind2], p2[teta.ind2],                   eta1=NULL, eta2=NULL, teta2, teta.st2, Cop2, VC$dof, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
  dHC01S <- copgHs(mm(p1[teta.ind2]-pdf1[teta.ind2], min.pr = VC$min.pr, max.pr = VC$max.pr), p2[teta.ind2],          eta1=NULL, eta2=NULL, teta2, teta.st2, Cop2, VC$dof, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
  dHC10S <- copgHs(p1[teta.ind2], mm(p2[teta.ind2]-pdf2[teta.ind2], min.pr = VC$min.pr, max.pr = VC$max.pr),          eta1=NULL, eta2=NULL, teta2, teta.st2, Cop2, VC$dof, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
  dHC00S <- copgHs(mm(p1[teta.ind2]-pdf1[teta.ind2], min.pr = VC$min.pr, max.pr = VC$max.pr), mm(p2[teta.ind2]-pdf2[teta.ind2], min.pr = VC$min.pr, max.pr = VC$max.pr), eta1=NULL, eta2=NULL, teta2, teta.st2, Cop2, VC$dof, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
 
} 
  
 
 
 
derC11.derp1 <- derC01.derp1 <- derC10.derp1 <- derC00.derp1 <- derC11.derp2 <- derC01.derp2 <- derC10.derp2 <- derC00.derp2 <- derC11.derthet <- derC01.derthet <- derC10.derthet <- derC00.derthet <- derteta.derteta.st <- NA 
 
 
 
if( length(teta1) != 0){ 
 
  derC11.derp1[teta.ind1]       <- dHC11F$c.copula.be1  
  derC01.derp1[teta.ind1]       <- dHC01F$c.copula.be1  
  derC10.derp1[teta.ind1]       <- dHC10F$c.copula.be1  
  derC00.derp1[teta.ind1]       <- dHC00F$c.copula.be1 
 
  derC11.derp2[teta.ind1]       <- dHC11F$c.copula.be2  
  derC01.derp2[teta.ind1]       <- dHC01F$c.copula.be2  
  derC10.derp2[teta.ind1]       <- dHC10F$c.copula.be2  
  derC00.derp2[teta.ind1]       <- dHC00F$c.copula.be2   
  
  derC11.derthet[teta.ind1]     <- dHC11F$c.copula.thet 
  derC01.derthet[teta.ind1]     <- dHC01F$c.copula.thet  
  derC10.derthet[teta.ind1]     <- dHC10F$c.copula.thet  
  derC00.derthet[teta.ind1]     <- dHC00F$c.copula.thet  
  
  derteta.derteta.st[teta.ind1] <- dHC11F$derteta.derteta.st
  
}  
  
  

if( length(teta2) != 0){ 
 
  derC11.derp1[teta.ind2]       <- dHC11S$c.copula.be1  
  derC01.derp1[teta.ind2]       <- dHC01S$c.copula.be1  
  derC10.derp1[teta.ind2]       <- dHC10S$c.copula.be1  
  derC00.derp1[teta.ind2]       <- dHC00S$c.copula.be1 
 
  derC11.derp2[teta.ind2]       <- dHC11S$c.copula.be2  
  derC01.derp2[teta.ind2]       <- dHC01S$c.copula.be2  
  derC10.derp2[teta.ind2]       <- dHC10S$c.copula.be2  
  derC00.derp2[teta.ind2]       <- dHC00S$c.copula.be2   
  
  derC11.derthet[teta.ind2]     <- dHC11S$c.copula.thet 
  derC01.derthet[teta.ind2]     <- dHC01S$c.copula.thet  
  derC10.derthet[teta.ind2]     <- dHC10S$c.copula.thet  
  derC00.derthet[teta.ind2]     <- dHC00S$c.copula.thet  
  
  derteta.derteta.st[teta.ind2] <- dHC11S$derteta.derteta.st
  
}  
   
  
  
  
  
  derpdf1.dereta1    <- dHs1$derpdf2.dereta2 
  derp1.dereta1      <- dHs1$derp2.dereta2
  derp1m1.dereta1    <- derp1.dereta1 - derpdf1.dereta1
    
  derpdf2.dereta2    <- dHs2$derpdf2.dereta2 
  derp2.dereta2      <- dHs2$derp2.dereta2
  derp2m1.dereta2    <- derp2.dereta2 - derpdf2.dereta2
  
  derpdf2.dersigma22.st  <- dHs2$derpdf2.dersigma2.st  
  derp2.dersigma22.st    <- dHs2$derp2.dersigma.st 
  derp2m1.dersigma22.st  <- derp2.dersigma22.st - derpdf2.dersigma22.st      

  fE1  <- (derC11.derp1 - derC10.derp1)*derp1.dereta1 - (derC01.derp1 - derC00.derp1)*derp1m1.dereta1    
  fE2  <- (derC11.derp2 - derC01.derp2)*derp2.dereta2 - (derC10.derp2 - derC00.derp2)*derp2m1.dereta2   
  fEt  <- (derC11.derthet - derC01.derthet - derC10.derthet + derC00.derthet)*derteta.derteta.st   
  fE2s <- (derC11.derp2 - derC01.derp2)*derp2.dersigma22.st - (derC10.derp2 - derC00.derp2)*derp2m1.dersigma22.st   

  dl.dbe1        <- VC$weights*fE1/E   
  dl.dbe2        <- VC$weights*fE2/E
  dl.dsigma22.st <- VC$weights*fE2s/E
  dl.dteta.st    <- VC$weights*fEt/E
  
  

##################


der2C11.derp1p1 <- der2C01.derp1p1 <- der2C10.derp1p1 <- der2C00.derp1p1 <- der2C11.derp2p2 <- der2C01.derp2p2 <- der2C10.derp2p2 <- der2C00.derp2p2 <- der2C11.derp1p2 <- der2C01.derp1p2 <- der2C10.derp1p2 <- der2C00.derp1p2 <- der2C11.derp1t <- der2C01.derp1t <- der2C10.derp1t <- der2C00.derp1t <- der2C11.derp2t <- der2C01.derp2t <- der2C10.derp2t <- der2C00.derp2t <- der2C11.derthet2 <- der2C01.derthet2 <- der2C10.derthet2 <- der2C00.derthet2 <- der2teta.derteta.stteta.st <- NA

if( length(teta1) != 0){ 

  der2C11.derp1p1[teta.ind1]            <- dHC11F$c.copula2.be1
  der2C01.derp1p1[teta.ind1]            <- dHC01F$c.copula2.be1  
  der2C10.derp1p1[teta.ind1]            <- dHC10F$c.copula2.be1  
  der2C00.derp1p1[teta.ind1]            <- dHC00F$c.copula2.be1 
  
  der2C11.derp2p2[teta.ind1]            <- dHC11F$c.copula2.be2
  der2C01.derp2p2[teta.ind1]            <- dHC01F$c.copula2.be2  
  der2C10.derp2p2[teta.ind1]            <- dHC10F$c.copula2.be2  
  der2C00.derp2p2[teta.ind1]            <- dHC00F$c.copula2.be2  
  
  der2C11.derp1p2[teta.ind1]            <- dHC11F$c.copula2.be1be2
  der2C01.derp1p2[teta.ind1]            <- dHC01F$c.copula2.be1be2 
  der2C10.derp1p2[teta.ind1]            <- dHC10F$c.copula2.be1be2  
  der2C00.derp1p2[teta.ind1]            <- dHC00F$c.copula2.be1be2
 
  der2C11.derp1t[teta.ind1]             <- dHC11F$c.copula2.be1t
  der2C01.derp1t[teta.ind1]             <- dHC01F$c.copula2.be1t 
  der2C10.derp1t[teta.ind1]             <- dHC10F$c.copula2.be1t  
  der2C00.derp1t[teta.ind1]             <- dHC00F$c.copula2.be1t 
   
  der2C11.derp2t[teta.ind1]             <- dHC11F$c.copula2.be2t
  der2C01.derp2t[teta.ind1]             <- dHC01F$c.copula2.be2t 
  der2C10.derp2t[teta.ind1]             <- dHC10F$c.copula2.be2t  
  der2C00.derp2t[teta.ind1]             <- dHC00F$c.copula2.be2t   
   
  der2C11.derthet2[teta.ind1]           <- dHC11F$bit1.th2ATE
  der2C01.derthet2[teta.ind1]           <- dHC01F$bit1.th2ATE  
  der2C10.derthet2[teta.ind1]           <- dHC10F$bit1.th2ATE  
  der2C00.derthet2[teta.ind1]           <- dHC00F$bit1.th2ATE   

  der2teta.derteta.stteta.st[teta.ind1] <- dHC11F$der2teta.derteta.stteta.st 
  
 } 
  



if( length(teta2) != 0){ 


  der2C11.derp1p1[teta.ind2]            <- dHC11S$c.copula2.be1
  der2C01.derp1p1[teta.ind2]            <- dHC01S$c.copula2.be1  
  der2C10.derp1p1[teta.ind2]            <- dHC10S$c.copula2.be1  
  der2C00.derp1p1[teta.ind2]            <- dHC00S$c.copula2.be1 
  
  der2C11.derp2p2[teta.ind2]            <- dHC11S$c.copula2.be2
  der2C01.derp2p2[teta.ind2]            <- dHC01S$c.copula2.be2  
  der2C10.derp2p2[teta.ind2]            <- dHC10S$c.copula2.be2  
  der2C00.derp2p2[teta.ind2]            <- dHC00S$c.copula2.be2  
  
  der2C11.derp1p2[teta.ind2]            <- dHC11S$c.copula2.be1be2
  der2C01.derp1p2[teta.ind2]            <- dHC01S$c.copula2.be1be2 
  der2C10.derp1p2[teta.ind2]            <- dHC10S$c.copula2.be1be2  
  der2C00.derp1p2[teta.ind2]            <- dHC00S$c.copula2.be1be2
 
  der2C11.derp1t[teta.ind2]             <- dHC11S$c.copula2.be1t
  der2C01.derp1t[teta.ind2]             <- dHC01S$c.copula2.be1t 
  der2C10.derp1t[teta.ind2]             <- dHC10S$c.copula2.be1t  
  der2C00.derp1t[teta.ind2]             <- dHC00S$c.copula2.be1t 
   
  der2C11.derp2t[teta.ind2]             <- dHC11S$c.copula2.be2t
  der2C01.derp2t[teta.ind2]             <- dHC01S$c.copula2.be2t 
  der2C10.derp2t[teta.ind2]             <- dHC10S$c.copula2.be2t  
  der2C00.derp2t[teta.ind2]             <- dHC00S$c.copula2.be2t   
   
  der2C11.derthet2[teta.ind2]           <- dHC11S$bit1.th2ATE
  der2C01.derthet2[teta.ind2]           <- dHC01S$bit1.th2ATE  
  der2C10.derthet2[teta.ind2]           <- dHC10S$bit1.th2ATE  
  der2C00.derthet2[teta.ind2]           <- dHC00S$bit1.th2ATE   

  der2teta.derteta.stteta.st[teta.ind2] <- dHC11S$der2teta.derteta.stteta.st 
  
 } 
  
    
  



  der2pdf1.dereta1     <- dHs1$der2pdf2.dereta2 
  der2p1.dereta1eta1   <- dHs1$der2p2.dereta2eta2
  der2p1m1.dereta1eta1 <- der2p1.dereta1eta1 - der2pdf1.dereta1

  der2pdf2.dereta2     <- dHs2$der2pdf2.dereta2
  der2p2.dereta2eta2   <- dHs2$der2p2.dereta2eta2
  der2p2m1.dereta2eta2 <- der2p2.dereta2eta2 - der2pdf2.dereta2 
  
  der2pdf2.dersigma22.st2 <- dHs2$der2pdf2.dersigma2.st2
  der2p2.dersigma22.st2   <- dHs2$der2p2.dersigma2.st2
  der2p2m1.dersigma22.st2 <- der2p2.dersigma22.st2 - der2pdf2.dersigma22.st2  
  
  der2pdf2.dereta2dersigma22.st <- dHs2$der2pdf2.dereta2dersigma2.st
  der2p2.dereta2dersigma22.st   <- dHs2$der2p2.dereta2dersigma2.st
  der2p2m1.dereta2dersigma22.st <- der2p2.dereta2dersigma22.st - der2pdf2.dereta2dersigma22.st

d2l.be1.be1 <- -VC$weights*( (E*(der2C11.derp1p1*derp1.dereta1^2 + derC11.derp1*der2p1.dereta1eta1 - 
(der2C01.derp1p1*derp1m1.dereta1^2 + derC01.derp1*der2p1m1.dereta1eta1) - 
(der2C10.derp1p1*derp1.dereta1^2 + derC10.derp1*der2p1.dereta1eta1) + 
(der2C00.derp1p1*derp1m1.dereta1^2 + derC00.derp1*der2p1m1.dereta1eta1)) - fE1^2)/E^2  )

d2l.be2.be2 <- -VC$weights*((E*(der2C11.derp2p2*derp2.dereta2^2 + derC11.derp2*der2p2.dereta2eta2 - 
(der2C01.derp2p2*derp2.dereta2^2 + derC01.derp2*der2p2.dereta2eta2) - 
(der2C10.derp2p2*derp2m1.dereta2^2 + derC10.derp2*der2p2m1.dereta2eta2) + 
(der2C00.derp2p2*derp2m1.dereta2^2 + derC00.derp2*der2p2m1.dereta2eta2)) - fE2^2)/E^2  )

d2l.sigma22.sigma22 <- -VC$weights*((E*(der2C11.derp2p2*derp2.dersigma22.st^2 + derC11.derp2*der2p2.dersigma22.st2 - 
(der2C01.derp2p2*derp2.dersigma22.st^2 + derC01.derp2*der2p2.dersigma22.st2) - 
(der2C10.derp2p2*derp2m1.dersigma22.st^2 + derC10.derp2*der2p2m1.dersigma22.st2) + 
(der2C00.derp2p2*derp2m1.dersigma22.st^2 + derC00.derp2*der2p2m1.dersigma22.st2)) - fE2s^2)/E^2  )

d2l.rho.rho <- -VC$weights*((E*(der2C11.derthet2*derteta.derteta.st^2 + derC11.derthet*der2teta.derteta.stteta.st - 
(der2C01.derthet2*derteta.derteta.st^2 + derC01.derthet*der2teta.derteta.stteta.st) - 
(der2C10.derthet2*derteta.derteta.st^2 + derC10.derthet*der2teta.derteta.stteta.st) + 
(der2C00.derthet2*derteta.derteta.st^2 + derC00.derthet*der2teta.derteta.stteta.st)) - fEt^2)/E^2  )                                   

d2l.be1.be2 <- -VC$weights*( (E*(der2C11.derp1p2*derp1.dereta1*derp2.dereta2  - 
der2C01.derp1p2*derp1m1.dereta1*derp2.dereta2 - 
der2C10.derp1p2*derp1.dereta1*derp2m1.dereta2 + 
der2C00.derp1p2*derp1m1.dereta1*derp2m1.dereta2) - fE1*fE2)/E^2  )

d2l.be1.sigma22 <- -VC$weights*((E*(der2C11.derp1p2*derp1.dereta1*derp2.dersigma22.st  - 
der2C01.derp1p2*derp1m1.dereta1*derp2.dersigma22.st - 
der2C10.derp1p2*derp1.dereta1*derp2m1.dersigma22.st + 
der2C00.derp1p2*derp1m1.dereta1*derp2m1.dersigma22.st) - fE1*fE2s)/E^2  )

d2l.be1.rho <- -VC$weights*((E*(der2C11.derp1t*derp1.dereta1*derteta.derteta.st  - 
der2C01.derp1t*derp1m1.dereta1*derteta.derteta.st - 
der2C10.derp1t*derp1.dereta1*derteta.derteta.st + 
der2C00.derp1t*derp1m1.dereta1*derteta.derteta.st) - fE1*fEt)/E^2  )


d2l.be2.sigma22 <- -VC$weights*((E*(der2C11.derp2p2*derp2.dereta2*derp2.dersigma22.st + derC11.derp2*der2p2.dereta2dersigma22.st - 
(der2C01.derp2p2*derp2.dereta2*derp2.dersigma22.st + derC01.derp2*der2p2.dereta2dersigma22.st) - 
(der2C10.derp2p2*derp2m1.dereta2*derp2m1.dersigma22.st + derC10.derp2*der2p2m1.dereta2dersigma22.st) + 
(der2C00.derp2p2*derp2m1.dereta2*derp2m1.dersigma22.st + derC00.derp2*der2p2m1.dereta2dersigma22.st)) - fE2*fE2s)/E^2
  )



d2l.be2.rho <- -VC$weights*((E*(der2C11.derp2t*derp2.dereta2*derteta.derteta.st  - 
der2C01.derp2t*derp2.dereta2*derteta.derteta.st - 
der2C10.derp2t*derp2m1.dereta2*derteta.derteta.st + 
der2C00.derp2t*derp2m1.dereta2*derteta.derteta.st) - fE2*fEt)/E^2  )

d2l.rho.sigma22 <- -VC$weights*((E*(der2C11.derp2t*derp2.dersigma22.st*derteta.derteta.st  - 
der2C01.derp2t*derp2.dersigma22.st*derteta.derteta.st - 
der2C10.derp2t*derp2m1.dersigma22.st*derteta.derteta.st + 
der2C00.derp2t*derp2m1.dersigma22.st*derteta.derteta.st) - fE2s*fEt)/E^2  )




    


if( is.null(VC$X3) ){



  G   <- -c(colSums( c(dl.dbe1)*VC$X1 ),
              colSums( c(dl.dbe2)*VC$X2 ),
              sum( dl.dsigma22.st ),
              sum( dl.dteta.st )            )
  
  
  
    be1.be1 <- crossprod(VC$X1*c(d2l.be1.be1),VC$X1)
    be2.be2 <- crossprod(VC$X2*c(d2l.be2.be2),VC$X2)
    be1.be2 <- crossprod(VC$X1*c(d2l.be1.be2),VC$X2)
    be1.rho <-   t(t(rowSums(t(VC$X1*c(d2l.be1.rho)))))
    be1.sigma22 <- t(t(rowSums(t(VC$X1*c(d2l.be1.sigma22))))) 
    be2.rho <-   t(t(rowSums(t(VC$X2*c(d2l.be2.rho)))))
    be2.sigma22 <- t(t(rowSums(t(VC$X2*c(d2l.be2.sigma22)))))
  
    H <- rbind( cbind( be1.be1    ,   be1.be2    ,       be1.sigma22,       be1.rho  ), 
                cbind( t(be1.be2) ,   be2.be2    ,       be2.sigma22,        be2.rho  ), 
                cbind( t(be1.sigma22) , t(be2.sigma22) , sum(d2l.sigma22.sigma22), sum(d2l.rho.sigma22)),
                cbind( t(be1.rho) ,   t(be2.rho) , sum(d2l.rho.sigma22), sum(d2l.rho.rho) ) 
                
              ) 

}



if( !is.null(VC$X3) ){



G   <- -c( colSums(       c(dl.dbe1)*VC$X1 ) ,
           colSums(       c(dl.dbe2)*VC$X2 ) ,
           colSums(c(dl.dsigma22.st)*VC$X3 ) ,
           colSums(   c(dl.dteta.st)*VC$X4 )  )  

                 
    be1.be1         <- crossprod(VC$X1*c(d2l.be1.be1),VC$X1)
    be2.be2         <- crossprod(VC$X2*c(d2l.be2.be2),VC$X2)
    be1.be2         <- crossprod(VC$X1*c(d2l.be1.be2),VC$X2)
    be1.rho         <- crossprod(VC$X1*c(d2l.be1.rho),VC$X4)    
    be2.rho         <- crossprod(VC$X2*c(d2l.be2.rho),VC$X4)   
    be1.sigma22     <- crossprod(VC$X1*c(d2l.be1.sigma22),VC$X3)      
    be2.sigma22     <- crossprod(VC$X2*c(d2l.be2.sigma22),VC$X3)   
    sigma22.sigma22 <- crossprod(VC$X3*c(d2l.sigma22.sigma22),VC$X3)
    rho.sigma22     <- crossprod(VC$X3*c(d2l.rho.sigma22),VC$X4)   
    rho.rho         <- crossprod(VC$X4*c(d2l.rho.rho),VC$X4)    
    
    
    H <- rbind( cbind( be1.be1        ,   be1.be2      ,       be1.sigma22,       be1.rho    ), 
                cbind( t(be1.be2)     ,   be2.be2      ,       be2.sigma22,       be2.rho    ), 
                cbind( t(be1.sigma22) , t(be2.sigma22) ,      sigma22.sigma22,   rho.sigma22),
                cbind( t(be1.rho)     ,   t(be2.rho)   ,   t(rho.sigma22),    rho.rho    ) ) 
                
                 
                 

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
  
  
  

  
  
  
  

         list(value=res, gradient=G, hessian=H, S.h=S.h, S.h1=S.h1, S.h2=S.h2, l=S.res, l.ln = l.ln, l.par=l.par, ps = ps, 
              eta1=eta1, eta2=eta2, etad=etad, etas1 = etas1, etas2 = etas2, 
              BivD=VC$BivD, p1 = p1, p2 = p2, pdf1 = pdf1, pdf2 = pdf2,          
              c.copula.be2 = c.copula.be2,
              c.copula.be1 = c.copula.be1,
              c.copula2.be1be2 = c.copula2.be1be2, 
              dl.dbe1          =dl.dbe1,       
              dl.dbe2          =dl.dbe2,       
              dl.dsigma22.st   =dl.dsigma22.st,
              dl.dteta.st      =dl.dteta.st,
                            teta.ind2 = teta.ind2, teta.ind1 = teta.ind1,
              Cop1 = Cop1, Cop2 = Cop2, teta1 = teta1, teta2 = teta2) 
              
              

  }
  

  
