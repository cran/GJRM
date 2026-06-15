bdiscrcontThetaOnly <- function(params, respvec, VC, ps, AT = FALSE){


p1 <- p2 <- pdf1 <- pdf2 <- c.copula.be2 <- c.copula.be1 <- c.copula2.be1be2 <- NA

Cop1 <- Cop2 <- nC1 <- nC2 <- teta.ind1 <- teta.ind2 <- teta1 <- teta2 <- NA


  if( is.null(VC$X3)) teta.st <- etad <- params
  if(!is.null(VC$X3)) teta.st <- etad <- VC$X5%*%params
  
##################
## Transformations
##################

    
resT    <- teta.tr(VC, teta.st)

teta.st1 <- teta.st2 <- teta.st <- resT$teta.st
teta1 <- teta2 <- teta <- resT$teta 
    
##################

Cop1 <- Cop2 <- VC$BivD 
nC1 <- nC2 <- VC$nC 

teta.ind1 <- as.logical(c(1,0,round(runif(VC$n-2))) ) # just random split here
teta.ind2 <- teta.ind1 == FALSE 


if(!(VC$BivD %in% VC$BivD2) && length(teta.st) > 1){ # not the best choice but implemented like this many years ago...for clarity

teta.st1 <- teta.st[teta.ind1]
teta.st2 <- teta.st[teta.ind2]

teta1 <- teta[teta.ind1]
teta2 <- teta[teta.ind2]

}


 
if(VC$BivD %in% VC$BivD2 && length(teta.st) > 1){

Cop1Cop2R <- Cop1Cop2(VC$BivD)
Cop1 <- Cop1Cop2R$Cop1
Cop2 <- Cop1Cop2R$Cop2
nC1  <- VC$ct[which(VC$ct[,1] == Cop1), 2] 
nC2  <- VC$ct[which(VC$ct[,1] == Cop2), 2]

teta.ind1 <- as.logical( etad >= 0 ) # for + coefs
teta.ind2 <- teta.ind1 == FALSE      # for - coefs


if(VC$BivD %in% VC$BivD2[c(1:4,13:16)]){
  teta1 <- teta.st1 <- etad[teta.ind1]
  teta2 <- teta.st2 <- etad[teta.ind2]
                                       }

if(VC$BivD %in% VC$BivD2[5:12]){
  teta1 <- teta.st1 <- sign(etad[teta.ind1]) * (abs(etad[teta.ind1]) + 1)
  teta2 <- teta.st2 <- sign(etad[teta.ind2]) * (abs(etad[teta.ind2]) + 1)
                                       }

} 

if(VC$BivD %in% VC$BivD2 && length(teta.st) == 1){


Cop1Cop2R <- Cop1Cop2(VC$BivD)
Cop1 <- Cop1Cop2R$Cop1
Cop2 <- Cop1Cop2R$Cop2
nC1  <- VC$ct[which(VC$ct[,1] == Cop1), 2] 
nC2  <- VC$ct[which(VC$ct[,1] == Cop2), 2]

teta.i1 <- as.logical( etad >= 0 ) # this becomes an index for + and - coefs

if(teta.i1 == TRUE){  Cop2 <- Cop1; nC2 <- nC1}
if(teta.i1 == FALSE){ Cop1 <- Cop2; nC1 <- nC2}   

if(VC$BivD %in% VC$BivD2[c(1:4,13:16)]) teta1 <- teta.st1 <- teta2 <- teta.st2 <- etad                                     
if(VC$BivD %in% VC$BivD2[5:12])         teta1 <- teta.st1 <- teta2 <- teta.st2 <- sign(etad) * (abs(etad) + 1)


teta.ind1 <- as.logical(c(1,0,round(runif(VC$n-2))) ) # again not the best choice but to make things run
teta.ind2 <- teta.ind1 == FALSE 

} 
 

#######################################################

  p1 <- respvec$y1
  p2 <- respvec$y2
  
  pdf1 <- VC$pdf1 
  pdf2 <- VC$pdf2
                               
  if( length(teta1) != 0) dH11 <- copgHsAT(p1[teta.ind1], p2[teta.ind1], teta1, Cop1, Ln = TRUE, par2 = VC$dof, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr, teta.st = teta.st1, T.BivD = VC$BivD)
  if( length(teta2) != 0) dH12 <- copgHsAT(p1[teta.ind2], p2[teta.ind2], teta2, Cop2, Ln = TRUE, par2 = VC$dof, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr, teta.st = teta.st2, T.BivD = VC$BivD)  
  h1  <- NA    
  if( length(teta1) != 0) h1[teta.ind1] <- dH11$c.copula.be2
  if( length(teta2) != 0) h1[teta.ind2] <- dH12$c.copula.be2      
    
  if( length(teta1) != 0) dH21 <- copgHsAT(mm(p1[teta.ind1]-pdf1[teta.ind1], min.pr = VC$min.pr, max.pr = VC$max.pr), p2[teta.ind1], teta1, Cop1, Ln = TRUE, par2 = VC$dof, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr, teta.st = teta.st1, T.BivD = VC$BivD)
  if( length(teta2) != 0) dH22 <- copgHsAT(mm(p1[teta.ind2]-pdf1[teta.ind2], min.pr = VC$min.pr, max.pr = VC$max.pr), p2[teta.ind2], teta2, Cop2, Ln = TRUE, par2 = VC$dof, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr, teta.st = teta.st2, T.BivD = VC$BivD)
  
  h2  <- NA    
  if( length(teta1) != 0) h2[teta.ind1] <- dH21$c.copula.be2
  if( length(teta2) != 0) h2[teta.ind2] <- dH22$c.copula.be2    

  c.copula2.be2thH1 <- NA
  if( length(teta1) != 0) c.copula2.be2thH1[teta.ind1] <- dH11$c.copula2.be2th
  if( length(teta2) != 0) c.copula2.be2thH1[teta.ind2] <- dH12$c.copula2.be2th
 
  c.copula2.be2thH2 <- NA
  if( length(teta1) != 0) c.copula2.be2thH2[teta.ind1] <- dH21$c.copula2.be2th
  if( length(teta2) != 0) c.copula2.be2thH2[teta.ind2] <- dH22$c.copula2.be2th  
 
   
  if( length(teta1) != 0) BITS.H11 <- copgHsContFM(p1[teta.ind1], p2[teta.ind1], teta1, teta.st1, Cop1, par2 = VC$dof, nu.st = log(VC$dof-2), T.BivD = VC$BivD)
  if( length(teta2) != 0) BITS.H12 <- copgHsContFM(p1[teta.ind2], p2[teta.ind2], teta2, teta.st2, Cop2, par2 = VC$dof, nu.st = log(VC$dof-2), T.BivD = VC$BivD)  
  
  if( length(teta1) != 0) BITS.H21 <- copgHsContFM(mm(p1[teta.ind1]-pdf1[teta.ind1], min.pr = VC$min.pr, max.pr = VC$max.pr), p2[teta.ind1], teta1, teta.st1, Cop1, par2 = VC$dof, nu.st = log(VC$dof-2), T.BivD = VC$BivD) 
  if( length(teta2) != 0) BITS.H22 <- copgHsContFM(mm(p1[teta.ind2]-pdf1[teta.ind2], min.pr = VC$min.pr, max.pr = VC$max.pr), p2[teta.ind2], teta2, teta.st2, Cop2, par2 = VC$dof, nu.st = log(VC$dof-2), T.BivD = VC$BivD) 
 
 
  der2h.derteta.teta.st.H1 <- NA
  if( length(teta1) != 0) der2h.derteta.teta.st.H1[teta.ind1] <- BITS.H11$der2h.derteta.teta.st
  if( length(teta2) != 0) der2h.derteta.teta.st.H1[teta.ind2] <- BITS.H12$der2h.derteta.teta.st

  derteta.derteta.st <- NA
  if( length(teta1) != 0) derteta.derteta.st[teta.ind1] <- BITS.H11$derteta.derteta.st  
  if( length(teta2) != 0) derteta.derteta.st[teta.ind2] <- BITS.H12$derteta.derteta.st 


  der2teta.derteta.stteta.st <- NA
  if( length(teta1) != 0) der2teta.derteta.stteta.st[teta.ind1] <- BITS.H11$der2teta.derteta.stteta.st  
  if( length(teta2) != 0) der2teta.derteta.stteta.st[teta.ind2] <- BITS.H12$der2teta.derteta.stteta.st 

  der2h.derteta.teta.st.H2 <- NA
  if( length(teta1) != 0) der2h.derteta.teta.st.H2[teta.ind1] <- BITS.H21$der2h.derteta.teta.st
  if( length(teta2) != 0) der2h.derteta.teta.st.H2[teta.ind2] <- BITS.H22$der2h.derteta.teta.st
   
#######################################################   
 

  diffh1.h2 <- h1 - h2
  diffh1.h2 <- mm(diffh1.h2, min.pr = VC$min.pr, max.pr = VC$max.pr)  

  l.par <- VC$weights*( log(diffh1.h2) )
  res   <- -sum(l.par)

  dl.dteta.st    <- VC$weights*( 1/diffh1.h2*(c.copula2.be2thH1 - c.copula2.be2thH2)  )
  
  d2l.rho.rho <- -VC$weights*( - (c.copula2.be2thH1 - c.copula2.be2thH2)^2/diffh1.h2^2 + diffh1.h2^-1*(der2h.derteta.teta.st.H1*derteta.derteta.st^2 - (c.copula2.be2thH1/derteta.derteta.st)*der2teta.derteta.stteta.st - (der2h.derteta.teta.st.H2*derteta.derteta.st^2 - (c.copula2.be2thH2/derteta.derteta.st)*der2teta.derteta.stteta.st)     ) )                                   
        
#######################################################   

if( is.null(VC$X3) ){

  G <- -sum( dl.dteta.st )
  H <- as.matrix(sum(d2l.rho.rho) ) 
                           
                    }

if( !is.null(VC$X3) ){

  G <- -c( colSums( c(dl.dteta.st)*VC$X5 ) )  
  H <- crossprod(VC$X5*c(d2l.rho.rho), VC$X5) 
                                        
                      }


 
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
  


 
######### ******************
C11 <- C01 <- NA

if (length(teta1) != 0){

C11[teta.ind1] <- BiCDF(p1[teta.ind1],                                                             p2[teta.ind1], family = nC1, par1 = teta1, par2 = VC$dof, test = TRUE)
C01[teta.ind1] <- BiCDF(mm(p1[teta.ind1]-pdf1[teta.ind1], min.pr = VC$min.pr, max.pr = VC$max.pr), p2[teta.ind1], family = nC1, par1 = teta1, par2 = VC$dof, test = TRUE)

                        }

if (length(teta2) != 0){

C11[teta.ind2] <- BiCDF(p1[teta.ind2],                                                             p2[teta.ind2], family = nC2, par1 = teta2, par2 = VC$dof, test = TRUE)
C01[teta.ind2] <- BiCDF(mm(p1[teta.ind2]-pdf1[teta.ind2], min.pr = VC$min.pr, max.pr = VC$max.pr), p2[teta.ind2], family = nC2, par1 = teta2, par2 = VC$dof, test = TRUE)

                       }

# for joint residuals
######### ****************** 
 




         list(value=res, gradient=G, hessian=H, S.h=S.h, S.h1=S.h1, S.h2=S.h2, l=S.res, l.par=l.par, ps = ps, 
              etad=etad,
              BivD=VC$BivD, p1 = p1, p2 = p2, pdf1 = pdf1, pdf2 = pdf2,          
              c.copula.be2 = c.copula.be2,
              c.copula.be1 = c.copula.be1,
              c.copula2.be1be2 = c.copula2.be1be2, 
              dl.dteta.st      =dl.dteta.st,
              teta.ind2 = teta.ind2, teta.ind1 = teta.ind1,
              Cop1 = Cop1, Cop2 = Cop2, nC1 = nC1, nC2 = nC2, teta1 = teta1, teta2 = teta2, C11 = C11, C01 = C01) 
              
              

  }
  
