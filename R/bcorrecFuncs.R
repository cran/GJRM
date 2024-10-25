




pPOtr <- function(y, mu, left.trunc){ (pPO(y, mu) - pPO(rep(left.trunc, length(mu)), mu = mu))/(1 - pPO(rep(left.trunc, length(mu)), mu)) } # PO does not recycle hence the use of rep, but ppois does it

dPOtr <- function(y, mu, left.trunc){ dPO(y, mu)/(1 - pPO(rep(left.trunc, length(mu)), mu)) }

rPtr  <- function(n, mu, left.trunc) qpois(runif(n, dpois(left.trunc, mu), 1), mu) # recycling rule ok here 




pNBItr <- function(y, mu, sigma, left.trunc){ (pNBI(y, mu, sigma) - pNBI(rep(left.trunc, length(mu)), mu, sigma))/(1 - pNBI(rep(left.trunc, length(mu)), mu, sigma)) } 

dNBItr <- function(y, mu, sigma, left.trunc){ dNBI(y, mu, sigma)/(1 - pNBI(rep(left.trunc, length(mu)), mu, sigma)) }

rNBItr <- function(n, mu, sigma, left.trunc){
               
              n  <- ceiling(n)
              p  <- runif(n)
              pp <- pNBI(rep(left.trunc, length(mu)), mu, sigma)
              p  <- (pp + p * (1 - pp))

              qnbinom(p, size = 1/sigma, mu = mu, lower.tail = TRUE, log.p = FALSE)

}



pNBIItr <- function(y, mu, sigma, left.trunc){ (pNBII(y, mu, sigma) - pNBII(rep(left.trunc, length(mu)), mu, sigma))/(1 - pNBII(rep(left.trunc, length(mu)), mu, sigma)) } 

dNBIItr <- function(y, mu, sigma, left.trunc){ dNBII(y, mu, sigma)/(1 - pNBII(rep(left.trunc, length(mu)), mu, sigma)) }

rNBIItr <- function(n, mu, sigma, left.trunc){
               
              n  <- ceiling(n)
              p  <- runif(n)
              pp <- pNBII(rep(left.trunc, length(mu)), mu, sigma)
              p  <- (pp + p * (1 - pp))

              qnbinom(p, size = mu/sigma, mu = mu, lower.tail = TRUE, log.p = FALSE)

}




pPIGtr <- function(y, mu, sigma, left.trunc){ (pPIG(y, mu, sigma) - pPIG(rep(left.trunc, length(mu)), mu, sigma))/(1 - pPIG(rep(left.trunc, length(mu)), mu, sigma)) } 

dPIGtr <- function(y, mu, sigma, left.trunc){ dPIG(y, mu, sigma)/(1 - pPIG(rep(left.trunc, length(mu)), mu, sigma)) }

rPIGtr <- function(n, mu, sigma, left.trunc){
               
              max.value <- 10000 
              n  <- ceiling(n)
              p  <- runif(n)
              pp <- pPIG(rep(left.trunc, length(mu)), mu, sigma)
              p  <- (pp + p * (1 - pp))

              
              ly <- length(p)
              QQQ <- rep(0, ly)
              nsigma <- rep(sigma, length = ly)
              nmu <- rep(mu, length = ly)
              for (i in seq(along = p)) {
              cumpro <- 0
              if (p[i] + 1e-09 >= 1) QQQ[i] <- Inf else{
                                 for(j in seq(from = 0, to = max.value)){
                                           cumpro <- pPIG(j, mu = nmu[i], sigma = nsigma[i], log.p = FALSE)
                                           QQQ[i] <- j
                                           if (p[i] <= cumpro) break
                                                                        }
                                                        }
                                         }
              QQQ
              
}












ordcont.m <- function(y, pk, pk1, eta.mu, sigma, nu, theta, margin, BivD, par2, min.dn, min.pr, max.pr, y1){           
     

p12.s <- 0

pd2  <- distrHsAT(y, eta2 = eta.mu, sigma2 = sigma, nu = NULL, margin2 = margin, min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)
p2   <- pd2$p2
pdf2 <- pd2$pdf2


if(y1 > 1){
            p12.f <- copgHs2(p1 = pk,  p2 = p2, teta = theta, teta.st = theta, BivD = BivD, par2 = par2, min.pr = min.pr, max.pr = max.pr, min.dn = min.dn)$c.copula.be2
            p12.s <- copgHs2(p1 = pk1, p2 = p2, teta = theta, teta.st = theta, BivD = BivD, par2 = par2, min.pr = min.pr, max.pr = max.pr, min.dn = min.dn)$c.copula.be2
            
          }

if(y1 == 1) p12.f <- copgHs2(p1 = pk,  p2 = p2, teta = theta, teta.st = theta, BivD = BivD, par2 = par2, min.pr = min.pr, max.pr = max.pr, min.dn = min.dn)$c.copula.be2


y*(p12.f - p12.s)*pdf2/(pk - pk1) 


 
}


ordcont.v <- function(y, pk, pk1, eta.mu, sigma, nu, theta, margin, BivD, par2, min.dn, min.pr, max.pr, y1){           
     

p12.s <- 0

pd2  <- distrHsAT(y, eta2 = eta.mu, sigma2 = sigma, nu = NULL, margin2 = margin, min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)
p2   <- pd2$p2
pdf2 <- pd2$pdf2


if(y1 > 1){
            p12.f <- copgHs2(p1 = pk,  p2 = p2, teta = theta, teta.st = theta, BivD = BivD, par2 = par2, min.pr = min.pr, max.pr = max.pr, min.dn = min.dn)$c.copula.be2
            p12.s <- copgHs2(p1 = pk1, p2 = p2, teta = theta, teta.st = theta, BivD = BivD, par2 = par2, min.pr = min.pr, max.pr = max.pr, min.dn = min.dn)$c.copula.be2
            
          }

if(y1 == 1) p12.f <- copgHs2(p1 = pk,  p2 = p2, teta = theta, teta.st = theta, BivD = BivD, par2 = par2, min.pr = min.pr, max.pr = max.pr, min.dn = min.dn)$c.copula.be2


y^2*(p12.f - p12.s)*pdf2/(pk - pk1) 


 
}







dicont.m <- function(y, y12, eta1, eta2, sig1, sig2, nu1 = NULL, nu2 = NULL, theta, margins, BivD, par2, min.dn, min.pr, max.pr, nC, left.trunc1){           
 
y1m <- y2m <- NA

p1pdf1 <- distrHsATDiscr(y,   eta2 = eta1, sigma2 = sig1, nu = 1, margin2 = margins[1], y2m = y1m, robust = FALSE, min.dn = min.pr, 
                         min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc1)
p1     <- as.numeric(p1pdf1$p2)
pmf1   <- as.numeric(p1pdf1$pdf2)


rpd2 <- distrHsAT(y12,   eta2 = eta2, sigma2 = sig2, nu = nu2, margin2 = margins[2], min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)     
p2   <- rpd2$p2
pdf2 <- rpd2$pdf2

                                             
  c.be2.h1 <- copgHs2(p1 = p1, p2 = p2, teta = theta, teta.st = theta, BivD = BivD, par2 = par2, min.pr = min.pr, max.pr = max.pr, min.dn = min.dn)$c.copula.be2
  c.be2.h2 <- copgHs2(p1 = mm(p1 - pmf1, min.pr = min.pr, max.pr = max.pr), p2 = p2, teta = theta, teta.st = theta, BivD = BivD, par2 = par2, min.pr = min.pr, max.pr = max.pr, min.dn = min.dn)$c.copula.be2
  

  
  diffh1.h2 <- c.be2.h1 - c.be2.h2  # diffh1.h2 <- mm(diffh1.h2, min.pr = min.pr, max.pr = max.pr)  

  
  y*diffh1.h2


}






dicont.v <- function(y, y12, eta1, eta2, sig1, sig2, nu1 = NULL, nu2 = NULL, theta, margins, BivD, par2, min.dn, min.pr, max.pr, nC, left.trunc1){           
 
y1m <- y2m <- NA

p1pdf1 <- distrHsATDiscr(y,   eta2 = eta1, sigma2 = sig1, nu = 1, margin2 = margins[1], y2m = y1m, robust = FALSE, min.dn = min.pr, 
                         min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc1)
p1     <- as.numeric(p1pdf1$p2)
pmf1   <- as.numeric(p1pdf1$pdf2)


rpd2 <- distrHsAT(y12,   eta2 = eta2, sigma2 = sig2, nu = nu2, margin2 = margins[2], min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)     
p2   <- rpd2$p2
pdf2 <- rpd2$pdf2

                                             
  c.be2.h1 <- copgHs2(p1 = p1, p2 = p2, teta = theta, teta.st = theta, BivD = BivD, par2 = par2, min.pr = min.pr, max.pr = max.pr, min.dn = min.dn)$c.copula.be2
  c.be2.h2 <- copgHs2(p1 = mm(p1 - pmf1, min.pr = min.pr, max.pr = max.pr), p2 = p2, teta = theta, teta.st = theta, BivD = BivD, par2 = par2, min.pr = min.pr, max.pr = max.pr, min.dn = min.dn)$c.copula.be2
  

  
  diffh1.h2 <- c.be2.h1 - c.be2.h2  # diffh1.h2 <- mm(diffh1.h2, min.pr = min.pr, max.pr = max.pr)  

  
  y^2*diffh1.h2


}









discrcont.m <- function(y, y12, eta1, eta2, sig1, sig2, nu1 = NULL, nu2 = NULL, theta, margins, BivD, par2, min.dn, min.pr, max.pr, left.trunc1){           
 
y1m <- y2m <- NA
 
   
  p1pdf1 <- distrHsATDiscr(y12, eta2 = eta1, sigma2 = sig1, nu = 1, margin2 = margins[1], y2m = y1m, robust = FALSE, 
                           min.dn = min.pr, min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc1)
  p1     <- as.numeric(p1pdf1$p2)
  pmf1   <- as.numeric(p1pdf1$pdf2)   
   

  rpd2 <- distrHsAT(y,   eta2 = eta2, sigma2 = sig2, nu = nu2, margin2 = margins[2], min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)     
  p2   <- rpd2$p2
  pdf2 <- rpd2$pdf2

  
  
  c.be2.h1 <- copgHs2(p1 = p1, p2 = p2, teta = theta, teta.st = theta, BivD = BivD, par2 = par2, min.pr = min.pr, max.pr = max.pr, min.dn = min.dn)$c.copula.be2
  c.be2.h2 <- copgHs2(p1 = mm(p1 - pmf1, min.pr = min.pr, max.pr = max.pr), p2 = p2, teta = theta, teta.st = theta, BivD = BivD, par2 = par2, min.pr = min.pr, max.pr = max.pr, min.dn = min.dn)$c.copula.be2
  

  
  diffh1.h2 <- c.be2.h1 - c.be2.h2  # diffh1.h2 <- mm(diffh1.h2, min.pr = min.pr, max.pr = max.pr)  

  
  y*pdf2*diffh1.h2/pmf1

}


discrcont.v <- function(y, y12, eta1, eta2, sig1, sig2, nu1 = NULL, nu2 = NULL, theta, margins, BivD, par2, min.dn, min.pr, max.pr, left.trunc1){           
 
y1m <- y2m <- NA
 
   
  p1pdf1 <- distrHsATDiscr(y12, eta2 = eta1, sigma2 = sig1, nu = 1, margin2 = margins[1], y2m = y1m, robust = FALSE, 
                           min.dn = min.pr, min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc1)
  p1     <- as.numeric(p1pdf1$p2)
  pmf1   <- as.numeric(p1pdf1$pdf2)   
   

  rpd2 <- distrHsAT(y,   eta2 = eta2, sigma2 = sig2, nu = nu2, margin2 = margins[2], min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)     
  p2   <- rpd2$p2
  pdf2 <- rpd2$pdf2

  
  
  c.be2.h1 <- copgHs2(p1 = p1, p2 = p2, teta = theta, teta.st = theta, BivD = BivD, par2 = par2, min.pr = min.pr, max.pr = max.pr, min.dn = min.dn)$c.copula.be2
  c.be2.h2 <- copgHs2(p1 = mm(p1 - pmf1, min.pr = min.pr, max.pr = max.pr), p2 = p2, teta = theta, teta.st = theta, BivD = BivD, par2 = par2, min.pr = min.pr, max.pr = max.pr, min.dn = min.dn)$c.copula.be2
  

  
  diffh1.h2 <- c.be2.h1 - c.be2.h2  # diffh1.h2 <- mm(diffh1.h2, min.pr = min.pr, max.pr = max.pr)  

  
  y^2*pdf2*diffh1.h2/pmf1

}







bindisc.m <- function(y1, y12, eta2, sig2, nu2 = NULL, theta, margins, BivD, par2, min.dn, min.pr, max.pr, nC, p0, left.trunc2){           
 
y1m <- y2m <- NA

    

p2pdf2 <- distrHsATDiscr(y12, eta2 = eta2, sigma2 = sig2, nu = 1, margin2 = margins[2], y2m = y2m, robust = FALSE, 
                         min.dn = min.pr, min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc2)
p2     <- as.numeric(p2pdf2$p2)
pmf2   <- as.numeric(p2pdf2$pdf2)



  C1 <- mm(BiCDF(p0, p2,                                              nC, theta, par2), min.pr = min.pr, max.pr = max.pr  )
  C2 <- mm(BiCDF(p0, mm(p2 - pmf2, min.pr = min.pr, max.pr = max.pr), nC, theta, par2), min.pr = min.pr, max.pr = max.pr  )
  
  A <- mm(C1 - C2,  min.pr = min.pr, max.pr = max.pr)     # for y1 = 0
  B <- mm(pmf2 - A, min.pr = min.pr, max.pr = max.pr)     # for y1 = 1
                                             

 if(y1 == 0) y12*A/p0 else y12*B/(1-p0) 
 


}


bindisc.v <- function(y1, y12, eta2, sig2, nu2 = NULL, theta, margins, BivD, par2, min.dn, min.pr, max.pr, nC, p0, left.trunc2){           
 
y1m <- y2m <- NA

    

p2pdf2 <- distrHsATDiscr(y12, eta2 = eta2, sigma2 = sig2, nu = 1, margin2 = margins[2], y2m = y2m, robust = FALSE, 
                         min.dn = min.pr, min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc2)
p2     <- as.numeric(p2pdf2$p2)
pmf2   <- as.numeric(p2pdf2$pdf2)



  C1 <- mm(BiCDF(p0, p2,                                              nC, theta, par2), min.pr = min.pr, max.pr = max.pr  )
  C2 <- mm(BiCDF(p0, mm(p2 - pmf2, min.pr = min.pr, max.pr = max.pr), nC, theta, par2), min.pr = min.pr, max.pr = max.pr  )
  
  A <- mm(C1 - C2,  min.pr = min.pr, max.pr = max.pr)     # for y1 = 0
  B <- mm(pmf2 - A, min.pr = min.pr, max.pr = max.pr)     # for y1 = 1
                                             

 if(y1 == 0) y12^2*A/p0 else y12^2*B/(1-p0) 
 


}









discdisc.m <- function(y, y12, eta1, eta2, sig1, sig2, nu1 = NULL, nu2 = NULL, theta, margins, BivD, par2, min.dn, min.pr, max.pr, cond, nC, left.trunc1, left.trunc2){           
 
y1m <- y2m <- NA

if(cond == 1){ 

p1pdf1 <- distrHsATDiscr(y12, eta2 = eta1, sigma2 = sig1, nu = 1, margin2 = margins[1], y2m = y1m, robust = FALSE, 
                         min.dn = min.pr, min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc1)
p1     <- as.numeric(p1pdf1$p2)
pmf1   <- as.numeric(p1pdf1$pdf2)


p2pdf2 <- distrHsATDiscr(y,   eta2 = eta2, sigma2 = sig2, nu = 1, margin2 = margins[2], y2m = y2m, robust = FALSE, 
                         min.dn = min.pr, min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc2)
p2     <- as.numeric(p2pdf2$p2)
pmf2   <- as.numeric(p2pdf2$pdf2)

                                             
              } 
 
 
if(cond == 2){ 

p1pdf1 <- distrHsATDiscr(y,   eta2 = eta1, sigma2 = sig1, nu = 1, margin2 = margins[1], y2m = y1m, robust = FALSE, 
                         min.dn = min.pr, min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc1)
p1     <- as.numeric(p1pdf1$p2)
pmf1   <- as.numeric(p1pdf1$pdf2)


p2pdf2 <- distrHsATDiscr(y12, eta2 = eta2, sigma2 = sig2, nu = 1, margin2 = margins[2], y2m = y2m, robust = FALSE, 
                         min.dn = min.pr, min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc2)
p2     <- as.numeric(p2pdf2$p2)
pmf2   <- as.numeric(p2pdf2$pdf2)

                                             
              } 
 

  C11 <- mm(BiCDF(   p1,                         p2,                            nC, theta, par2), min.pr = min.pr, max.pr = max.pr  )
  C01 <- mm(BiCDF(mm(p1 - pmf1, min.pr, max.pr), p2,                            nC, theta, par2), min.pr = min.pr, max.pr = max.pr  )
  C10 <- mm(BiCDF(   p1,                         mm(p2 - pmf2, min.pr, max.pr), nC, theta, par2), min.pr = min.pr, max.pr = max.pr  )
  C00 <- mm(BiCDF(mm(p1 - pmf1, min.pr, max.pr), mm(p2 - pmf2, min.pr, max.pr), nC, theta, par2), min.pr = min.pr, max.pr = max.pr  )

  p12 <- C11 - C01 - C10 + C00
  p12 <- ifelse(p12 < min.pr, min.pr, p12)



 if(cond == 1) y*p12/pmf1 else y*p12/pmf2 
 
 #if(cond == 1) y*pmf2 else y*pmf1 # test 


}



discdisc.v <- function(y, y12, eta1, eta2, sig1, sig2, nu1 = NULL, nu2 = NULL, theta, margins, BivD, par2, min.dn, min.pr, max.pr, cond, nC, left.trunc1, left.trunc2){           
 
y1m <- y2m <- NA

  
if(cond == 1){ 

p1pdf1 <- distrHsATDiscr(y12, eta2 = eta1, sigma2 = sig1, nu = 1, margin2 = margins[1], y2m = y1m, robust = FALSE, 
                         min.dn = min.pr, min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc1)
p1     <- as.numeric(p1pdf1$p2)
pmf1   <- as.numeric(p1pdf1$pdf2)


p2pdf2 <- distrHsATDiscr(y,   eta2 = eta2, sigma2 = sig2, nu = 1, margin2 = margins[2], y2m = y2m, robust = FALSE, min.dn = min.pr, 
                         min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc2)
p2     <- as.numeric(p2pdf2$p2)
pmf2   <- as.numeric(p2pdf2$pdf2)

                                             
              } 
 
 
if(cond == 2){ 

p1pdf1 <- distrHsATDiscr(y,   eta2 = eta1, sigma2 = sig1, nu = 1, margin2 = margins[1], y2m = y1m, robust = FALSE, 
                         min.dn = min.pr, min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc1)
p1     <- as.numeric(p1pdf1$p2)
pmf1   <- as.numeric(p1pdf1$pdf2)


p2pdf2 <- distrHsATDiscr(y12, eta2 = eta2, sigma2 = sig2, nu = 1, margin2 = margins[2], y2m = y2m, robust = FALSE, 
                         min.dn = min.pr, min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc2)
p2     <- as.numeric(p2pdf2$p2)
pmf2   <- as.numeric(p2pdf2$pdf2)

                                             
              } 
 

  C11 <- mm(BiCDF(   p1,                         p2,                            nC, theta, par2), min.pr = min.pr, max.pr = max.pr  )
  C01 <- mm(BiCDF(mm(p1 - pmf1, min.pr, max.pr), p2,                            nC, theta, par2), min.pr = min.pr, max.pr = max.pr  )
  C10 <- mm(BiCDF(   p1,                         mm(p2 - pmf2, min.pr, max.pr), nC, theta, par2), min.pr = min.pr, max.pr = max.pr  )
  C00 <- mm(BiCDF(mm(p1 - pmf1, min.pr, max.pr), mm(p2 - pmf2, min.pr, max.pr), nC, theta, par2), min.pr = min.pr, max.pr = max.pr  )

  p12 <- C11 - C01 - C10 + C00
  p12 <- ifelse(p12 < min.pr, min.pr, p12)



 if(cond == 1) y^2*p12/pmf1 else y^2*p12/pmf2 
 


}
















contcont.m <- function(y, y12, eta1, eta2, sig1, sig2, nu1 = NULL, nu2 = NULL, theta, margins, BivD, par2, min.dn, min.pr, max.pr, cond){           
     
     
 if(cond == 1){    
     
 rpd1 <- distrHsAT(y12, eta2 = eta1, sigma2 = sig1, nu = nu1, margin2 = margins[1], min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)
 rpd2 <- distrHsAT(y,   eta2 = eta2, sigma2 = sig2, nu = nu2, margin2 = margins[2], min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)     

              }
              
 if(cond == 2){    
     
 rpd1 <- distrHsAT(y,   eta2 = eta1, sigma2 = sig1, nu = nu1, margin2 = margins[1], min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)
 rpd2 <- distrHsAT(y12, eta2 = eta2, sigma2 = sig2, nu = nu2, margin2 = margins[2], min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)     
 
              }             

 pdf1 <- rpd1$pdf2
 p1   <- rpd1$p2
 
 pdf2 <- rpd2$pdf2
 p2   <- rpd2$p2 
 
 c.be1be2 <- copgHs2(p1 = p1, p2 = p2, teta = theta, teta.st = theta, BivD = BivD, par2 = par2, min.pr = min.pr, max.pr = max.pr, min.dn = min.dn)$c.copula2.be1be2

 if(cond == 1) y*c.be1be2*pdf2 else y*c.be1be2*pdf1 

}


contcont.v <- function(y, y12, eta1, eta2, sig1, sig2, nu1 = NULL, nu2 = NULL, theta, margins, BivD, par2, min.dn, min.pr, max.pr, cond){           
     
     
 if(cond == 1){    
     
 rpd1 <- distrHsAT(y12, eta2 = eta1, sigma2 = sig1, nu = nu1, margin2 = margins[1], min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)
 rpd2 <- distrHsAT(y,   eta2 = eta2, sigma2 = sig2, nu = nu2, margin2 = margins[2], min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)     

              }
              
 if(cond == 2){    
     
 rpd1 <- distrHsAT(y,   eta2 = eta1, sigma2 = sig1, nu = nu1, margin2 = margins[1], min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)
 rpd2 <- distrHsAT(y12, eta2 = eta2, sigma2 = sig2, nu = nu2, margin2 = margins[2], min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)     
 
              }             

 pdf1 <- rpd1$pdf2
 p1   <- rpd1$p2
 
 pdf2 <- rpd2$pdf2
 p2   <- rpd2$p2 
 
 c.be1be2 <- copgHs2(p1 = p1, p2 = p2, teta = theta, teta.st = theta, BivD = BivD, par2 = par2, min.pr = min.pr, max.pr = max.pr, min.dn = min.dn)$c.copula2.be1be2

 if(cond == 1) y^2*c.be1be2*pdf2 else y^2*c.be1be2*pdf1 

}









bincont.m <- function(y, p0, eta.mu, sigma, nu, theta, margin, BivD, par2, min.dn, min.pr, max.pr, y1){           
     
 rpd <- distrHsAT(y, eta2 = eta.mu, sigma2 = sigma, nu = nu, margin2 = margin, min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)     
     
 p2 <- rpd$p2; pdf2 <- rpd$pdf2
 
 p0 <- rep(p0, length(p2))
 c.be2 <- copgHs2(p1 = p0, p2 = p2, teta = theta, teta.st = theta, BivD = BivD, par2 = par2, min.pr = min.pr, max.pr = max.pr, min.dn = min.dn)$c.copula.be2

 if(y1 == 1) y*(1 - c.be2)*pdf2/(1-p0) else y*(c.be2*pdf2/p0) 
 
}







bincont.v <- function(y, p0, eta.mu, sigma, nu, theta, margin, BivD, par2, min.dn, min.pr, max.pr, y1){           
     
 rpd <- distrHsAT(y, eta2 = eta.mu, sigma2 = sigma, nu = nu, margin2 = margin, min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)     
     
 p2 <- rpd$p2; pdf2 <- rpd$pdf2
 
 p0 <- rep(p0, length(p2))
 c.be2 <- copgHs2(p1 = p0, p2 = p2, teta = theta, teta.st = theta, BivD = BivD, par2 = par2, min.pr = min.pr, max.pr = max.pr, min.dn = min.dn)$c.copula.be2

 if(y1 == 1) y^2*(1 - c.be2)*pdf2/(1-p0) else y^2*(c.be2*pdf2/p0) 
 
}



ifef <- function(dv){ # this will only work for ill-defined models that will result in convergence failure anyway

 dv <- ifelse(is.na(dv)  , .Machine$double.eps, dv ) 
 dv <- ifelse(dv == Inf  ,        8.218407e+20, dv )
 dv <- ifelse(dv == -Inf ,       -8.218407e+20, dv )
 dv

}




intB <- function(y, eta, sigma2, sigma2.st, nu, nu.st, margin, rc, min.dn, min.pr, max.pr, discr = FALSE, ym = NULL, left.trunc = 0){           
                                       
 if(discr == FALSE) pdf <- distrHsAT(y, eta, sigma2, nu, margin, min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)$pdf2
 if(discr == TRUE)  pdf <- distrHsDiscr(y, eta, sigma2, 1, 1, 1, margin, naive = TRUE, ym, min.dn = min.dn, 
                                        min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc)$pdf2

 log( 1 + exp( log( pdf ) + rc ) ) 
   
}















gradBbit1 <- function(y, eta, sigma2, sigma2.st, nu, nu.st, margin, rc, min.dn, min.pr, max.pr, discr = FALSE, ym = NULL, left.trunc = 0){ 

 

       if(discr == FALSE) dHs <- distrHs(y, eta, sigma2, sigma2.st, nu, nu.st, margin, naive = TRUE, min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)
       if(discr == TRUE)  dHs <- distrHsDiscr(y, eta, sigma2, sigma2.st, 1, 1, margin, naive = TRUE, ym, min.dn = min.dn, 
                                              min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc)

       pdf2            <- dHs$pdf2
       derpdf2.dereta2 <- dHs$derpdf2.dereta2 

       comp1 <- 1 + exp(log( pdf2 ) + rc) 
       comp2 <- pdf2/comp1

       dl.dbe <- derpdf2.dereta2/pdf2

      comp2*dl.dbe 
     
}




gradBbit2 <- function(y, eta, sigma2, sigma2.st, nu, nu.st, margin, rc, min.dn, min.pr, max.pr, discr = FALSE, ym = NULL, left.trunc = 0){ 



       if(discr == FALSE) dHs <- distrHs(y, eta, sigma2, sigma2.st, nu, nu.st, margin, naive = TRUE, min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)
       if(discr == TRUE)  dHs <- distrHsDiscr(y, eta, sigma2, sigma2.st, 1, 1, margin, naive = TRUE, ym, min.dn = min.dn, 
                                              min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc)

       pdf2                 <- dHs$pdf2
       derpdf2.dersigma2.st <- dHs$derpdf2.dersigma2.st  

       comp1 <- 1 + exp(log( pdf2 ) + rc) 
       comp2 <- pdf2/comp1

       dl.dsigma.st <- derpdf2.dersigma2.st/pdf2

      comp2*dl.dsigma.st 
       
}




gradBbit3 <- function(y, eta, sigma2, sigma2.st, nu, nu.st, margin, rc, min.dn, min.pr, max.pr, discr = FALSE, ym = NULL){ 

 

       dHs <- distrHs(y, eta, sigma2, sigma2.st, nu, nu.st, margin, naive = TRUE, min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)

       pdf2                 <- dHs$pdf2
       derpdf2.dernu.st     <- dHs$derpdf2.dernu.st  

       comp1 <- 1 + exp(log( pdf2 ) + rc) 
       comp2 <- pdf2/comp1

       dl.dnu.st <- derpdf2.dernu.st/pdf2

       comp2*dl.dnu.st 
       
       }
       
       
       
       
       
hessBbit1 <- function(y, eta, sigma2, sigma2.st, nu, nu.st, margin, rc, min.dn, min.pr, max.pr, discr = FALSE, ym = NULL, left.trunc = 0){ 



       if(discr == FALSE) dHs <- distrHs(y, eta, sigma2, sigma2.st, nu, nu.st, margin, naive = TRUE, min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)
       if(discr == TRUE)  dHs <- distrHsDiscr(y, eta, sigma2, sigma2.st, 1, 1, margin, naive = TRUE, ym, min.dn = min.dn, 
                                              min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc)

        pdf2                         <- dHs$pdf2
        derpdf2.dereta2              <- dHs$derpdf2.dereta2
        der2pdf2.dereta2             <- dHs$der2pdf2.dereta2
        
        comp1 <- 1 + exp(log(pdf2) + rc)
        comp2 <- pdf2/comp1
        comp3 <- pdf2/comp1^2
        
        d2l.be.be <- (der2pdf2.dereta2 * pdf2 - (derpdf2.dereta2)^2)/pdf2^2
        dl.dbe    <- derpdf2.dereta2/pdf2
        
        comp2*d2l.be.be + dl.dbe^2*comp3 

} 




hessBbit2 <- function(y, eta, sigma2, sigma2.st, nu, nu.st, margin, rc, min.dn, min.pr, max.pr, discr = FALSE, ym = NULL, left.trunc = 0){ 



       if(discr == FALSE) dHs <- distrHs(y, eta, sigma2, sigma2.st, nu, nu.st, margin, naive = TRUE, min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)
       if(discr == TRUE)  dHs <- distrHsDiscr(y, eta, sigma2, sigma2.st, 1, 1, margin, naive = TRUE, ym, min.dn = min.dn, 
                                              min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc)

        pdf2                         <- dHs$pdf2
        derpdf2.dersigma2.st         <- dHs$derpdf2.dersigma2.st
        der2pdf2.dersigma2.st2       <- dHs$der2pdf2.dersigma2.st2
        
        comp1 <- 1 + exp(log(pdf2) + rc)
        comp2 <- pdf2/comp1
        comp3 <- pdf2/comp1^2
        
        d2l.sigma.sigma <- (der2pdf2.dersigma2.st2 * pdf2 - (derpdf2.dersigma2.st)^2)/pdf2^2
        dl.dsigma.st    <- derpdf2.dersigma2.st/pdf2

        comp2*d2l.sigma.sigma + dl.dsigma.st^2*comp3 
        
   
}




hessBbit3 <- function(y, eta, sigma2, sigma2.st, nu, nu.st, margin, rc, min.dn, min.pr, max.pr, discr = FALSE, ym = NULL, left.trunc = 0){ 



       if(discr == FALSE) dHs <- distrHs(y, eta, sigma2, sigma2.st, nu, nu.st, margin, naive = TRUE, min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)
       if(discr == TRUE)  dHs <- distrHsDiscr(y, eta, sigma2, sigma2.st, 1, 1, margin, naive = TRUE, ym, min.dn = min.dn, 
                                              min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc)

        pdf2                         <- dHs$pdf2
        derpdf2.dereta2              <- dHs$derpdf2.dereta2
        derpdf2.dersigma2.st         <- dHs$derpdf2.dersigma2.st
        der2pdf2.dereta2dersigma2.st <- dHs$der2pdf2.dereta2dersigma2.st
        
        comp1 <- 1 + exp(log(pdf2) + rc)
        comp2 <- pdf2/comp1
        comp3 <- pdf2/comp1^2
        
        d2l.be.sigma <- (der2pdf2.dereta2dersigma2.st * pdf2 - derpdf2.dereta2 * derpdf2.dersigma2.st)/pdf2^2
        
        dl.dbe       <- derpdf2.dereta2/pdf2
        dl.dsigma.st <- derpdf2.dersigma2.st/pdf2

       comp2*d2l.be.sigma + dl.dbe*dl.dsigma.st*comp3 
    
}




hessBbit4 <- function(y, eta, sigma2, sigma2.st, nu, nu.st, margin, rc, min.dn, min.pr, max.pr, discr = FALSE, ym = NULL){ 

 

        dHs <- distrHs(y, eta, sigma2, sigma2.st, nu, nu.st, margin, naive = TRUE, min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)

        pdf2               <- dHs$pdf2
        derpdf2.dernu.st   <- dHs$derpdf2.dernu.st           
        der2pdf2.dernu.st2 <- dHs$der2pdf2.dernu.st2           
        
        comp1 <- 1 + exp(log(pdf2) + rc)
        comp2 <- pdf2/comp1
        comp3 <- pdf2/comp1^2
        
        d2l.nu.nu <- (der2pdf2.dernu.st2*pdf2-(derpdf2.dernu.st)^2)/pdf2^2
        dl.dnu.st <- derpdf2.dernu.st/pdf2

        comp2*d2l.nu.nu + dl.dnu.st^2*comp3 
 
}




hessBbit5 <- function(y, eta, sigma2, sigma2.st, nu, nu.st, margin, rc, min.dn, min.pr, max.pr, discr = FALSE, ym = NULL){ 

 

        dHs <- distrHs(y, eta, sigma2, sigma2.st, nu, nu.st, margin, naive = TRUE, min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)

        pdf2                     <- dHs$pdf2
        derpdf2.dereta2          <- dHs$derpdf2.dereta2
        der2pdf2.dereta2dernu.st <- dHs$der2pdf2.dereta2dernu.st   
        derpdf2.dernu.st         <- dHs$derpdf2.dernu.st           
        
        comp1 <- 1 + exp(log(pdf2) + rc)
        comp2 <- pdf2/comp1
        comp3 <- pdf2/comp1^2
        
        d2l.be.nu  <- (der2pdf2.dereta2dernu.st*pdf2 - derpdf2.dereta2*derpdf2.dernu.st)/pdf2^2 
        dl.dbe     <- derpdf2.dereta2/pdf2
        dl.dnu.st  <- derpdf2.dernu.st/pdf2

        comp2*d2l.be.nu + dl.dbe*dl.dnu.st*comp3 

}




hessBbit6 <- function(y, eta, sigma2, sigma2.st, nu, nu.st, margin, rc, min.dn, min.pr, max.pr, discr = FALSE, ym = NULL){ 



        dHs <- distrHs(y, eta, sigma2, sigma2.st, nu, nu.st, margin, naive = TRUE, min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)

        pdf2                          <- dHs$pdf2
        derpdf2.dersigma2.st          <- dHs$derpdf2.dersigma2.st
        der2pdf2.dersigma2.stdernu.st <- dHs$der2pdf2.sigma2.st2dernu.st
        derpdf2.dernu.st              <- dHs$derpdf2.dernu.st           
        
        comp1 <- 1 + exp(log(pdf2) + rc)
        comp2 <- pdf2/comp1
        comp3 <- pdf2/comp1^2
        
        d2l.sigma.nu  <- (der2pdf2.dersigma2.stdernu.st*pdf2-(derpdf2.dersigma2.st*derpdf2.dernu.st))/pdf2^2 

        dl.dsigma.st  <- derpdf2.dersigma2.st/pdf2
        dl.dnu.st     <- derpdf2.dernu.st/pdf2

        comp2*d2l.sigma.nu + dl.dsigma.st*dl.dnu.st*comp3
       
}


###################################################################################################
###################################################################################################

int1f <- function(y, eta, sigma2, nu, margin, rc, min.dn, min.pr, max.pr){ 



   pdf <- distrHsAT(y, eta, sigma2, nu, margin, min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)$pdf2
   
   log( 1 + exp( log( pdf ) + rc ) )
   
}




d.bpsi <- function(y, X1, X2, X3, eta, sigma2, sigma2.st, nu, nu.st, margin, rc, j, min.dn, min.pr, max.pr){ 


       dHs <- distrHs(y, eta, sigma2, sigma2.st, nu, nu.st, margin, naive = TRUE, min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)

       pdf2                 <- dHs$pdf2
       derpdf2.dereta2      <- dHs$derpdf2.dereta2 
       derpdf2.dersigma2.st <- dHs$derpdf2.dersigma2.st  
       derpdf2.dernu.st     <- dHs$derpdf2.dernu.st  

       comp1 <- 1 + exp(log( pdf2 ) + rc) 
       comp2 <- pdf2/comp1

       dl.dbe       <- derpdf2.dereta2/pdf2
       dl.dsigma.st <- derpdf2.dersigma2.st/pdf2
       dl.dnu.st    <- derpdf2.dernu.st/pdf2


       if( margin %in% c("DAGUM","SM","TW") ) res <- cbind( comp2*as.numeric(dl.dbe)%*%t(X1), comp2*as.numeric(dl.dsigma.st)%*%t(X2), comp2*as.numeric(dl.dnu.st)%*%t(X3) ) else
                                              res <- cbind( comp2*as.numeric(dl.dbe)%*%t(X1), comp2*as.numeric(dl.dsigma.st)%*%t(X2) )
      
       res[, j]
}




gradF <- function(params, n, VC, margin, lB, uB, rc, min.dn, min.pr, max.pr){


  G <- matrix(NA, n, length(params))

for(i in 1:n){

  X1 <- VC$X1[i,]
  X2 <- VC$X2[i,]
  X3 <- VC$X3[i,]
  nu <- nu.st <- 1
  
  eta       <- X1%*%params[1:VC$X1.d2]
  sigma2.st <- X2%*%params[(1+VC$X1.d2):(VC$X1.d2+VC$X2.d2)]
  
     if( margin %in% c("DAGUM","SM","TW") ){ 
  
       nu.st <- X3%*%params[(1+VC$X1.d2+VC$X2.d2):(VC$X1.d2+VC$X2.d2+VC$X3.d2)] 
 
       ss <- enu.tr(nu.st, margin)  
       nu.st <- ss$vrb.st 
       nu    <- ss$vrb   
                                      }
 
  ss        <- esp.tr(sigma2.st, margin)
  sigma2.st <- ss$vrb.st
  sigma2    <- ss$vrb
                                    
  for(j in 1:length(params)) G[i, j] <- integrate(d.bpsi, lB, uB, X1, X2, X3, eta, sigma2, sigma2.st, nu, nu.st, margin, rc, j, min.dn, min.pr, max.pr)$value

             }
             
  colSums(G)
  
}


