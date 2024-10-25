PE <- function(x1, idx, n.sim = 100, prob.lev = 0.05, 
                  plot = FALSE, 
                  main = "Histogram of Simulated Average Effects", 
                  xlab = "Simulated Average Effects", ...){


etap.noi <- X.int <- X.noi <- eti1 <- eti0 <- etno <- indS <- bs <- ind.excl <- p.int1 <- p.int0 <- d.int1 <- d.int0 <- p.etn <- d.etn <- ass.p <- ass.pst <- C.11 <- C.10 <- sig2 <- peti1s <- peti0s <- sigma2.st <- sigma2s <- eti1s <- eti0s <- d0 <- d1 <- p.etns <- etnos <- etds <- ass.ps <- 1
diffEf <- fy1.y2 <- est.ATso <- y2 <- CIF <- Pr <- Effects <- C.11 <- C.10 <- C.11s <- C.10s <- v0s <- v1s <- ATs <- AUXs <- NULL
nu  <- 1 
nus <- rep(1, n.sim)



if(!(x1$margins[1] %in% x1$bl) ) stop("Error in first margin value. It should be one of:\nprobit, logit, cloglog.")
if(!(x1$margins[2] %in% x1$bl) ) stop("Error in second margin value. It should be one of:\nprobit, logit, cloglog.")



if( missing(idx) ) stop("You must provide a value for idx.")  



m2       <- x1$VC$m2 
m3       <- x1$VC$m3 
bin.link <- x1$VC$bl  

end <- 0
est.ATb <- NA
indD <- list()


########################################################
# model x1 - binary binary
########################################################

ind1  <- 1:x1$X1.d2 
ind2  <- x1$X1.d2 + (1:x1$X2.d2)




X1 <- as.matrix(x1$X1[idx,])
X2 <- as.matrix(x1$X2[idx,])

eta1 <- t(X1)%*%x1$coefficients[ind1] 
eta2 <- t(X2)%*%x1$coefficients[ind2]  


p.1 <- probm(eta1, x1$margins[1], min.dn = x1$VC$min.dn, min.pr = x1$VC$min.pr, max.pr = x1$VC$max.pr)$pr 
p.2 <- probm(eta2, x1$margins[1], min.dn = x1$VC$min.dn, min.pr = x1$VC$min.pr, max.pr = x1$VC$max.pr)$pr

if( is.null(x1$X3)  )   ass.p <- x1$theta    

if( !is.null(x1$X3) ) { X3    <- as.matrix( x1$X3[idx,] )
                        etd   <- t(X3)%*%x1$coefficients[(x1$X1.d2+x1$X2.d2+1):(x1$X1.d2+x1$X2.d2+x1$X3.d2)] 
                        ass.p <- teta.tr(x1$VC, etd)$teta
                       }                                                        #  ass.p <- (x1$theta[index1])[idx]  



AUX   <- mm( BiCDF(p.1, p.2, x1$nC, ass.p), min.pr = x1$VC$min.pr, max.pr = x1$VC$max.pr  )  # this is p11
C.11  <- AUX / p.1
C.10  <- mm( p.2 - AUX, min.pr = x1$VC$min.pr, max.pr = x1$VC$max.pr )   / mm(1 - p.1, min.pr = x1$VC$min.pr, max.pr = x1$VC$max.pr)   # this is actually C01, this has been corrected
              
              
######
# CIs
######
              
bs <- rMVN(n.sim, mean = x1$coefficients, sigma=x1$Vb)

eta1s <- t(X1)%*%t(bs[,ind1]) 
eta2s <- t(X2)%*%t(bs[,ind2])  

p.1s <- probm(eta1s, x1$margins[1], min.dn = x1$VC$min.dn, min.pr = x1$VC$min.pr, max.pr = x1$VC$max.pr)$pr 
p.2s <- probm(eta2s, x1$margins[1], min.dn = x1$VC$min.dn, min.pr = x1$VC$min.pr, max.pr = x1$VC$max.pr)$pr

if( !is.null(x1$X3) ) etds <- t(X3)%*%t(bs[,(x1$X1.d2+x1$X2.d2+1):(x1$X1.d2+x1$X2.d2+x1$X3.d2)])
if(  is.null(x1$X3) ) etds <- bs[, length(x1$coefficients)]
   

resT    <- teta.tr(x1$VC, etds)
ass.ps  <- resT$teta 



if(x1$BivD == "N") {

	for(i in 1:n.sim){ 
	
	         AUXs[i]  <- mm( BiCDF(p.1s[i], p.2s[i], x1$nC, ass.ps[i], test = FALSE), min.pr = x1$VC$min.pr, max.pr = x1$VC$max.pr  ) 
		 C.11s[i] <- AUXs[i]  / p.1s[i]
		 C.10s[i] <- mm( p.2s[i] - AUXs[i], min.pr = x1$VC$min.pr, max.pr = x1$VC$max.pr ) / mm(1 - p.1s[i], min.pr = x1$VC$min.pr, max.pr = x1$VC$max.pr)  
	                  }                 
}

if(x1$BivD != "N") {


 AUXs  <- mm( BiCDF(p.1s, p.2s, x1$nC, ass.ps, test = FALSE), min.pr = x1$VC$min.pr, max.pr = x1$VC$max.pr  ) 
 C.11s <- AUXs / p.1s 
 C.10s <- mm( p.2s - AUXs, min.pr = x1$VC$min.pr, max.pr = x1$VC$max.pr ) / mm(1 - p.1s, min.pr = x1$VC$min.pr, max.pr = x1$VC$max.pr)
             
}






AT <- C.11 - C.10  

   
   
##################



ATs <- C.11s - C.10s



 if(plot == TRUE){
  
  hist(ATs, freq = FALSE, main = main, xlab = xlab, 
       ylim = c(0, max(density(ATs)$y, hist(ATs, plot = FALSE)$density)), ...)
  lines(density(ATs))
 
 }




CIs <- as.numeric(quantile(ATs, c(prob.lev/2, 1 - prob.lev/2), na.rm = TRUE))

res <- c(CIs[1], AT, CIs[2])

out <- list(res = res, prob.lev = prob.lev) # , lims = c(lil, uil) ) 

class(out) <- "PE"

out

}

