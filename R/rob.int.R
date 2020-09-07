rob.int <- function(x, rc, l.grid = 1000, tol = 1e-4, var.range = NULL){


lo <- l.grid

margin <- x$VC$margins[2]
params <- x$coefficients


min.dn <- 1e-160 # x$VC$min.dn
min.pr <- x$VC$min.pr
max.pr <- x$VC$max.pr 



j <- 1; inds <- posi <- NULL 

if(is.null(var.range)){

rlo <- (max(x$VC$y1) - min(x$VC$y1))/lo
if(rlo > 1) lo <- round(lo*rlo) 


if( margin %in% c("N","N2","GU","rGU","LO","LN") )                            seq.y <- seq(min(x$VC$y1) - ((max(x$VC$y1) - min(x$VC$y1))/2), max(x$VC$y1) + ((max(x$VC$y1) - min(x$VC$y1))/2), length.out = lo)
if( margin %in% c("WEI","iG","GA","DAGUM","SM","FISK","GP","GPII","GPo")  )   seq.y <- seq(1e-12, max(x$VC$y1) + ((max(x$VC$y1) - min(x$VC$y1))/2), length.out = lo)
if( margin %in% c("TW")  )                                                    seq.y <- seq(0, max(x$VC$y1) + ((max(x$VC$y1) - min(x$VC$y1))/2), length.out = lo)
if( margin %in% c("BE")  )                                                    seq.y <- seq(1e-12, 0.999999,  length.out = lo)      

}else{

rlo <- (var.range[2] - var.range[1])/lo
if(rlo > 1) lo <- round(lo*rlo) 
seq.y <- seq(var.range[1], var.range[2], length.out = lo)  


}



###################################################################################################
###################################################################################################

n <- x$VC$n
if(is.null(x$VC$X2)){x$VC$X2 <- matrix(1, n, 1); x$VC$X2.d2 <- 1} 
if(is.null(x$VC$X3)){x$VC$X3 <- matrix(1, n, 1); x$VC$X3.d2 <- 1} 


eta <- eta.tr(x$VC$X1%*%params[1:x$VC$X1.d2], margin)
ss  <- esp.tr(x$VC$X2%*%params[(1+x$VC$X1.d2):(x$VC$X1.d2+x$VC$X2.d2)], margin)

sigma2    <- ss$vrb
sigma2.st <- ss$vrb.st

if( margin %in% c("DAGUM","SM","TW") ){

            nus   <- enu.tr(x$VC$X3%*%params[(1+x$VC$X1.d2+x$VC$X2.d2):(x$VC$X1.d2+x$VC$X2.d2+x$VC$X3.d2)], margin)
            nu    <- nus$vrb
            nu.st <- nus$vrb.st
            
} else nu <- nu.st <- 1


###################################################################################################
###################################################################################################




for(i in 1:lo){ 

  ires <- intB(seq.y[i], eta, sigma2, sigma2.st, nu, nu.st, margin, rc, min.dn, min.pr, max.pr) 
  
  if(min(ires) == max(ires) && min(ires) < tol) inds[i] <- TRUE else inds[i] <- FALSE

  
  if(i > 1 && inds[i] != inds[i-1]) { if(j==1){ posi[j] <- i-1; j <- j + 1} else posi[2] <- i}  

}

if((is.null(posi) || length(posi) < 2) && margin %in% c("N","N2","GU","rGU","LO","LN")) stop("Increase the tolerance value or try a different range.")





if((is.null(posi) || length(posi) < 2) && margin %in% c("WEI","iG","GA","DAGUM","TW","SM","FISK","GP","GPII","GPo") ){


posi <- NULL; j <- 1

  for(i in 1:lo){ 

    ires <- intB(seq.y[i], eta, sigma2, sigma2.st, nu, nu.st, margin, rc, min.dn, min.pr, max.pr) 
  
    if(min(ires) == max(ires) && min(ires) < tol) inds[i] <- TRUE else inds[i] <- FALSE 
    
    if(i > 2 && inds[i] != inds[i-1]) { if(j == 1) posi[2] <- i; j <- j + 1}
    
    if(!is.null(posi[2])) posi[1] <- 1 

    }

  if(is.null(posi) || length(posi) < 2) stop("Increase the tolerance value or try a different range.")


}



bs <- c(seq.y[posi[1]], seq.y[posi[2]])
names(bs) <- c("lB", "uB")
bs

}
