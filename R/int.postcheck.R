int.postcheck <- function(mo, margin, n.rep = 50, prob.lev = 0.05, y2m, eq = 1){

cont <- c("N", "GU", "rGU", "LO", "LN", "WEI","iG", "GA", "DAGUM", "SM", "BE", "FISK","GP","GPII","GPo","TW")
disc <- c("NBI", "NBII", "PIG", "PO","DGP0" ,"ZTP","DGP","DGPII") 

n   <- mo$n

if(mo$Cont == "NO"){

y2     <- mo$y2
eta2   <- mo$eta2
sigma2 <- mo$sigma2 
nu     <- mo$nu


if(mo$VC$ccss == "yes"){
   eta2 <- eta2[mo$inde]
   if(!is.null(mo$X3) && !(mo$VC$margins[2] %in% mo$VC$m1d) ){
   sigma2 <- sigma2[mo$inde]
   nu     <- nu[mo$inde]
                                                             }
n <- mo$n.sel                                                             
                      }


}

if(mo$univar.gamlss == TRUE){

y2     <- mo$y1
eta2   <- mo$eta1
sigma2 <- mo$sigma21 <- mo$sigma2 
nu     <- mo$nu

}else{



if(mo$Cont == "YES"){

if(eq == 1){

y2     <- mo$y1
eta2   <- mo$eta1
sigma2 <- mo$sigma21 
nu     <- mo$nu1

}

if(eq == 2){

y2     <- mo$y2
eta2   <- mo$eta2
sigma2 <- mo$sigma22
nu     <- mo$nu2

}


}


}



qrs <- matrix(0, n, n.rep)



for(i in 1:n.rep){


y2s <- sim.resp(margin, n, eta2, sigma2, nu, setseed = FALSE)





if(margin %in% cont){


p22s <- distrHsAT(y2s, eta2, sigma2, nu, margin, min.dn = mo$VC$min.dn, min.pr = mo$VC$min.pr, max.pr = mo$VC$max.pr)$p2

if(margin %in% c("TW")) p22s[y2s == 0] <- runif(sum(y2s == 0), min = 0, max = p22s[y2s == 0]) 

qrs[,i] <- sort(  qnorm( mm(p22s, min.pr = 1e-300, max.pr = 0.9999999999999999)  )  )     



}


if(margin %in% disc){

if(margin %in% c("ZTP","DGP","DGPII","DGP0")){
    ly2 <- length(y2s)
    y2ms <- list()
    my2s <- max(y2s)
    
    
    
    if(margin %in% c("DGP","DGPII","DGP0")) for(j in 1:ly2){ y2ms[[j]] <- seq(0, y2s[j]); length(y2ms[[j]]) <- my2s+1} 
    if(margin %in% c("ZTP"))  for(j in 1:ly2){ y2ms[[j]] <- seq(1, y2s[j]); length(y2ms[[j]]) <- my2s} 
    
    y2ms <- do.call(rbind, y2ms) 
    
    #if(margins != "ZTP") for(i in 1:ly1){ y1m[[i]] <- seq(0, y1[i]); length(y1m[[i]]) <- my1+1} 
    #if(margins == "ZTP") for(i in 1:ly1){ y1m[[i]] <- seq(1, y1[i]); length(y1m[[i]]) <- my1}     
    
                         }

tmp  <- distrHsATDiscr(y2s, eta2, sigma2, nu, margin, y2m = y2ms, min.dn = mo$VC$min.dn, min.pr = mo$VC$min.pr, max.pr = mo$VC$max.pr)
tmpp <- tmp$p2
tmpd <- tmp$pdf2   
qrs[,i] <- sort( qnorm( mm(runif(y2s, tmpp - tmpd, tmpp), min.pr = 1e-300, max.pr = 0.9999999999999999) )  )
                    }
                    
  
  
  
                    
if(mo$VC$ccss == "yes") qrs[,i] <- qrs[,i] - mean(qrs[,i]) 
                    
}

Qq <- quantile(as.numeric(qrs), (1:n - 0.5)/n)
n  <- length(Qq)

lim <- apply(qrs, 1, FUN = quantile, p = c(prob.lev/2, 1 - prob.lev/2))




if(margin %in% cont){


p22 <- distrHsAT(y2, eta2, sigma2, nu, margin, min.dn = mo$VC$min.dn, min.pr = mo$VC$min.pr, max.pr = mo$VC$max.pr)$p2

if(margin %in% c("TW")) p22[y2 == 0] <- runif(sum(y2 == 0), min = 0, max = p22[y2 == 0]) 


qr <- qnorm( p22 ) 



}


if(margin %in% disc){

tmp  <- distrHsATDiscr(y2, eta2, sigma2, nu, margin, y2m = y2m, min.dn = mo$VC$min.dn, min.pr = mo$VC$min.pr, max.pr = mo$VC$max.pr)
tmpp <- tmp$p2
tmpd <- tmp$pdf2   
set.seed(100)
qr <- qnorm( runif(y2, tmpp - tmpd, tmpp) )  

                    }


if(mo$VC$ccss == "yes") qr <- qr - mean(qr) 

qqplot(Qq, qr, pch = ".", ylim = range(c(lim, qr)), xlab = "Theoretical Quantiles",  ylab = "Quantile Residuals")
polygon(c(Qq, Qq[n:1], Qq[1]), c(lim[1, ], lim[2, n:1], lim[1, 1]), col = "grey80", border = NA)
abline(0, 1, col = "red")
points(Qq, sort(qr), pch = ".")


}
