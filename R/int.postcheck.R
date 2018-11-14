int.postcheck <- function(mo, margin, n.rep = 50, prob.lev = 0.05, y2m, eq = 1){

cont <- c("N", "N2", "GU", "rGU", "LO", "LN", "WEI","iG", "GA", "DAGUM", "SM", "BE", "FISK")
disc <- c("NBI", "NBII", "PIG", "PO", "ZTP") 



n   <- mo$n

if(margin == "ZTP") rZTP <- function(n, mu) qpois(runif(n, dpois(0, mu), 1), mu)


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

}



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






qrs <- matrix(0, n, n.rep)



for(i in 1:n.rep){

if(margin == "N")     y2s <- rNO(   n,    mu = eta2,         sigma = sqrt(sigma2)) 
if(margin == "N2")    y2s <- rNO(   n,    mu = eta2,         sigma = sigma2) 
if(margin == "GU")    y2s <- rGU(   n,    mu = eta2,         sigma = sqrt(sigma2)) 
if(margin == "rGU")   y2s <- rRG(   n,    mu = eta2,         sigma = sqrt(sigma2)) 
if(margin == "LO")    y2s <- rLO(   n,    mu = eta2,         sigma = sqrt(sigma2)) 
if(margin == "LN")    y2s <- rLOGNO(n,    mu = eta2,         sigma = sqrt(sigma2)) 
if(margin == "WEI")   y2s <- rWEI(  n,    mu = exp(eta2),    sigma = sqrt(sigma2)) 
#if(margin == "GO")    y2s <- rgompertz(n, shape = exp(eta2), rate = sqrt(sigma2))
if(margin == "iG")    y2s <- rIG(   n,    mu = exp(eta2),    sigma = sqrt(sigma2)) 
if(margin == "GA")    y2s <- rGA(   n,    mu = exp(eta2),    sigma = sqrt(sigma2)) 
if(margin == "GAi")   y2s <- rGA(   n,    mu = eta2,         sigma = sqrt(sigma2)) 

#if(margin == "GA2")   y2s <- rgamma(n, shape = sqrt(sigma2), rate = exp(eta2) ) 
#if(margin == "GGA")   y2s <- exp( log(rgamma(n, shape = nu))/sqrt(sigma2) + log(exp(eta2)) )
 
if(margin == "DAGUM") y2s <- rGB2(  n,    mu = exp(eta2),    sigma = sqrt(sigma2), nu = nu, tau = 1) 
if(margin == "SM")    y2s <- rGB2(  n,    mu = exp(eta2),    sigma = sqrt(sigma2), nu = 1,  tau = nu) 
if(margin == "BE")    y2s <- rBE(   n,    mu = plogis(eta2), sigma = sqrt(sigma2))
if(margin == "FISK")  y2s <- rGB2(  n,    mu = exp(eta2),    sigma = sqrt(sigma2), nu = 1,  tau = 1)
if(margin == "NBI")   y2s <- rNBI(  n,    mu = exp(eta2),    sigma = sqrt(sigma2)) 
if(margin == "NBII")  y2s <- rNBII( n,    mu = exp(eta2),    sigma = sqrt(sigma2)) 
if(margin == "PIG")   y2s <- rPIG(  n,    mu = exp(eta2),    sigma = sqrt(sigma2)) 
if(margin == "PO")    y2s <- rPO(   n,    mu = exp(eta2)) 
if(margin == "ZTP")   y2s <- rZTP(  n,    mu = exp(eta2))   

if(margin %in% cont) qrs[,i] <- sort(  qnorm(  distrHsAT(y2s, eta2, sigma2, nu, margin)$p2 )  )


if(margin %in% disc){

if(margin %in% c("ZTP")){
    ly2 <- length(y2s)
    y2ms <- list()
    my2s <- max(y2s)
    for(j in 1:ly2){ y2ms[[j]] <- seq(0, y2s[j]); length(y2ms[[j]]) <- my2s+1} 
    y2ms <- do.call(rbind, y2ms)          
                         }

tmp  <- distrHsATDiscr(y2s, eta2, sigma2, nu, margin, y2m = y2ms)
tmpp <- tmp$p2
tmpd <- tmp$pdf2   
qrs[,i] <- sort( qnorm( runif(y2s, tmpp - tmpd, tmpp) )  )
                    }
                    
                    
                    
if(mo$VC$ccss == "yes") qrs[,i] <- qrs[,i] - mean(qrs[,i]) 
                    
}

Qq <- quantile(as.numeric(qrs), (1:n - 0.5)/n)
n  <- length(Qq)

lim <- apply(qrs, 1, FUN = quantile, p = c(prob.lev/2, 1 - prob.lev/2))




if(margin %in% cont) qr <- qnorm( distrHsAT(y2, eta2, sigma2, nu, margin)$p2 ) 


if(margin %in% disc){

tmp  <- distrHsATDiscr(y2, eta2, sigma2, nu, margin, y2m = y2m)
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
