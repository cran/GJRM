pred.mvt <- function(x, eq, fun = "mean", n.sim = 100, prob.lev = 0.05, ...){


 if( !(fun %in% c("mean", "variance", "tau")) ) stop("The value for fun can either be mean, variance or tau.")
 
 
 if(fun != "tau"){
 
 if(missing(eq)) stop("You must provide the equation number.")
 if(!(eq %in% c(1, 2))) stop("The value for eq can be either 1 or 2.")

 if(!(x$margins[1] %in% c(x$VC$m2, x$VC$m2d)) && !(x$margins[2] %in% c(x$VC$m2, x$VC$m2d)) ) stop("This is currently implemented for margins with two parameters.\nGet in touch to check progress on the other cases.") 

 if(eq == 1) mar <- x$margins[1]
 if(eq == 2) mar <- x$margins[2] 
 
 bs <- rMVN(n.sim, mean = x$coefficients, sigma = x$Vb)
 
 
 
if( x$margins[1] %in% c(x$VC$m2, x$VC$m2d) && x$margins[2] %in% c(x$VC$m2, x$VC$m2d) ){ ###############
 
  
 if(eq == 1){ind1 <- (1:x$X1.d2); ind2 <- (x$X1.d2 + x$X2.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2) } 
 if(eq == 2){ind1 <- (x$X1.d2 + 1):(x$X1.d2 + x$X2.d2); ind2 <- (x$X1.d2 + x$X2.d2 + x$X3.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2) }  

if(eq == 1){
 Xfit1 <- predict(x, eq = 1, type = "lpmatrix", ...)
 Xfit2 <- predict(x, eq = 3, type = "lpmatrix", ...) 
 fit1  <- predict(x, eq = 1, type = "link", ...)
 fit2  <- predict(x, eq = 3, type = "link", ...) 
           }

if(eq == 2){
 Xfit1 <- predict(x, eq = 2, type = "lpmatrix", ...)
 Xfit2 <- predict(x, eq = 4, type = "lpmatrix", ...) 
 fit1  <- predict(x, eq = 2, type = "link", ...)
 fit2  <- predict(x, eq = 4, type = "link", ...)  
           }

fit1Sim <- Xfit1%*%t(bs[, ind1])
fit2Sim <- Xfit2%*%t(bs[, ind2])


}############




if( x$margins[1] %in% c(x$VC$bl, x$VC$m1d) && x$margins[2] %in% c(x$VC$m2, x$VC$m2d) ){ #############
 
  
 if(eq == 1){ind1 <- (1:x$X1.d2) } 
 if(eq == 2){ind1 <- (x$X1.d2 + 1):(x$X1.d2 + x$X2.d2); ind2 <- (x$X1.d2 + x$X2.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2) }  

if(eq == 1){
 Xfit1 <- predict(x, eq = 1, type = "lpmatrix", ...)
 fit1  <- predict(x, eq = 1, type = "link", ...)
           }

if(eq == 2){
 Xfit1 <- predict(x, eq = 2, type = "lpmatrix", ...)
 Xfit2 <- predict(x, eq = 3, type = "lpmatrix", ...) 
 fit1  <- predict(x, eq = 2, type = "link", ...)
 fit2  <- predict(x, eq = 3, type = "link", ...)  
           }

fit1Sim <- Xfit1%*%t(bs[, ind1])
if(eq == 2) fit2Sim <- Xfit2%*%t(bs[, ind2])


}##################















 if(mar %in% c("LO")){
 
    if(fun == "mean")     {fit <- fit1;             fitSim <- fit1Sim            } 
    if(fun == "variance") {fit <- pi^2*exp(fit2)/3; fitSim <- pi^2*exp(fit2Sim)/3}
                                     
 } 
 

 if(mar %in% c("LN")){
 
    if(fun == "mean")     {fit <- exp(fit1)*sqrt(exp(exp(fit2)));                    fitSim <- exp(fit1Sim)*sqrt(exp(exp(fit2Sim)))                       } 
    if(fun == "variance") {fit <- exp(exp(fit2))*( exp(exp(fit2)) - 1 )*exp(2*fit1); fitSim <- exp(exp(fit2Sim))*( exp(exp(fit2Sim)) - 1 )*exp(2*fit1Sim) }
                                    
 } 
 

 if(mar %in% c("iG")){
 
    if(fun == "mean")     {fit <- exp(fit1)^3;           fitSim <- exp(fit1Sim)^3              } 
    if(fun == "variance") {fit <- exp(fit1)^3*exp(fit2); fitSim <- exp(fit1Sim)^3*exp(fit2Sim) }
                                    
 }  
 
 
 
  
 
CIpred <- rowQuantiles(fitSim, probs = c(prob.lev/2, 1-prob.lev/2), na.rm = TRUE)


}



if(fun == "tau"){


res <- jc.probs(x, mean(x$y1), mean(x$y2), intervals = TRUE, n.sim = n.sim, prob.lev = prob.lev, ...)

fit    <- res[, 6]
CIpred <- res[, 7:8]


}


list(pred = fit, CIpred = CIpred)

}


