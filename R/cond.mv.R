cond.mv <- function(x, eq, y1 = NULL, y2 = NULL, newdata, fun = "mean", n.sim = 100, prob.lev = 0.05){# , bin.model = NULL){

nu <- res.mCI <- res.mCITemp <- res.vCI <- res.vCITemp <- nu1 <- nu2 <- dof <- dofS <- thetS <- nu1S <- nu2S <- sig1S <- sig2S <- sig1 <- sig2 <- sigma <- NULL
pb <- pbS <- 1


if(class(x)[1] == "gamlss")              stop("Function not suitable for univariate models.")

if(missing(eq))                          stop("An equation number must be provided.")
if(!missing(eq) && !(eq %in% c(1, 2)))   stop("The equation number can be either 1 or 2.")


if(missing(newdata))                               stop("You must provide a data frame.")
if(!missing(newdata)){ 
                       if(!is.data.frame(newdata)) stop("You must provide a data frame.")
                       if(dim(newdata)[1] > 1)     stop("This calculation is supported for single row data frames.")
                     }

if( !(fun %in% c("mean", "variance")) ) stop("The value for fun can be either mean or variance.")


m2  <- x$VC$m2
m3  <- x$VC$m3
m1d <- c("P", "tP")
m2d <- c("NBI", "NBII", "PIG","tNBI", "tNBII", "tPIG")

lt1 <- x$VC$left.trunc1
lt2 <- x$VC$left.trunc2



margin <- x$margins[eq]

if(margin %in% x$VC$bl)                                    stop("This function is of no use for a binary margin.")
if(margin %in% c("GP","GPII","GPo","DGP","DGPII","DGP0") ) stop("Function not ready yet for the chosen distribution(s).")


if( margin %in% c("N","GU","rGU","LO") )                            { lB <- -Inf;                      uB <- Inf}
if( margin %in% c("LN","WEI","IG","GA","DAGUM","SM","FISK","TW")  ) { lB <- sqrt(.Machine$double.eps); uB <- Inf} 
if( margin %in% c("BE")  )                                          { lB <- sqrt(.Machine$double.eps); uB <- 0.999999}



 #if(!is.null(bin.model)){
 #
 #if(class(bin.model)[1] != "gam")           stop("The binary model has to be fitted using gam() from mgcv.")
 #if(family(bin.model)$family != "binomial") stop("A binary model is required here.")
 #
 #
 #Xb  <- predict(bin.model, newdata = newdata, type = "lpmatrix")
 #bss <- rMVN(n.sim, mean = bin.model$coefficients, sigma = bin.model$Vp)
 # 
 #pb  <- predict(bin.model, newdata = newdata, type = "response")
 #
 #if( family(bin.model)$link == "probit")  pbS <- pnorm(Xb%*%t(bss)) 
 #if( family(bin.model)$link == "logit")   pbS <- plogis(Xb%*%t(bss)) 
 #if( family(bin.model)$link == "cloglog") pbS <- 1-exp(-exp(Xb%*%t(bss)))
 #
 #} 
 


########################################################################################################################################


if( !is.null(x$VC$K1) && x$margins[2] %in% m2){ # ord-cont

   if( is.null(y1) )                                   stop("A value for y1 must be provided.")
   if(!(y1 >= range(x$y1)[1] && y1 <= range(x$y1)[2])) stop("y1 is not within its observed range.") 


pk1 <- 0

pkt <- predict(x, eq = 1, newdata = newdata, type = "response")$p1.cum 
pk  <- pkt[, y1]
if(y1 > 1) pk1 <- pkt[, y1 - 1]
                                                    
eta2 <- c(eta.tr(predict(x, eq = 2, newdata = newdata, type = "link"), x$margins[2]))


if( !is.null(x$X3) ){

	sigma <- c(esp.tr(predict(x, eq = 3, newdata = newdata), x$margins[2])$vrb)
	theta <- c(teta.tr(x$VC, predict(x, eq = 4, newdata = newdata))$teta)

                     }

if( is.null(x$X3) ){

   	sigma <- x$sigma2
   	theta <- x$theta 

                   }


	res.m <- distrExIntegrate(ordcont.m, lower = lB, upper = uB, pk = pk, pk1 = pk1, eta.mu = eta2, sigma = sigma, nu = NULL, theta = theta, 
	                          margin = margin, BivD = x$BivD, par2 = x$dof, min.dn = x$VC$min.dn, 
	                          min.pr = x$VC$min.pr, max.pr = x$VC$max.pr, y1 = y1)            

if(fun == "variance"){

        res.v <- distrExIntegrate(ordcont.v, lower = lB, upper = uB, pk = pk, pk1 = pk1, eta.mu = eta2, sigma = sigma, nu = NULL, theta = theta, 
	                          margin = margin, BivD = x$BivD, par2 = x$dof, min.dn = x$VC$min.dn, 
	                          min.pr = x$VC$min.pr, max.pr = x$VC$max.pr, y1 = y1)
	                          
        res.v <- res.v - res.m^2 
        
                      }



######## intervals 




 
for(ii in 1:n.sim){

xx <- x # must re-use always original estimated coef vector

cut.sim <- xx$coefficients[1 : (xx$VC$K1 - 1)]
	cut.sim.ti <- rep(0, xx$VC$K1 - 1)
	cut.sim.ti[1] <- cut.sim[1] ; for(i in 2 : (xx$VC$K1 - 1)) {cut.sim.ti[i] <- sqrt(cut.sim[i] - cut.sim[i - 1])}

coefficients.s <- xx$coefficients
	coefficients.s[1 : (xx$VC$K1 - 1)] <- cut.sim.ti


bs <- rMVN(1, mean = coefficients.s, sigma = xx$Vb)  


# ... and then transformed back

cut.sim.ti <- bs[, 1 : (xx$VC$K1 - 1)]
	cut.sim <- matrix(nrow = 1, ncol = xx$VC$K1 - 1, 0)	
	cut.sim[, 1] <- cut.sim.ti[1] ; for (i in 2 : (xx$VC$K1 - 1)) cut.sim[, i] <- cut.sim[, i - 1] + cut.sim.ti[i]^2

bs[, 1 : (xx$VC$K1 - 1)] <- cut.sim


xx$coefficients <- coef.s <- c(bs) 

pk1 <- 0

pkt <- predict(xx, eq = 1, newdata = newdata, type = "response")$p1.cum 
pk  <- pkt[, y1]
if(y1 > 1) pk1 <- pkt[, y1 - 1]
                                                    
eta2 <- c(eta.tr(predict(xx, eq = 2, newdata = newdata, type = "link"), xx$margins[2]))


if( !is.null(xx$X3) ){

	sigma <- c(esp.tr(predict(xx, eq = 3, newdata = newdata), xx$margins[2])$vrb)
	theta <- c(teta.tr(xx$VC, predict(xx, eq = 4, newdata = newdata))$teta)

                     }

if( is.null(xx$X3) ){

   	sigma <- c(esp.tr(coef.s[length(coef.s)-1], xx$margins[2])$vrb)
   	theta <- c(teta.tr(xx$VC, coef.s[length(coef.s)])$teta)

                   }



res.mCITemp <- suppressMessages(try(distrExIntegrate(ordcont.m, lower = lB, upper = uB, pk = pk, pk1 = pk1, eta.mu = eta2, sigma = sigma, nu = NULL, theta = theta, 
	                          margin = margin, BivD = x$BivD, par2 = x$dof, min.dn = x$VC$min.dn, 
	                          min.pr = x$VC$min.pr, max.pr = x$VC$max.pr, y1 = y1), silent = TRUE))

   if( is.numeric(res.mCITemp) == FALSE ) next
   res.mCI[ii] <- as.numeric(res.mCITemp)
                                     
if(fun == "variance"){

res.vCITemp <- suppressMessages(try(distrExIntegrate(ordcont.v, lower = lB, upper = uB, pk = pk, pk1 = pk1, eta.mu = eta2, sigma = sigma, nu = NULL, theta = theta, 
	                          margin = margin, BivD = x$BivD, par2 = x$dof, min.dn = x$VC$min.dn, 
	                          min.pr = x$VC$min.pr, max.pr = x$VC$max.pr, y1 = y1), silent = TRUE) )   
                                     
   if( is.numeric(res.vCITemp) == FALSE ) next
   res.vCI[ii] <- as.numeric(res.vCITemp)

res.vCI[ii] <- res.vCI[ii] - res.mCI[ii]^2

                    }

} 







rm(xx)


} # ord-cont


#######


if( x$margins[1] %in% c(m1d, m2d) && x$margins[2] %in% c(m2, m3) ){ ### discr 1/2 - cont 2/3 ### 


if(eq == 1){ 
    if( is.null(y2) )                                     stop("A value for y2 must be provided.")
    if( !(y2 >= range(x$y2)[1] && y2 <= range(x$y2)[2]) ) stop("A value for y2, in its observed range, must be provided.")
    y12 <- y2 
    } 

if(eq == 2){ 
    if( is.null(y1) )                                     stop("A value for y1 must be provided.")
    if( !(y1 >= range(x$y1)[1] && y1 <= range(x$y1)[2]) ) stop("A value for y1, in its observed range, must be provided.")
    y12 <- y1  
    }




 fit1  <- c(eta.tr(predict(x, eq = 1, newdata = newdata, type = "link"), x$margins[1]))    
 fit2  <- c(eta.tr(predict(x, eq = 2, newdata = newdata, type = "link"), x$margins[2]))
 
 Xfit1 <- predict(x, eq = 1, newdata = newdata, type = "lpmatrix")
 Xfit2 <- predict(x, eq = 2, newdata = newdata, type = "lpmatrix")
 
 
 bs <- rMVN(n.sim, mean = x$coefficients, sigma = x$Vb)
 
 fit1S <- Xfit1%*%t(bs[, 1:x$X1.d2]) 
 fit2S <- Xfit2%*%t(bs[, (x$X1.d2 + 1):(x$X1.d2 + x$X2.d2)]) 
 

#######################


 if(!is.null(x$X3)){   
 
 
         if( x$margins[1] %in% m1d && x$margins[2] %in% m2){
  
                      fit3  <- predict(x, eq = 3, newdata = newdata, type = "link")
                      fit4  <- predict(x, eq = 4, newdata = newdata, type = "link")
                      Xfit3 <- predict(x, eq = 3, newdata = newdata, type = "lpmatrix")
                      Xfit4 <- predict(x, eq = 4, newdata = newdata, type = "lpmatrix")  
                      
                      fit3S <- Xfit3%*%t(bs[, (x$X1.d2 + x$X2.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2)])
                      fit4S <- Xfit4%*%t(bs[, (x$X1.d2 + x$X2.d2 + x$X3.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2)])
                      
                      sig2  <- c(esp.tr(fit3, margin)$vrb)                        
                      thet  <- c(teta.tr(x$VC, fit4)$teta)
 
                      sig2S  <- c(esp.tr(fit3S, margin)$vrb)                        
                      thetS  <- c(teta.tr(x$VC, fit4S)$teta)                     
 
 
                                                           }
 
 
        if( x$margins[1] %in% m1d && x$margins[2] %in% m3){ 
 
                     fit3  <- predict(x, eq = 3, newdata = newdata, type = "link")
                     fit4  <- predict(x, eq = 4, newdata = newdata, type = "link")
                     fit5  <- predict(x, eq = 5, newdata = newdata, type = "link")
                     Xfit3 <- predict(x, eq = 3, newdata = newdata, type = "lpmatrix")
                     Xfit4 <- predict(x, eq = 4, newdata = newdata, type = "lpmatrix")  
                     Xfit5 <- predict(x, eq = 5, newdata = newdata, type = "lpmatrix")    

                     fit3S <- Xfit3%*%t(bs[, (x$X1.d2 + x$X2.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2)])
                     fit4S <- Xfit4%*%t(bs[, (x$X1.d2 + x$X2.d2 + x$X3.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2)])
                     fit5S <- Xfit5%*%t(bs[, (x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2)])

                     sig2  <- c(esp.tr(fit3, margin)$vrb)                        
                     nu2   <- c(enu.tr(fit4, margin)$vrb)
                     thet  <- c(teta.tr(x$VC, fit5)$teta)
                     
                     sig2S  <- c(esp.tr(fit3S, margin)$vrb)                        
                     nu2S   <- c(enu.tr(fit4S, margin)$vrb)
                     thetS  <- c(teta.tr(x$VC, fit5S)$teta)                     
                                                             }  




                                                                                            

        if( x$margins[1] %in% m2d && x$margins[2] %in% m2){
 
                     fit3  <- predict(x, eq = 3, newdata = newdata, type = "link")
                     fit4  <- predict(x, eq = 4, newdata = newdata, type = "link")
                     fit5  <- predict(x, eq = 5, newdata = newdata, type = "link")
                     Xfit3 <- predict(x, eq = 3, newdata = newdata, type = "lpmatrix")
                     Xfit4 <- predict(x, eq = 4, newdata = newdata, type = "lpmatrix")  
                     Xfit5 <- predict(x, eq = 5, newdata = newdata, type = "lpmatrix")   
                     
                     fit3S <- Xfit3%*%t(bs[, (x$X1.d2 + x$X2.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2)])
                     fit4S <- Xfit4%*%t(bs[, (x$X1.d2 + x$X2.d2 + x$X3.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2)])
                     fit5S <- Xfit5%*%t(bs[, (x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2)])
                     
                     sig1  <- c(esp.tr(fit3, margin)$vrb)
                     sig2  <- c(esp.tr(fit4, margin)$vrb)                        
                     thet  <- c(teta.tr(x$VC, fit5)$teta)

                     sig1S  <- c(esp.tr(fit3S, margin)$vrb)
                     sig2S  <- c(esp.tr(fit4S, margin)$vrb)                        
                     thetS  <- c(teta.tr(x$VC, fit5S)$teta)                     


                                                           }
 
                                                                                            
        if( x$margins[1] %in% m2d && x$margins[2] %in% m3){ 
 
                     fit3  <- predict(x, eq = 3, newdata = newdata, type = "link")
                     fit4  <- predict(x, eq = 4, newdata = newdata, type = "link")
                     fit5  <- predict(x, eq = 5, newdata = newdata, type = "link")
                     fit6  <- predict(x, eq = 6, newdata = newdata, type = "link")
                     Xfit3 <- predict(x, eq = 3, newdata = newdata, type = "lpmatrix")
                     Xfit4 <- predict(x, eq = 4, newdata = newdata, type = "lpmatrix")  
                     Xfit5 <- predict(x, eq = 5, newdata = newdata, type = "lpmatrix")    
                     Xfit6 <- predict(x, eq = 6, newdata = newdata, type = "lpmatrix")   

                     fit3S <- Xfit3%*%t(bs[, (x$X1.d2 + x$X2.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2)])
                     fit4S <- Xfit4%*%t(bs[, (x$X1.d2 + x$X2.d2 + x$X3.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2)])
                     fit5S <- Xfit5%*%t(bs[, (x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2)])
                     fit6S <- Xfit6%*%t(bs[, (x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2 + x$X6.d2)])

                     sig1  <- c(esp.tr(fit3, margin)$vrb)
                     sig2  <- c(esp.tr(fit4, margin)$vrb)                        
                     nu2   <- c(enu.tr(fit5, margin)$vrb)
                     thet  <- c(teta.tr(x$VC, fit6)$teta)
                     
                     sig1S  <- c(esp.tr(fit3S, margin)$vrb)
                     sig2S  <- c(esp.tr(fit4S, margin)$vrb)                        
                     nu2S   <- c(enu.tr(fit5S, margin)$vrb)
                     thetS  <- c(teta.tr(x$VC, fit6S)$teta)                     
                                                             }                                                                                            
 
 



 
 
                   }  

 
 
 if(is.null(x$X3)){ 
 
                    sig1 <- x$sigma1
                    sig2 <- x$sigma2
                    nu1  <- x$nu1
                    nu2  <- x$nu2                    
                    dof  <- x$dof
                    thet <- x$theta   



        if( x$margins[1] %in% m1d && x$margins[2] %in% m2){ 
 
                     fit3S <- bs[, x$X1.d2 + x$X2.d2 + 1]
                     fit4S <- bs[, x$X1.d2 + x$X2.d2 + 1 + 1]
                     
                     sig2S  <- c(esp.tr(fit3S, margin)$vrb)                        
                     thetS  <- c(teta.tr(x$VC, fit4S)$teta) 
                     
                                                           }                     



        if( x$margins[1] %in% m1d && x$margins[2] %in% m3){ 
 
                     fit3S <- bs[, x$X1.d2 + x$X2.d2 + 1]
                     fit4S <- bs[, x$X1.d2 + x$X2.d2 + 1 + 1]
                     fit5S <- bs[, x$X1.d2 + x$X2.d2 + 1 + 1 + 1]
                     
                     sig2S  <- c(esp.tr(fit3S, margin)$vrb)                        
                     nu2S   <- c(enu.tr(fit4S, margin)$vrb)
                     thetS  <- c(teta.tr(x$VC, fit5S)$teta) 
                     
                                                           } 


        if( x$margins[1] %in% m2d && x$margins[2] %in% m2){ 
 
                     fit3S <- bs[, x$X1.d2 + x$X2.d2 + 1]
                     fit4S <- bs[, x$X1.d2 + x$X2.d2 + 1 + 1]
                     fit5S <- bs[, x$X1.d2 + x$X2.d2 + 1 + 1 + 1]
                     
                     sig1S  <- c(esp.tr(fit3S, margin)$vrb)
                     sig2S  <- c(esp.tr(fit4S, margin)$vrb)                        
                     thetS  <- c(teta.tr(x$VC, fit5S)$teta) 
                     
                                                           }                                                                                             

                                                                                            
 
                                                                                            
        if( x$margins[1] %in% m2d && x$margins[2] %in% m3){ 
 
                     fit3S <- bs[, x$X1.d2 + x$X2.d2 + 1]
                     fit4S <- bs[, x$X1.d2 + x$X2.d2 + 1 + 1]
                     fit5S <- bs[, x$X1.d2 + x$X2.d2 + 1 + 1 + 1]
                     fit6S <- bs[, x$X1.d2 + x$X2.d2 + 1 + 1 + 1 + 1]
                     
                     sig1S  <- c(esp.tr(fit3S, margin)$vrb)
                     sig2S  <- c(esp.tr(fit4S, margin)$vrb)                        
                     nu2S   <- c(enu.tr(fit5S, margin)$vrb)
                     thetS  <- c(teta.tr(x$VC, fit6S)$teta) 
                     
                                                           }                                                                                            
 
 

 




                   } 


#######################


if(eq == 2){ ##### equation of interest is no. 2 so contin. variable ####


	res.m <- distrExIntegrate(discrcont.m, lower = lB, upper = uB, y12 = y12, eta1 = fit1, eta2 = fit2, sig1 = sig1, sig2 = sig2, nu1 = nu1, 
	                                nu2 = nu2, theta = thet, margins = x$margins, BivD = x$BivD, par2 = dof, min.dn = x$VC$min.dn, 
	                                min.pr = x$VC$min.pr, max.pr = x$VC$max.pr, left.trunc1 = lt1)            

if(fun == "variance"){

        res.v <- distrExIntegrate(discrcont.v, lower = lB, upper = uB, y12 = y12, eta1 = fit1, eta2 = fit2, sig1 = sig1, sig2 = sig2, nu1 = nu1, 
	                                nu2 = nu2, theta = thet, margins = x$margins, BivD = x$BivD, par2 = dof, min.dn = x$VC$min.dn, 
	                                min.pr = x$VC$min.pr, max.pr = x$VC$max.pr, left.trunc1 = lt1) 
        res.v <- res.v - res.m^2 
        
                      }



######## intervals 

for(i in 1:n.sim){

res.mCITemp <- suppressMessages(try(distrExIntegrate(discrcont.m, lower = lB, upper = uB, y12 = y12, eta1 = fit1S[i], eta2 = fit2S[i], sig1 = sig1S[i], sig2 = sig2S[i], nu1 = nu1S[i], 
	                                nu2 = nu2S[i], theta = thetS[i], margins = x$margins, BivD = x$BivD, par2 = dof, min.dn = x$VC$min.dn, 
	                                min.pr = x$VC$min.pr, max.pr = x$VC$max.pr, left.trunc1 = lt1), silent = TRUE))

   if( is.numeric(res.mCITemp) == FALSE ) next 
   res.mCI[i] <- as.numeric(res.mCITemp)


                                     
if(fun == "variance"){

res.vCITemp <- suppressMessages(try(distrExIntegrate(discrcont.v, lower = lB, upper = uB, y12 = y12, eta1 = fit1S[i], eta2 = fit2S[i], sig1 = sig1S[i], sig2 = sig2S[i], nu1 = nu1S[i], 
	                                nu2 = nu2S[i], theta = thetS[i], margins = x$margins, BivD = x$BivD, par2 = dof, min.dn = x$VC$min.dn, 
	                                min.pr = x$VC$min.pr, max.pr = x$VC$max.pr, left.trunc1 = lt1), silent = TRUE))


   if( is.numeric(res.vCITemp) == FALSE ) next
   res.vCI[i] <- as.numeric(res.vCITemp)



res.vCI[i] <- res.vCI[i] - res.mCI[i]^2

                    }

                 }


} #### end eq = 2 contin. var  #####







##########################################################################


if(eq == 1){ ##### EQUAT of interest is no. 1 so count variable ####



eps <- 1e-05                  
m.y <- max(c(x$y1, x$y2))*5 # arbitrary, generous enough but might be changed in future


res.m <- Stemp <- NA
y <- i <- tl <- 1 

while ( tl > 0  ){

res.m[i] <- tl <- dicont.m(y, y12 = y12, eta1 = fit1, eta2 = fit2, sig1 = sig1, sig2 = sig2, nu1 = nu1, 
	                     nu2 = nu2, theta = thet, margins = x$margins, BivD = x$BivD, par2 = dof, min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, 
	                     max.pr = x$VC$max.pr, nC = nC, left.trunc1 = lt1)           
Stemp[i] <- sum(res.m)
if(i > 1){ if( (Stemp[i] - Stemp[i - 1])/Stemp[i - 1]*100 < eps ) break }
if(y > m.y) break
y <- y + 1; i <- i + 1

                   }

res.m <- sum(res.m)


           

if(fun == "variance"){

res.v <- Stemp <- NA
y <- i <- tl <- 1

	while ( tl > 0  ){
		res.v[i] <- tl <- dicont.v(y, y12 = y12, eta1 = fit1, eta2 = fit2, sig1 = sig1, sig2 = sig2, nu1 = nu1, 
	                     nu2 = nu2, theta = thet, margins = x$margins, BivD = x$BivD, par2 = dof, min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, 
	                     max.pr = x$VC$max.pr, nC = nC, left.trunc1 = lt1)           

                Stemp[i] <- sum(res.v)
                if(i > 1){ if( (Stemp[i] - Stemp[i - 1])/Stemp[i - 1]*100 < eps ) break}		
		if(y > m.y) break
		y <- y + 1; i <- i + 1		
                           }
                    
        res.v <- sum(res.v)
        res.v <- res.v - res.m^2 
        
}



######## intervals 

for(i in 1:n.sim){ # INTERVALS


res.mCITemp <- Stemp <- NA
y <- j <- tl <- 1 



while ( tl > 0 ){

res.mCITemp[j] <- tl <- dicont.m(y, y12 = y12, eta1 = fit1S[i], eta2 = fit2S[i], sig1 = sig1S[i], sig2 = sig2S[i], nu1 = nu1S[i], 
	                     nu2 = nu2S[i], theta = thetS[i], margins = x$margins, BivD = x$BivD, par2 = dof, min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, 
	                     max.pr = x$VC$max.pr, nC = nC, left.trunc1 = lt1)           

Stemp[j] <- sum(res.mCITemp)
if(j > 1){ if( (Stemp[j] - Stemp[j - 1])/Stemp[j - 1]*100 < eps ) break}
if(y > m.y) break
y <- y + 1; j <- j + 1

                    }

res.mCI[i] <- sum(res.mCITemp)



 

if(fun == "variance"){

res.vCITemp <- Stemp <- NA
y <- j <- tl <- 1 


while ( tl > eps  ){

res.vCITemp[j] <- tl <- dicont.v(y, y12 = y12, eta1 = fit1S[i], eta2 = fit2S[i], sig1 = sig1S[i], sig2 = sig2S[i], nu1 = nu1S[i], 
	                     nu2 = nu2S[i], theta = thetS[i], margins = x$margins, BivD = x$BivD, par2 = dof, min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, 
	                     max.pr = x$VC$max.pr, nC = nC, left.trunc1 = lt1)          

Stemp[j] <- sum(res.vCITemp)
if(j > 1){ if( (Stemp[j] - Stemp[j - 1])/Stemp[j - 1]*100 < eps ) break}
if(y > m.y) break
y <- y + 1; j <- j + 1
                    }

res.vCI[i] <- sum(res.vCITemp)

res.vCI[i] <- res.vCI[i] - res.mCI[i]^2


       
}



} # INTERVALS








} #### end eq = 1 cont var  #####























} ### DISCR - cont ### 







########################################################################################################################################





if( x$margins[1] %in% c(x$VC$bl) && x$margins[2] %in% c(x$VC$m1d, x$VC$m2d) ){ ### bin - DISCR ### 

 if( is.null(y1) )       stop("A value (0 or 1) for y1 must be provided.")
 if( !(y1 %in% c(0,1)) ) stop("y1 can be either 0 or 1.")

 eq <- 2
 nC <- x$nC
 eps <- 1e-05                  
 m.y <- max(x$y2)*5 # arbitrary, generous enough but might be changed in future
 

 fit0 <- predict(x, eq = eq - 1, newdata = newdata, type = "link")           # binary eq.
 fit1 <- predict(x, eq = eq,     newdata = newdata, type = "link") 
 
 
 if(!is.null(x$X3)){                     
                                          fit2 <- predict(x, eq = eq + 1, newdata = newdata, type = "link")
                     if( margin %in% m2d) fit3 <- predict(x, eq = eq + 2, newdata = newdata, type = "link")

                     if( margin %in% m1d) theta <- c(teta.tr(x$VC, fit2)$teta)
                     if( margin %in% m2d){sigma <- c(esp.tr(fit2, margin)$vrb); theta <- c(teta.tr(x$VC, fit3)$teta) }
 
                   }  


 if( is.null(x$X3)){ sigma <- x$sigma; theta <- x$theta } 


 eta.mu <- c(eta.tr(fit1, margin)) 
 p0 <- 1 - probm(fit0, margin = x$margins[eq - 1], only.pr = TRUE, bc = FALSE, min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, max.pr = x$VC$max.pr)$pr 
 


res.m <- Stemp <- NA
y <- i <- tl <- 1 

while ( tl > 0  ){

res.m[i] <- tl <- bindisc.m(y1 = y1, y12 = y, eta2 = fit1, sig2 = sigma, nu2 = NULL, theta = theta, margins = x$margins, BivD = x$BivD, par2 = x$dof, 
                            min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, max.pr = x$VC$max.pr, nC = nC, p0 = p0, left.trunc2 = lt2)
         
Stemp[i] <- sum(res.m)
if(i > 1){ if( (Stemp[i] - Stemp[i - 1])/Stemp[i - 1]*100 < eps ) break }
if(y > m.y) break
y <- y + 1; i <- i + 1

                   }

res.m <- sum(res.m)


           

if(fun == "variance"){

res.v <- Stemp <- NA
y <- i <- tl <- 1

	while ( tl > 0  ){
		res.v[i] <- tl <- bindisc.v(y1 = y1, y12 = y, eta2 = fit1, sig2 = sigma, nu2 = NULL, theta = theta, margins = x$margins, BivD = x$BivD, par2 = x$dof, 
                                            min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, max.pr = x$VC$max.pr, nC = nC, p0 = p0, left.trunc2 = lt2)         

                Stemp[i] <- sum(res.v)
                if(i > 1){ if( (Stemp[i] - Stemp[i - 1])/Stemp[i - 1]*100 < eps ) break}		
		if(y > m.y) break
		y <- y + 1; i <- i + 1		
                           }
                    
        res.v <- sum(res.v)
        res.v <- res.v - res.m^2 
        
}




######## intervals 

 bs <- rMVN(n.sim, mean = x$coefficients, sigma = x$Vb)

 			     Xlp0 <- predict(x, eq = eq - 1, newdata = newdata, type = "lpmatrix"); eta0   <- Xlp0 %*% t(bs[, 1:x$X1.d2]) 
 	                     Xlp1 <- predict(x, eq = eq,     newdata = newdata, type = "lpmatrix"); eta.mu <- Xlp1 %*% t(bs[, (x$X1.d2 + 1):(x$X1.d2 + x$X2.d2)]) 

 
 
 if(!is.null(x$X3)){ 
 	                     Xlp2 <- predict(x, eq = eq + 1, newdata = newdata, type = "lpmatrix"); eta.2 <- Xlp2 %*% t(bs[, (x$X1.d2 + x$X2.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2)]) 
      if( margin %in% m2d){   Xlp3 <- predict(x, eq = eq + 2, newdata = newdata, type = "lpmatrix"); eta.3 <- Xlp3 %*% t(bs[, (x$X1.d2 + x$X2.d2 + x$X3.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2)])}
        
                   }  
                         
 if( is.null(x$X3)){ 
 	                     eta.2 <- bs[, x$X1.d2 + x$X2.d2 + 1] 
       if( margin %in% m2d)  eta.3 <- bs[, x$X1.d2 + x$X2.d2 + 1 + 1]
        
                   }                     



p0     <- 1 - probm(eta0, margin = x$margins[eq - 1], only.pr = TRUE, bc = FALSE, min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, max.pr = x$VC$max.pr)$pr
eta.mu <- eta.tr(eta.mu, margin)


if( margin %in% m1d) theta <- teta.tr(x$VC, eta.2)$teta 
if( margin %in% m2d){sigma <- c(esp.tr(eta.2, margin)$vrb); theta <- c(teta.tr(x$VC, eta.3)$teta) }


 
 
 
 for(i in 1:n.sim){ # INTERVALS
 
 
 res.mCITemp <- Stemp <- NA
 y <- j <- tl <- 1 
 
 
 
 while ( tl > 0 ){
 
 res.mCITemp[j] <- tl <- bindisc.m(y1 = y1, y12 = y, eta2 = eta.mu[i], sig2 = sigma[i], nu2 = NULL, theta = theta[i], margins = x$margins, BivD = x$BivD, par2 = x$dof, 
                                   min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, max.pr = x$VC$max.pr, nC = nC, p0 = p0[i], left.trunc2 = lt2)           
 
 Stemp[j] <- sum(res.mCITemp)
 if(j > 1){ if( (Stemp[j] - Stemp[j - 1])/Stemp[j - 1]*100 < eps ) break}
 if(y > m.y) break
 y <- y + 1; j <- j + 1
 
                     }
 
 res.mCI[i] <- sum(res.mCITemp)
 
 
 
  
 
 if(fun == "variance"){
 
 res.vCITemp <- Stemp <- NA
 y <- j <- tl <- 1 
 
 
 while ( tl > eps  ){
 
 res.vCITemp[j] <- tl <- bindisc.v(y1 = y1, y12 = y, eta2 = eta.mu[i], sig2 = sigma[i], nu2 = NULL, theta = theta[i], margins = x$margins, BivD = x$BivD, par2 = x$dof, 
                                   min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, max.pr = x$VC$max.pr, nC = nC, p0 = p0[i], left.trunc2 = lt2)         
 
 Stemp[j] <- sum(res.vCITemp)
 if(j > 1){ if( (Stemp[j] - Stemp[j - 1])/Stemp[j - 1]*100 < eps ) break}
 if(y > m.y) break
 y <- y + 1; j <- j + 1
                     }
 
 res.vCI[i] <- sum(res.vCITemp)
 
 res.vCI[i] <- res.vCI[i] - res.mCI[i]^2
 
 
        
 }
 
 
 
 } # INTERVALS
 

 
 

} ### bin - DISCR ### 


##############################################################





if(is.null(x$VC$K1) && x$margins[1] %in% c(x$VC$bl) && x$margins[2] %in% c(x$VC$m2, x$VC$m3) ){ ### bin - cont ### 

 if( is.null(y1) )       stop("A value (0 or 1) for y1 must be provided.")
 if( !(y1 %in% c(0,1)) ) stop("y1 can be either 0 or 1.")


 eq <- 2

 fit0 <- predict(x, eq = eq - 1, newdata = newdata, type = "link")           # binary eq.
 fit1 <- predict(x, eq = eq,     newdata = newdata, type = "link") 
 
 
 if(!is.null(x$X3)){                     
                                          fit2 <- predict(x, eq = eq + 1, newdata = newdata, type = "link")
                                          fit3 <- predict(x, eq = eq + 2, newdata = newdata, type = "link")
                     if( margin %in% m3 ) fit4 <- predict(x, eq = eq + 3, newdata = newdata, type = "link")
                     
                                          sigma <- c(esp.tr(fit2, margin)$vrb)
                     if( margin %in% m3 ) nu    <- c(enu.tr(fit3, margin)$vrb)
                     if( margin %in% m2 ) theta <- c(teta.tr(x$VC, fit3)$teta) else theta <- c(teta.tr(x$VC, fit4)$teta)
 
                   }  


 if( is.null(x$X3)){     
                                          sigma <- x$sigma
                                          nu <- x$nu
                                          theta <- x$theta 
                    } 


 eta.mu <- c(eta.tr(fit1, margin)) # no need for original scale here as distrHsAT will do it
 p0 <- 1 - probm(fit0, margin = x$margins[eq - 1], only.pr = TRUE, bc = FALSE, min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, max.pr = x$VC$max.pr)$pr 
 

	res.m <- distrExIntegrate(bincont.m, lower = lB, upper = uB, p0 = p0, eta.mu = eta.mu, sigma = sigma, nu = nu, theta = theta, 
	                                margin = margin, BivD = x$BivD, par2 = x$dof, min.dn = x$VC$min.dn, 
	                                min.pr = x$VC$min.pr, max.pr = x$VC$max.pr, y1 = y1)            

if(fun == "variance"){

        res.v <- distrExIntegrate(bincont.v, lower = lB, upper = uB, p0 = p0, eta.mu = eta.mu, sigma = sigma, nu = nu, theta = theta, 
                                        margin = margin, BivD = x$BivD, par2 = x$dof, min.dn = x$VC$min.dn, 
                                        min.pr = x$VC$min.pr, max.pr = x$VC$max.pr, y1 = y1) 
        res.v <- res.v - res.m^2 
        
                      }



######## intervals 

 bs <- rMVN(n.sim, mean = x$coefficients, sigma = x$Vb)

 			     Xlp0 <- predict(x, eq = eq - 1, newdata = newdata, type = "lpmatrix"); eta0   <- Xlp0 %*% t(bs[, 1:x$X1.d2]) 
 	                     Xlp1 <- predict(x, eq = eq,     newdata = newdata, type = "lpmatrix"); eta.mu <- Xlp1 %*% t(bs[, (x$X1.d2 + 1):(x$X1.d2 + x$X2.d2)]) 

 if(!is.null(x$X3)){ 
 	                     Xlp2 <- predict(x, eq = eq + 1, newdata = newdata, type = "lpmatrix"); eta.si <- Xlp2 %*% t(bs[, (x$X1.d2 + x$X2.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2)]) 
                             Xlp3 <- predict(x, eq = eq + 2, newdata = newdata, type = "lpmatrix"); eta.3  <- Xlp3 %*% t(bs[, (x$X1.d2 + x$X2.d2 + x$X3.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2)])
       if( margin %in% m3 ){ Xlp4 <- predict(x, eq = eq + 3, newdata = newdata, type = "lpmatrix"); eta.4  <- Xlp4 %*% t(bs[, (x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2)])}
        
                   }  
                   
 if( is.null(x$X3)){ 
 	                     eta.si <- bs[, x$X1.d2 + x$X2.d2 + 1] 
                             eta.3  <- bs[, x$X1.d2 + x$X2.d2 + 1 + 1]
       if( margin %in% m3 )  eta.4  <- bs[, x$X1.d2 + x$X2.d2 + 1 + 1 + 1]
        
                   }                     



p0     <- 1 - probm(eta0, margin = x$margins[eq - 1], only.pr = TRUE, bc = FALSE, min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, max.pr = x$VC$max.pr)$pr
eta.mu <- eta.tr(eta.mu, margin)


                     sigma  <- esp.tr(eta.si, margin)$vrb
if( margin %in% m3 ) nu     <- enu.tr(eta.3, margin)$vrb
if( margin %in% m2 ) theta  <- teta.tr(x$VC, eta.3)$teta else theta <- teta.tr(x$VC, eta.4)$teta

 
 
for(i in 1:n.sim){

res.mCITemp <- suppressMessages(try(distrExIntegrate(bincont.m, lower = lB, upper = uB, p0 = p0[i], eta.mu = eta.mu[i], sigma = sigma[i], nu = nu[i], theta = theta[i], 
                                     margin = margin, BivD = x$BivD, par2 = x$dof, min.dn = x$VC$min.dn, 
                                     min.pr = x$VC$min.pr, max.pr = x$VC$max.pr, y1 = y1), silent = TRUE))

   if( is.numeric(res.mCITemp) == FALSE ) next
   res.mCI[i] <- as.numeric(res.mCITemp)
                                     
if(fun == "variance"){

res.vCITemp <- suppressMessages(try(distrExIntegrate(bincont.v, lower = lB, upper = uB, p0 = p0[i], eta.mu = eta.mu[i], sigma = sigma[i], nu = nu[i], theta = theta[i], 
                                     margin = margin, BivD = x$BivD, par2 = x$dof, min.dn = x$VC$min.dn, 
                                     min.pr = x$VC$min.pr, max.pr = x$VC$max.pr, y1 = y1), silent = TRUE) )   
                                     
   if( is.numeric(res.vCITemp) == FALSE ) next
   res.vCI[i] <- as.numeric(res.vCITemp)

res.vCI[i] <- res.vCI[i] - res.mCI[i]^2

                    }

} 


} ### bin - cont ### 


##############################################################
##############################################################

if( x$margins[1] %in% c(x$VC$m2, x$VC$m3) && x$margins[2] %in% c(x$VC$m2, x$VC$m3) ){ ### cont2/3 - cont2/3 ### 


if(eq == 1){ 
    cond <- 2
    if( is.null(y2) )                                     stop("A value for y2 must be provided.")
    if( !(y2 >= range(x$y2)[1] && y2 <= range(x$y2)[2]) ) stop("A value for y2, in its observed range, must be provided.")
    y12 <- y2 
    } 

if(eq == 2){ 
    cond <- 1
    if( is.null(y1) )                                     stop("A value for y1 must be provided.")
    if( !(y1 >= range(x$y1)[1] && y1 <= range(x$y1)[2]) ) stop("A value for y1, in its observed range, must be provided.")
    y12 <- y1  
    }




 fit1  <- c(eta.tr(predict(x, eq = 1, newdata = newdata, type = "link"), x$margins[1]))    
 fit2  <- c(eta.tr(predict(x, eq = 2, newdata = newdata, type = "link"), x$margins[2]))
 
 Xfit1 <- predict(x, eq = 1, newdata = newdata, type = "lpmatrix")
 Xfit2 <- predict(x, eq = 2, newdata = newdata, type = "lpmatrix")
 
 
 bs <- rMVN(n.sim, mean = x$coefficients, sigma = x$Vb)
 
 fit1S <- Xfit1%*%t(bs[, 1:x$X1.d2]) 
 fit2S <- Xfit2%*%t(bs[, (x$X1.d2 + 1):(x$X1.d2 + x$X2.d2)]) 
 

#######################


 if(!is.null(x$X3)){   
 
 
        if( x$margins[1] %in% c(x$VC$m2) && x$margins[2] %in% c(x$VC$m2) &&  x$BivD != "T"){
 
                     fit3  <- predict(x, eq = 3, newdata = newdata, type = "link")
                     fit4  <- predict(x, eq = 4, newdata = newdata, type = "link")
                     fit5  <- predict(x, eq = 5, newdata = newdata, type = "link")
                     Xfit3 <- predict(x, eq = 3, newdata = newdata, type = "lpmatrix")
                     Xfit4 <- predict(x, eq = 4, newdata = newdata, type = "lpmatrix")  
                     Xfit5 <- predict(x, eq = 5, newdata = newdata, type = "lpmatrix")   
                     
                     fit3S <- Xfit3%*%t(bs[, (x$X1.d2 + x$X2.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2)])
                     fit4S <- Xfit4%*%t(bs[, (x$X1.d2 + x$X2.d2 + x$X3.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2)])
                     fit5S <- Xfit5%*%t(bs[, (x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2)])
                     
                     sig1  <- c(esp.tr(fit3, margin)$vrb)
                     sig2  <- c(esp.tr(fit4, margin)$vrb)                        
                     thet  <- c(teta.tr(x$VC, fit5)$teta)

                     sig1S  <- c(esp.tr(fit3S, margin)$vrb)
                     sig2S  <- c(esp.tr(fit4S, margin)$vrb)                        
                     thetS  <- c(teta.tr(x$VC, fit5S)$teta)                     


                                                                                            }

        if( x$margins[1] %in% c(x$VC$m2) && x$margins[2] %in% c(x$VC$m2) &&  x$BivD == "T"){
 
                     fit3  <- predict(x, eq = 3, newdata = newdata, type = "link")
                     fit4  <- predict(x, eq = 4, newdata = newdata, type = "link")
                     fit5  <- predict(x, eq = 5, newdata = newdata, type = "link")
                     fit6  <- predict(x, eq = 6, newdata = newdata, type = "link")
                     Xfit3 <- predict(x, eq = 3, newdata = newdata, type = "lpmatrix")
                     Xfit4 <- predict(x, eq = 4, newdata = newdata, type = "lpmatrix")  
                     Xfit5 <- predict(x, eq = 5, newdata = newdata, type = "lpmatrix")                     
                     Xfit6 <- predict(x, eq = 6, newdata = newdata, type = "lpmatrix")  

                     fit3S <- Xfit3%*%t(bs[, (x$X1.d2 + x$X2.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2)])
                     fit4S <- Xfit4%*%t(bs[, (x$X1.d2 + x$X2.d2 + x$X3.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2)])
                     fit5S <- Xfit5%*%t(bs[, (x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2)])
                     fit6S <- Xfit6%*%t(bs[, (x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2 + x$X6.d2)])

                     sig1  <- c(esp.tr(fit3, margin)$vrb)
                     sig2  <- c(esp.tr(fit4, margin)$vrb)                        
                     dof   <- c(dof.tr(fit5)$vao)
                     thet  <- c(teta.tr(x$VC, fit6)$teta)
                     
                     sig1S  <- c(esp.tr(fit3S, margin)$vrb)
                     sig2S  <- c(esp.tr(fit4S, margin)$vrb)                        
                     dofS   <- c(dof.tr(fit5S)$vao)
                     thetS  <- c(teta.tr(x$VC, fit6S)$teta)                     
                     
                                                                                            } 
                                                                                            
        if( x$margins[1] %in% c(x$VC$m3) && x$margins[2] %in% c(x$VC$m2) &&  x$BivD != "T"){ # as above
 
                     fit3  <- predict(x, eq = 3, newdata = newdata, type = "link")
                     fit4  <- predict(x, eq = 4, newdata = newdata, type = "link")
                     fit5  <- predict(x, eq = 5, newdata = newdata, type = "link")
                     fit6  <- predict(x, eq = 6, newdata = newdata, type = "link")
                     Xfit3 <- predict(x, eq = 3, newdata = newdata, type = "lpmatrix")
                     Xfit4 <- predict(x, eq = 4, newdata = newdata, type = "lpmatrix")  
                     Xfit5 <- predict(x, eq = 5, newdata = newdata, type = "lpmatrix")    
                     Xfit6 <- predict(x, eq = 6, newdata = newdata, type = "lpmatrix")  

                     fit3S <- Xfit3%*%t(bs[, (x$X1.d2 + x$X2.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2)])
                     fit4S <- Xfit4%*%t(bs[, (x$X1.d2 + x$X2.d2 + x$X3.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2)])
                     fit5S <- Xfit5%*%t(bs[, (x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2)])
                     fit6S <- Xfit6%*%t(bs[, (x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2 + x$X6.d2)])

                     sig1  <- c(esp.tr(fit3, margin)$vrb)
                     sig2  <- c(esp.tr(fit4, margin)$vrb)                        
                     nu1   <- c(enu.tr(fit5, margin)$vrb)
                     thet  <- c(teta.tr(x$VC, fit6)$teta)
                     
                     sig1S  <- c(esp.tr(fit3S, margin)$vrb)
                     sig2S  <- c(esp.tr(fit4S, margin)$vrb)                        
                     nu1S   <- c(enu.tr(fit5S, margin)$vrb)
                     thetS  <- c(teta.tr(x$VC, fit6S)$teta)                     
                     
                                                                                            } 
                                                                                            
        if( x$margins[1] %in% c(x$VC$m3) && x$margins[2] %in% c(x$VC$m2) &&  x$BivD == "T"){
 
                     fit3  <- predict(x, eq = 3, newdata = newdata, type = "link")
                     fit4  <- predict(x, eq = 4, newdata = newdata, type = "link")
                     fit5  <- predict(x, eq = 5, newdata = newdata, type = "link")
                     fit6  <- predict(x, eq = 6, newdata = newdata, type = "link")
                     fit7  <- predict(x, eq = 7, newdata = newdata, type = "link")
                     
                     Xfit3 <- predict(x, eq = 3, newdata = newdata, type = "lpmatrix")
                     Xfit4 <- predict(x, eq = 4, newdata = newdata, type = "lpmatrix")  
                     Xfit5 <- predict(x, eq = 5, newdata = newdata, type = "lpmatrix")    
                     Xfit6 <- predict(x, eq = 6, newdata = newdata, type = "lpmatrix")  
                     Xfit7 <- predict(x, eq = 7, newdata = newdata, type = "lpmatrix")   

                     fit3S <- Xfit3%*%t(bs[, (x$X1.d2 + x$X2.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2)])
                     fit4S <- Xfit4%*%t(bs[, (x$X1.d2 + x$X2.d2 + x$X3.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2)])
                     fit5S <- Xfit5%*%t(bs[, (x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2)])
                     fit6S <- Xfit6%*%t(bs[, (x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2 + x$X6.d2)])
                     fit7S <- Xfit7%*%t(bs[, (x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2 + x$X6.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2 + x$X6.d2 + x$X7.d2)])

                     sig1  <- c(esp.tr(fit3, margin)$vrb)
                     sig2  <- c(esp.tr(fit4, margin)$vrb)                        
                     nu1   <- c(enu.tr(fit5, margin)$vrb)
                     dof   <- c(dof.tr(fit6)$vao)
                     thet  <- c(teta.tr(x$VC, fit7)$teta)
                     
                     sig1S  <- c(esp.tr(fit3S, margin)$vrb)
                     sig2S  <- c(esp.tr(fit4S, margin)$vrb)                        
                     nu1S   <- c(enu.tr(fit5S, margin)$vrb)
                     dofS   <- c(dof.tr(fit6S)$vao)
                     thetS  <- c(teta.tr(x$VC, fit7S)$teta)                     
                     
                                                                                            }  
                                                                                            
        if( x$margins[1] %in% c(x$VC$m2) && x$margins[2] %in% c(x$VC$m3) &&  x$BivD != "T"){ 
 
                     fit3  <- predict(x, eq = 3, newdata = newdata, type = "link")
                     fit4  <- predict(x, eq = 4, newdata = newdata, type = "link")
                     fit5  <- predict(x, eq = 5, newdata = newdata, type = "link")
                     fit6  <- predict(x, eq = 6, newdata = newdata, type = "link")
                     Xfit3 <- predict(x, eq = 3, newdata = newdata, type = "lpmatrix")
                     Xfit4 <- predict(x, eq = 4, newdata = newdata, type = "lpmatrix")  
                     Xfit5 <- predict(x, eq = 5, newdata = newdata, type = "lpmatrix")    
                     Xfit6 <- predict(x, eq = 6, newdata = newdata, type = "lpmatrix")   

                     fit3S <- Xfit3%*%t(bs[, (x$X1.d2 + x$X2.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2)])
                     fit4S <- Xfit4%*%t(bs[, (x$X1.d2 + x$X2.d2 + x$X3.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2)])
                     fit5S <- Xfit5%*%t(bs[, (x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2)])
                     fit6S <- Xfit6%*%t(bs[, (x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2 + x$X6.d2)])

                     sig1  <- c(esp.tr(fit3, margin)$vrb)
                     sig2  <- c(esp.tr(fit4, margin)$vrb)                        
                     nu2   <- c(enu.tr(fit5, margin)$vrb)
                     thet  <- c(teta.tr(x$VC, fit6)$teta)
                     
                     sig1S  <- c(esp.tr(fit3S, margin)$vrb)
                     sig2S  <- c(esp.tr(fit4S, margin)$vrb)                        
                     nu2S   <- c(enu.tr(fit5S, margin)$vrb)
                     thetS  <- c(teta.tr(x$VC, fit6S)$teta)                     
                                                                                            }                                                                                            
 
        if( x$margins[1] %in% c(x$VC$m2) && x$margins[2] %in% c(x$VC$m3) &&  x$BivD == "T"){
 
                     fit3  <- predict(x, eq = 3, newdata = newdata, type = "link")
                     fit4  <- predict(x, eq = 4, newdata = newdata, type = "link")
                     fit5  <- predict(x, eq = 5, newdata = newdata, type = "link")
                     fit6  <- predict(x, eq = 6, newdata = newdata, type = "link")
                     fit7  <- predict(x, eq = 7, newdata = newdata, type = "link")
                     
                     Xfit3 <- predict(x, eq = 3, newdata = newdata, type = "lpmatrix")
                     Xfit4 <- predict(x, eq = 4, newdata = newdata, type = "lpmatrix")  
                     Xfit5 <- predict(x, eq = 5, newdata = newdata, type = "lpmatrix")    
                     Xfit6 <- predict(x, eq = 6, newdata = newdata, type = "lpmatrix")  
                     Xfit7 <- predict(x, eq = 7, newdata = newdata, type = "lpmatrix")   
 
                     fit3S <- Xfit3%*%t(bs[, (x$X1.d2 + x$X2.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2)])
                     fit4S <- Xfit4%*%t(bs[, (x$X1.d2 + x$X2.d2 + x$X3.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2)])
                     fit5S <- Xfit5%*%t(bs[, (x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2)])
                     fit6S <- Xfit6%*%t(bs[, (x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2 + x$X6.d2)])
                     fit7S <- Xfit7%*%t(bs[, (x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2 + x$X6.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2 + x$X6.d2 + x$X7.d2)])


                     sig1  <- c(esp.tr(fit3, margin)$vrb)
                     sig2  <- c(esp.tr(fit4, margin)$vrb)    
                     nu2   <- c(enu.tr(fit5, margin)$vrb)
                     dof   <- c(dof.tr(fit6)$vao)
                     thet  <- c(teta.tr(x$VC, fit7)$teta)
                     
                     sig1S  <- c(esp.tr(fit3S, margin)$vrb)
                     sig2S  <- c(esp.tr(fit4S, margin)$vrb)    
                     nu2S   <- c(enu.tr(fit5S, margin)$vrb)
                     dofS   <- c(dof.tr(fit6S)$vao)
                     thetS  <- c(teta.tr(x$VC, fit7S)$teta)                     
                     
                                                                                            }  

        if( x$margins[1] %in% c(x$VC$m3) && x$margins[2] %in% c(x$VC$m3) &&  x$BivD != "T"){
 
                     fit3  <- predict(x, eq = 3, newdata = newdata, type = "link")
                     fit4  <- predict(x, eq = 4, newdata = newdata, type = "link")
                     fit5  <- predict(x, eq = 5, newdata = newdata, type = "link")
                     fit6  <- predict(x, eq = 6, newdata = newdata, type = "link")
                     fit7  <- predict(x, eq = 7, newdata = newdata, type = "link")
                     
                     Xfit3 <- predict(x, eq = 3, newdata = newdata, type = "lpmatrix")
                     Xfit4 <- predict(x, eq = 4, newdata = newdata, type = "lpmatrix")  
                     Xfit5 <- predict(x, eq = 5, newdata = newdata, type = "lpmatrix")    
                     Xfit6 <- predict(x, eq = 6, newdata = newdata, type = "lpmatrix")  
                     Xfit7 <- predict(x, eq = 7, newdata = newdata, type = "lpmatrix")  

                     fit3S <- Xfit3%*%t(bs[, (x$X1.d2 + x$X2.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2)])
                     fit4S <- Xfit4%*%t(bs[, (x$X1.d2 + x$X2.d2 + x$X3.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2)])
                     fit5S <- Xfit5%*%t(bs[, (x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2)])
                     fit6S <- Xfit6%*%t(bs[, (x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2 + x$X6.d2)])
                     fit7S <- Xfit7%*%t(bs[, (x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2 + x$X6.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2 + x$X6.d2 + x$X7.d2)])

                     sig1  <- c(esp.tr(fit3, margin)$vrb)
                     sig2  <- c(esp.tr(fit4, margin)$vrb)                        
                     nu1   <- c(enu.tr(fit5, margin)$vrb)
                     nu2   <- c(enu.tr(fit6, margin)$vrb)
                     thet  <- c(teta.tr(x$VC, fit7)$teta)
                     
                     sig1S  <- c(esp.tr(fit3S, margin)$vrb)
                     sig2S  <- c(esp.tr(fit4S, margin)$vrb)                        
                     nu1S   <- c(enu.tr(fit5S, margin)$vrb)
                     nu2S   <- c(enu.tr(fit6S, margin)$vrb)
                     thetS  <- c(teta.tr(x$VC, fit7S)$teta)                     
                     
                                                                                            } 

        if( x$margins[1] %in% c(x$VC$m3) && x$margins[2] %in% c(x$VC$m3) &&  x$BivD == "T"){
 
                     fit3  <- predict(x, eq = 3, newdata = newdata, type = "link")
                     fit4  <- predict(x, eq = 4, newdata = newdata, type = "link")
                     fit5  <- predict(x, eq = 5, newdata = newdata, type = "link")
                     fit6  <- predict(x, eq = 6, newdata = newdata, type = "link")
                     fit7  <- predict(x, eq = 7, newdata = newdata, type = "link")
                     fit8  <- predict(x, eq = 8, newdata = newdata, type = "link")
                     
                     Xfit3 <- predict(x, eq = 3, newdata = newdata, type = "lpmatrix")
                     Xfit4 <- predict(x, eq = 4, newdata = newdata, type = "lpmatrix")  
                     Xfit5 <- predict(x, eq = 5, newdata = newdata, type = "lpmatrix")    
                     Xfit6 <- predict(x, eq = 6, newdata = newdata, type = "lpmatrix")  
                     Xfit7 <- predict(x, eq = 7, newdata = newdata, type = "lpmatrix") 
                     Xfit8 <- predict(x, eq = 8, newdata = newdata, type = "lpmatrix")    

                     fit3S <- Xfit3%*%t(bs[, (x$X1.d2 + x$X2.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2)])
                     fit4S <- Xfit4%*%t(bs[, (x$X1.d2 + x$X2.d2 + x$X3.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2)])
                     fit5S <- Xfit5%*%t(bs[, (x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2)])
                     fit6S <- Xfit6%*%t(bs[, (x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2 + x$X6.d2)])
                     fit7S <- Xfit7%*%t(bs[, (x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2 + x$X6.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2 + x$X6.d2 + x$X7.d2)])
                     fit8S <- Xfit8%*%t(bs[, (x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2 + x$X6.d2 + x$X7.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2 + x$X6.d2 + x$X7.d2 + x$X8.d2)])


                     sig1  <- c(esp.tr(fit3, margin)$vrb)
                     sig2  <- c(esp.tr(fit4, margin)$vrb)                      
                     nu1   <- c(enu.tr(fit5, margin)$vrb)
                     nu2   <- c(enu.tr(fit6, margin)$vrb)
                     dof   <- c(dof.tr(fit7)$vao)
                     thet  <- c(teta.tr(x$VC, fit8)$teta)
                     
                     sig1S  <- c(esp.tr(fit3S, margin)$vrb)
                     sig2S  <- c(esp.tr(fit4S, margin)$vrb)                      
                     nu1S   <- c(enu.tr(fit5S, margin)$vrb)
                     nu2S   <- c(enu.tr(fit6S, margin)$vrb)
                     dofS   <- c(dof.tr(fit7S)$vao)
                     thetS  <- c(teta.tr(x$VC, fit8S)$teta)                     
                     
                     
                                                                                            }  
 
                   }  

 
 
 if(is.null(x$X3)){ 
 
                    sig1 <- x$sigma1
                    sig2 <- x$sigma2
                    nu1  <- x$nu1
                    nu2  <- x$nu2                    
                    dof  <- x$dof
                    thet <- x$theta   
                    
# 1    2     3     4    5    6         7        8
#mu1, mu2, sig1, sig2, nu1, nu2, dof (for T), theta


        if( x$margins[1] %in% c(x$VC$m2) && x$margins[2] %in% c(x$VC$m2) &&  x$BivD != "T"){
 
                     fit3S <- bs[, x$X1.d2 + x$X2.d2 + 1]
                     fit4S <- bs[, x$X1.d2 + x$X2.d2 + 1 + 1]
                     fit5S <- bs[, x$X1.d2 + x$X2.d2 + 1 + 1 + 1]
                     
                     sig1S  <- c(esp.tr(fit3S, margin)$vrb)
                     sig2S  <- c(esp.tr(fit4S, margin)$vrb)                        
                     thetS  <- c(teta.tr(x$VC, fit5S)$teta)                     

                                                                                            }

        if( x$margins[1] %in% c(x$VC$m2) && x$margins[2] %in% c(x$VC$m2) &&  x$BivD == "T"){
 
                     fit3S <- bs[, x$X1.d2 + x$X2.d2 + 1]
                     fit4S <- bs[, x$X1.d2 + x$X2.d2 + 1 + 1]
                     fit5S <- bs[, x$X1.d2 + x$X2.d2 + 1 + 1 + 1]
                     fit6S <- bs[, x$X1.d2 + x$X2.d2 + 1 + 1 + 1 + 1]
                     
                     sig1S  <- c(esp.tr(fit3S, margin)$vrb)
                     sig2S  <- c(esp.tr(fit4S, margin)$vrb)                        
                     dofS   <- c(dof.tr(fit5S)$vao)
                     thetS  <- c(teta.tr(x$VC, fit6S)$teta)                     
                     
                                                                                            } 
                                                                                            
        if( x$margins[1] %in% c(x$VC$m3) && x$margins[2] %in% c(x$VC$m2) &&  x$BivD != "T"){ 
 
                     fit3S <- bs[, x$X1.d2 + x$X2.d2 + 1]
                     fit4S <- bs[, x$X1.d2 + x$X2.d2 + 1 + 1]
                     fit5S <- bs[, x$X1.d2 + x$X2.d2 + 1 + 1 + 1]
                     fit6S <- bs[, x$X1.d2 + x$X2.d2 + 1 + 1 + 1 + 1]
                     
                     sig1S  <- c(esp.tr(fit3S, margin)$vrb)
                     sig2S  <- c(esp.tr(fit4S, margin)$vrb)                        
                     nu1S   <- c(enu.tr(fit5S, margin)$vrb)
                     thetS  <- c(teta.tr(x$VC, fit6S)$teta)                     
                     
                                                                                            } 
                                                                                            
        if( x$margins[1] %in% c(x$VC$m3) && x$margins[2] %in% c(x$VC$m2) &&  x$BivD == "T"){
 
                     fit3S <- bs[, x$X1.d2 + x$X2.d2 + 1]
                     fit4S <- bs[, x$X1.d2 + x$X2.d2 + 1 + 1]
                     fit5S <- bs[, x$X1.d2 + x$X2.d2 + 1 + 1 + 1]
                     fit6S <- bs[, x$X1.d2 + x$X2.d2 + 1 + 1 + 1 + 1]
                     fit7S <- bs[, x$X1.d2 + x$X2.d2 + 1 + 1 + 1 + 1 + 1]
                     
                     
                     sig1S  <- c(esp.tr(fit3S, margin)$vrb)
                     sig2S  <- c(esp.tr(fit4S, margin)$vrb)                        
                     nu1S   <- c(enu.tr(fit5S, margin)$vrb)
                     dofS   <- c(dof.tr(fit6S)$vao)
                     thetS  <- c(teta.tr(x$VC, fit7S)$teta)                     
                     
                                                                                            }  
                                                                                            
        if( x$margins[1] %in% c(x$VC$m2) && x$margins[2] %in% c(x$VC$m3) &&  x$BivD != "T"){ 
 
                     fit3S <- bs[, x$X1.d2 + x$X2.d2 + 1]
                     fit4S <- bs[, x$X1.d2 + x$X2.d2 + 1 + 1]
                     fit5S <- bs[, x$X1.d2 + x$X2.d2 + 1 + 1 + 1]
                     fit6S <- bs[, x$X1.d2 + x$X2.d2 + 1 + 1 + 1 + 1]
                     
                     sig1S  <- c(esp.tr(fit3S, margin)$vrb)
                     sig2S  <- c(esp.tr(fit4S, margin)$vrb)                        
                     nu2S   <- c(enu.tr(fit5S, margin)$vrb)
                     thetS  <- c(teta.tr(x$VC, fit6S)$teta)                     
                                                                                            }                                                                                            
 
        if( x$margins[1] %in% c(x$VC$m2) && x$margins[2] %in% c(x$VC$m3) &&  x$BivD == "T"){
 
                     fit3S <- bs[, x$X1.d2 + x$X2.d2 + 1]
                     fit4S <- bs[, x$X1.d2 + x$X2.d2 + 1 + 1]
                     fit5S <- bs[, x$X1.d2 + x$X2.d2 + 1 + 1 + 1]
                     fit6S <- bs[, x$X1.d2 + x$X2.d2 + 1 + 1 + 1 + 1]
                     fit7S <- bs[, x$X1.d2 + x$X2.d2 + 1 + 1 + 1 + 1 + 1]
                     
                     sig1S  <- c(esp.tr(fit3S, margin)$vrb)
                     sig2S  <- c(esp.tr(fit4S, margin)$vrb)    
                     nu2S   <- c(enu.tr(fit5S, margin)$vrb)
                     dofS   <- c(dof.tr(fit6S)$vao)
                     thetS  <- c(teta.tr(x$VC, fit7S)$teta)                     
                     
                                                                                            }  

        if( x$margins[1] %in% c(x$VC$m3) && x$margins[2] %in% c(x$VC$m3) &&  x$BivD != "T"){
 
                     fit3S <- bs[, x$X1.d2 + x$X2.d2 + 1]
                     fit4S <- bs[, x$X1.d2 + x$X2.d2 + 1 + 1]
                     fit5S <- bs[, x$X1.d2 + x$X2.d2 + 1 + 1 + 1]
                     fit6S <- bs[, x$X1.d2 + x$X2.d2 + 1 + 1 + 1 + 1]
                     fit7S <- bs[, x$X1.d2 + x$X2.d2 + 1 + 1 + 1 + 1 + 1]
                     
                     sig1S  <- c(esp.tr(fit3S, margin)$vrb)
                     sig2S  <- c(esp.tr(fit4S, margin)$vrb)                        
                     nu1S   <- c(enu.tr(fit5S, margin)$vrb)
                     nu2S   <- c(enu.tr(fit6S, margin)$vrb)
                     thetS  <- c(teta.tr(x$VC, fit7S)$teta)                     
                     
                                                                                            } 

        if( x$margins[1] %in% c(x$VC$m3) && x$margins[2] %in% c(x$VC$m3) &&  x$BivD == "T"){
 

                     fit3S <- bs[, x$X1.d2 + x$X2.d2 + 1]
                     fit4S <- bs[, x$X1.d2 + x$X2.d2 + 1 + 1]
                     fit5S <- bs[, x$X1.d2 + x$X2.d2 + 1 + 1 + 1]
                     fit6S <- bs[, x$X1.d2 + x$X2.d2 + 1 + 1 + 1 + 1]
                     fit7S <- bs[, x$X1.d2 + x$X2.d2 + 1 + 1 + 1 + 1 + 1]
                     fit8S <- bs[, x$X1.d2 + x$X2.d2 + 1 + 1 + 1 + 1 + 1 + 1]
                     
                     sig1S  <- c(esp.tr(fit3S, margin)$vrb)
                     sig2S  <- c(esp.tr(fit4S, margin)$vrb)                      
                     nu1S   <- c(enu.tr(fit5S, margin)$vrb)
                     nu2S   <- c(enu.tr(fit6S, margin)$vrb)
                     dofS   <- c(dof.tr(fit7S)$vao)
                     thetS  <- c(teta.tr(x$VC, fit8S)$teta)                                         
                     
                                                                                            }  


                   } 


#######################




	res.m <- distrExIntegrate(contcont.m, lower = lB, upper = uB, y12 = y12, eta1 = fit1, eta2 = fit2, sig1 = sig1, sig2 = sig2, nu1 = nu1, 
	                                nu2 = nu2, theta = thet, margins = x$margins, BivD = x$BivD, par2 = dof, min.dn = x$VC$min.dn, 
	                                min.pr = x$VC$min.pr, max.pr = x$VC$max.pr, cond = cond)            

if(fun == "variance"){

        res.v <- distrExIntegrate(contcont.v, lower = lB, upper = uB, y12 = y12, eta1 = fit1, eta2 = fit2, sig1 = sig1, sig2 = sig2, nu1 = nu1, 
	                                nu2 = nu2, theta = thet, margins = x$margins, BivD = x$BivD, par2 = dof, min.dn = x$VC$min.dn, 
	                                min.pr = x$VC$min.pr, max.pr = x$VC$max.pr, cond = cond) 
        res.v <- res.v - res.m^2 
        
                      }



######## intervals 

for(i in 1:n.sim){

res.mCITemp <- suppressMessages(try(distrExIntegrate(contcont.m, lower = lB, upper = uB, y12 = y12, eta1 = fit1S[i], eta2 = fit2S[i], sig1 = sig1S[i], sig2 = sig2S[i], nu1 = nu1S[i], 
	                                nu2 = nu2S[i], theta = thetS[i], margins = x$margins, BivD = x$BivD, par2 = dofS[i], min.dn = x$VC$min.dn, 
	                                min.pr = x$VC$min.pr, max.pr = x$VC$max.pr, cond = cond), silent = TRUE))

   if( is.numeric(res.mCITemp) == FALSE ) next 
   res.mCI[i] <- as.numeric(res.mCITemp)


                                     
if(fun == "variance"){

res.vCITemp <- suppressMessages(try(distrExIntegrate(contcont.v, lower = lB, upper = uB, y12 = y12, eta1 = fit1S[i], eta2 = fit2S[i], sig1 = sig1S[i], sig2 = sig2S[i], nu1 = nu1S[i], 
	                                nu2 = nu2S[i], theta = thetS[i], margins = x$margins, BivD = x$BivD, par2 = dofS[i], min.dn = x$VC$min.dn, 
	                                min.pr = x$VC$min.pr, max.pr = x$VC$max.pr, cond = cond), silent = TRUE))


   if( is.numeric(res.vCITemp) == FALSE ) next
   res.vCI[i] <- as.numeric(res.vCITemp)



res.vCI[i] <- res.vCI[i] - res.mCI[i]^2

                    }

                 }






} ### cont - cont ### 









if( x$margins[1] %in% c(m1d, m2d) && x$margins[2] %in% c(m1d, m2d) ){ ### DISCR - DISCR ### 


if(eq == 1){ 
    cond <- 2
    if( is.null(y2) )                                     stop("A value for y2 must be provided.")
    if( !(y2 >= range(x$y2)[1] && y2 <= range(x$y2)[2]) ) stop("A value for y2, in its observed range, must be provided.")
    y12 <- y2 
    
    } 

if(eq == 2){ 
    cond <- 1
    if( is.null(y1) )                                     stop("A value for y1 must be provided.")
    if( !(y1 >= range(x$y1)[1] && y1 <= range(x$y1)[2]) ) stop("A value for y1, in its observed range, must be provided.")
    y12 <- y1  
    }


 fit1  <- c(eta.tr(predict(x, eq = 1, newdata = newdata, type = "link"), x$margins[1]))    
 fit2  <- c(eta.tr(predict(x, eq = 2, newdata = newdata, type = "link"), x$margins[2]))
 
 Xfit1 <- predict(x, eq = 1, newdata = newdata, type = "lpmatrix")
 Xfit2 <- predict(x, eq = 2, newdata = newdata, type = "lpmatrix")
 
 
 bs <- rMVN(n.sim, mean = x$coefficients, sigma = x$Vb)
 
 fit1S <- Xfit1%*%t(bs[, 1:x$X1.d2]) 
 fit2S <- Xfit2%*%t(bs[, (x$X1.d2 + 1):(x$X1.d2 + x$X2.d2)]) 
 
 nC <- x$nC
 

#######################


 if(!is.null(x$X3)){   
 
 
        if( x$margins[1] %in% m2d && x$margins[2] %in% m2d){
 
                     fit3  <- predict(x, eq = 3, newdata = newdata, type = "link")
                     fit4  <- predict(x, eq = 4, newdata = newdata, type = "link")
                     fit5  <- predict(x, eq = 5, newdata = newdata, type = "link")
                     Xfit3 <- predict(x, eq = 3, newdata = newdata, type = "lpmatrix")
                     Xfit4 <- predict(x, eq = 4, newdata = newdata, type = "lpmatrix")  
                     Xfit5 <- predict(x, eq = 5, newdata = newdata, type = "lpmatrix")   
                     
                     fit3S <- Xfit3%*%t(bs[, (x$X1.d2 + x$X2.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2)])
                     fit4S <- Xfit4%*%t(bs[, (x$X1.d2 + x$X2.d2 + x$X3.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2)])
                     fit5S <- Xfit5%*%t(bs[, (x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2)])
                     
                     sig1  <- c(esp.tr(fit3, margin)$vrb)
                     sig2  <- c(esp.tr(fit4, margin)$vrb)                        
                     thet  <- c(teta.tr(x$VC, fit5)$teta)

                     sig1S  <- c(esp.tr(fit3S, margin)$vrb)
                     sig2S  <- c(esp.tr(fit4S, margin)$vrb)                        
                     thetS  <- c(teta.tr(x$VC, fit5S)$teta)                                       
                                                            }
       
       
        if( x$margins[1] %in% m1d && x$margins[2] %in% m2d){
 
                     fit3  <- predict(x, eq = 3, newdata = newdata, type = "link")
                     fit4  <- predict(x, eq = 4, newdata = newdata, type = "link")
                     Xfit3 <- predict(x, eq = 3, newdata = newdata, type = "lpmatrix")
                     Xfit4 <- predict(x, eq = 4, newdata = newdata, type = "lpmatrix")  
                     
                     fit3S <- Xfit3%*%t(bs[, (x$X1.d2 + x$X2.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2)])
                     fit4S <- Xfit4%*%t(bs[, (x$X1.d2 + x$X2.d2 + x$X3.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2)])
                     
                     sig2  <- c(esp.tr(fit3, margin)$vrb)                        
                     thet  <- c(teta.tr(x$VC, fit4)$teta)

                     sig2S  <- c(esp.tr(fit3S, margin)$vrb)                        
                     thetS  <- c(teta.tr(x$VC, fit4S)$teta)                                        
                                                            }  
          
          
        if( x$margins[1] %in% m1d && x$margins[2] %in% m1d){
 
                     fit3  <- predict(x, eq = 3, newdata = newdata, type = "link")
                     Xfit3 <- predict(x, eq = 3, newdata = newdata, type = "lpmatrix")
                     
                     fit3S <- Xfit3%*%t(bs[, (x$X1.d2 + x$X2.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2)])
                     
                     thet  <- c(teta.tr(x$VC, fit3)$teta)

                     thetS  <- c(teta.tr(x$VC, fit3S)$teta)                                       
                                                             }                                                                               


  }  # !is.nullX3

 
 
 if(is.null(x$X3)){ 
 
                    sig1 <- x$sigma1
                    sig2 <- x$sigma2                 
                    dof  <- x$dof
                    thet <- x$theta   
                    
        if( x$margins[1] %in% m2d && x$margins[2] %in% m2d){
 
                     fit3S <- bs[, x$X1.d2 + x$X2.d2 + 1]
                     fit4S <- bs[, x$X1.d2 + x$X2.d2 + 1 + 1]
                     fit5S <- bs[, x$X1.d2 + x$X2.d2 + 1 + 1 + 1]
                     
                     sig1S <- c(esp.tr(fit3S, margin)$vrb)
                     sig2S <- c(esp.tr(fit4S, margin)$vrb)                        
                     thetS <- c(teta.tr(x$VC, fit5S)$teta)                     
                                                            }


                                                                                            
        if( x$margins[1] %in% m1d && x$margins[2] %in% m2d){
 
                     fit3S <- bs[, x$X1.d2 + x$X2.d2 + 1]
                     fit4S <- bs[, x$X1.d2 + x$X2.d2 + 1 + 1]
                     
                     sig2S <- c(esp.tr(fit3S, margin)$vrb)                        
                     thetS <- c(teta.tr(x$VC, fit4S)$teta)                     
                                                           } 
                                                                                            
 
                                                                                            
        if( x$margins[1] %in% m1d && x$margins[2] %in% m1d){
 
                     fit3S <- bs[, x$X1.d2 + x$X2.d2 + 1]                    
                     thetS <- c(teta.tr(x$VC, fit3S)$teta)                     
                                                            }                                                                                            
 
 
 } 


#######################
#thet <- 0 # used for testing
#thetS <- rep(0, n.sim)



eps <- 1e-05                  
m.y <- max(c(x$y1, x$y2))*5 # arbitrary, generous enough but might be changed in future


res.m <- Stemp <- NA
y <- i <- tl <- 1 

while ( tl > 0  ){

res.m[i] <- tl <- discdisc.m(y, y12 = y12, eta1 = fit1, eta2 = fit2, sig1 = sig1, sig2 = sig2, nu1 = NULL, 
	                     nu2 = NULL, theta = thet, margins = x$margins, BivD = x$BivD, par2 = dof, min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, 
	                     max.pr = x$VC$max.pr, cond = cond, nC = nC, left.trunc1 = lt1, left.trunc2 = lt2)           
Stemp[i] <- sum(res.m)
if(i > 1){ if( (Stemp[i] - Stemp[i - 1])/Stemp[i - 1]*100 < eps ) break }
if(y > m.y) break
y <- y + 1; i <- i + 1

                   }

res.m <- sum(res.m)


           

if(fun == "variance"){

res.v <- Stemp <- NA
y <- i <- tl <- 1

	while ( tl > 0  ){
		res.v[i] <- tl <- discdisc.v(y, y12 = y12, eta1 = fit1, eta2 = fit2, sig1 = sig1, sig2 = sig2, nu1 = NULL, 
			                     nu2 = NULL, theta = thet, margins = x$margins, BivD = x$BivD, par2 = dof, min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, 
			                     max.pr = x$VC$max.pr, cond = cond, nC = nC, left.trunc1 = lt1, left.trunc2 = lt2)           

                Stemp[i] <- sum(res.v)
                if(i > 1){ if( (Stemp[i] - Stemp[i - 1])/Stemp[i - 1]*100 < eps ) break}		
		if(y > m.y) break
		y <- y + 1; i <- i + 1		
                           }
                    
        res.v <- sum(res.v)
        res.v <- res.v - res.m^2 
        
}



######## intervals 

for(i in 1:n.sim){ # INTERVALS


res.mCITemp <- Stemp <- NA
y <- j <- tl <- 1 



while ( tl > 0 ){

res.mCITemp[j] <- tl <- discdisc.m(y, y12 = y12, eta1 = fit1S[i], eta2 = fit2S[i], sig1 = sig1S[i], sig2 = sig2S[i], nu1 = NULL, 
	                     nu2 = NULL, theta = thetS[i], margins = x$margins, BivD = x$BivD, par2 = dof, min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, 
	                     max.pr = x$VC$max.pr, cond = cond, nC = nC, left.trunc1 = lt1, left.trunc2 = lt2)           

Stemp[j] <- sum(res.mCITemp)
if(j > 1){ if( (Stemp[j] - Stemp[j - 1])/Stemp[j - 1]*100 < eps ) break}
if(y > m.y) break
y <- y + 1; j <- j + 1

                    }

res.mCI[i] <- sum(res.mCITemp)



 

if(fun == "variance"){

res.vCITemp <- Stemp <- NA
y <- j <- tl <- 1 


while ( tl > eps  ){

res.vCITemp[j] <- tl <- discdisc.v(y, y12 = y12, eta1 = fit1S[i], eta2 = fit2S[i], sig1 = sig1S[i], sig2 = sig2S[i], nu1 = NULL, 
	                     nu2 = NULL, theta = thetS[i], margins = x$margins, BivD = x$BivD, par2 = dof, min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, 
	                     max.pr = x$VC$max.pr, cond = cond, nC = nC, left.trunc1 = lt1, left.trunc2 = lt2)           

Stemp[j] <- sum(res.vCITemp)
if(j > 1){ if( (Stemp[j] - Stemp[j - 1])/Stemp[j - 1]*100 < eps ) break}
if(y > m.y) break
y <- y + 1; j <- j + 1
                    }

res.vCI[i] <- sum(res.vCITemp)

res.vCI[i] <- res.vCI[i] - res.mCI[i]^2


       
}



} # INTERVALS






} ### discr - discr ### 






lbn <- paste(prob.lev/2*100, "%", sep = "")
ubn <- paste((1-(prob.lev/2))*100, "%", sep = "")



###########################################################################################################################################




if(fun == "mean") res.mvCIs <- res.mCI else res.mvCIs <- res.vCI 


res.mvCIss <- as.numeric(quantile(pbS*res.mvCIs, c(prob.lev/2, 1 - prob.lev/2), na.rm = TRUE))

  if(fun == "mean")     res <- c(res.mvCIss[1], pb*res.m, res.mvCIss[2])
  if(fun == "variance") res <- c(res.mvCIss[1], pb*res.v, res.mvCIss[2])
  
  if(fun == "mean")     names(res) <- c(lbn, "mean",     ubn)
  if(fun == "variance") names(res) <- c(lbn, "variance", ubn)



out <- list(res = res, sim.mv = pbS*res.mvCIs, prob.lev = prob.lev, margins = x$margins, fun = fun)
 							 
class(out) <- "cond.mv"

out


}    # end of function
