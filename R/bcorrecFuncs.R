  testStat  <- getFromNamespace("testStat", "mgcv")
  liu2      <- getFromNamespace("liu2", "mgcv")
  
  ldTweedie    <- getFromNamespace("ldTweedie", "mgcv")
  gamlss.etamu <- getFromNamespace("gamlss.etamu", "mgcv")
  gamlss.gH    <- getFromNamespace("gamlss.gH", "mgcv")



transform_bb <- function(Y) {

         m  <- rowSums(Y)
         p  <- (Y[,1] + 0.5) / (m + 1)      # rescaled proportion
         z  <- qlogis(p)                    # logit transform
         return(z)
                                  }


mc.copula.prob <- function(obj, y.vec, newdata, event.fun, joint = TRUE,
                           n.draw = 50000, n.sim = 100, prob.lev = 0.05) {
 
lbn <- paste(prob.lev/2*100, "%", sep = "")
ubn <- paste((1-(prob.lev/2))*100, "%", sep = "") 
 
  simY0 <- rpcc(obj, newdata, n.draw = n.draw, sim.par = FALSE, independent = !joint)
  p.hat <- mean(event.fun(simY0, y.vec))
  
  p.vals <- numeric(n.sim)
  
  for (b in seq_len(n.sim)) {
    simY <- rpcc(obj, newdata, n.draw, sim.par = TRUE, independent = !joint)
    p.vals[b] <- mean(event.fun(simY, y.vec))
  }
  
  ci <- as.numeric(quantile(p.vals, c(prob.lev/2, 1 - prob.lev/2)))
  
  res.p <- data.frame(pr = p.hat,
           lbn = ci[1],
           ubn = ci[2])
     
  names(res.p) <- c("pr", paste(lbn), paste(ubn)) 
  dimnames(res.p)[[1]] <- ""  
  
   
  attr(res.p, "p.vals") <- list(p.vals = p.vals)
  res.p
}





rpcc <- function(obj, newdata, n.draw = 1000, sim.par = FALSE, independent = FALSE){

if(missing(newdata))                               stop("You must provide a data frame.")
if(!missing(newdata)){ 
                       if(!is.data.frame(newdata)) stop("You must provide a data frame.")
                       if(dim(newdata)[1] > 1)     stop("This calculation is supported for single row data frames.")
                     }


if(n.draw < 2) stop("n.draw has to be a value greater than 2.")


eta.nu1 <- eta.nu2 <- eta.nu3 <- 0

#####*************####
#### margin stuff ####
#####*************####

  marg1 <- obj$m1$margins[1]
  marg2 <- obj$m2$margins[1]
  marg3 <- obj$m3$margins[1]
  
  m1.coef <- obj$m1$coefficients
  m2.coef <- obj$m2$coefficients
  m3.coef <- obj$m3$coefficients
  
  if(sim.par == TRUE){
  
  
  m1.coef <-  c(rMVN(1, mean = m1.coef, sigma=obj$m1$Vb)) 
  m2.coef <-  c(rMVN(1, mean = m2.coef, sigma=obj$m2$Vb)) 
  m3.coef <-  c(rMVN(1, mean = m3.coef, sigma=obj$m3$Vb)) 

  }
  
 
  eta.mu1 <- predict.SemiParBIV(obj$m1, eq = 1, newdata = newdata, type = "lpmatrix")%*%m1.coef[1:obj$m1$X1.d2]
  eta.mu2 <- predict.SemiParBIV(obj$m2, eq = 1, newdata = newdata, type = "lpmatrix")%*%m2.coef[1:obj$m2$X1.d2]
  eta.mu3 <- predict.SemiParBIV(obj$m3, eq = 1, newdata = newdata, type = "lpmatrix")%*%m3.coef[1:obj$m3$X1.d2]

  if(is.null(obj$m1$X2)) eta.sigma1 <- m1.coef[obj$m1$X1.d2 + 1] else eta.sigma1 <- predict.SemiParBIV(obj$m1, eq = 2, newdata = newdata, type = "lpmatrix")%*%m1.coef[(obj$m1$X1.d2 + 1):(obj$m1$X1.d2 + obj$m1$X2.d2)] 
  if(is.null(obj$m2$X2)) eta.sigma2 <- m2.coef[obj$m2$X1.d2 + 1] else eta.sigma2 <- predict.SemiParBIV(obj$m2, eq = 2, newdata = newdata, type = "lpmatrix")%*%m2.coef[(obj$m2$X1.d2 + 1):(obj$m2$X1.d2 + obj$m2$X2.d2)] 
  if(is.null(obj$m3$X2)) eta.sigma3 <- m3.coef[obj$m3$X1.d2 + 1] else eta.sigma3 <- predict.SemiParBIV(obj$m3, eq = 2, newdata = newdata, type = "lpmatrix")%*%m3.coef[(obj$m3$X1.d2 + 1):(obj$m3$X1.d2 + obj$m3$X2.d2)] 
 
  
  if(marg1 %in% c(obj$m1$VC$m3)){
  
  	if(is.null(obj$m1$X3)) eta.nu1 <- m1.coef[obj$m1$X1.d2 + 2] else eta.nu1 <- predict.SemiParBIV(obj$m1, eq = 3, newdata = newdata, type = "lpmatrix")%*%m1.coef[(obj$m1$X1.d2 + obj$m1$X2.d2 + 1):(obj$m1$X1.d2 + obj$m1$X2.d2 + obj$m1$X3.d2)] 

                                }
  

  if(marg2 %in% c(obj$m2$VC$m3)){
  
  	if(is.null(obj$m2$X3)) eta.nu2 <- m2.coef[obj$m2$X1.d2 + 2] else eta.nu2 <- predict.SemiParBIV(obj$m2, eq = 3, newdata = newdata, type = "lpmatrix")%*%m2.coef[(obj$m2$X1.d2 + obj$m2$X2.d2 + 1):(obj$m2$X1.d2 + obj$m2$X2.d2 + obj$m2$X3.d2)] 

                                 }
  

  if(marg3 %in% c(obj$m3$VC$m3)){
  
  	if(is.null(obj$m3$X3)) eta.nu3 <- m3.coef[obj$m3$X1.d2 + 2] else eta.nu3 <- predict.SemiParBIV(obj$m3, eq = 3, newdata = newdata, type = "lpmatrix")%*%m3.coef[(obj$m3$X1.d2 + obj$m3$X2.d2 + 1):(obj$m3$X1.d2 + obj$m3$X2.d2 + obj$m3$X3.d2)] 

                                }
  
# part below and several parts above inefficient but clear to work with while developing it   
  
eta.mu1 <- as.numeric(eta.mu1) 
eta.mu2 <- as.numeric(eta.mu2)  
eta.mu3 <- as.numeric(eta.mu3)  

eta.sigma1 <- as.numeric(eta.sigma1)
eta.sigma2 <- as.numeric(eta.sigma2)
eta.sigma3 <- as.numeric(eta.sigma3)

eta.nu1 <- as.numeric(eta.nu1)
eta.nu2 <- as.numeric(eta.nu2)
eta.nu3 <- as.numeric(eta.nu3)

#####*************####
#### margin stuff ####
#####*************####

#if(obj$simplified == FALSE){
#
#y1 <- y.vec[1]; y2 <- y.vec[2]; y3 <- y.vec[3]
#
#
#if(marg1 == "N")     p1 <- pNO(   y1, mu =     eta.mu1,     sigma = esp.tr(eta.sigma1, marg1)$vrb) 
#if(marg1 == "GU")    p1 <- pGU(   y1, mu =     eta.mu1,     sigma = esp.tr(eta.sigma1, marg1)$vrb) 
#if(marg1 == "rGU")   p1 <- pRG(   y1, mu =     eta.mu1,     sigma = esp.tr(eta.sigma1, marg1)$vrb) 
#if(marg1 == "LO")    p1 <- pLO(   y1, mu =     eta.mu1,     sigma = esp.tr(eta.sigma1, marg1)$vrb) 
#if(marg1 == "LN")    p1 <- pLOGNO(y1, mu =     eta.mu1,     sigma = esp.tr(eta.sigma1, marg1)$vrb) 
#if(marg1 == "WEI")   p1 <- pWEI(  y1, mu = exp(eta.mu1),    sigma = esp.tr(eta.sigma1, marg1)$vrb) 
#if(marg1 == "IG")    p1 <- pIG(   y1, mu = exp(eta.mu1),    sigma = esp.tr(eta.sigma1, marg1)$vrb) 
#if(marg1 == "GA")    p1 <- pGA(   y1, mu = exp(eta.mu1),    sigma = esp.tr(eta.sigma1, marg1)$vrb) 
#if(marg1 == "DAGUM") p1 <- pGB2(  y1, mu = exp(eta.mu1),    sigma = esp.tr(eta.sigma1, marg1)$vrb, nu = exp(eta.nu1), tau = 1) 
#if(marg1 == "SM")    p1 <- pGB2(  y1, mu = exp(eta.mu1),    sigma = esp.tr(eta.sigma1, marg1)$vrb, nu = 1 , tau = exp(eta.nu1)) 
#if(marg1 == "BE")    p1 <- pBE(   y1, mu = plogis(eta.mu1), sigma = esp.tr(eta.sigma1, marg1)$vrb) 
#if(marg1 == "FISK")  p1 <- pGB2(  y1, mu = exp(eta.mu1),    sigma = esp.tr(eta.sigma1, marg1)$vrb, nu = 1 , tau = 1)
#
#
#if(marg2 == "N")     p2 <- pNO(   y2, mu =     eta.mu2,     sigma = esp.tr(eta.sigma2, marg2)$vrb) 
#if(marg2 == "GU")    p2 <- pGU(   y2, mu =     eta.mu2,     sigma = esp.tr(eta.sigma2, marg2)$vrb) 
#if(marg2 == "rGU")   p2 <- pRG(   y2, mu =     eta.mu2,     sigma = esp.tr(eta.sigma2, marg2)$vrb) 
#if(marg2 == "LO")    p2 <- pLO(   y2, mu =     eta.mu2,     sigma = esp.tr(eta.sigma2, marg2)$vrb) 
#if(marg2 == "LN")    p2 <- pLOGNO(y2, mu =     eta.mu2,     sigma = esp.tr(eta.sigma2, marg2)$vrb) 
#if(marg2 == "WEI")   p2 <- pWEI(  y2, mu = exp(eta.mu2),    sigma = esp.tr(eta.sigma2, marg2)$vrb) 
#if(marg2 == "IG")    p2 <- pIG(   y2, mu = exp(eta.mu2),    sigma = esp.tr(eta.sigma2, marg2)$vrb) 
#if(marg2 == "GA")    p2 <- pGA(   y2, mu = exp(eta.mu2),    sigma = esp.tr(eta.sigma2, marg2)$vrb) 
#if(marg2 == "DAGUM") p2 <- pGB2(  y2, mu = exp(eta.mu2),    sigma = esp.tr(eta.sigma2, marg2)$vrb, nu = exp(eta.nu2), tau = 1) 
#if(marg2 == "SM")    p2 <- pGB2(  y2, mu = exp(eta.mu2),    sigma = esp.tr(eta.sigma2, marg2)$vrb, nu = 1 , tau = exp(eta.nu2)) 
#if(marg2 == "BE")    p2 <- pBE(   y2, mu = plogis(eta.mu2), sigma = esp.tr(eta.sigma2, marg2)$vrb) 
#if(marg2 == "FISK")  p2 <- pGB2(  y2, mu = exp(eta.mu2),    sigma = esp.tr(eta.sigma2, marg2)$vrb, nu = 1 , tau = 1)
#
#
#if(marg3 == "N")     p3 <- pNO(   y3, mu =     eta.mu3,     sigma = esp.tr(eta.sigma3, marg3)$vrb) 
#if(marg3 == "GU")    p3 <- pGU(   y3, mu =     eta.mu3,     sigma = esp.tr(eta.sigma3, marg3)$vrb) 
#if(marg3 == "rGU")   p3 <- pRG(   y3, mu =     eta.mu3,     sigma = esp.tr(eta.sigma3, marg3)$vrb) 
#if(marg3 == "LO")    p3 <- pLO(   y3, mu =     eta.mu3,     sigma = esp.tr(eta.sigma3, marg3)$vrb) 
#if(marg3 == "LN")    p3 <- pLOGNO(y3, mu =     eta.mu3,     sigma = esp.tr(eta.sigma3, marg3)$vrb) 
#if(marg3 == "WEI")   p3 <- pWEI(  y3, mu = exp(eta.mu3),    sigma = esp.tr(eta.sigma3, marg3)$vrb) 
#if(marg3 == "IG")    p3 <- pIG(   y3, mu = exp(eta.mu3),    sigma = esp.tr(eta.sigma3, marg3)$vrb) 
#if(marg3 == "GA")    p3 <- pGA(   y3, mu = exp(eta.mu3),    sigma = esp.tr(eta.sigma3, marg3)$vrb) 
#if(marg3 == "DAGUM") p3 <- pGB2(  y3, mu = exp(eta.mu3),    sigma = esp.tr(eta.sigma3, marg3)$vrb, nu = exp(eta.nu3), tau = 1) 
#if(marg3 == "SM")    p3 <- pGB2(  y3, mu = exp(eta.mu3),    sigma = esp.tr(eta.sigma3, marg3)$vrb, nu = 1 , tau = exp(eta.nu3)) 
#if(marg3 == "BE")    p3 <- pBE(   y3, mu = plogis(eta.mu3), sigma = esp.tr(eta.sigma3, marg3)$vrb) 
#if(marg3 == "FISK")  p3 <- pGB2(  y3, mu = exp(eta.mu3),    sigma = esp.tr(eta.sigma3, marg3)$vrb, nu = 1 , tau = 1)

#}


if(independent == FALSE){


newdata <- cbind(newdata, .gjrm_p1 = 0.5, .gjrm_p2 = 0.5, .gjrm_p3 = 0.5)


  cop.map <- list(N = 1, T = 2, F = 5, C0 = 3, C180 = 13, C90 = 23, C270 = 33, G0 = 4, G180 = 14, G90 = 24, G270 = 34, J0 = 6, J180 = 16, J90 = 26, J270 = 36)
  cop.sel <- obj$copula.sel
  cop.fam <- NA
  dof.vec <- obj$dof.vec
  fampcc <- parpcc <- matrix(0, 3, 3)
  par2pcc <- matrix(3, 3, 3)
 
 
  if(obj$cond.var == 1){
  
     t1.coef <- obj$b12$coefficients
     t2.coef <- obj$b13$coefficients
     t3.coef <- obj$b23.1$coefficients
     
     if(sim.par == TRUE){
     
     t1.coef <-  c(rMVN(1, mean = t1.coef, sigma=obj$b12$Vb))   
     t2.coef <-  c(rMVN(1, mean = t2.coef, sigma=obj$b13$Vb))   
     t3.coef <-  c(rMVN(1, mean = t3.coef, sigma=obj$b23.1$Vb)) 
   
     }
  
  

    theta12   <- teta.tr(obj$b12$VC,   predict(obj$b12,   newdata = newdata, type = "lpmatrix") %*% t1.coef)$teta
    theta13   <- teta.tr(obj$b13$VC,   predict(obj$b13,   newdata = newdata, type = "lpmatrix") %*% t2.coef)$teta
    theta23.1 <- teta.tr(obj$b23.1$VC, predict(obj$b23.1, newdata = newdata, type = "lpmatrix") %*% t3.coef)$teta

    Theta <- c(theta12, theta13, theta23.1)  
  
                       }
  

  if(obj$cond.var == 2){
  
     t1.coef <- obj$b12$coefficients
     t2.coef <- obj$b23$coefficients
     t3.coef <- obj$b13.2$coefficients
     
     if(sim.par == TRUE){
     
     t1.coef <-  c(rMVN(1, mean = t1.coef, sigma=obj$b12$Vb))   
     t2.coef <-  c(rMVN(1, mean = t2.coef, sigma=obj$b23$Vb))   
     t3.coef <-  c(rMVN(1, mean = t3.coef, sigma=obj$b13.2$Vb)) 
   
     }  
  
    theta12   <- teta.tr(obj$b12$VC,   predict(obj$b12,   newdata = newdata, type = "lpmatrix") %*% t1.coef)$teta
    theta23   <- teta.tr(obj$b23$VC,   predict(obj$b23,   newdata = newdata, type = "lpmatrix") %*% t2.coef)$teta
    theta13.2 <- teta.tr(obj$b13.2$VC, predict(obj$b13.2, newdata = newdata, type = "lpmatrix") %*% t3.coef)$teta


    Theta <- c(theta12, theta23, theta13.2)
  
                        }
  
  
  if(obj$cond.var == 3){
  
       t1.coef <- obj$b13$coefficients
       t2.coef <- obj$b23$coefficients
       t3.coef <- obj$b12.3$coefficients
       
       if(sim.par == TRUE){
       
       t1.coef <-  c(rMVN(1, mean = t1.coef, sigma=obj$b13$Vb))   
       t2.coef <-  c(rMVN(1, mean = t2.coef, sigma=obj$b23$Vb))   
       t3.coef <-  c(rMVN(1, mean = t3.coef, sigma=obj$b12.3$Vb)) 
     
     } 
  
  
    theta13   <- teta.tr(obj$b13$VC,   predict(obj$b13,   newdata = newdata, type = "lpmatrix") %*% t1.coef)$teta
    theta23   <- teta.tr(obj$b23$VC,   predict(obj$b23,   newdata = newdata, type = "lpmatrix") %*% t2.coef)$teta
    theta12.3 <- teta.tr(obj$b12.3$VC, predict(obj$b12.3, newdata = newdata, type = "lpmatrix") %*% t3.coef)$teta

    Theta <- c(theta13, theta23, theta12.3)
  
                       }  
  
  


  for(i in 1:3){

    cop.name <- cop.sel[i]
     if(grepl("C|G|J", cop.name) & grepl("\\d.*\\d", cop.name)) {

         parts <- unlist(regmatches(cop.name, gregexpr("[CGJ]\\d{1,3}", cop.name)))
         
         if(Theta[i] >= 0) cop.name <- parts[1] else cop.name <- parts[2]  
                                                                 }                                                          
    cop.fam[i] <- cop.map[[cop.name]]
    
                }  
  
  
################################  


  if(obj$cond.var == 1){    # b12, b13, b23.1

    matpcc <- matrix(c(2, 0, 0,
                       3, 3, 0,
                       1, 1, 1                                           
                               ), 3, 3, byrow = TRUE)
 
                       } 

 if(obj$cond.var == 2){

    matpcc <- matrix(c(1, 0, 0,
                       3, 3, 0,
                       2, 2, 2
                               ), 3, 3, byrow = TRUE)

                        }

if(obj$cond.var == 3){

    matpcc <- matrix(c(1, 0, 0,
                       2, 2, 0,
                       3, 3, 3
                               ), 3, 3, byrow = TRUE)

                     }


    fampcc <- matrix(c(0,                   0, 0,                         # as an example
                       cop.fam[3],          0, 0,                         # b13.2
                       cop.fam[1], cop.fam[2], 0                          # b12, b23
                               ), 3, 3, byrow = TRUE)

    parpcc <- matrix(c(0,               0, 0,
                       Theta[3],        0, 0,                             # b13.2
                       Theta[1], Theta[2], 0                              #  b12, b23
                               ), 3, 3, byrow = TRUE)

   par2pcc <- matrix(c(0,                   0, 0,
                       dof.vec[3],          0, 0,                         # b13.2
                       dof.vec[1], dof.vec[2], 0                          #  b12, b23
                               ), 3, 3, byrow = TRUE)


  RVM.pcc <- RVineMatrix(matpcc, fampcc, parpcc, par2 = par2pcc)
  
  simU <- RVineSim(n.draw, RVM.pcc)
 
}

if(independent == TRUE) simU <- cbind(runif(n.draw), runif(n.draw), runif(n.draw)) 

  
if(marg1 == "N")     Y1 <- qNO(   simU[, 1], mu =     eta.mu1,     sigma = esp.tr(eta.sigma1, marg1)$vrb) 
if(marg1 == "GU")    Y1 <- qGU(   simU[, 1], mu =     eta.mu1,     sigma = esp.tr(eta.sigma1, marg1)$vrb) 
if(marg1 == "rGU")   Y1 <- qRG(   simU[, 1], mu =     eta.mu1,     sigma = esp.tr(eta.sigma1, marg1)$vrb) 
if(marg1 == "LO")    Y1 <- qLO(   simU[, 1], mu =     eta.mu1,     sigma = esp.tr(eta.sigma1, marg1)$vrb) 
if(marg1 == "LN")    Y1 <- qLOGNO(simU[, 1], mu =     eta.mu1,     sigma = esp.tr(eta.sigma1, marg1)$vrb) 
if(marg1 == "WEI")   Y1 <- qWEI(  simU[, 1], mu = exp(eta.mu1),    sigma = esp.tr(eta.sigma1, marg1)$vrb) 
if(marg1 == "IG")    Y1 <- qIG(   simU[, 1], mu = exp(eta.mu1),    sigma = esp.tr(eta.sigma1, marg1)$vrb) 
if(marg1 == "GA")    Y1 <- qGA(   simU[, 1], mu = exp(eta.mu1),    sigma = esp.tr(eta.sigma1, marg1)$vrb) 
if(marg1 == "DAGUM") Y1 <- qGB2(  simU[, 1], mu = exp(eta.mu1),    sigma = esp.tr(eta.sigma1, marg1)$vrb, nu = exp(eta.nu1), tau = 1) 
if(marg1 == "SM")    Y1 <- qGB2(  simU[, 1], mu = exp(eta.mu1),    sigma = esp.tr(eta.sigma1, marg1)$vrb, nu = 1 , tau = exp(eta.nu1)) 
if(marg1 == "BE")    Y1 <- qBE(   simU[, 1], mu = plogis(eta.mu1), sigma = esp.tr(eta.sigma1, marg1)$vrb) 
if(marg1 == "FISK")  Y1 <- qGB2(  simU[, 1], mu = exp(eta.mu1),    sigma = esp.tr(eta.sigma1, marg1)$vrb, nu = 1 , tau = 1)


if(marg2 == "N")     Y2 <- qNO(   simU[, 2], mu =     eta.mu2,     sigma = esp.tr(eta.sigma2, marg2)$vrb) 
if(marg2 == "GU")    Y2 <- qGU(   simU[, 2], mu =     eta.mu2,     sigma = esp.tr(eta.sigma2, marg2)$vrb) 
if(marg2 == "rGU")   Y2 <- qRG(   simU[, 2], mu =     eta.mu2,     sigma = esp.tr(eta.sigma2, marg2)$vrb) 
if(marg2 == "LO")    Y2 <- qLO(   simU[, 2], mu =     eta.mu2,     sigma = esp.tr(eta.sigma2, marg2)$vrb) 
if(marg2 == "LN")    Y2 <- qLOGNO(simU[, 2], mu =     eta.mu2,     sigma = esp.tr(eta.sigma2, marg2)$vrb) 
if(marg2 == "WEI")   Y2 <- qWEI(  simU[, 2], mu = exp(eta.mu2),    sigma = esp.tr(eta.sigma2, marg2)$vrb) 
if(marg2 == "IG")    Y2 <- qIG(   simU[, 2], mu = exp(eta.mu2),    sigma = esp.tr(eta.sigma2, marg2)$vrb) 
if(marg2 == "GA")    Y2 <- qGA(   simU[, 2], mu = exp(eta.mu2),    sigma = esp.tr(eta.sigma2, marg2)$vrb) 
if(marg2 == "DAGUM") Y2 <- qGB2(  simU[, 2], mu = exp(eta.mu2),    sigma = esp.tr(eta.sigma2, marg2)$vrb, nu = exp(eta.nu2), tau = 1) 
if(marg2 == "SM")    Y2 <- qGB2(  simU[, 2], mu = exp(eta.mu2),    sigma = esp.tr(eta.sigma2, marg2)$vrb, nu = 1 , tau = exp(eta.nu2)) 
if(marg2 == "BE")    Y2 <- qBE(   simU[, 2], mu = plogis(eta.mu2), sigma = esp.tr(eta.sigma2, marg2)$vrb) 
if(marg2 == "FISK")  Y2 <- qGB2(  simU[, 2], mu = exp(eta.mu2),    sigma = esp.tr(eta.sigma2, marg2)$vrb, nu = 1 , tau = 1)


if(marg3 == "N")     Y3 <- qNO(   simU[, 3], mu =     eta.mu3,     sigma = esp.tr(eta.sigma3, marg3)$vrb) 
if(marg3 == "GU")    Y3 <- qGU(   simU[, 3], mu =     eta.mu3,     sigma = esp.tr(eta.sigma3, marg3)$vrb) 
if(marg3 == "rGU")   Y3 <- qRG(   simU[, 3], mu =     eta.mu3,     sigma = esp.tr(eta.sigma3, marg3)$vrb) 
if(marg3 == "LO")    Y3 <- qLO(   simU[, 3], mu =     eta.mu3,     sigma = esp.tr(eta.sigma3, marg3)$vrb) 
if(marg3 == "LN")    Y3 <- qLOGNO(simU[, 3], mu =     eta.mu3,     sigma = esp.tr(eta.sigma3, marg3)$vrb) 
if(marg3 == "WEI")   Y3 <- qWEI(  simU[, 3], mu = exp(eta.mu3),    sigma = esp.tr(eta.sigma3, marg3)$vrb) 
if(marg3 == "IG")    Y3 <- qIG(   simU[, 3], mu = exp(eta.mu3),    sigma = esp.tr(eta.sigma3, marg3)$vrb) 
if(marg3 == "GA")    Y3 <- qGA(   simU[, 3], mu = exp(eta.mu3),    sigma = esp.tr(eta.sigma3, marg3)$vrb) 
if(marg3 == "DAGUM") Y3 <- qGB2(  simU[, 3], mu = exp(eta.mu3),    sigma = esp.tr(eta.sigma3, marg3)$vrb, nu = exp(eta.nu3), tau = 1) 
if(marg3 == "SM")    Y3 <- qGB2(  simU[, 3], mu = exp(eta.mu3),    sigma = esp.tr(eta.sigma3, marg3)$vrb, nu = 1 , tau = exp(eta.nu3)) 
if(marg3 == "BE")    Y3 <- qBE(   simU[, 3], mu = plogis(eta.mu3), sigma = esp.tr(eta.sigma3, marg3)$vrb) 
if(marg3 == "FISK")  Y3 <- qGB2(  simU[, 3], mu = exp(eta.mu3),    sigma = esp.tr(eta.sigma3, marg3)$vrb, nu = 1 , tau = 1)

 
return(data.frame(Y1 = Y1, Y2 = Y2, Y3 = Y3))


}






not_empty <- function(x) !is.null(x) && length(x) > 0 && any(!is.na(x))


jointCDF <- function(obj, p1, p2, newdata, tail = c("lower","lower")){

if(missing(newdata))                               stop("You must provide a data frame.")
if(!missing(newdata)){ 
                       if(!is.data.frame(newdata)) stop("You must provide a data frame.")
                       if(dim(newdata)[1] > 1)     stop("This calculation is supported for single row data frames.")
                     }
 
nC   <- obj$nC
BivD <- obj$BivD

theta <- teta.tr(obj$VC, predict.SemiParBIV(obj, newdata = newdata))$teta

if(obj$BivD %in% obj$BivD2){ 

            if(theta >= 0){ nC <- obj$fit$nC1; BivD <- obj$fit$Cop1 } else{ nC <- obj$fit$nC2; BivD <- obj$fit$Cop2 } 

                            }
                     
    if(all(tail == c("lower","lower"))) p12 <- mm(BiCDF(p1, p2, nC, theta, obj$VC$dof), 
                                                  min.pr = obj$VC$min.pr, max.pr = obj$VC$max.pr)

    if(all(tail == c("upper","upper"))) p12 <- mm(1 - p1 - p2 + BiCDF(p1, p2, nC, theta, obj$VC$dof), 
                                                  min.pr = obj$VC$min.pr, max.pr = obj$VC$max.pr)

    if(all(tail == c("upper","lower"))) p12 <- mm(p2 - BiCDF(p1, p2, nC, theta, obj$VC$dof), 
                                                  min.pr = obj$VC$min.pr, max.pr = obj$VC$max.pr)

    if(all(tail == c("lower","upper"))) p12 <- mm(p1 - BiCDF(p1, p2, nC, theta, obj$VC$dof), 
                                                  min.pr = obj$VC$min.pr, max.pr = obj$VC$max.pr)
 
p12
 
}



jointPDF <- function(obj, y1, y2, newdata){

if(missing(newdata))                               stop("You must provide a data frame.")
if(!missing(newdata)){ 
                       if(!is.data.frame(newdata)) stop("You must provide a data frame.")
                       if(dim(newdata)[1] > 1)     stop("This calculation is supported for single row data frames.")
                     }

    VC     <- obj$VC


# Allowed combinations:
# 1) both continuous (m2)
# 2) both discrete (m1d, m2d)
# 3) one discrete (m1d/m2d) and one continuous (m2/m3)

ok1 <- (VC$margins[1] %in% VC$m2   & VC$margins[2] %in% VC$m2)
ok2 <- (VC$margins[1] %in% c(VC$m1d, VC$m2d) & VC$margins[2] %in% c(VC$m1d, VC$m2d))
ok3 <- (VC$margins[1] %in% c(VC$m1d, VC$m2d) & VC$margins[2] %in% c(VC$m2, VC$m3))

if (!(ok1 | ok2 | ok3)) {
    stop("The chosen combination of margins is not implemented in this version.\n",
         "Please get in touch if you need support for this case.")
}
   


    params <- obj$coefficients
    


    X1mat <- predict.SemiParBIV(obj, eq = 1, newdata = newdata, type = "lpmatrix")
    X2mat <- predict.SemiParBIV(obj, eq = 2, newdata = newdata, type = "lpmatrix")

    eta1 <- eta.tr( X1mat%*%params[1:VC$X1.d2],                           VC$margins[1])
    eta2 <- eta.tr( X2mat%*%params[(VC$X1.d2 + 1):(VC$X1.d2 + VC$X2.d2)], VC$margins[2])
    etad <- etas1 <- etas2 <- l.ln <- NULL 
  
 
##################**********########## 

if(VC$margins[1] %in% c(VC$m2,VC$m2d) & VC$margins[2] %in% c(VC$m2,VC$m2d)){
 
  if(is.null(VC$X3)){  
    sigma21.st <- etas1 <- params[(VC$X1.d2 + VC$X2.d2 + 1)]
    sigma22.st <- etas2 <- params[(VC$X1.d2 + VC$X2.d2 + 2)]
    teta.st    <- etad  <- params[(VC$X1.d2 + VC$X2.d2 + 3)]
                    } 
  
  
  if(!is.null(VC$X3)){
  
    X3mat <- predict.SemiParBIV(obj, eq = 3, newdata = newdata, type = "lpmatrix")
    X4mat <- predict.SemiParBIV(obj, eq = 4, newdata = newdata, type = "lpmatrix") 
    X5mat <- predict.SemiParBIV(obj, eq = 5, newdata = newdata, type = "lpmatrix") 
  
  
    sigma21.st <- etas1 <- X3mat%*%params[(VC$X1.d2 + VC$X2.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2)]
    sigma22.st <- etas2 <- X4mat%*%params[(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2)]
    teta.st    <- etad  <- X5mat%*%params[(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2 + VC$X5.d2)]
                     }  
  
    sstr1 <- esp.tr(sigma21.st, VC$margins[1])  
    sstr2 <- esp.tr(sigma22.st, VC$margins[2])  
  
    sigma21.st <- sstr1$vrb.st 
    sigma22.st <- sstr2$vrb.st 
    
    sigma21    <- sstr1$vrb 
    sigma22    <- sstr2$vrb  
 
 nu2 <- nu2.st <- 0

}  

##################**********########## 

  
if(VC$margins[1] %in% VC$m1d & VC$margins[2] %in% c(VC$m2,VC$m2d)){
 


  if(is.null(VC$X3)){  
    sigma21.st <- etas1 <- 0
    sigma22.st <- etas2 <- params[(VC$X1.d2 + VC$X2.d2 + 1)]
    teta.st    <- etad  <- params[(VC$X1.d2 + VC$X2.d2 + 2)]
                     } 
  
  
  if(!is.null(VC$X3)){ 
  
    X3mat <- predict.SemiParBIV(obj, eq = 3, newdata = newdata, type = "lpmatrix")
    X4mat <- predict.SemiParBIV(obj, eq = 4, newdata = newdata, type = "lpmatrix")   
  
    sigma21.st <- etas1 <- 0
    sigma22.st <- etas2 <- X3mat%*%params[(VC$X1.d2 + VC$X2.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2)]
    teta.st    <- etad  <- X4mat%*%params[(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2)]
                      }  
  

    sstr1 <- 0  
    sstr2 <- esp.tr(sigma22.st, VC$margins[2])  
  
    sigma21.st <- 0
    sigma22.st <- sstr2$vrb.st 
    
    sigma21    <- 0
    sigma22    <- sstr2$vrb  
 
nu2 <- nu2.st <- 0
 
}  

##################**********########## 

if(VC$margins[1] %in% VC$m1d & VC$margins[2] %in% c(VC$m1d)){
 


  if(is.null(VC$X3)){  
    sigma21.st <- etas1 <- 0
    sigma22.st <- etas2 <- 0
    teta.st    <- etad  <- params[(VC$X1.d2 + VC$X2.d2 + 1)]
                     } 
  
  
  if(!is.null(VC$X3)){ 
  
    X3mat <- predict.SemiParBIV(obj, eq = 3, newdata = newdata, type = "lpmatrix")
  
    sigma21.st <- etas1 <- 0
    sigma22.st <- etas2 <- 0
    teta.st    <- etad  <- X3mat%*%params[(VC$X1.d2 + VC$X2.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2)]
                      }  
  

    sstr1 <- 0  
    sstr2 <- 0  
  
    sigma21.st <- 0
    sigma22.st <- 0 
    
    sigma21    <- 0
    sigma22    <- 0  
 
nu2 <- nu2.st <- 0
 
}

##################**********########## 



if(VC$margins[1] %in% VC$m1d & VC$margins[2] %in% VC$m3){
 

 if(is.null(VC$X3)){  
    sigma21.st <- etas1 <- 0
    sigma22.st <- etas2 <- params[(VC$X1.d2 + VC$X2.d2 + 1)]
    nu2.st     <- etan2 <- params[(VC$X1.d2 + VC$X2.d2 + 2)]    
    teta.st    <- etad  <- params[(VC$X1.d2 + VC$X2.d2 + 3)]
  } 
  
  
 if(!is.null(VC$X3)){  
 
 
    X3mat <- predict.SemiParBIV(obj, eq = 3, newdata = newdata, type = "lpmatrix")
    X4mat <- predict.SemiParBIV(obj, eq = 4, newdata = newdata, type = "lpmatrix") 
    X5mat <- predict.SemiParBIV(obj, eq = 5, newdata = newdata, type = "lpmatrix") 
 
    sigma21.st <- etas1 <- 0
    sigma22.st <- etas2 <- X3mat%*%params[(VC$X1.d2 + VC$X2.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2)]
    nu2.st     <- etan2 <- X4mat%*%params[(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2)]
    teta.st    <- etad  <- X5mat%*%params[(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2 + VC$X5.d2)]
  }    
  
  

  
sstr1 <- 0 
sstr2 <- esp.tr(sigma22.st, VC$margins[2])  
  
sigma21.st <- 0 
sigma22.st <- sstr2$vrb.st 
    
sigma21    <- 0 
sigma22    <- sstr2$vrb 
    

sstr2 <- enu.tr(nu2.st, VC$margins[2])  
nu2.st <- sstr2$vrb.st 
nu2    <- sstr2$vrb 
    

}  



##################**********########## 


if(VC$margins[1] %in% c(VC$m2,VC$m2d) & VC$margins[2] %in% VC$m3){
 

    
  if(is.null(VC$X3)){  
    sigma21.st <- etas1 <- params[(VC$X1.d2 + VC$X2.d2 + 1)]
    sigma22.st <- etas2 <- params[(VC$X1.d2 + VC$X2.d2 + 2)]
    nu2.st     <- etan2 <- params[(VC$X1.d2 + VC$X2.d2 + 3)]    
    teta.st    <- etad  <- params[(VC$X1.d2 + VC$X2.d2 + 4)]
  } 
  
  
 if(!is.null(VC$X3)){  
 
    X3mat <- predict.SemiParBIV(obj, eq = 3, newdata = newdata, type = "lpmatrix")
    X4mat <- predict.SemiParBIV(obj, eq = 4, newdata = newdata, type = "lpmatrix") 
    X5mat <- predict.SemiParBIV(obj, eq = 5, newdata = newdata, type = "lpmatrix") 
    X6mat <- predict.SemiParBIV(obj, eq = 6, newdata = newdata, type = "lpmatrix") 
    
    
    sigma21.st <- etas1 <- X3mat%*%params[(VC$X1.d2 + VC$X2.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2)]
    sigma22.st <- etas2 <- X4mat%*%params[(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2)]
    nu2.st     <- etan2 <- X5mat%*%params[(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2 + VC$X5.d2)]
    teta.st    <- etad  <- X6mat%*%params[(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2 + VC$X5.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2 + VC$X5.d2 + VC$X6.d2)]
  }    
  
  
    sstr1 <- esp.tr(sigma21.st, VC$margins[1])  
    sstr2 <- esp.tr(sigma22.st, VC$margins[2])  
  
    sigma21.st <- sstr1$vrb.st 
    sigma22.st <- sstr2$vrb.st 
    
    sigma21    <- sstr1$vrb 
    sigma22    <- sstr2$vrb 
    

sstr2 <- enu.tr(nu2.st, VC$margins[2])  
nu2.st <- sstr2$vrb.st 
nu2    <- sstr2$vrb 
    
    
}  
    
    
  
  
##################
## Transformations
##################
  

resT    <- teta.tr(VC, teta.st)

teta.st1 <- teta.st2 <- teta.st <- resT$teta.st
teta1 <- teta2 <- teta <- resT$teta 
    
##################

Cop1 <- Cop2 <- VC$BivD 
nC1 <- nC2 <- VC$nC 


if(VC$BivD %in% VC$BivD2 && length(teta.st) == 1){


Cop1Cop2R <- Cop1Cop2(VC$BivD)
Cop1 <- Cop1Cop2R$Cop1
Cop2 <- Cop1Cop2R$Cop2
nC1  <- VC$ct[which(VC$ct[,1] == Cop1), 2] 
nC2  <- VC$ct[which(VC$ct[,1] == Cop2), 2]

teta.i1 <- as.logical( etad >= 0 ) # this becomes an index for + and - coefs

if(teta.i1 == TRUE){  Cop2 <- Cop1; nC2 <- nC1}
if(teta.i1 == FALSE){ Cop1 <- Cop2; nC1 <- nC2}   

if(VC$BivD %in% VC$BivD2[c(1:4,13:16)]) teta <- etad                                     
if(VC$BivD %in% VC$BivD2[5:12])         teta <- sign(etad) * (abs(etad) + 1)


} 

 
#####*************############# 
#####*************############# 

 
if(VC$margins[1] %in% c(VC$m1d, VC$m2d) & VC$margins[2] %in% c(VC$m2,VC$m3)){
 

##################

  dHs1 <- distrHsATDiscr(y1, c(eta1), c(sigma21), nu = c(1), margin2=VC$margins[1], y2m = VC$y1m, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr, left.trunc = VC$left.trunc1)                 
  dHs2 <-      distrHsAT(y2, c(eta2), c(sigma22), c(nu2),    margin2=VC$margins[2],               min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr, left.trunc = VC$left.trunc2)

  pdf1 <- dHs1$pdf2
  pdf2 <- dHs2$pdf2

  p1 <- dHs1$p2
  p2 <- dHs2$p2
    
  dH11 <- copgHsCond(c(p1), c(p2), c(teta), VC$dof, Cop1, min.pr = VC$min.pr, max.pr = VC$max.pr)
  h1 <- dH11$c.copula.be2
     

  dH21 <- copgHsCond(c(mm(p1-pdf1, min.pr = VC$min.pr, max.pr = VC$max.pr)), c(p2), c(teta), VC$dof, Cop1, min.pr = VC$min.pr, max.pr = VC$max.pr)
  h2 <- dH21$c.copula.be2
  
  diffh1.h2 <- h1 - h2
  diffh1.h2 <- mm(diffh1.h2, min.pr = VC$min.pr, max.pr = VC$max.pr)  

  res <- c(pdf2*diffh1.h2) 
  



  
}



if(VC$margins[1] %in% c(VC$m1d, VC$m2d) & VC$margins[2] %in% c(VC$m1d,VC$m2d)){



  dHs1 <- distrHsATDiscr(y1, c(eta1), c(sigma21), nu = c(1), margin2=VC$margins[1], y2m = VC$y1m, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr, left.trunc = VC$left.trunc1)
  dHs2 <- distrHsATDiscr(y2, c(eta2), c(sigma22), nu = c(1), margin2=VC$margins[2], y2m = VC$y2m, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr, left.trunc = VC$left.trunc2)

  pdf1 <- dHs1$pdf2
  pdf2 <- dHs2$pdf2

  p1 <- dHs1$p2
  p2 <- dHs2$p2
  

  
  C11 <- mm(BiCDF(p1, p2, nC1, teta, VC$dof), min.pr = VC$min.pr, max.pr = VC$max.pr  )
  C01 <- mm(BiCDF(mm(p1-pdf1, min.pr = VC$min.pr, max.pr = VC$max.pr), p2, nC1, teta, VC$dof), min.pr = VC$min.pr, max.pr = VC$max.pr  )
  C10 <- mm(BiCDF(p1, mm(p2-pdf2, min.pr = VC$min.pr, max.pr = VC$max.pr), nC1, teta, VC$dof), min.pr = VC$min.pr, max.pr = VC$max.pr  )
  C00 <- mm(BiCDF(mm(p1-pdf1, min.pr = VC$min.pr, max.pr = VC$max.pr), mm(p2-pdf2, min.pr = VC$min.pr, max.pr = VC$max.pr), nC1, teta, VC$dof), min.pr = VC$min.pr, max.pr = VC$max.pr  )
  

  
  res <- c( mm(C11 - C01 - C10 + C00, min.pr = VC$min.pr, max.pr = VC$max.pr) ) 

}




 
if(VC$margins[1] %in% c(VC$m2) & VC$margins[2] %in% c(VC$m2)){
 

##################


  dHs1 <- distrHsAT(y1, c(eta1), c(sigma21), c(1), margin2=VC$margins[1], min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr, left.trunc = VC$left.trunc1)
  dHs2 <- distrHsAT(y2, c(eta2), c(sigma22), c(1), margin2=VC$margins[2], min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr, left.trunc = VC$left.trunc2)


  pdf1 <- dHs1$pdf2
  pdf2 <- dHs2$pdf2

  p1 <- dHs1$p2
  p2 <- dHs2$p2
  
  
  dH1 <- copgHsAT(p1, p2, teta, Cop1, Ln = TRUE, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
  c.copula2.be1be2 <- dH1$c.copula2.be1be2
  

  res <- c( pdf1*pdf2*c.copula2.be1be2 )
  
}

  
as.numeric(res)  
  

}





.qqenv.chisq2 <- function(res, n.sim, prob.lev, df = 2){

  n <- length(res)
  sims <- replicate(n.sim, sort(rchisq(n, df = df)))
  theo <- qchisq((1:n - 0.5)/n, df = df)
  xlab <- "Theoretical Quantiles"
  ylab <- "Sample Quantiles"

  lim <- apply(sims, 1, quantile, c(prob.lev/2, 1 - prob.lev/2))

  qqplot(theo, res, pch = ".", ylim = range(c(lim, res)), xlab = xlab, ylab = ylab)

  polygon(c(theo, theo[n:1], theo[1]), c(lim[1, ], lim[2, n:1], lim[1, 1]), col = "grey80", border = NA)

  abline(0, 1)
  points(theo, sort(res), pch = ".")
}




.qqenv <- function(res, n.sim, prob.lev){

  n <- length(res)

    sims <- replicate(n.sim, sort(rnorm(n)))
    theo <- qnorm((1:n - 0.5)/n)
    xlab <- "Theoretical Quantiles" 
    ylab <- "Sample Quantiles"

  lim <- apply(sims, 1, quantile, c(prob.lev/2, 1 - prob.lev/2))

  qqplot(theo, res, pch = ".", ylim = range(c(lim, res)), xlab = xlab, ylab = ylab)

  polygon(c(theo, theo[n:1], theo[1]), c(lim[1, ], lim[2, n:1], lim[1, 1]), col = "grey80", border = NA)

  abline(0, 1)
  points(theo, sort(res), pch = ".")
}





rqres_ordinal <- function(Fm, Fp, min.pr = 1e-16, max.pr = 1 - 1e-16) {
  u <- runif(length(Fm), Fm, Fp)
  u <- pmin(pmax(u, min.pr), max.pr)
  qnorm(u)
}





rqres_nrsp_F <- function(F1, cens, F2 = NULL, seed = NULL, min.pr = 1e-16, max.pr = 1 - 1e-16){

  if(!is.null(seed)) set.seed(seed)

  n <- length(F1)
  u <- numeric(n)


  if(is.numeric(cens)){

    event <- cens == 1

    # uncensored
    u[event] <- F1[event]

    # right censored
    u[!event] <- runif(sum(!event), F1[!event], 1)
  }



  if(is.factor(cens)){

    cens <- as.character(cens)

    # U = uncensored
    indU <- cens == "U"
    u[indU] <- F1[indU]

    # R = right censored
    indR <- cens == "R"
    u[indR] <- runif(sum(indR), F1[indR], 1)

    # L = left censored
    indL <- cens == "L"
    u[indL] <- runif(sum(indL), 0, F1[indL])

    # I = interval censored
    indI <- cens == "I"

    if(any(indI)) u[indI] <- runif(sum(indI), F1[indI], F2[indI])
     
  }


  u <- pmin(pmax(u, min.pr), max.pr)
  r <- qnorm(u)

  return(r)
}





fix.fenv <- function(x, env = .GlobalEnv) {
  x$formula <- lapply(x$formula, function(f) {
    if (inherits(f, "formula"))
      environment(f) <- env
    f
  })
  x
}


## ****** this needs to be extended for count case
cop.mod.sel <- function(eq, data, sig.lev, crit, dof.grid, infl.fac) {

  #all.cop <- c("N",
  #  "C0C90","C0C270","C180C90","C180C270",
  #  "G0G90","G0G270","G180G90","G180G270",
  #  "J0J90","J0J270","J180J90","J180J270",
  #  "F","AMH","FGM","PL","T")
    
  all.cop <- c("N",
    "C0C90","C0C270","C180C90","C180C270",
    "G0G90","G0G270","G180G90","G180G270",
    "J0J90","J0J270","J180J90","J180J270",
    "F","T")

  model.list <- list()
  idx <- 1
  n <- NULL

  ## --------------------------------------------------
  ## Fit all copulas
  ## --------------------------------------------------
  for (cop in all.cop) {

    if (cop == "T") {

      for (dof in dof.grid) {

        fit <- suppressWarnings(
          tryCatch(
            gjrm(eq, fixed.margins = TRUE, data = data,
                 copula = "T", dof = dof, infl.fac = infl.fac),
            error = function(e) NULL
          )
        )
        
if (!is.null(fit)) {
  if (!is.null(fit$fit$gradient) &&
      max(abs(fit$fit$gradient)) > 0.01) {
    fit <- NULL
  }
}
        if (!is.null(fit)) {
          model.list[[idx]] <- fit
          idx <- idx + 1
          if (is.null(n)) n <- fit$n
        }
      }

    } else {

      fit <- suppressWarnings(
        tryCatch(
          gjrm(eq, fixed.margins = TRUE, data = data,
               copula = cop, infl.fac = infl.fac),
          error = function(e) NULL
        )
      )
      
if (!is.null(fit)) {
  if (!is.null(fit$fit$gradient) &&
      max(abs(fit$fit$gradient)) > 0.01) {
    fit <- NULL
  }
}
      if (!is.null(fit)) {
        model.list[[idx]] <- fit
        idx <- idx + 1
        if (is.null(n)) n <- fit$n
      }
    }
  }

  ## --------------------------------------------------
  ## Remove failed fits
  ## --------------------------------------------------
  model.list <- Filter(Negate(is.null), model.list)

  l.c <- length(model.list)

  if (l.c == 0) stop("No models were successfully fitted.")

  preference.count <- integer(l.c)

  ## ==================================================
  ## MODEL SELECTION
  ## ==================================================

  if (crit == "vuong") {

    for (i in 1:(l.c - 1)) {
      for (j in (i + 1):l.c) {

        vtr <- invisible(
          capture.output(
            tmp <- vuong.test(model.list[[i]], model.list[[j]])
          )
        )

        vtr <- tmp
        if (is.null(vtr)) next

        vuong.stat <- vtr[[1]]
        thr <- qnorm(1 - sig.lev / 2)

        if (abs(vuong.stat) <= thr) next

        if (vuong.stat >  thr)
          preference.count[i] <- preference.count[i] + 1

        if (vuong.stat < -thr)
          preference.count[j] <- preference.count[j] + 1
      }
    }

    fm <- model.list[[which.max(preference.count)]]
  }

  ## --------------------------------------------------

  if (crit == "clarke") {

    for (i in 1:(l.c - 1)) {
      for (j in (i + 1):l.c) {

        b <- invisible(
          capture.output(
            tmp <- clarke.test(model.list[[i]], model.list[[j]])
          )
        )

        b <- tmp[1]
        if (is.null(b)) next

        if (b >= n / 2) {
          pvalue <- 2 * (1 - pbinom(b - 1, n, 0.5))
          if (pvalue <= sig.lev)
            preference.count[i] <- preference.count[i] + 1
        }

        if (b < n / 2) {
          pvalue <- 2 * pbinom(b, n, 0.5)
          if (pvalue <= sig.lev)
            preference.count[j] <- preference.count[j] + 1
        }
      }
    }

    fm <- model.list[[which.max(preference.count)]]
  }

  ## --------------------------------------------------

  if (crit == "aic") {
    aic.values <- sapply(model.list, AIC)
    fm <- model.list[[which.min(aic.values)]]
  }

  if (crit == "bic") {
    bic.values <- sapply(model.list, BIC)
    fm <- model.list[[which.min(bic.values)]]
  }

  ## --------------------------------------------------
  ## Fix formula environments
  ## --------------------------------------------------
  fm$formula <- lapply(fm$formula, function(f) {
    environment(f) <- parent.env(environment(eq[[1]]))
    f
  })

  return(fm)
}




cdf.pcc <- function(y, y1 = NULL, y2 = NULL, y3 = NULL, eta1, eta2, eta3, sig1, sig2, sig3, nu1 = NULL, nu2 = NULL, nu3 = NULL, 
                       theta1, theta2, theta3 = NULL, beta.theta3 = NULL, 
                       margins, cop.sel, dof.vec, obj.c, simplified, eq.theta, newdata, 
                       min.dn, min.pr, max.pr, eq, left.trunc.v, cond.var, nC, sign.theta = NULL){           
 
ys <- list(y1, y2, y3)
ys[[cond.var]] <- y
l.cv <- length(y)

rpd1 <- distrHsAT(ys[[1]], eta2 = eta1, sigma2 = sig1, nu = nu1, margin2 = margins[1], min.dn = min.dn, min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc.v[1])
rpd2 <- distrHsAT(ys[[2]], eta2 = eta2, sigma2 = sig2, nu = nu2, margin2 = margins[2], min.dn = min.dn, min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc.v[2])
rpd3 <- distrHsAT(ys[[3]], eta2 = eta3, sigma2 = sig3, nu = nu3, margin2 = margins[3], min.dn = min.dn, min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc.v[3])
                  
p1   <- rpd1$p2 
pdf1 <- rpd1$pdf2 
p2   <- rpd2$p2 
pdf2 <- rpd2$pdf2 
p3   <- rpd3$p2 
pdf3 <- rpd3$pdf2 
 
pdf.c <- list(pdf1, pdf2, pdf3)[[cond.var]]


if(cond.var == 1){
 
   p1 <- if (length(p1) == 1) rep(p1, l.cv) else p1 # done for T case
   p2 <- if (length(p2) == 1) rep(p2, l.cv) else p2
 
   res.F <- copgHsAT(p1 = p1, p2 = p2, teta = theta1, BivD = cop.sel[1], Ln = TRUE, par2 = dof.vec[1], min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)
   p1.h  <- res.F$c.copula.be1

   p1 <- if (length(p1) == 1) rep(p1, l.cv) else p1
   p3 <- if (length(p3) == 1) rep(p3, l.cv) else p3

   res.S <- copgHsAT(p1 = p1, p2 = p3, teta = theta2, BivD = cop.sel[2], Ln = TRUE, par2 = dof.vec[2], min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)
   p2.h  <- res.S$c.copula.be1
   
                }  
                
if(cond.var == 2){

   p1 <- if (length(p1) == 1) rep(p1, l.cv) else p1
   p2 <- if (length(p2) == 1) rep(p2, l.cv) else p2
 
   res.F <- copgHsAT(p1 = p1, p2 = p2, teta = theta1, BivD = cop.sel[1], Ln = TRUE, par2 = dof.vec[1], min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)
   p1.h  <- res.F$c.copula.be2

   p2 <- if (length(p2) == 1) rep(p2, l.cv) else p2
   p3 <- if (length(p3) == 1) rep(p3, l.cv) else p3

   res.S <- copgHsAT(p1 = p2, p2 = p3, teta = theta2, BivD = cop.sel[2], Ln = TRUE, par2 = dof.vec[2], min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)
   p2.h  <- res.S$c.copula.be1
   
                }
                
if(cond.var == 3){

   p1 <- if (length(p1) == 1) rep(p1, l.cv) else p1
   p3 <- if (length(p3) == 1) rep(p3, l.cv) else p3
        
   res.F <- copgHsAT(p1 = p1, p2 = p3, teta = theta1, BivD = cop.sel[1], Ln = TRUE, par2 = dof.vec[1], min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)
   p1.h  <- res.F$c.copula.be2
   
   p2 <- if (length(p2) == 1) rep(p2, l.cv) else p2
   p3 <- if (length(p3) == 1) rep(p3, l.cv) else p3   

   res.S <- copgHsAT(p1 = p2, p2 = p3, teta = theta2, BivD = cop.sel[2], Ln = TRUE, par2 = dof.vec[2], min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)
   p2.h  <- res.S$c.copula.be2
   
                }                
                             

   if(!(simplified == TRUE && is.null(eq.theta))){
   
          if(simplified == FALSE){ 
                     if(cond.var == 1) Xpr.t <- predict(obj.c, newdata = merge(data.frame(.gjrm_p1 = p1), newdata), type = "lpmatrix" ) 
                     if(cond.var == 2) Xpr.t <- predict(obj.c, newdata = merge(data.frame(.gjrm_p2 = p2), newdata), type = "lpmatrix" ) 
                     if(cond.var == 3) Xpr.t <- predict(obj.c, newdata = merge(data.frame(.gjrm_p3 = p3), newdata), type = "lpmatrix" ) 
                                 }
          
          if(simplified == TRUE)  Xpr.t <- predict(obj.c, newdata = newdata, type = "lpmatrix" ) # a bit inefficient here but kept for clarity
          
          eta.theta3 <- Xpr.t %*% beta.theta3
          theta3    <- c(teta.tr(obj.c$VC, eta.theta3)$teta)     


    if(obj.c$copula %in% obj.c$BivD2){
  
      if(sign.theta ==  1) nC <- obj.c$fit$nC1; theta3 <-    abs(theta3)
      if(sign.theta == -1) nC <- obj.c$fit$nC2; theta3 <- -1*abs(theta3)
      
                                      }
                                                  }

c.T <- mm( BiCDF(p1.h, p2.h, nC, theta3, dof.vec[3]), min.pr = min.pr, max.pr = max.pr)

c.T * pdf.c

}






mar.cond.var <- function(y, y1 = NULL, y2 = NULL, y3 = NULL, eta1, eta2, eta3, sig1, sig2, sig3, nu1 = NULL, nu2 = NULL, nu3 = NULL, theta1, theta2, 
                         theta3 = NULL, beta.theta3 = NULL, 
                         margins, cop.sel, dof.vec, obj.c, simplified, eq.theta, newdata, 
                         min.dn, min.pr, max.pr, eq, left.trunc.v, cond.var, sI = NULL){           
 
ys <- list(y1, y2, y3)
ys[[cond.var]] <- y
l.cv <- length(y)

rpd1 <- distrHsAT(ys[[1]], eta2 = eta1, sigma2 = sig1, nu = nu1, margin2 = margins[1], min.dn = min.dn, min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc.v[1])
rpd2 <- distrHsAT(ys[[2]], eta2 = eta2, sigma2 = sig2, nu = nu2, margin2 = margins[2], min.dn = min.dn, min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc.v[2])
rpd3 <- distrHsAT(ys[[3]], eta2 = eta3, sigma2 = sig3, nu = nu3, margin2 = margins[3], min.dn = min.dn, min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc.v[3])
                  
p1   <- rpd1$p2 
pdf1 <- rpd1$pdf2 
p2   <- rpd2$p2 
pdf2 <- rpd2$pdf2 
p3   <- rpd3$p2 
pdf3 <- rpd3$pdf2 
 
pdf.c <- list(pdf1, pdf2, pdf3)[[cond.var]]


if(cond.var == 1){
 
   p1 <- if (length(p1) == 1) rep(p1, l.cv) else p1 # done for T case
   p2 <- if (length(p2) == 1) rep(p2, l.cv) else p2
 
   res.F <- copgHsAT(p1 = p1, p2 = p2, teta = theta1, BivD = cop.sel[1], Ln = TRUE, par2 = dof.vec[1], min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)
   c.F   <- res.F$c.copula2.be1be2
   p1.h  <- res.F$c.copula.be1

   p1 <- if (length(p1) == 1) rep(p1, l.cv) else p1
   p3 <- if (length(p3) == 1) rep(p3, l.cv) else p3

   res.S <- copgHsAT(p1 = p1, p2 = p3, teta = theta2, BivD = cop.sel[2], Ln = TRUE, par2 = dof.vec[2], min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)
   c.S   <- res.S$c.copula2.be1be2
   p2.h  <- res.S$c.copula.be1
   
                }  
                
if(cond.var == 2){

   p1 <- if (length(p1) == 1) rep(p1, l.cv) else p1
   p2 <- if (length(p2) == 1) rep(p2, l.cv) else p2
 
   res.F <- copgHsAT(p1 = p1, p2 = p2, teta = theta1, BivD = cop.sel[1], Ln = TRUE, par2 = dof.vec[1], min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)
   c.F   <- res.F$c.copula2.be1be2
   p1.h  <- res.F$c.copula.be2

   p2 <- if (length(p2) == 1) rep(p2, l.cv) else p2
   p3 <- if (length(p3) == 1) rep(p3, l.cv) else p3

   res.S <- copgHsAT(p1 = p2, p2 = p3, teta = theta2, BivD = cop.sel[2], Ln = TRUE, par2 = dof.vec[2], min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)
   c.S   <- res.S$c.copula2.be1be2
   p2.h  <- res.S$c.copula.be1
   
                }
                
if(cond.var == 3){

   p1 <- if (length(p1) == 1) rep(p1, l.cv) else p1
   p3 <- if (length(p3) == 1) rep(p3, l.cv) else p3
        
   res.F <- copgHsAT(p1 = p1, p2 = p3, teta = theta1, BivD = cop.sel[1], Ln = TRUE, par2 = dof.vec[1], min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)
   c.F   <- res.F$c.copula2.be1be2
   p1.h  <- res.F$c.copula.be2
   
   p2 <- if (length(p2) == 1) rep(p2, l.cv) else p2
   p3 <- if (length(p3) == 1) rep(p3, l.cv) else p3   

   res.S <- copgHsAT(p1 = p2, p2 = p3, teta = theta2, BivD = cop.sel[2], Ln = TRUE, par2 = dof.vec[2], min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)
   c.S   <- res.S$c.copula2.be1be2
   p2.h  <- res.S$c.copula.be2
   
                }                
                             

   if(!(simplified == TRUE && is.null(eq.theta))){
   
          if(simplified == FALSE){ 
                     if(cond.var == 1) Xpr.t <- predict(obj.c, newdata = merge(data.frame(.gjrm_p1 = p1), newdata), type = "lpmatrix" ) 
                     if(cond.var == 2) Xpr.t <- predict(obj.c, newdata = merge(data.frame(.gjrm_p2 = p2), newdata), type = "lpmatrix" ) 
                     if(cond.var == 3) Xpr.t <- predict(obj.c, newdata = merge(data.frame(.gjrm_p3 = p3), newdata), type = "lpmatrix" ) 
                                 }
          
          if(simplified == TRUE)  Xpr.t <- predict(obj.c, newdata = newdata, type = "lpmatrix" ) 
          
          eta.theta3 <- Xpr.t %*% beta.theta3
          theta3    <- c(teta.tr(obj.c$VC, eta.theta3)$teta)     
    
    if(obj.c$copula %in% obj.c$BivD2){
  
      if(sI ==  1) nC <- obj.c$fit$nC1; theta3 <-    abs(theta3)
      if(sI == -1) nC <- obj.c$fit$nC2; theta3 <- -1*abs(theta3)
      
                                      }    
    
                           }


  c.T <- copgHsAT(p1 = p1.h, p2 = p2.h, teta = theta3, BivD = cop.sel[3], Ln = TRUE, par2 = dof.vec[3], 
                  min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)$c.copula2.be1be2


  c.F * c.S * c.T * pdf.c

}




cond.var.i <- function(y, y1 = NULL, y2 = NULL, y3 = NULL, eta1, eta2, eta3, sig1, sig2, sig3, nu1 = NULL, nu2 = NULL, nu3 = NULL, theta1, theta2, 
                       theta3 = NULL, beta.theta3 = NULL, 
                       margins, cop.sel, dof.vec, obj.c, simplified, eq.theta, newdata, 
                       min.dn, min.pr, max.pr, eq, left.trunc.v, nor.const = NULL, var.y = FALSE, cond.var, sI = NULL){           
 
if(var.y == FALSE) y.a <- y else y.a <- y^2

ys <- list(y1, y2, y3)

ys[[eq]] <- y 
l.cv <- length(y)

rpd1 <- distrHsAT(c(ys[[1]]), eta2 = eta1, sigma2 = sig1, nu = nu1, margin2 = margins[1], min.dn = min.dn, min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc.v[1])
rpd2 <- distrHsAT(c(ys[[2]]), eta2 = eta2, sigma2 = sig2, nu = nu2, margin2 = margins[2], min.dn = min.dn, min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc.v[2])   
rpd3 <- distrHsAT(c(ys[[3]]), eta2 = eta3, sigma2 = sig3, nu = nu3, margin2 = margins[3], min.dn = min.dn, min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc.v[3])   

pdf1 <- rpd1$pdf2
p1   <- rpd1$p2 
pdf2 <- rpd2$pdf2
p2   <- rpd2$p2 
pdf3 <- rpd3$pdf2
p3   <- rpd3$p2  

pdf.c <- list(pdf1, pdf2, pdf3)[[eq]]


 if(cond.var == 1){

   p1 <- if (length(p1) == 1) rep(p1, l.cv) else p1 # done for T case
   p2 <- if (length(p2) == 1) rep(p2, l.cv) else p2
  
    res.F <- copgHsAT(p1 = p1, p2 = p2, teta = theta1, BivD = cop.sel[1], Ln = TRUE, par2 = dof.vec[1], min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)
    c.F   <- res.F$c.copula2.be1be2
    p1.h  <- res.F$c.copula.be1

   p1 <- if (length(p1) == 1) rep(p1, l.cv) else p1
   p3 <- if (length(p3) == 1) rep(p3, l.cv) else p3

 
    res.S <- copgHsAT(p1 = p1, p2 = p3, teta = theta2, BivD = cop.sel[2], Ln = TRUE, par2 = dof.vec[2], min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)
    c.S   <- res.S$c.copula2.be1be2
    p2.h  <- res.S$c.copula.be1
    
                 }  
                 
 if(cond.var == 2){

   p1 <- if (length(p1) == 1) rep(p1, l.cv) else p1
   p2 <- if (length(p2) == 1) rep(p2, l.cv) else p2
  
    res.F <- copgHsAT(p1 = p1, p2 = p2, teta = theta1, BivD = cop.sel[1], Ln = TRUE, par2 = dof.vec[1], min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)
    c.F   <- res.F$c.copula2.be1be2
    p1.h  <- res.F$c.copula.be2
 
   p2 <- if (length(p2) == 1) rep(p2, l.cv) else p2
   p3 <- if (length(p3) == 1) rep(p3, l.cv) else p3
 
    res.S <- copgHsAT(p1 = p2, p2 = p3, teta = theta2, BivD = cop.sel[2], Ln = TRUE, par2 = dof.vec[2], min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)
    c.S   <- res.S$c.copula2.be1be2
    p2.h  <- res.S$c.copula.be1
    
                 }
                 
 if(cond.var == 3){

   p1 <- if (length(p1) == 1) rep(p1, l.cv) else p1
   p3 <- if (length(p3) == 1) rep(p3, l.cv) else p3
         
    res.F <- copgHsAT(p1 = p1, p2 = p3, teta = theta1, BivD = cop.sel[1], Ln = TRUE, par2 = dof.vec[1], min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)
    c.F   <- res.F$c.copula2.be1be2
    p1.h  <- res.F$c.copula.be2

   p2 <- if (length(p2) == 1) rep(p2, l.cv) else p2
   p3 <- if (length(p3) == 1) rep(p3, l.cv) else p3 
 
    res.S <- copgHsAT(p1 = p2, p2 = p3, teta = theta2, BivD = cop.sel[2], Ln = TRUE, par2 = dof.vec[2], min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)
    c.S   <- res.S$c.copula2.be1be2
    p2.h  <- res.S$c.copula.be2
    
                 }                
                              
 
    if(!(simplified == TRUE && is.null(eq.theta))){
    
           if(simplified == FALSE){ 
                      if(cond.var == 1) Xpr.t <- predict(obj.c, newdata = merge(data.frame(.gjrm_p1 = p1), newdata), type = "lpmatrix" ) 
                      if(cond.var == 2) Xpr.t <- predict(obj.c, newdata = merge(data.frame(.gjrm_p2 = p2), newdata), type = "lpmatrix" ) 
                      if(cond.var == 3) Xpr.t <- predict(obj.c, newdata = merge(data.frame(.gjrm_p3 = p3), newdata), type = "lpmatrix" ) 
                                  }
           
           if(simplified == TRUE)  Xpr.t <- predict(obj.c, newdata = newdata, type = "lpmatrix" ) 
           
           eta.theta3 <- Xpr.t %*% beta.theta3
           theta3    <- c(teta.tr(obj.c$VC, eta.theta3)$teta)     
     
    if(obj.c$copula %in% obj.c$BivD2){
  
      if(sI ==  1) nC <- obj.c$fit$nC1; theta3 <-    abs(theta3)
      if(sI == -1) nC <- obj.c$fit$nC2; theta3 <- -1*abs(theta3)
      
                                      }       
     
     
     
                            }
 
   c.T <- copgHsAT(p1 = p1.h, p2 = p2.h, teta = theta3, BivD = cop.sel[3], Ln = TRUE, par2 = dof.vec[3], 
                  min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)$c.copula2.be1be2



#         if(eq == 1 && cond.var == 1) y.a * c.F * c.S * c.T * pdf.c / nor.const
#         if(eq == 2 && cond.var == 1) y.a * c.F * c.T * pdf.c
#         if(eq == 3 && cond.var == 1) y.a * c.S * c.T * pdf.c
#         
#         if(eq == 1 && cond.var == 2) y.a * c.F * c.T * pdf.c
#         if(eq == 2 && cond.var == 2) y.a * c.F * c.S * c.T * pdf.c  / nor.const
#         if(eq == 3 && cond.var == 2) y.a * c.S * c.T * pdf.c
#         
#         if(eq == 1 && cond.var == 3) y.a * c.F * c.T * pdf.c
#         if(eq == 2 && cond.var == 3) y.a * c.S * c.T * pdf.c
#         if(eq == 3 && cond.var == 3) y.a * c.F * c.S * c.T * pdf.c  / nor.const


if (eq == 1 && cond.var == 1) {
  y.a * c.F * c.S * c.T * pdf.c / nor.const
} else if (eq == 2 && cond.var == 1) {
  y.a * c.F * c.T * pdf.c
} else if (eq == 3 && cond.var == 1) {
  y.a * c.S * c.T * pdf.c
} else if (eq == 1 && cond.var == 2) {
  y.a * c.F * c.T * pdf.c
} else if (eq == 2 && cond.var == 2) {
  y.a * c.F * c.S * c.T * pdf.c / nor.const
} else if (eq == 3 && cond.var == 2) {
  y.a * c.S * c.T * pdf.c
} else if (eq == 1 && cond.var == 3) {
  y.a * c.F * c.T * pdf.c
} else if (eq == 2 && cond.var == 3) {
  y.a * c.S * c.T * pdf.c
} else if (eq == 3 && cond.var == 3) {
  y.a * c.F * c.S * c.T * pdf.c / nor.const
}


                 

}












############################





rNtr <- function(n, mu, sigma, left.trunc){
               
              n  <- ceiling(n)
              p  <- runif(n)
              pp <- pnorm(rep(left.trunc, length(mu)), mu, sigma)
              p  <- (pp + p * (1 - pp))
              qnorm(p, mean = mu, sd = sigma, lower.tail = TRUE)

}




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












ordcont.m <- function(y, pk, pk1, eta.mu, sigma, nu, theta, margin, BivD, par2, min.dn, min.pr, max.pr, y1, left.trunc2){           
     

p12.s <- 0

pd2  <- distrHsAT(y, eta2 = eta.mu, sigma2 = sigma, nu = NULL, margin2 = margin, min.dn = min.dn, min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc2)
p2   <- pd2$p2
pdf2 <- pd2$pdf2


if(y1 > 1){
            p12.f <- copgHs2(p1 = pk,  p2 = p2, teta = theta, teta.st = theta, BivD = BivD, par2 = par2, min.pr = min.pr, max.pr = max.pr, min.dn = min.dn)$c.copula.be2
            p12.s <- copgHs2(p1 = pk1, p2 = p2, teta = theta, teta.st = theta, BivD = BivD, par2 = par2, min.pr = min.pr, max.pr = max.pr, min.dn = min.dn)$c.copula.be2
            
          }

if(y1 == 1) p12.f <- copgHs2(p1 = pk,  p2 = p2, teta = theta, teta.st = theta, BivD = BivD, par2 = par2, min.pr = min.pr, max.pr = max.pr, min.dn = min.dn)$c.copula.be2


y*(p12.f - p12.s)*pdf2/(pk - pk1) 


 
}


ordcont.v <- function(y, pk, pk1, eta.mu, sigma, nu, theta, margin, BivD, par2, min.dn, min.pr, max.pr, y1, left.trunc2){           
     

p12.s <- 0

pd2  <- distrHsAT(y, eta2 = eta.mu, sigma2 = sigma, nu = NULL, margin2 = margin, min.dn = min.dn, min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc2)
p2   <- pd2$p2
pdf2 <- pd2$pdf2


if(y1 > 1){
            p12.f <- copgHs2(p1 = pk,  p2 = p2, teta = theta, teta.st = theta, BivD = BivD, par2 = par2, min.pr = min.pr, max.pr = max.pr, min.dn = min.dn)$c.copula.be2
            p12.s <- copgHs2(p1 = pk1, p2 = p2, teta = theta, teta.st = theta, BivD = BivD, par2 = par2, min.pr = min.pr, max.pr = max.pr, min.dn = min.dn)$c.copula.be2
            
          }

if(y1 == 1) p12.f <- copgHs2(p1 = pk,  p2 = p2, teta = theta, teta.st = theta, BivD = BivD, par2 = par2, min.pr = min.pr, max.pr = max.pr, min.dn = min.dn)$c.copula.be2


y^2*(p12.f - p12.s)*pdf2/(pk - pk1) 


 
}







dicont.m <- function(y, y12, eta1, eta2, sig1, sig2, nu1 = NULL, nu2 = NULL, theta, margins, BivD, par2, min.dn, min.pr, max.pr, nC, left.trunc1, left.trunc2){           
 
y1m <- y2m <- NA

p1pdf1 <- distrHsATDiscr(y,   eta2 = eta1, sigma2 = sig1, nu = 1, margin2 = margins[1], y2m = y1m, robust = FALSE, min.dn = min.pr, 
                         min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc1)
p1     <- as.numeric(p1pdf1$p2)
pmf1   <- as.numeric(p1pdf1$pdf2)


rpd2 <- distrHsAT(y12,   eta2 = eta2, sigma2 = sig2, nu = nu2, margin2 = margins[2], min.dn = min.dn, min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc2)     
p2   <- rpd2$p2
pdf2 <- rpd2$pdf2

                                             
  c.be2.h1 <- copgHs2(p1 = p1, p2 = p2, teta = theta, teta.st = theta, BivD = BivD, par2 = par2, min.pr = min.pr, max.pr = max.pr, min.dn = min.dn)$c.copula.be2
  c.be2.h2 <- copgHs2(p1 = mm(p1 - pmf1, min.pr = min.pr, max.pr = max.pr), p2 = p2, teta = theta, teta.st = theta, BivD = BivD, par2 = par2, min.pr = min.pr, max.pr = max.pr, min.dn = min.dn)$c.copula.be2
  

  
  diffh1.h2 <- c.be2.h1 - c.be2.h2  # diffh1.h2 <- mm(diffh1.h2, min.pr = min.pr, max.pr = max.pr)  

  
  if(length(y) == 2) y <- y[1]
  
  y*diffh1.h2


}






dicont.v <- function(y, y12, eta1, eta2, sig1, sig2, nu1 = NULL, nu2 = NULL, theta, margins, BivD, par2, min.dn, min.pr, max.pr, nC, left.trunc1, left.trunc2){           
 
y1m <- y2m <- NA

p1pdf1 <- distrHsATDiscr(y,   eta2 = eta1, sigma2 = sig1, nu = 1, margin2 = margins[1], y2m = y1m, robust = FALSE, min.dn = min.pr, 
                         min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc1)
p1     <- as.numeric(p1pdf1$p2)
pmf1   <- as.numeric(p1pdf1$pdf2)


rpd2 <- distrHsAT(y12,   eta2 = eta2, sigma2 = sig2, nu = nu2, margin2 = margins[2], min.dn = min.dn, min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc2)     
p2   <- rpd2$p2
pdf2 <- rpd2$pdf2

                                             
  c.be2.h1 <- copgHs2(p1 = p1, p2 = p2, teta = theta, teta.st = theta, BivD = BivD, par2 = par2, min.pr = min.pr, max.pr = max.pr, min.dn = min.dn)$c.copula.be2
  c.be2.h2 <- copgHs2(p1 = mm(p1 - pmf1, min.pr = min.pr, max.pr = max.pr), p2 = p2, teta = theta, teta.st = theta, BivD = BivD, par2 = par2, min.pr = min.pr, max.pr = max.pr, min.dn = min.dn)$c.copula.be2
  

  
  diffh1.h2 <- c.be2.h1 - c.be2.h2  # diffh1.h2 <- mm(diffh1.h2, min.pr = min.pr, max.pr = max.pr)  

  
  if(length(y) == 2) y <- y[1]
  
  
  y^2*diffh1.h2


}









discrcont.m <- function(y, y12, eta1, eta2, sig1, sig2, nu1 = NULL, nu2 = NULL, theta, margins, BivD, par2, min.dn, min.pr, max.pr, left.trunc1, left.trunc2){           
 
y1m <- y2m <- NA
 
   
  p1pdf1 <- distrHsATDiscr(y12, eta2 = eta1, sigma2 = sig1, nu = 1, margin2 = margins[1], y2m = y1m, robust = FALSE, 
                           min.dn = min.pr, min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc1)
  p1     <- as.numeric(p1pdf1$p2)
  pmf1   <- as.numeric(p1pdf1$pdf2)   
   

  rpd2 <- distrHsAT(y,   eta2 = eta2, sigma2 = sig2, nu = nu2, margin2 = margins[2], min.dn = min.dn, min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc2)     
  p2   <- rpd2$p2
  pdf2 <- rpd2$pdf2

  
  
  c.be2.h1 <- copgHs2(p1 = p1, p2 = p2, teta = theta, teta.st = theta, BivD = BivD, par2 = par2, min.pr = min.pr, max.pr = max.pr, min.dn = min.dn)$c.copula.be2
  c.be2.h2 <- copgHs2(p1 = mm(p1 - pmf1, min.pr = min.pr, max.pr = max.pr), p2 = p2, teta = theta, teta.st = theta, BivD = BivD, par2 = par2, min.pr = min.pr, max.pr = max.pr, min.dn = min.dn)$c.copula.be2
  

  
  diffh1.h2 <- c.be2.h1 - c.be2.h2  # diffh1.h2 <- mm(diffh1.h2, min.pr = min.pr, max.pr = max.pr)  

  
  y*pdf2*diffh1.h2/pmf1

}


discrcont.v <- function(y, y12, eta1, eta2, sig1, sig2, nu1 = NULL, nu2 = NULL, theta, margins, BivD, par2, min.dn, min.pr, max.pr, left.trunc1, left.trunc2){           
 
y1m <- y2m <- NA
 
   
  p1pdf1 <- distrHsATDiscr(y12, eta2 = eta1, sigma2 = sig1, nu = 1, margin2 = margins[1], y2m = y1m, robust = FALSE, 
                           min.dn = min.pr, min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc1)
  p1     <- as.numeric(p1pdf1$p2)
  pmf1   <- as.numeric(p1pdf1$pdf2)   
   

  rpd2 <- distrHsAT(y,   eta2 = eta2, sigma2 = sig2, nu = nu2, margin2 = margins[2], min.dn = min.dn, min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc2)     
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
                                             


if( length(y12) == 2 ) y12 <- y12[1]

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


if( length(y12) == 2 ) y12 <- y12[1]
                                             
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


if(margins[2] == "BB") y <- y[1]

                                             
              } 
 
 
if(cond == 2){ 

p1pdf1 <- distrHsATDiscr(y,   eta2 = eta1, sigma2 = sig1, nu = 1, margin2 = margins[1], y2m = y1m, robust = FALSE, 
                         min.dn = min.pr, min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc1)
p1     <- as.numeric(p1pdf1$p2)
pmf1   <- as.numeric(p1pdf1$pdf2)

if(margins[1] == "BB") y <- y[1]

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


contcont.m <- function(y, y12, eta1, eta2, sig1, sig2, nu1 = NULL, nu2 = NULL, theta, margins, BivD, par2, min.dn, min.pr, max.pr, cond, left.trunc1, left.trunc2){           
     
     
 if(cond == 1){    
     
 rpd1 <- distrHsAT(y12, eta2 = eta1, sigma2 = sig1, nu = nu1, margin2 = margins[1], min.dn = min.dn, min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc1)
 rpd2 <- distrHsAT(y,   eta2 = eta2, sigma2 = sig2, nu = nu2, margin2 = margins[2], min.dn = min.dn, min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc2)     

              }
              
 if(cond == 2){    
     
 rpd1 <- distrHsAT(y,   eta2 = eta1, sigma2 = sig1, nu = nu1, margin2 = margins[1], min.dn = min.dn, min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc1)
 rpd2 <- distrHsAT(y12, eta2 = eta2, sigma2 = sig2, nu = nu2, margin2 = margins[2], min.dn = min.dn, min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc2)     
 
              }             

 pdf1 <- rpd1$pdf2
 p1   <- rpd1$p2
 
 pdf2 <- rpd2$pdf2
 p2   <- rpd2$p2 
 
 c.be1be2 <- copgHs2(p1 = p1, p2 = p2, teta = theta, teta.st = theta, BivD = BivD, par2 = par2, min.pr = min.pr, max.pr = max.pr, min.dn = min.dn)$c.copula2.be1be2

 if(cond == 1) y*c.be1be2*pdf2 else y*c.be1be2*pdf1 

}


contcont.v <- function(y, y12, eta1, eta2, sig1, sig2, nu1 = NULL, nu2 = NULL, theta, margins, BivD, par2, min.dn, min.pr, max.pr, cond, left.trunc1, left.trunc2){           
     
     
 if(cond == 1){    
     
 rpd1 <- distrHsAT(y12, eta2 = eta1, sigma2 = sig1, nu = nu1, margin2 = margins[1], min.dn = min.dn, min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc1)
 rpd2 <- distrHsAT(y,   eta2 = eta2, sigma2 = sig2, nu = nu2, margin2 = margins[2], min.dn = min.dn, min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc2)     

              }
              
 if(cond == 2){    
     
 rpd1 <- distrHsAT(y,   eta2 = eta1, sigma2 = sig1, nu = nu1, margin2 = margins[1], min.dn = min.dn, min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc1)
 rpd2 <- distrHsAT(y12, eta2 = eta2, sigma2 = sig2, nu = nu2, margin2 = margins[2], min.dn = min.dn, min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc2)     
 
              }             

 pdf1 <- rpd1$pdf2
 p1   <- rpd1$p2
 
 pdf2 <- rpd2$pdf2
 p2   <- rpd2$p2 
 
 c.be1be2 <- copgHs2(p1 = p1, p2 = p2, teta = theta, teta.st = theta, BivD = BivD, par2 = par2, min.pr = min.pr, max.pr = max.pr, min.dn = min.dn)$c.copula2.be1be2

 if(cond == 1) y^2*c.be1be2*pdf2 else y^2*c.be1be2*pdf1 

}









bincont.m <- function(y, p0, eta.mu, sigma, nu, theta, margin, BivD, par2, min.dn, min.pr, max.pr, y1, left.trunc2){           
     
 rpd <- distrHsAT(y, eta2 = eta.mu, sigma2 = sigma, nu = nu, margin2 = margin, min.dn = min.dn, min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc2)     
     
 p2 <- rpd$p2; pdf2 <- rpd$pdf2
 
 p0 <- rep(p0, length(p2))
 c.be2 <- copgHs2(p1 = p0, p2 = p2, teta = theta, teta.st = theta, BivD = BivD, par2 = par2, min.pr = min.pr, max.pr = max.pr, min.dn = min.dn)$c.copula.be2

 if(y1 == 1) y*(1 - c.be2)*pdf2/(1-p0) else y*(c.be2*pdf2/p0) 
 
}







bincont.v <- function(y, p0, eta.mu, sigma, nu, theta, margin, BivD, par2, min.dn, min.pr, max.pr, y1, left.trunc2){           
     
 rpd <- distrHsAT(y, eta2 = eta.mu, sigma2 = sigma, nu = nu, margin2 = margin, min.dn = min.dn, min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc2)     
     
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
                                       
 if(discr == FALSE) pdf <- distrHsAT(y, eta, sigma2, nu, margin, min.dn = min.dn, min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc)$pdf2
 if(discr == TRUE)  pdf <- distrHsDiscr(y, eta, sigma2, 1, 1, 1, margin, naive = TRUE, ym, min.dn = min.dn, 
                                        min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc)$pdf2

 log( 1 + exp( log( pdf ) + rc ) ) 
   
}










##################################################################




gradBbit1 <- function(y, eta, sigma2, sigma2.st, nu, nu.st, margin, rc, min.dn, min.pr, max.pr, discr = FALSE, ym = NULL, left.trunc = 0){ 

 

       if(discr == FALSE) dHs <- distrHs(y, eta, sigma2, sigma2.st, nu, nu.st, margin, naive = TRUE, min.dn = min.dn, min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc)
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



       if(discr == FALSE) dHs <- distrHs(y, eta, sigma2, sigma2.st, nu, nu.st, margin, naive = TRUE, min.dn = min.dn, min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc)
       if(discr == TRUE)  dHs <- distrHsDiscr(y, eta, sigma2, sigma2.st, 1, 1, margin, naive = TRUE, ym, min.dn = min.dn, 
                                              min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc)

       pdf2                 <- dHs$pdf2
       derpdf2.dersigma2.st <- dHs$derpdf2.dersigma2.st  

       comp1 <- 1 + exp(log( pdf2 ) + rc) 
       comp2 <- pdf2/comp1

       dl.dsigma.st <- derpdf2.dersigma2.st/pdf2

      comp2*dl.dsigma.st 
       
}




gradBbit3 <- function(y, eta, sigma2, sigma2.st, nu, nu.st, margin, rc, min.dn, min.pr, max.pr, discr = FALSE, ym = NULL, left.trunc = 0){ 

 

       dHs <- distrHs(y, eta, sigma2, sigma2.st, nu, nu.st, margin, naive = TRUE, min.dn = min.dn, min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc)

       pdf2                 <- dHs$pdf2
       derpdf2.dernu.st     <- dHs$derpdf2.dernu.st  

       comp1 <- 1 + exp(log( pdf2 ) + rc) 
       comp2 <- pdf2/comp1

       dl.dnu.st <- derpdf2.dernu.st/pdf2

       comp2*dl.dnu.st 
       
       }
       
       
       
       
       
hessBbit1 <- function(y, eta, sigma2, sigma2.st, nu, nu.st, margin, rc, min.dn, min.pr, max.pr, discr = FALSE, ym = NULL, left.trunc = 0){ 



       if(discr == FALSE) dHs <- distrHs(y, eta, sigma2, sigma2.st, nu, nu.st, margin, naive = TRUE, min.dn = min.dn, min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc)
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



       if(discr == FALSE) dHs <- distrHs(y, eta, sigma2, sigma2.st, nu, nu.st, margin, naive = TRUE, min.dn = min.dn, min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc)
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



       if(discr == FALSE) dHs <- distrHs(y, eta, sigma2, sigma2.st, nu, nu.st, margin, naive = TRUE, min.dn = min.dn, min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc)
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




hessBbit4 <- function(y, eta, sigma2, sigma2.st, nu, nu.st, margin, rc, min.dn, min.pr, max.pr, discr = FALSE, ym = NULL, left.trunc = 0){ 

 

        dHs <- distrHs(y, eta, sigma2, sigma2.st, nu, nu.st, margin, naive = TRUE, min.dn = min.dn, min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc)

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




hessBbit5 <- function(y, eta, sigma2, sigma2.st, nu, nu.st, margin, rc, min.dn, min.pr, max.pr, discr = FALSE, ym = NULL, left.trunc = 0){ 

 

        dHs <- distrHs(y, eta, sigma2, sigma2.st, nu, nu.st, margin, naive = TRUE, min.dn = min.dn, min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc)

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




hessBbit6 <- function(y, eta, sigma2, sigma2.st, nu, nu.st, margin, rc, min.dn, min.pr, max.pr, discr = FALSE, ym = NULL, left.trunc = 0){ 



        dHs <- distrHs(y, eta, sigma2, sigma2.st, nu, nu.st, margin, naive = TRUE, min.dn = min.dn, min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc)

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

int1f <- function(y, eta, sigma2, nu, margin, rc, min.dn, min.pr, max.pr, left.trunc = 0){ 



   pdf <- distrHsAT(y, eta, sigma2, nu, margin, min.dn = min.dn, min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc)$pdf2
   
   log( 1 + exp( log( pdf ) + rc ) )
   
}




d.bpsi <- function(y, X1, X2, X3, eta, sigma2, sigma2.st, nu, nu.st, margin, rc, j, min.dn, min.pr, max.pr, left.trunc = 0){ 


       dHs <- distrHs(y, eta, sigma2, sigma2.st, nu, nu.st, margin, naive = TRUE, min.dn = min.dn, min.pr = min.pr, max.pr = max.pr, left.trunc = left.trunc)

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
                                    
  for(j in 1:length(params)) G[i, j] <- integrate(d.bpsi, lB, uB, X1, X2, X3, eta, sigma2, sigma2.st, nu, nu.st, margin, rc, j, min.dn, min.pr, max.pr, left.trunc = VC$left.trunc)$value

             }
             
  colSums(G)
  
}


