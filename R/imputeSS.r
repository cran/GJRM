imputeSS <- function(x, m = 10){

if(x$VC$Cont != "NO" && x$VC$ccss != "yes") stop("This function can only be used for a selection model.")

epsilon <- sqrt(.Machine$double.eps)

posit    <- c(which(x$y1 == 0))
vec.mi   <- list()
yst.aver <- mean(x$y2) 
margin   <- x$margins[2]

if(margin %in% c("TW") ) stop("Tweedie not tested. Get in touch for details.") 


#if(margin %in% c("PO","ZTP","NBI","NBII","NBIa","NBIIa","PIG","probit","logit","cloglog","DGP","DGPII") ) stop("Check next release for the tested versions of these develoments.")


# tw not tested, sqrt would be better

if(margin %in% c("LN","WEI","iG","GA","GAi","DAGUM","SM","FISK","GP","GPII","GPo","TW") ) {yst.aver <- log(yst.aver); if(yst.aver == "-Inf") yst.aver <- log(1e-14) } 
if(margin %in% c("BE") )                                                {yst.aver <- qlogis(mm(yst.aver)) }


#########################
#########################

sim.beta.y <- function(x, s){ 

sigma2 <- nu <- NULL

l1 <- length(x$gam1$coefficients)
l2 <- length(x$gam2$coefficients)
l3 <- length(x$gam3$coefficients)
l4 <- length(x$gam4$coefficients)
l5 <- length(x$gam5$coefficients)

betahatSim <- rMVN(1, x$coefficients, x$Vb)
 
p1 <- as.numeric( 1 - probm(x$X1[s, ]%*%betahatSim[1:l1], x$margins[1])$pr )

eta2 <- eta.tr( as.numeric( x$X2s[s, ]%*%betahatSim[(l1 + 1):(l1 + l2)] ), x$margins[2])



if(is.null(x$X3)){


if(x$margins[2] %in% c(x$VC$m1d, x$VC$bl)) teta <- teta.tr(x$VC, betahatSim[l1 + l2 + 1])$teta 


if(x$margins[2] %in% c(x$VC$m2,x$VC$m2d)){ 

sigma2 <- esp.tr(betahatSim[l1 + l2 + 1], x$margins[2])$vrb 
teta   <- teta.tr(x$VC, betahatSim[l1 + l2 + 2])$teta 

                                         }
                                         
if(x$margins[2] %in% c(x$VC$m3)){

sigma2 <- esp.tr(betahatSim[l1 + l2 + 1], x$margins[2])$vrb 
nu     <- enu.tr(betahatSim[l1 + l2 + 2], x$margins[2])$vrb 
teta   <- teta.tr(x$VC, betahatSim[l1 + l2 + 3])$teta 

                                }
                                
}


if(!is.null(x$X3)){


if(x$margins[2] %in% c(x$VC$m1d, x$VC$bl)) teta <- teta.tr(x$VC, x$X3s[s, ]%*%betahatSim[(l1 + l2 + 1):(l1 + l2 + l3)] )$teta


if(x$margins[2] %in% c(x$VC$m2,x$VC$m2d)){ 

sigma2 <- esp.tr(x$X3s[s, ]%*%betahatSim[(l1 + l2 + 1):(l1 + l2 + l3)], x$margins[2])$vrb 
teta   <- teta.tr(x$VC, x$X4s[s, ]%*%betahatSim[(l1 + l2 + l3 + 1):(l1 + l2 + l3 + l4)])$teta 

                                         }
                                         
if(x$margins[2] %in% c(x$VC$m3)){

sigma2 <- esp.tr(x$X3s[s, ]%*%betahatSim[(l1 + l2 + 1):(l1 + l2 + l3)], x$margins[2])$vrb 
nu     <- enu.tr(x$X4s[s, ]%*%betahatSim[(l1 + l2 + l3 + 1):(l1 + l2 + l3 + l4)], x$margins[2])$vrb 
teta   <- teta.tr(x$VC, x$X5s[s, ]%*%betahatSim[(l1 + l2 + l3 + l4 + 1):(l1 + l2 + l3 + l4 + l5)])$teta 

                                }
                                
}


y2 <- sim.resp(x$margins[2], 1, eta2, sigma2, nu, setseed = FALSE)

list(p1 = p1, eta2 = eta2, sigma2 = sigma2, nu = nu, teta = teta, y2 = y2, margin = x$margins[2], BivD = x$BivD, sigma = sigma2)


}


#########################

obj.grad.hess <- function(y2.st, out.beta){ 

p1     <- out.beta$p1
eta2   <- out.beta$eta2
sigma2 <- out.beta$sigma2 
nu     <- out.beta$nu 
teta   <- out.beta$teta  

margin <- out.beta$margin
BivD   <- out.beta$BivD

y2 <- y2.st 


# tw not tested, sqrt would be better

if(margin %in% c("LN","WEI","iG","GA","GAi","DAGUM","SM","FISK","GP","GPII","GPo","TW") ) y2 <- esp.tr(y2.st, "LN")$vrb 
if(margin %in% c("BE") )                                           y2 <- esp.tr(y2.st, "BE")$vrb

ppdf <- distrHsAT(y2, eta2, sigma2, nu, margin)

p2          <- ppdf$p2
derp2.dery2 <- ppdf$pdf2 

c.copula.be2 <- copgHsAT(p1, p2, teta, BivD, Ln = FALSE, par2 = x$dof)$c.copula.be2

cc12          <- copgHs3(p1, p2, eta1 = NULL, eta2 = NULL, teta, teta.st = NULL, BivD, par2 = x$dof)
c.copula2.be2 <- cc12$c.copula2.be2 
der2h.derp2p2 <- cc12$der2h.derp2p2


ddd <- distrHsAT1(y2.st, eta2, sigma2, nu, margin)


dery.dery.st   <- ddd$dery.dery.st    
der2y.dery.st2 <- ddd$der2y.dery.st2  
der2p2.dery22  <- ddd$der2p2.dery22



f.g  <- c.copula.be2/p1 
gf.g <- c.copula2.be2*derp2.dery2*dery.dery.st/p1 
hf.g <- (der2h.derp2p2 * derp2.dery2^2 + c.copula2.be2 * der2p2.dery22)/p1 * dery.dery.st^2 + c.copula2.be2*derp2.dery2/p1*der2y.dery.st2

list(value = f.g, gradient = gf.g, hessian = as.matrix(hf.g))

}


#########################

obj.Discr <- function(y2.st, out.beta){ 

p1     <- out.beta$p1
eta2   <- out.beta$eta2
sigma2 <- out.beta$sigma2 
nu     <- out.beta$nu 
teta   <- out.beta$teta  

margin <- out.beta$margin
BivD   <- out.beta$BivD

y2  <- y2.st 


if( margin == "ZTP"){

  mu2 <- c(exp(eta2))
     
     if( y2 > 170){ 
                    prec <- pmax(53, getPrec(mu2), getPrec(y2))         
                    mu2 <- mpfr(mu2, prec)
                    y2  <- mpfr( y2, prec) 
                   }
                   
  y2m <- seq(1, y2)                 
  pdf2FUNC <- function(y2, mu2) mu2^y2/(exp(mu2)-1)*1/factorial(y2)  

  p2 <- sum( as.numeric( pdf2FUNC(y2m, mu2) ) ) 

}


if( margin %in% c("DGP","DGPII")){


if( margin == "DGP")   mu2    <- c(eta2)
if( margin == "DGPII") mu2    <- c(eta2^2)


sigma2 <- c(sigma2)
y2m <- seq(1, y2)                 


# sigma2 is sigma
                                         
  pdf2FUNC2 <- function(y2, mu2, sigma2) (1 + mu2*y2/sigma2)^(-1/mu2)    - (1 + mu2*(1+y2)/sigma2)^(-1/mu2)  

  p2  <- sum( as.numeric( pdf2FUNC2(y2m, mu2, sigma2) )  ) 


}


ppd <- distrHsATDiscr(y2, eta2, sigma2, nu, margin, y2m = NULL, robust = FALSE) 
if( margin %in% c("ZTP","DGP","DGPII")) ppd$p2 <- p2

# test if it works for p2 when margin is ZTP


p2 <- ppd$p2
pdf2 <- ppd$pdf2

C1 <- BiCDF(p1, p2, x$VC$nC, teta, x$VC$dof)
C2 <- BiCDF(p1, mm(p2 - pdf2), x$VC$nC, teta, x$VC$dof)
  
A <- ifelse( (C1 - C2) < epsilon, epsilon, C1 - C2)

f.g <- A/(p1*pdf2) 

list(value = f.g)

}


#########################

obj.Bin <- function(y2.st, out.beta){ 

p1     <- out.beta$p1
eta2   <- out.beta$eta2
sigma2 <- out.beta$sigma2 
nu     <- out.beta$nu 
teta   <- out.beta$teta  

margin <- out.beta$margin
BivD   <- out.beta$BivD

y2 <- y2.st 


if(y2 == 1){

# recall that here 1-p1=P(Y_1=1) and that p1=P(Y_1=0)

p2 <- probm(eta2, margin)$pr 

p11 <- BiCDF(1-p1, p2, x$VC$nC, teta, x$VC$dof)
p01 <- pmax(p2 - p11, epsilon)

f.g <- p01/(p1*p2)

            }


if(y2 == 0){

p2 <- 1 - probm(eta2, margin)$pr 

p11 <- BiCDF(1-p1, 1-p2, x$VC$nC, teta, x$VC$dof)
p01 <- pmax((1-p2) - p11, epsilon)
p10 <- pmax((1-p1) - p11, epsilon)
p00 <- pmax(1 - p11 - p10 - p01, epsilon)

f.g <- p00/(p1*p2)

            }

list(value = f.g)


}



#########################
#########################






for(i in 1:m){ ###

y2.imp <- NA; j <- 1


      for(s in posit){ ## check trust - next, skip s?
      
repeat{ #
   out.beta <- sim.beta.y(x, s)
   y2.imp[j] <- y2.st <- out.beta$y2 
   

   if(margin %in% c("LN","WEI","iG","GA","GAi","DAGUM","SM","FISK","GP","GPII","GPo","TW") ) {y2.st <- log(y2.imp[j]); if(y2.st == "-Inf") y2.st <- log(1e-14)} 
   if(margin %in% c("BE") )                                           {y2.st <- qlogis(mm(y2.imp[j]))}
 


if( margin %in% c("probit","logit","cloglog") ){
    
   f.g     <- obj.Bin(y2.st, out.beta)$value  
   M.value <- max(c(obj.Bin(0, out.beta)$value , obj.Bin(1, out.beta)$value) )    
   
                                                } 
 

if( margin %in% c("PO","ZTP","NBI","NBII","NBIa","NBIIa","PIG","DGP","DGPII") ){
   
   
   test.oD <- NA
   max.y   <- max(x$y2)  +  ( max(x$y2) - min(x$y2) )/2
   seq.y   <- 0:max.y
   
   f.g     <- obj.Discr(y2.st, out.beta)$value
  
   for(jj in 1:length(seq.y)) test.oD[jj] <- obj.Discr(seq.y[jj], out.beta)$value  
  
   M.value <- max(test.oD)    
   
                                                                   }



if(!(margin %in% c("PO","ZTP","NBI","NBII","NBIa","NBIIa","PIG","probit","logit","cloglog","DGP","DGPII")) ){
   
   f.g      <- obj.grad.hess(y2.st, out.beta)$value
   M.value  <- try(trust(obj.grad.hess, parinit = yst.aver, rinit = 1, rmax = 100, minimize = F, out.beta = out.beta)$value); if ('try-error' %in% class(M.value)) next
   
                                                                }



   if(runif(1) < f.g/M.value){ j <- j + 1; break}
   
      } #
      
                     } ##


                                                                                                                           
vec.mi[[i]] <- y2.imp
print(i)

             } #### 


vec.mi

} # end
