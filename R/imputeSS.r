imputeSS <- function(x, m){


if(x$VC$Cont != "NO" && x$VC$ccss != "yes") stop("This function can only be used with a copulaSampleSel object.")

if(x$margins[2] %in% c("NBI", "NBII", "PIG", "PO", "ZTP")) stop("This function is not currently suitable for discrete margins.")

sim.beta.y <- function(x, s){

nu <- NULL

l1 <- length(coef(x$gam1))
l2 <- length(coef(x$gam2))
l3 <- length(coef(x$gam3))
l4 <- length(coef(x$gam4))
l5 <- length(coef(x$gam5))

betahatSim <- rMVN(1, x$coefficients, x$Vb)
 
p1 <- as.numeric( 1 - probm(x$X1[s, ]%*%betahatSim[1:l1], x$margins[1])$pr )

eta2 <- eta.tr( as.numeric( x$X2s[s, ]%*%betahatSim[(l1 + 1):(l1 + l2)] ), x$margins[2])


if(is.null(x$X3)){

if(x$margins[2] %in% c(x$VC$m1d,x$VC$m2,x$VC$m2d)){
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

if(x$margins[2] %in% c(x$VC$m1d,x$VC$m2,x$VC$m2d)){
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

list(p1 = p1, eta2 = eta2, sigma2 = sigma2, nu = nu, teta = teta, y2 = y2, margin = x$margins[2], BivD = x$BivD)


}




obj.grad.hess <- function(y2.st, out.beta){ 

p1     <- out.beta$p1
eta2   <- out.beta$eta2
sigma2 <- out.beta$sigma2 
nu     <- out.beta$nu 
teta   <- out.beta$teta  

margin <- out.beta$margin
BivD   <- out.beta$BivD

y2 <- y2.st 

if(margin %in% c("LN","WEI","iG","GA","GAi","DAGUM","SM","FISK") ) y2 <- esp.tr(y2.st, "LN")$vrb # it can be any distr
if(margin %in% c("BE") )                                                            y2 <- esp.tr(y2.st, "BE")$vrb


ppdf <- distrHsAT(y2, eta2, sigma2, nu, margin)
p2   <- ppdf$p2
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





posit <- c(which(x$y1 == 0))

vec.mi <- list()

yst.aver <- mean(x$y2) 

margin <- x$margins[2]

if(margin %in% c("LN","WEI","iG","GA","GAi","DAGUM","SM","FISK") ) {yst.aver <- log(yst.aver); if(yst.aver == "-Inf") yst.aver <- log(1e-14) } 
if(margin %in% c("BE") )                                           {yst.aver <- qlogis(mm(yst.aver)) }




for(i in 1:m){

y2.imp <- NA; j <- 1


      for(s in posit){ # check trust - next, skip s?
      
repeat{
   out.beta <- sim.beta.y(x, s)
   y2.imp[j] <- y2.st <- out.beta$y2 
   

   if(margin %in% c("LN","WEI","iG","GA","GAi","DAGUM","SM","FISK") ) {y2.st <- log(y2.imp[j]); if(y2.st == "-Inf") y2.st <- log(1e-14)} 
   if(margin %in% c("BE") )                                           {y2.st <- qlogis(mm(y2.imp[j]))}
 
   f.g      <- obj.grad.hess(y2.st, out.beta)$value
   M.value  <- try(trust(obj.grad.hess, parinit = yst.aver, rinit = 1, rmax = 100, minimize = F, out.beta = out.beta)$value); if ('try-error' %in% class(M.value)) next
   if(runif(1) < f.g/M.value){ j <- j + 1; break}
   
      }
      
                     }
                                                                                                                           
vec.mi[[i]] <- y2.imp
print(i)

             } # end for


vec.mi

} # end
