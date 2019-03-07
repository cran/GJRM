imputeCounter <- function(x, m = 10, nm.end){

# if(!(x$VC$Cont != "NO" && x$VC$ccss != "yes")) stop("This function can only be used for a selection model.")
# put appropriate checks

if(missing(nm.end)) stop("You must provide the name of the endogenous variable.")
if(!(x$margins[1] %in% x$bl) ) stop("First equation must be for the binary response.")


epsilon <- 1e-07

posit0 <- c(which(x$y1 == 0))
posit1 <- c(which(x$y1 == 1))

vec.mi0 <- vec.mi1 <- Res <- list()

imp.vec <- rep(list(NA),m)


yst.aver <- mean(x$y2) 
margin   <- x$margins[2]
y1.c <- AT <- NA

if(margin %in% c("LN","WEI","iG","GA","GAi","DAGUM","SM","FISK","GP") ) {yst.aver <- log(yst.aver); if(yst.aver == "-Inf") yst.aver <- log(1e-14) } 
if(margin %in% c("BE") )                                                 {yst.aver <- qlogis(mm(yst.aver)) }


#########################
#########################

sim.beta.y <- function(x, s, zo){ 

p1 <- p0 <- eta2 <- eta2.0 <- eta2.1 <- sigma2 <- sigma2.0 <- sigma2.1 <- nu <- nu.0 <- nu.1 <- y2.0 <- y2.1 <- NULL

l1 <- length(x$gam1$coefficients)
l2 <- length(x$gam2$coefficients)
l3 <- length(x$gam3$coefficients)
l4 <- length(x$gam4$coefficients)
l5 <- length(x$gam5$coefficients)


betahatSim <- rMVN(1, x$coefficients, x$Vb)
 
 
 
p1 <- as.numeric( probm(x$X1[s, ]%*%betahatSim[1:l1], x$margins[1])$pr )
p0 <- as.numeric( 1 - p1 )




X2 <- as.matrix(x$X2)

if( nm.end %in% all.vars(x$formula[[2]][[3]]) ){

if(zo == 0) X2[, nm.end] <- 0
if(zo == 1) X2[, nm.end] <- 1

}

eta2 <- eta.tr( as.numeric( X2[s, ]%*%betahatSim[(l1 + 1):(l1 + l2)] ), x$margins[2])





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


if(x$margins[2] %in% c(x$VC$m1d, x$VC$bl)) teta <- teta.tr(x$VC, x$X3[s, ]%*%betahatSim[(l1 + l2 + 1):(l1 + l2 + l3)] )$teta




if(x$margins[2] %in% c(x$VC$m2,x$VC$m2d)){ 



X3 <- as.matrix(x$X3)

if( nm.end %in% all.vars(x$formula[[3]][[3]]) ){

  if(zo == 0) X3[, nm.end] <- 0
  if(zo == 1) X3[, nm.end] <- 1

                                               }

sigma2 <- esp.tr(X3[s, ]%*%betahatSim[(l1 + l2 + 1):(l1 + l2 + l3)], x$margins[2])$vrb 

teta   <- teta.tr(x$VC, x$X4[s, ]%*%betahatSim[(l1 + l2 + l3 + 1):(l1 + l2 + l3 + l4)])$teta 

                                         }
                                         
                                         
                                         
                                         
                                         
if(x$margins[2] %in% c(x$VC$m3)){


X3 <- as.matrix(x$X3)


if( nm.end %in% all.vars(x$formula[[3]][[3]]) ){

if(zo == 0) X3[, nm.end] <- 0
if(zo == 1) X3[, nm.end] <- 1
                                               }

sigma2 <- esp.tr(X3[s, ]%*%betahatSim[(l1 + l2 + 1):(l1 + l2 + l3)], x$margins[2])$vrb 




X4 <- as.matrix(x$X4)

if( nm.end %in% all.vars(x$formula[[4]][[3]]) ){

if(zo == 0) X4[, nm.end] <- 0
if(zo == 1) X4[, nm.end] <- 1

                                               }

nu   <- enu.tr(X4[s, ]%*%betahatSim[(l1 + l2 + l3 + 1):(l1 + l2 + l3 + l4)], x$margins[2])$vrb 



teta <- teta.tr(x$VC, x$X5[s, ]%*%betahatSim[(l1 + l2 + l3 + l4 + 1):(l1 + l2 + l3 + l4 + l5)])$teta 

                                }
                                
}


y2 <- sim.resp(x$margins[2], 1, eta2, sigma2, nu, setseed = FALSE)




list(p0 = p0, p1 = p1, eta2 = eta2, sigma2 = sigma2, nu = nu, 
     teta = teta, y2 = y2, margin = x$margins[2], BivD = x$BivD)






}


#########################

obj.grad.hess <- function(y2.st, out.beta, zo){ 

eta2   <- out.beta$eta2
sigma2 <- out.beta$sigma2
nu     <- out.beta$nu
teta   <- out.beta$teta  


if(zo == 0) ppr <- out.beta$p0 
if(zo == 1) ppr <- out.beta$p1

margin <- out.beta$margin
BivD   <- out.beta$BivD

y2 <- y2.st 


if(margin %in% c("LN","WEI","iG","GA","GAi","DAGUM","SM","FISK","GP") ) y2 <- esp.tr(y2.st, "LN")$vrb 
if(margin %in% c("BE") )                                           y2 <- esp.tr(y2.st, "BE")$vrb

ppdf <- distrHsAT(y2, eta2, sigma2, nu, margin)

p2          <- ppdf$p2
derp2.dery2 <- ppdf$pdf2 

c.copula.be2 <- copgHsAT(ppr, p2, teta, BivD, Ln = FALSE, par2 = x$dof)$c.copula.be2

cc12          <- copgHs3(ppr, p2, eta1 = NULL, eta2 = NULL, teta, teta.st = NULL, BivD, par2 = x$dof)
c.copula2.be2 <- cc12$c.copula2.be2 
der2h.derp2p2 <- cc12$der2h.derp2p2


ddd <- distrHsAT1(y2.st, eta2, sigma2, nu, margin)


dery.dery.st   <- ddd$dery.dery.st    
der2y.dery.st2 <- ddd$der2y.dery.st2  
der2p2.dery22  <- ddd$der2p2.dery22



f.g  <- c.copula.be2/ppr 
gf.g <- c.copula2.be2*derp2.dery2*dery.dery.st/ppr 
hf.g <- (der2h.derp2p2 * derp2.dery2^2 + c.copula2.be2 * der2p2.dery22)/ppr * dery.dery.st^2 + c.copula2.be2*derp2.dery2/ppr*der2y.dery.st2


list(value = f.g, gradient = gf.g, hessian = as.matrix(hf.g))

}


#########################

obj.Discr <- function(y2.st, out.beta, zo){ 

eta2   <- out.beta$eta2
sigma2 <- out.beta$sigma2 
nu     <- out.beta$nu 
teta   <- out.beta$teta  

if(zo == 0) ppr <- out.beta$p0 
if(zo == 1) ppr <- out.beta$p1


margin <- out.beta$margin
BivD   <- out.beta$BivD

y2 <- y2.st 

ppd <- distrHsATDiscr(y2, eta2, sigma2, nu, margin, y2m = NULL, robust = FALSE) 
# p2 not available for ZTP unless we supply y2m and we set robust = TRUE

p2 <- ppd$p2
pdf2 <- ppd$pdf2

C1 <- BiCDF(ppr, p2, x$VC$nC, teta, x$VC$dof)
C2 <- BiCDF(ppr, mm(p2 - pdf2), x$VC$nC, teta, x$VC$dof)
  
A <- ifelse( (C1 - C2) < epsilon, epsilon, C1 - C2)

f.g <- A/(ppr*pdf2)

list(value = f.g)

}


#########################

obj.Bin <- function(y2.st, out.beta, zo){ 

if(zo == 0) ppr <- out.beta$p0 
if(zo == 1) ppr <- out.beta$p1

eta2   <- out.beta$eta2
sigma2 <- out.beta$sigma2 
nu     <- out.beta$nu 
teta   <- out.beta$teta  

margin <- out.beta$margin
BivD   <- out.beta$BivD

y2 <- y2.st 


if(zo == 0 && y2 == 1){


p2 <- probm(eta2, margin)$pr 

p11 <- BiCDF(1 - ppr, p2, x$VC$nC, teta, x$VC$dof)
p01 <- pmax(p2 - p11, epsilon)

f.g <- p01/(ppr*p2)


}



if(zo == 0 && y2 == 0){ 

p2 <- 1 - probm(eta2, margin)$pr 

p11 <- BiCDF(1-ppr, 1-p2, x$VC$nC, teta, x$VC$dof)
p01 <- pmax((1-p2) - p11, epsilon)
p10 <- pmax((1 - ppr) - p11, epsilon)
p00 <- pmax(1 - p11 - p10 - p01, epsilon)

f.g <- p00/(ppr*p2)

}




if(zo == 1 && y2 == 1){


p2 <- probm(eta2, margin)$pr 

p11 <- BiCDF(ppr, p2, x$VC$nC, teta, x$VC$dof)

f.g <- p11/(ppr*p2)


}



if(zo == 1 && y2 == 0){ 

p2 <- 1 - probm(eta2, margin)$pr 

p11 <- BiCDF(ppr, 1-p2, x$VC$nC, teta, x$VC$dof)
p10 <- pmax(ppr - p11, epsilon)

f.g <- p10/(ppr*p2)

}




list(value = f.g)


}



#########################
#########################






for(i in 1:m){ ####

y2.imp <- NA; j <- 1


      for(s in posit0){ ## check trust - next, skip s?
      
repeat{ #
   out.beta <- sim.beta.y(x, s, zo = 1)
   y2.imp[j] <- y2.st <- out.beta$y2 
   

   if(margin %in% c("LN","WEI","iG","GA","GAi","DAGUM","SM","FISK","GP") ) {y2.st <- log(y2.imp[j]); if(y2.st == "-Inf") y2.st <- log(1e-14)} 
   if(margin %in% c("BE") )                                           {y2.st <- qlogis(mm(y2.imp[j]))}
 


if( margin %in% c("probit","logit","cloglog") ){
    
   f.g     <- obj.Bin(y2.st, out.beta, zo = 1)$value  
   M.value <- max(c(obj.Bin(0, out.beta, zo = 1)$value , obj.Bin(1, out.beta, zo = 1)$value) )    
   
                                                } 
 


if( margin %in% c("PO","ZTP","NBI","NBII","NBIa","NBIIa","PIG","DGP") ){
   
   
   test.oD <- NA
   max.y   <- max(x$y2)  +  ( max(x$y2) - min(x$y2) )/2
   seq.y   <- 0:max.y
   
   f.g     <- obj.Discr(y2.st, out.beta, zo = 1)$value
  
   for(jj in 1:length(seq.y)) test.oD[jj] <- obj.Discr(seq.y[jj], out.beta, zo = 1)$value  
  
   M.value <- max(test.oD)    
   
                                                                   }



if(!(margin %in% c("DGP","PO","ZTP","NBI","NBII","NBIa","NBIIa","PIG","probit","logit","cloglog")) ){
   
   f.g      <- obj.grad.hess(y2.st, out.beta, zo = 1)$value
   M.value  <- try(trust(obj.grad.hess, parinit = yst.aver, rinit = 1, rmax = 100, minimize = F, out.beta = out.beta, zo = 1)$value); if ('try-error' %in% class(M.value)) next
   
                                                                }



   if(runif(1) < f.g/M.value){ j <- j + 1; break}
   
      } #
      
                     } ##


                                                                                                                           
vec.mi0[[i]] <- y2.imp
print(i)

             } #### 












for(i in 1:m){ ####

y2.imp <- NA; j <- 1


      for(s in posit1){ ## check trust - next, skip s?
      
repeat{ #
   out.beta <- sim.beta.y(x, s, zo = 0)
   y2.imp[j] <- y2.st <- out.beta$y2 
   

   if(margin %in% c("LN","WEI","iG","GA","GAi","DAGUM","SM","FISK","GP") ) {y2.st <- log(y2.imp[j]); if(y2.st == "-Inf") y2.st <- log(1e-14)} 
   if(margin %in% c("BE") )                                           {y2.st <- qlogis(mm(y2.imp[j]))}
 


if( margin %in% c("probit","logit","cloglog") ){
    
   f.g     <- obj.Bin(y2.st, out.beta, zo = 0)$value  
   M.value <- max(c(obj.Bin(0, out.beta, zo = 0)$value , obj.Bin(1, out.beta, zo = 0)$value) )    
   
                                                } 
 


if( margin %in% c("PO","ZTP","NBI","NBII","NBIa","NBIIa","PIG","DGP") ){
   
   
   test.oD <- NA
   max.y   <- max(x$y2)  +  ( max(x$y2) - min(x$y2) )/2
   seq.y   <- 0:max.y
   
   f.g     <- obj.Discr(y2.st, out.beta, zo = 0)$value
  
   for(jj in 1:length(seq.y)) test.oD[jj] <- obj.Discr(seq.y[jj], out.beta, zo = 0)$value  
  
   M.value <- max(test.oD)    
   
                                                                   }



if(!(margin %in% c("PO","ZTP","NBI","NBII","NBIa","NBIIa","PIG","probit","logit","cloglog","DGP")) ){
   
   f.g      <- obj.grad.hess(y2.st, out.beta, zo = 0)$value
   M.value  <- try(trust(obj.grad.hess, parinit = yst.aver, rinit = 1, rmax = 100, minimize = F, out.beta = out.beta, zo = 0)$value); if ('try-error' %in% class(M.value)) next
   
                                                                }



   if(runif(1) < f.g/M.value){ j <- j + 1; break}
   
      } #
      
                     } ##


                                                                                                                           
vec.mi1[[i]] <- y2.imp
print(i)

             } #### 
             
             
             
# assemble results






for(i in 1:m){ imp.vec[[i]][posit0] <- vec.mi0[[i]]; imp.vec[[i]][posit1] <- vec.mi1[[i]] }  


y1 <- x$y1
y2 <- x$y2


for(i in 1:m) Res[[i]] <- data.frame( y2 = y2, y2.c = imp.vec[[i]] )



for(j in 1:m){

  for(i in 1:x$n){

   if(y1[i] == 1) Res[[j]][i, c(1,2)] <- Res[[j]][i, c(2,1)] 

                 }
                 
                 
   AT[j] <- mean(Res[[j]][,2] - Res[[j]][,1])               

}





list(Res = Res, AT = AT, imp.vec = imp.vec, posit0 = posit0, posit1 = posit1, vec.mi0 = vec.mi0, vec.mi1 = vec.mi1)





} # end
