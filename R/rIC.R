rIC <- function(obj){


if(obj$robust == FALSE) stop("This function allows only for robust model objects.")


fit <- obj
m1d  <- c("PO", "ZTP")
m2d  <- c("NBI", "NBII","NBIa", "NBIIa","PIG","DGP","DGPII")

n      <- fit$VC$n
rc     <- fit$VC$rc
margin <- fit$VC$margins[2]
weights <- fit$VC$weights


#################################################################
#################################################################
#################################################################


if(margin %in% c(m1d, m2d) ){


###############################
# fixed quantities
###############################

ym     <- fit$VC$y1m
y      <- fit$VC$y
ygrid  <- fit$VC$ygrid
chs    <- fit$VC$chunk.size



no.splits <- ceiling(length(ygrid)*n/chs)
if(no.splits > 1) num.ind <- cut(1:n, no.splits, labels = FALSE)

###################################################################################################

if(is.null(fit$VC$X2)){fit$VC$X2 <- fit$VC$X3 <- matrix(1, n, 1); fit$VC$X2.d2 <- fit$VC$X3.d2 <- 1} 

eta <- eta.tr(fit$VC$X1%*%fit$coefficients[1:fit$VC$X1.d2], margin)

if( !(margin %in% c(m1d)) ){

  ss <- esp.tr(fit$VC$X2%*%fit$coefficients[(1+fit$VC$X1.d2):(fit$VC$X1.d2+fit$VC$X2.d2)], margin)
  sigma2    <- ss$vrb
  sigma2.st <- ss$vrb.st

} else sigma2 <- sigma2.st <- 1 

if( length(sigma2)==1){ sigma2 <- rep(sigma2, n); sigma2.st <- rep(sigma2.st, n)} 

###################################################################################################


if(no.splits > 1){

sumb <- 0
gradbit1b <- rep(0, dim(fit$VC$X1)[2]) # 
gradbit2b <- rep(0, dim(fit$VC$X2)[2]) # 

iter.split <- 1

for(i in 1:no.splits){

n.red      <- length(eta[num.ind==i])
etaL       <- rep(eta[num.ind==i], each = length(ygrid))
weightsL   <- rep(weights[num.ind==i], each = length(ygrid))
sigma2L    <- rep(sigma2[num.ind==i], each = length(ygrid))
sigma2.stL <- rep(sigma2.st[num.ind==i], each = length(ygrid))
ygridL     <- rep(ygrid, n.red)

                      X1L <- apply(as.matrix(fit$VC$X1[num.ind==i,]), 2, rep, each = length(ygrid))
if( margin %in% m2d ) X2L <- apply(as.matrix(fit$VC$X2[num.ind==i,]), 2, rep, each = length(ygrid))

####


 if( margin %in% m1d ){
 
   g1b <- c( weightsL*gradBbit1(ygridL, etaL, sigma2L, sigma2.stL, 1, 1, margin, rc, discr = TRUE, ym) )*X1L 
   G1B <- t(g1b)%*%g1b

     if(iter.split == 1) gradbit1b <- matrix(0, dim(G1B), dim(G1B))

   gradbit1b <- gradbit1b - G1B # sign minus should be ok


                     }


  if( margin %in% m2d ){
  
     
     g1b <- c( weightsL*gradBbit1(ygridL, etaL, sigma2L, sigma2.stL, 1, 1, margin, rc, discr = TRUE, ym) )*X1L 
     g2b <- c( weightsL*gradBbit2(ygridL, etaL, sigma2L, sigma2.stL, 1, 1, margin, rc, discr = TRUE, ym) )*X2L
     gb  <- cbind(g1b, g2b)
     GB  <- t(gb)%*%gb
     
       if(iter.split == 1) gradbit2b <- matrix(0, dim(GB), dim(GB))

     gradbit2b <- gradbit2b - GB
                
                       }
  
  iter.split <- iter.split + 1
  
  } 
                     
                     
}else{

###########################


etaL       <- rep(eta, each = length(ygrid))
weightsL   <- rep(weights, each = length(ygrid))
sigma2L    <- rep(sigma2, each = length(ygrid))
sigma2.stL <- rep(sigma2.st, each = length(ygrid))
ygridL     <- rep(ygrid, n)

                      X1L <- apply(fit$VC$X1, 2, rep, each = length(ygrid))
if( margin %in% m2d ) X2L <- apply(fit$VC$X2, 2, rep, each = length(ygrid))

####

 if( margin %in% m1d ){

g1b <- c( weightsL*gradBbit1(ygridL, etaL, sigma2L, sigma2.stL, 1, 1, margin, rc, discr = TRUE, ym) )*X1L

gradbit1b <- -t(g1b)%*%g1b

}


if( margin %in% m2d ){

   g1b <- c( weightsL*gradBbit1(ygridL, etaL, sigma2L, sigma2.stL, 1, 1, margin, rc, discr = TRUE, ym) )*X1L

   g2b <- c( weightsL*gradBbit2(ygridL, etaL, sigma2L, sigma2.stL, 1, 1, margin, rc, discr = TRUE, ym) )*X2L

   gb <- cbind(g1b, g2b)

   gradbit2b <- -t(gb)%*%gb
                
                     }

}




if(margin %in% fit$VC$m2d) bs <- gradbit2b else bs <- gradbit1b  


###################################################################################################

if(margin %in% m1d) g <- c(fit$fit$dl.dbe)*fit$VC$X1  

if(margin %in% m2d) g <- cbind(  c(fit$fit$dl.dbe)*fit$VC$X1, c(fit$fit$dl.dsigma.st)*fit$VC$X2 )




}


















#####################################################################
#####################################################################
#####################################################################



if(!(margin %in% c(m1d, m2d)) ){


m2sel <- c("WEI","iG","GA","BE","FISK","DAGUM","SM","TW","GP","GPII","GPo")

###############################
# fixed quantities
###############################

if(is.null(fit$VC$X2)){fit$VC$X2 <- fit$VC$X3 <- matrix(1, n, 1); fit$VC$X2.d2 <- fit$VC$X3.d2 <- 1} 

###################################################################################################
###################################################################################################

if(is.null(fit$VC$lB) && is.null(fit$VC$uB)){

if( margin %in% c("N","N2","GU","rGU","LO","LN") )             { lB <- -Inf;      uB <- Inf}
if( margin %in% c("WEI","iG","GA","DAGUM","SM","FISK","GP","GPII","GPo","TW")  ) { lB <- sqrt(.Machine$double.eps); uB <- Inf} # tw not test, 0 included?
if( margin %in% c("BE")  )                                     { lB <- sqrt(.Machine$double.eps); uB <- 0.9999999}

}else{lB <- fit$VC$lB; uB <- fit$VC$uB}


###################################################################################################
###################################################################################################

eta <- eta.tr(fit$VC$X1%*%fit$coefficients[1:fit$VC$X1.d2], margin)
ss  <- esp.tr(fit$VC$X2%*%fit$coefficients[(1+fit$VC$X1.d2):(fit$VC$X1.d2+fit$VC$X2.d2)], margin)

sigma2    <- ss$vrb
sigma2.st <- ss$vrb.st

if( margin %in% c("DAGUM","SM","TW") ){

            nus   <- enu.tr(fit$VC$X3%*%fit$coefficients[(1+fit$VC$X1.d2+fit$VC$X2.d2):(fit$VC$X1.d2+fit$VC$X2.d2+fit$VC$X3.d2)], margin)
            nu    <- nus$vrb
            nu.st <- nus$vrb.st
            
} else nu <- nu.st <- 1

###################################################################################################
###################################################################################################

if( margin %in% m2sel ){

intgrad   <- function(eta, sigma2, sigma2.st, nu, nu.st, margin, rc, lB, uB) distrExIntegrate(gradBbit1, lower = lB, upper = uB, eta = eta, sigma2 = sigma2, sigma2.st = sigma2.st, nu = nu, nu.st = nu.st, margin = margin, rc = rc)              
v.intgrad <- Vectorize(intgrad) 
intgrad   <- v.intgrad(eta = eta, sigma2 = sigma2, sigma2.st = sigma2.st, nu = nu, nu.st = nu.st, margin = margin, rc = rc, lB = lB, uB = uB)

gradbit1  <- -c(weights*intgrad)*fit$VC$X1 

}

intgrad   <- function(eta, sigma2, sigma2.st, nu, nu.st, margin, rc, lB, uB) distrExIntegrate(gradBbit2, lower = lB, upper = uB, eta = eta, sigma2 = sigma2, sigma2.st = sigma2.st, nu = nu, nu.st = nu.st, margin = margin, rc = rc)              
v.intgrad <- Vectorize(intgrad) 
intgrad   <- v.intgrad(eta = eta, sigma2 = sigma2, sigma2.st = sigma2.st, nu = nu, nu.st = nu.st, margin = margin, rc = rc, lB = lB, uB = uB)


gradbit2  <- -c(weights*intgrad)*fit$VC$X2 


if( margin %in% c("DAGUM","SM","TW") ){

intgrad   <- function(eta, sigma2, sigma2.st, nu, nu.st, margin, rc, lB, uB) distrExIntegrate(gradBbit3, lower = lB, upper = uB, eta = eta, sigma2 = sigma2, sigma2.st = sigma2.st, nu = nu, nu.st = nu.st, margin = margin, rc = rc)              
v.intgrad <- Vectorize(intgrad) 
intgrad   <- v.intgrad(eta = eta, sigma2 = sigma2, sigma2.st = sigma2.st, nu = nu, nu.st = nu.st, margin = margin, rc = rc, lB = lB, uB = uB)

gradbit3  <- -c(weights*intgrad)*fit$VC$X3

}



###################################################################################################
###################################################################################################

if(margin %in% fit$VC$m2){


   if(margin %in% m2sel){
   
   bp <- cbind(gradbit1, gradbit2) 
   
   #bs <- matrix(0, 21, 21) #for(i in 1:n)  bs <- bs + t(t(bp[i,]))%*%bp[i,]
   
   
   
                        }
   
   
   if(!(margin %in% m2sel)){
   
   gradbit1 <- matrix(0, nrow = n, ncol = length(fit$coefficients[1:fit$VC$X1.d2]))
   bp <- cbind(gradbit1, gradbit2) 

                          }
      
}





if(margin %in% fit$VC$m3) bp <- cbind(gradbit1, gradbit2, gradbit3) 




   bs <- t(bp)%*%bp



#############


if(!(margin %in% fit$VC$m3)) g <- cbind(  c(fit$fit$dl.dbe)*fit$VC$X1, c(fit$fit$dl.dsigma.st)*fit$VC$X2 ) else  g <- cbind(  c(fit$fit$dl.dbe)*fit$VC$X1, c(fit$fit$dl.dsigma.st)*fit$VC$X2, c(fit$fit$dl.dnu.st)*fit$VC$X3 )





}





###################################################################################################
###################################################################################################


h <- t(g)%*%g

hbs <- PDef(h - bs)$res

F <- fit$Vb1%*%hbs

edf <- sum(diag(F))
                                                                
list(bs = bs, h = h, hbs = hbs, rAIC = 2*edf - 2*-fit$fit$l, rBIC = log(n)*edf - 2*-fit$fit$l, edf = edf, F = F)



}    

