AT <- function(x, nm.end, eq = NULL, E = TRUE, treat = TRUE, type = "simultaneous", ind = NULL, 
   n.sim = 100, prob.lev = 0.05, length.out = NULL, hd.plot = FALSE, te.plot = FALSE, 
   main = "Histogram and Kernel Density of Simulated Average Effects", 
   xlab = "Simulated Average Effects", ...){


if(x$Cont == "YES") stop("This function is not suitable for bivariate models with continuous/discrete margins.")
if(x$triv == TRUE && x$Model == "TSS") stop("This function is not suitable for trivariate probit models with double sample selection.")
if(x$Cont == "NO" && x$VC$ccss == "yes" && x$margins[2] != "GA") stop("This function is not suitable for selection models.")

# TESS case? not done yet

delta <- FALSE

CIs <- est.AT <- NULL

##################################################
##################################################
##################################################

if(x$triv == TRUE){

if( is.null(eq)  ) stop("You need to provide the number of the equation containing the endogenous variable.")
if(!is.null(ind) ) stop("This option is not currently allowed for. Get in touch to check progress.")
if(type == "naive" ) stop("This type is not currently implemented. Get in touch to check progress.")
if(missing(nm.end)) stop("You must provide the name of the endogenous variable.")



if(eq==1){ X.int <- as.matrix(x$X1); ind.int <- 1:x$X1.d2                       } 
if(eq==2){ X.int <- as.matrix(x$X2); ind.int <- (1:x$X2.d2) + x$X1.d2           }  
if(eq==3){ X.int <- as.matrix(x$X3); ind.int <- (1:x$X3.d2) + x$X1.d2 + x$X2.d2 }  

if(type == "simultaneous") coef.int <- as.numeric(coef(x)[ind.int])

if(type == "univariate"){

	if(eq==1) ngam <- x$gam1 
	if(eq==2) ngam <- x$gam2 
	if(eq==3) ngam <- x$gam3 
	
	coef.int <- coef(ngam)  

                        }

d0 <- d1 <- X.int
d0[, nm.end] <- 0
d1[, nm.end] <- 1

eti1 <- d1%*%coef.int 
eti0 <- d0%*%coef.int 

p.int1 <- probm(eti1, x$margins[eq])$pr 
p.int0 <- probm(eti0, x$margins[eq])$pr

est.AT <- mean(p.int1, na.rm = TRUE) - mean(p.int0, na.rm = TRUE) 



if(type == "univariate") {bs <- rMVN(n.sim, mean = coef.int, sigma=ngam$Vp)
                          eti1s <- d1%*%t(bs)
                          eti0s <- d0%*%t(bs) }

if(type == "simultaneous")  {bs <- rMVN(n.sim, mean = coef(x), sigma=x$Vb)
                          eti1s <- d1%*%t(bs[,ind.int])
                          eti0s <- d0%*%t(bs[,ind.int]) } 

 peti1s  <- probm(eti1s, x$margins[eq])$pr 
 peti0s  <- probm(eti0s, x$margins[eq])$pr 
 
 est.ATb <- colMeans(peti1s, na.rm = TRUE) - colMeans(peti0s, na.rm = TRUE) 
 CIs     <- as.numeric(quantile(est.ATb, c(prob.lev/2, 1 - prob.lev/2), na.rm = TRUE))

                  
  if(hd.plot == TRUE){
  
  mult <- 100
  
  hist(est.ATb*mult, freq = FALSE, main=main, 
       xlab=xlab, 
       ylim=c(0,max(density(est.ATb*mult)$y,hist(est.ATb*mult, plot = FALSE)$density)), ...)
  lines(density(est.ATb*mult))

                     }


res <- c(CIs[1], est.AT, CIs[2])


out <- list(res=res, prob.lev=prob.lev, sim.AT=est.ATb, type = type, eq = eq, triv = x$triv)


}


##################################################
##################################################
##################################################



if(x$triv == FALSE){



etap.noi <- X.int <- X.noi <- eti1 <- eti0 <- etno <- indS <- bs <- ind.excl <- p.int1 <- p.int0 <- d.int1 <- d.int0 <- p.etn <- d.etn <- ass.p <- ass.pst <- C.11 <- C.10 <- sig2 <- peti1s <- peti0s <- sigma2.st <- sigma2s <- eti1s <- eti0s <- d0 <- d1 <- p.etns <- etnos <- etds <- ass.ps <- 1
diffEf <- fy1.y2 <- est.ATso <- y2 <- CIF <- Pr <- Effects <- C.11 <- C.10 <- NULL

m1d  <- x$VC$m1d 
m2d  <- x$VC$m2d 
m2  <- x$VC$m2 
m3  <- x$VC$m3 
bin.link <- x$VC$bl  
mat <- c("SM","DAGUM","GU","rGU","LO","LN","WEI","iG","GA","BE","FISK") # excludes "N"


end <- 0
epsilon <- 0.0000001 # 0.9999999 sqrt(.Machine$double.eps)
max.p   <- 0.9999999
est.ATb <- NA
indD <- list()


if(x$v1[1] %in% x$v2[-1]) {end <- 1; eq <- 2} 
if(x$v2[1] %in% x$v1[-1]) {end <- 2; eq <- 1}


if(x$margins[2] == "DAGUM" && eq == 2) { if( min(sqrt(x$sigma2)) <= 1) stop("sigma parameter has value(s) smaller than 1, hence the mean is indeterminate.")}
if(x$margins[2] == "SM"    && eq == 2) { if( min(sqrt(x$sigma2)*x$nu) <= 1) stop("sigma*nu has value(s) smaller than 1, hence the mean is indeterminate.")}
if(x$margins[2] == "FISK"  && eq == 2) { if( min(sqrt(x$sigma2)) <= 1) stop("sigma parameter has value(s) smaller than 1, hence the mean is indeterminate.")}
# not sure about FISK2 

if( !( type %in% c("naive","univariate","simultaneous") ) ) stop("Error in parameter type value. It should be one of: naive, univariate or simultaneous.")


if(missing(nm.end)) stop("You must provide the name of the endogenous variable.")


if(x$VC$ccss == "yes" && x$margins[2] != "GA"){
if(x$Model=="BSS" || x$Model=="BPO" || x$Model=="BPO0" || end==0) stop("Calculation of this average treatment effect is valid for recursive models only.")
}


if(type == "univariate" && x$margins[2] %in% c("N","N2") && eq == 2 && x$gamlssfit == FALSE) stop("You need to fit the univariate model to obtain the AT. Refit the model and set gamlssfit = TRUE.")



if(x$VC$ccss == "yes" && x$margins[2] != "GA"){
if(x$margins[2] %in% mat && eq == 2) stop("AT currently available for Gaussian outcome margin only.")
}


if(is.character(nm.end)==FALSE) stop("nm.end is not a character!")
if( !is.null(ind) && E == FALSE) stop("ind is not designed to be used when some observations are excluded from the AT's calculation.")
if( type == "naive" && E == FALSE) stop("It does not make sense to calculate the naive estimate from the treated only.")

if( !is.null(ind) ){ 

    if(is.logical(ind) == FALSE) stop("ind must be a logical variable.")
    if(length(table(ind))!=2 )   stop("ind must be a logical binary variable.")
    if( length(ind) != x$n )     stop("ind must have the same length as the number of observations used in fitting.")   

}


######################################################################


if(is.null(ind)) ind <- 1:x$n


if(E == FALSE ) {

 if( !(x$margins[2] %in% bin.link)) ind <- 1:x$n  
 
 if(eq==1) X.int <- as.matrix(x$X1[ind,])
 if(eq==2) X.int <- as.matrix(x$X2[ind,]) 

    if(treat == TRUE)  ind <- as.logical(X.int[, nm.end]) 
    if(treat == FALSE) ind <- as.logical(X.int[, nm.end])!=TRUE
                                              
}




######################################################################

if(type == "naive" && !(x$margins[2] %in% bin.link) ) stop("Please fit a bivariate model with intercept and endogenous variable only and then use AT with the univariate type option.")

######################################################################


if(type == "naive" && x$margins[2] %in% bin.link){ ## binary binary case with eq = 1 or eq = 2

if(eq==2){
y1 <- x$y1[ind] 
y2 <- x$y2[ind]
}

if(eq==1){
y1 <- x$y2[ind] 
y2 <- x$y1[ind]
}

tab2 <- table(y1, y2)

pY1cT1 <- prop.table(tab2,1)[4] 
pY1cT0 <- prop.table(tab2,1)[3] 

est.AT <- (pY1cT1 - pY1cT0)

sv <- qnorm(prob.lev/2,lower.tail = FALSE) * sqrt( (pY1cT1*(1-pY1cT1))/x$n + (pY1cT0*(1-pY1cT0))/x$n )

CIs <- c(est.AT - sv, est.AT + sv)

est.ATb <- est.ATso <- NULL

}

######################################################################
######################################################################

if(type != "naive" && x$margins[2] %in% bin.link){ ## binary binary case with eq = 1 or eq = 2

#############

if(type == "simultaneous"){
	indD[[1]] <- 1:x$X1.d2 
	indD[[2]] <- x$X1.d2+(1:x$X2.d2)
                       }

if(eq==1){ X.int <- as.matrix(x$X1[ind,])
    if(type == "simultaneous") ind.int <- indD[[1]]
         }

if(eq==2){ X.int <- as.matrix(x$X2[ind,])
    if(type == "simultaneous") ind.int <- indD[[2]] 
         }

if(type == "simultaneous") coef.int <- as.numeric(coef(x)[ind.int])
	   
             
d0 <- d1 <- X.int
d0[,nm.end] <- 0
d1[,nm.end] <- 1


if(type == "simultaneous"){
	eti1 <- d1%*%coef.int 
	eti0 <- d0%*%coef.int 
                       }

if(type == "univariate"){
	if(eq==1) ngam <- x$gam1
	if(eq==2) ngam <- x$gam2

	eti1 <- d1%*%coef(ngam) 
	eti0 <- d0%*%coef(ngam) 
                         }

#############

p.int1 <- probm(eti1, x$margins[eq])$pr 
p.int0 <- probm(eti0, x$margins[eq])$pr

est.AT <- mean(p.int1, na.rm = TRUE) - mean(p.int0, na.rm = TRUE) 


#############


if(delta == FALSE){

 if(type == "univariate") {bs <- rMVN(n.sim, mean = coef(ngam), sigma=ngam$Vp); eti1s <- d1%*%t(bs);           eti0s <- d0%*%t(bs) }
 if(type == "simultaneous")  {bs <- rMVN(n.sim, mean = coef(x), sigma=x$Vb);       eti1s <- d1%*%t(bs[,ind.int]); eti0s <- d0%*%t(bs[,ind.int]) } 

 peti1s  <- probm(eti1s, x$margins[eq])$pr 
 peti0s  <- probm(eti0s, x$margins[eq])$pr 
 est.ATb <- colMeans(peti1s, na.rm = TRUE) - colMeans(peti0s, na.rm = TRUE) 
 
 CIs <- as.numeric(quantile(est.ATb, c(prob.lev/2, 1 - prob.lev/2), na.rm = TRUE))

                   }



  if(hd.plot == TRUE && delta == FALSE){
  
  if(x$margins[2] %in% c(m2, m3) && eq == 2) mult <- 1 else mult <- 100
  
  hist(est.ATb*mult, freq = FALSE, main=main, 
       xlab=xlab, 
       ylim=c(0,max(density(est.ATb*mult)$y,hist(est.ATb*mult, plot = FALSE)$density)), ...)
  lines(density(est.ATb*mult))

                     }





}

######################################################################
######################################################################

if(type != "naive" && x$margins[2] %in% c("GA") && x$VC$ccss == "yes"){ ## binary binary case with eq = 1 or eq = 2

#############

eq <- 2

if(type == "simultaneous"){
	indD[[1]] <- 1:x$X1.d2 
	indD[[2]] <- x$X1.d2+(1:x$X2.d2)
                       }

if(eq==2){ X.int <- as.matrix(x$X2s[ind,])
    if(type == "simultaneous") ind.int <- indD[[2]] 
         }

if(type == "simultaneous") coef.int <- as.numeric(coef(x)[ind.int])
	   
             
d0 <- d1 <- X.int
d0[,nm.end] <- 0
d1[,nm.end] <- 1


if(type == "simultaneous"){
	eti1 <- d1%*%coef.int 
	eti0 <- d0%*%coef.int 
                       }

if(type == "univariate"){

	ngam <- x$gamlss 

	eti1 <- d1%*%coef(ngam)[1:x$X2.d2] 
	eti0 <- d0%*%coef(ngam)[1:x$X2.d2] 
                         }

#############

p.int1 <- exp(eti1) 
p.int0 <- exp(eti0)

est.AT <- mean(p.int1, na.rm = TRUE) - mean(p.int0, na.rm = TRUE) 


#############


if(delta == FALSE){

 if(type == "univariate")    {bs <- rMVN(n.sim, mean = coef(ngam), sigma=ngam$Vb); eti1s <- d1%*%t(bs[,1:x$X2.d2]); eti0s <- d0%*%t(bs[,1:x$X2.d2]) }
 if(type == "simultaneous")  {bs <- rMVN(n.sim, mean = coef(x),    sigma=x$Vb);    eti1s <- d1%*%t(bs[,ind.int]);   eti0s <- d0%*%t(bs[,ind.int]) } 

 peti1s  <- exp(eti1s) 
 peti0s  <- exp(eti0s) 
 est.ATb <- colMeans(peti1s, na.rm = TRUE) - colMeans(peti0s, na.rm = TRUE) 
 
 CIs <- as.numeric(quantile(est.ATb, c(prob.lev/2, 1 - prob.lev/2), na.rm = TRUE))

                   }



  if(hd.plot == TRUE && delta == FALSE){
  
  if(x$margins[2] %in% c(m2, m3) && eq == 2) mult <- 1 else mult <- 100
  
  hist(est.ATb*mult, freq = FALSE, main=main, 
       xlab=xlab, 
       ylim=c(0,max(density(est.ATb*mult)$y,hist(est.ATb*mult, plot = FALSE)$density)), ...)
  lines(density(est.ATb*mult))

                     }





}









######################################################################
######################################################################






if(type != "naive" && !(x$margins[2] %in% bin.link) && eq == 1){


if(is.null(length.out)) length.out <- length( seq( min(ceiling(x$y2)) , max(floor(x$y2)) ) ) 
y2   <- round( seq( min(ceiling(x$y2)) , max(floor(x$y2)), length.out = length.out  ), 2 ) 
 
 ly2  <- length(y2)
 data <- x$dataset[ind,]
 
 if(type == "simultaneous")  {
                           ind.int <- 1:x$X1.d2
                           bs <- rMVN(n.sim, mean = coef(x), sigma = x$Vb) 
                           coefe  <- x$coef[ind.int] 
                           coefes <- t(bs[, ind.int]) 
 
                          }
 
 if(type == "univariate") {bs <- rMVN(n.sim, mean = coef(x$gam1), sigma = x$gam1$Vp) 
                           coefe  <- x$gam1$coefficient 
                           coefes <- t(bs) 
                          }
 
 
 
 
 
 
sratio <- function(x1, x2) x1 - x2  
fy1.y2 <- fy1.y2S <- list()
diffE  <- NA 

diffES <- list()
diffEfSquant <- as.data.frame(matrix(NA, ly2 - 1, 2))


for(i in 1:ly2) {

data[, 2]   <- y2[i]
lpm    <- as.matrix( predict.gam(x$gam1, newdata = data, type = "lpmatrix") )
eta1   <- lpm%*%coefe
etins  <- lpm%*%coefes

fy1.y2[[i]]  <- mean(probm(eta1, x$margins[eq])$pr )
fy1.y2S[[i]] <- colMeans( probm(etins, x$margins[eq])$pr  )

}




for(i in 1:(ly2-1)) {

  diffE[i]          <- sratio(fy1.y2[[i+1]] , fy1.y2[[i]])
  diffES[[i]]       <- sratio(fy1.y2S[[i+1]], fy1.y2S[[i]])      
  diffEfSquant[i, ] <- quantile(diffES[[i]], probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE) 
                            } 




Effects <- data.frame(Effects = diffE, diffEfSquant)  
names(Effects)[2:3] <- names(quantile(c(1,1), probs = c(prob.lev/2,1-prob.lev/2)))
dimnames(Effects)[[1]] <- y2[2:ly2]


if(te.plot == TRUE){

plot(y2[2:ly2], diffE, ylab = "TE", xlab = "Treatment", pch = 16, ylim = c(min(diffEfSquant[,1]),max(diffEfSquant[,2])))
lines(y2[2:ly2], diffE, type = "l")
for (i in 1:(ly2-1)) lines( y = c(diffEfSquant[i,1], diffEfSquant[i,2]), x = c(y2[i+1],y2[i+1]))

}







}






if(type != "naive" && x$margins[2] %in% c("N","N2") && eq == 2){

 if(type == "univariate") {bs <- rMVN(n.sim, mean = coef(x$gamlss), sigma=x$gamlss$Vb)
                           est.AT  <- est.ATso <- coef(x$gamlss)[nm.end] 
                           est.ATb <- bs[, which(names(coef(x$gamlss))==nm.end) ]
                           } 
                           
 if(type == "simultaneous")  {bs <- rMVN(n.sim, mean = coef(x), sigma=x$Vb)
                           est.AT  <- est.ATso <- x$coefficients[nm.end]
                           est.ATb <- bs[, nm.end]        
                           }
                           
 CIs <- as.numeric(quantile(est.ATb, c(prob.lev/2, 1 - prob.lev/2), na.rm = TRUE))
                           

}
  


######################################################################
######################################################################



rm(etap.noi, X.int, X.noi, eti1, eti0, etno, indS, bs, ind.excl, p.int1, p.int0, d.int1, d.int0,
   p.etn, d.etn, ass.p, ass.pst, C.11, C.10, sig2, peti1s, peti0s, sigma2.st, sigma2s, eti1s, eti0s, d0, d1,
   p.etns, etnos, etds, ass.ps) 


res <- c(CIs[1], est.AT, CIs[2])


out <- list(res=res, prob.lev=prob.lev, sim.AT=est.ATb, mar2=x$margins[2], type = type, 
            Effects = Effects, treat = y2, eq = eq, bl = x$VC$bl, triv = x$triv)
 							 
}#### triv   




 
class(out) <- "AT"

out





}


