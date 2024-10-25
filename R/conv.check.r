conv.check <- function(x, blather = FALSE){


#if(!is.null(x$conv.sp)){
#if(x$iter.sp >= x$iterlimsp) cat("Smoothing algorithm reached the max. number of loops allowed.\n\n")
#}


est.probs <- c(x$fit$p1, x$fit$p2, x$fit$p3, x$fit$c.copula.be2, x$fit$c.copula.be1)
est.dens  <- c(x$fit$pdf1, x$fit$pdf2, x$fit$pdf3, x$fit$c.copula2.be1be2)

e.v <- eigen(x$fit$hessian, symmetric = TRUE, only.values = TRUE)$values

cat("\nMaximum absolute gradient value:",max(abs(x$fit$gradient)))

if(x$hess==TRUE) mv <- "Observed" else mv <- "Expected" 

if (min(e.v) > 0) cat("\n",mv," information matrix is positive definite\n",sep="") else cat("\n",mv," information matrix is not positive definite\n",sep="")

if(blather == TRUE){

cat("Eigenvalue range: [",min(e.v),",",max(e.v),"]\n", sep = "")

if( (x$l.sp1!=0 || x$l.sp2!=0 || x$l.sp3!=0 || x$l.sp4!=0 || x$l.sp5!=0 || x$l.sp6!=0 || x$l.sp7!=0 || x$l.sp8!=0 || x$l.sp9!=0) && x$fp==FALSE ){

cat("\nTrust region iterations before smoothing parameter estimation:",x$iter.if)
cat("\nLoops for smoothing parameter estimation:",x$iter.sp) 
cat("\nTrust region iterations within smoothing loops:",x$iter.inner)

} else cat("\nTrust region iterations:",x$iter.if)

if( any( is.na(est.probs) == FALSE ) ) cat("\nOverall probability range:",range(est.probs, na.rm = TRUE)) 
if( any( is.na(est.dens)  == FALSE ) ) cat("\nOverall density range:",range(est.dens, na.rm = TRUE)) 

}

cat("\n\n")


}