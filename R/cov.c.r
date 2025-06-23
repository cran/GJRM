cov.c <- function(SemiParFit){

#e.v   <- round(min(eigen(SemiParFit$fit$hessian, symmetric = TRUE, only.values = TRUE)$values), 6)
#gradi <- round(max(abs(SemiParFit$fit$gradient)), 1)
#if(SemiParFit$conv.sp == FALSE) {cat("Check convergence using conv.check().\n")}

#e.v   <- min(eigen(SemiParFit$fit$hessian, symmetric = TRUE, only.values = TRUE)$values)
e.v   <- min(eigen(SemiParFit$fit$hessian)$values)

gradi <- max(abs(SemiParFit$fit$gradient))

me1 <- "Maximum absolute gradient value is not close to 0."
me2 <- "Information matrix is not positive definite."
me3 <- "Read the WARNINGS section in ?gjrm."

if(gradi > 10 && e.v < 0){ warning(me1, call. = FALSE); warning(paste(me2,"\n",me3), call. = FALSE)} 
if(gradi > 10 && e.v > 0)  warning(paste(me1,"\n",me3), call. = FALSE)

if(gradi < 10 && e.v < 0)  warning(paste(me2,"\n",me3), call. = FALSE)


}

