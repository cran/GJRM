plot.SemiParBIV <- function(x, eq, ...){

 if(missing(eq)) stop("You must provide the equation number.")

 if(eq > x$l.flist) stop("The fitted model has a smaller number of equations.") 

if (!is.null(x$VC$K1)) {
    
        K1  <- x$VC$K1
	CLM.shift  <- K1 - 2
	CLM.shift2 <- CLM.shift + 1 # This is needed because in CopulaCLM the intercept has been already removed from X1.d2
} else {
	CLM.shift <- CLM.shift2 <- 0
}  


# if(!is.null(x$VC$K1) && eq == 1 && seWithMean == TRUE) stop("Option seWithMean not allowed for the case of ordinal response.")
# plot.gam(ss.plot, ...)


 
 if(eq==1){ ss.plot <- x$gam1
            ind <- (CLM.shift + 1) : (x$X1.d2 + CLM.shift2) # intercept removed from gam1 in ordinal case for convenience and here we pretend to have it again but value is not correct
            } 
 if(eq==2){ ss.plot <- x$gam2
            ind <- (x$X1.d2 + 1 + CLM.shift2):(x$X1.d2 + x$X2.d2 + CLM.shift2) 
            }
 if(eq==3){ ss.plot <- x$gam3
            ind <- (x$X1.d2 + x$X2.d2 + 1 + CLM.shift2):(x$X1.d2 + x$X2.d2 + x$X3.d2 + CLM.shift2) }
            
 if(eq==4){ ss.plot <- x$gam4
            ind <- (x$X1.d2 + x$X2.d2 + x$X3.d2 + 1 + CLM.shift2):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + CLM.shift2) }
            
 if(eq==5){ ss.plot <- x$gam5
            ind <- (x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + 1 + CLM.shift2):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2 + CLM.shift2) }
            
 if(eq==6){ ss.plot <- x$gam6
            ind <- (x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2 + 1 + CLM.shift2):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2 + x$X6.d2 + CLM.shift2) }   
            
 if(eq==7){ ss.plot <- x$gam7
            ind <- (x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2 + x$X6.d2 + 1 + CLM.shift2):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2 + x$X6.d2 + x$X7.d2 + CLM.shift2) }             
               
 if(eq==8){ ss.plot <- x$gam8
            ind <- (x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2 + x$X6.d2 + x$X7.d2 + 1 + CLM.shift2):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2 + x$X6.d2 + x$X7.d2 + x$X8.d2 + CLM.shift2) }             
 
 
 if(any(names(ss.plot$coefficients) == "(Intercept)")  &&  length(ss.plot$coefficients) == 1) stop("There is nothing to plot.")
 
  
           ss.plot$coefficients <- x$coefficients[ind]
           ss.plot$Vp <- x$Vb[ind,ind]
           ss.plot$Vp.t <- x$Vb.t[ind,ind]
           ss.plot$sig2 <- 1
           ss.plot$edf <- diag(x$F)[ind]
           ss.plot$scale.estimated <- FALSE 
           ss.plot$call$data <- x$call$data

  plot(ss.plot, ...) 
  
       
}

