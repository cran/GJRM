predict2.CLM <- function(object, eq, type, newdata = NULL, ...){

# CopulaCLM: type added to the function's inputs


if(missing(eq)) stop("You must provide the equation number.")

if(eq > object$l.flist) stop("The fitted model has a smaller number of equations.") 


# CLM.shift accounts for the presence of the cut-points and the removed intercept in the ordinal model

if (!is.null(object$VC$K1)) {
        K1 <- object$VC$K1  
	CLM.shift  <- K1 - 2
	CLM.shift2 <- CLM.shift + 1 # This is needed because in CopulaCLM the intercept has been already removed from X1.d2
} else {
	CLM.shift <- 0 ; CLM.shift2 <- 0
}

###
                     
 if(eq==1){ ss.pred <- object$gam1
            ind <- (1:object$X1.d2) + CLM.shift2 
            } 
 if(eq==2){ ss.pred <- object$gam2
            
            #if(object$VC$informative == "yes") ind <- object$X1.d2 + (1:dim(object$VC$fgam$X)[2])
            #else 
            
            ind <- ((object$X1.d2+1):(object$X1.d2+object$X2.d2)) + CLM.shift2
            
            }
 if(eq==3){ ss.pred <- object$gam3
            ind <- ((object$X1.d2+object$X2.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2)) + CLM.shift2 }   

 if(eq==4){ ss.pred <- object$gam4
            ind <- ((object$X1.d2+object$X2.d2+object$X3.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2)) + CLM.shift2 }             

 if(eq==5){ ss.pred <- object$gam5
            ind <- ((object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2+object$X5.d2)) + CLM.shift2 }             
  
 if(eq==6){ ss.pred <- object$gam6
            ind <- ((object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2+object$X5.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2+object$X5.d2+object$X6.d2)) + CLM.shift2 }             

 if(eq==7){ ss.pred <- object$gam7
            ind <- ((object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2+object$X5.d2+object$X6.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2+object$X5.d2+object$X6.d2+object$X7.d2)) + CLM.shift2 }             

 if(eq==8){ ss.pred <- object$gam8
            ind <- ((object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2+object$X5.d2+object$X6.d2+object$X7.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2+object$X5.d2+object$X6.d2+object$X7.d2+object$X8.d2)) + CLM.shift2 } 
                                             
                                   
           ss.pred$coefficients   <- object$coefficients[ind]
           ss.pred$coefficients.t <- object$coef.t[ind]
           ss.pred$Vp   <- object$Vb[ind,ind]
           ss.pred$Vp.t <- object$Vb.t[ind,ind]
           ss.pred$sig2 <- 1
           ss.pred$scale.estimated <- FALSE 


# Adjustment for the ordinal-continuous model:

	   if (eq == 1) {
		if (length(ss.pred$smooth) != 0) {
			for (i in 1 : length(ss.pred$smooth)) {
				ss.pred$smooth[[i]]$first.para <- ss.pred$smooth[[i]]$first.para - 1
				ss.pred$smooth[[i]]$last.para  <- ss.pred$smooth[[i]]$last.para  - 1
			}
		}
		prediction.CLM(ss.pred, object, newdata, type)
	   } else {

###

  #predict.gam(ss.pred, ...)
  predict(ss.pred, type, ...)
	   }

}
