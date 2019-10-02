predict.CopulaCLM <- function(object, eq, type = "link", ...) {

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


# Adjustment for the ordinal-continuous model. The intercept is "fictionally" added (and set arbitrarily to 1) to apply predict():
# NOTE: ss.pred$coefficients.t and ss.pred$Vp.t are not adjusted for now because they result to be NULL for eq == 1.

	if (eq == 1) {
		ss.pred$coefficients <- c(1, ss.pred$coefficients)
			names(ss.pred$coefficients)[1] <- "(Intercept)"
		eq1_c_Vp <- rep(0, length(ind))
		eq1_r_Vp <- c(0, eq1_c_Vp)
		
		ss.pred$Vp <- cbind(eq1_c_Vp, ss.pred$Vp)
		ss.pred$Vp <- rbind(eq1_r_Vp, ss.pred$Vp)
			colnames(ss.pred$Vp)[1] <- rownames(ss.pred$Vp)[1] <- "(Intercept)"
		
		# Predictions are computed and the fictional intercept removed from predict():
		
		if (type == "link") { res <- predict(ss.pred, type = "link", ...) - 1 } # This works as long as the fooling intercept is set to 1.
		if (type == "terms") { res <- predict(ss.pred, type = "terms", ...) }
		if (type == "iterms") { stop("Type equal to iterms is not currently allowed as option.") }
		
		if (type == "response") {
			eta1 <- predict(ss.pred, type = "link", ...) - 1 ; n <- length(eta1)
				eta1 <- matrix(nrow = n, ncol = K1 - 1, eta1)
			cut  <- t(matrix(nrow = K1 - 1, ncol = n, object$coefficients[1 : (K1 - 1)]))
			lp1  <- cbind(cut - eta1)

			pk <- probm(lp1, object$VC$margins[1], only.pr = FALSE, bc = TRUE)$pr	
			p1 <- matrix(nrow = n, ncol = K1, 0)
				p1[,  1] <- pk[, 1]
				p1[, K1] <- 1 - pk[, K1 - 1]

			if (K1 > 2) {
				for (k in 2 : (K1 - 1)) p1[, k] <- pk[, k] - pk[, k - 1]
			}

			pk <- cbind(pk, rep(1))
			res <- list(p1.mass = p1,
                                    p1.cum  = pk )
		}
		if (type == "lpmatrix") { res <- predict(ss.pred, type = "lpmatrix", ...) ; res <- res[, -1] }
		
	} else {
 		 res <- predict(ss.pred, type = type, ...)
	}

return(res)

}
