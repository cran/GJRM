susu <- function(object, SE, Vb, informative = "no", K1 = NULL){

  tableN <- table <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)

  testStat <- getFromNamespace("testStat", "mgcv")
  liu2   <- getFromNamespace("liu2", "mgcv") 

  index <- 1:2
  
  if(informative == "yes") index <- 1
  
  
if (!is.null(K1)) {
  
	CLM.shift  <- K1 - 2
	CLM.shift2 <- CLM.shift + 1 # This is needed because in CopulaCLM the intercept has been already removed from X1.d2
} else {
	CLM.shift <- 0 ; CLM.shift2 <- 0
}  
  
  
  ind1 <- 1:object$gp1
  ind2 <- object$X1.d2 + (1:object$gp2) + CLM.shift2
  ind3 <- ind4 <- ind5 <- ind6 <- ind7 <- ind8 <- NULL 
  
  if(informative == "yes"){ index <- 1; ind2 <- NULL }
  
  
  
  if(!is.null(object$X3) ) {
  
       ind3 <- object$X1.d2 + object$X2.d2 + (1:object$gp3) + CLM.shift2
       index <- 1:3
       
       if(!is.null(object$X4) ) {
       ind4 <- object$X1.d2 + object$X2.d2 + object$X3.d2 + (1:object$gp4) + CLM.shift2
       index <- 1:4
       }
                                
       if(!is.null(object$X5) ) {
       ind5 <- object$X1.d2 + object$X2.d2 + object$X3.d2 + object$X4.d2 + (1:object$gp5) + CLM.shift2
       index <- 1:5  
       }     
                                
       if(!is.null(object$X6) ) {
       ind6 <- object$X1.d2 + object$X2.d2 + object$X3.d2 + object$X4.d2 + object$X5.d2 + (1:object$gp6) + CLM.shift2
       index <- 1:6  
       }  
       
       if(!is.null(object$X7) ) {
       ind7 <- object$X1.d2 + object$X2.d2 + object$X3.d2 + object$X4.d2 + object$X5.d2 + object$X6.d2 + (1:object$gp7) + CLM.shift2
       index <- 1:7  
       }
       
       if(!is.null(object$X8) ) {
       ind8 <- object$X1.d2 + object$X2.d2 + object$X3.d2 + object$X4.d2 + object$X5.d2 + object$X6.d2 + object$X7.d2 + (1:object$gp8) + CLM.shift2
       index <- 1:8  
       }         
                            
  }
                            
  ind <- list( ind1 = ind1,
               ind2 = ind2,
               ind3 = ind3, 
               ind4 = ind4,
               ind5 = ind5,
               ind6 = ind6,
               ind7 = ind7,
               ind8 = ind8)
                

  for(i in index){
  estimate <- object$coefficients[ind[[i]]]
  se       <- SE[ind[[i]]]
  ratio    <- estimate/se
  pv       <- 2*pnorm(abs(ratio), lower.tail = FALSE)
  table[[i]] <- cbind(estimate,se,ratio,pv)
  dimnames(table[[i]])[[2]] <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  }
  

  
  if( object$l.sp1!=0 || object$l.sp2!=0 || object$l.sp3!=0 || object$l.sp4!=0 || object$l.sp5!=0 || object$l.sp6!=0 || object$l.sp7!=0 || object$l.sp8!=0){

  	pTerms.df <- pTerms.chi.sq <- pTerms.pv <- tableN <- list(0, 0, 0, 0, 0, 0, 0, 0)
        XX <- object$R # this should be ok for monotonic effects as well but I may come back to this
        
           for(i in index){

             if(i==1) {mmm <- object$VC$lsgam1; if(mmm==0) next}
             if(i==2) {mmm <- object$VC$lsgam2; if(mmm==0) next} 
             if(i==3) {mmm <- object$VC$lsgam3; if(mmm==0) next} 
             if(i==4) {mmm <- object$VC$lsgam4; if(mmm==0) next} 
             if(i==5) {mmm <- object$VC$lsgam5; if(mmm==0) next} 
             if(i==6) {mmm <- object$VC$lsgam6; if(mmm==0) next} 
             if(i==7) {mmm <- object$VC$lsgam7; if(mmm==0) next}              
             if(i==8) {mmm <- object$VC$lsgam8; if(mmm==0) break} 
  
		for(k in 1:mmm){

                        if(i==1){ gam <- object$gam1; ind <-  gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para  + CLM.shift                               } 
                        if(i==2){ gam <- object$gam2; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para) + CLM.shift2 + object$X1.d2                } 
                        if(i==3){ gam <- object$gam3; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para) + CLM.shift2 + object$X1.d2 + object$X2.d2 }
                        if(i==4){ gam <- object$gam4; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para) + CLM.shift2 + object$X1.d2 + object$X2.d2 + object$X3.d2 }
                        if(i==5){ gam <- object$gam5; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para) + CLM.shift2 + object$X1.d2 + object$X2.d2 + object$X3.d2 + object$X4.d2 }
                        if(i==6){ gam <- object$gam6; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para) + CLM.shift2 + object$X1.d2 + object$X2.d2 + object$X3.d2 + object$X4.d2 + object$X5.d2 }
                        if(i==7){ gam <- object$gam7; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para) + CLM.shift2 + object$X1.d2 + object$X2.d2 + object$X3.d2 + object$X4.d2 + object$X5.d2 + object$X6.d2 }
                        if(i==8){ gam <- object$gam8; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para) + CLM.shift2 + object$X1.d2 + object$X2.d2 + object$X3.d2 + object$X4.d2 + object$X5.d2 + object$X6.d2 + object$X7.d2 }
                          
                          
                        gam$sig2            <- 1
                        gam$scale.estimated <- FALSE                          
                          
                        if(max(gam$smooth[[k]]$null.space.dim) == 0){
                        
                        LRB <- rbind(XX, t(mroot(object$fit$S.h)))
			LRB <- cbind(LRB[, -ind], LRB[, ind])
			ind1 <- (ncol(LRB) - length(ind) + 1):ncol(LRB)
			Rm <- qr.R(qr(LRB, tol = 0, LAPACK = FALSE))[ind1, ind1]
                        B <- mroot(object$Ve[ind, ind, drop = FALSE])
                          
			b.hat <- object$coefficients[ind]
			d <- Rm %*% b.hat
			stat <- sum(d^2)
			ev <- eigen(crossprod(Rm %*% B), symmetric = TRUE, only.values = TRUE)$values
			ev[ev < 0] <- 0
			rank <- sum(ev > max(ev) * .Machine$double.eps^0.8)
			pval <- liu2(stat, ev)                          
                        Tp <- list(stat = stat, pval = pval, rank = rank)  
                          
                        }
                          
			if(max(gam$smooth[[k]]$null.space.dim) != 0){
			
			if(object$VC$surv.flex == FALSE) b <- object$coefficients[ind] else b <- object$coef.t[ind] 
			
			V <- Vb[ind,ind, drop = FALSE]
			
			Xt <- XX[, ind, drop = FALSE] 
			pTerms.df[[i]][k] <- min(ncol(Xt), object$edf11[[i]][k])
			Tp <- testStat(b, Xt, V, pTerms.df[[i]][k], type = 0, res.df = -1)
			
			} 
			
			
			pTerms.chi.sq[[i]][k] <- Tp$stat 
			pTerms.df[[i]][k] <- Tp$rank
                        pTerms.pv[[i]][k] <- Tp$pval
			                 
                }
                
              tableN[[i]] <- cbind(object$edf[[i]], pTerms.df[[i]], pTerms.chi.sq[[i]], pTerms.pv[[i]])
              dimnames(tableN[[i]])[[2]] <- c("edf", "Ref.df", "Chi.sq", "p-value")
              
            }

  if(object$VC$gc.l == TRUE) gc()

  }
  

list(tableN = tableN, table = table)


}

