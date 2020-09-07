summary.SemiParTRIV <- function(object, n.sim = 100, prob.lev = 0.05, ...){

n <- object$n
  
if(is.null(object$X4)) epd12s <- epd13s <- epd23s <- NA  
if(!is.null(object$X4)) epd12s <- epd13s <- epd23s <- matrix(NA, n, n.sim)  
  
  
  lf <- length(object$coefficients)
  Vb <- object$Vb 
  SE <- sqrt(diag(Vb)) 
  bs <- rMVN(n.sim, mean = object$coefficients, sigma=Vb)  

if(!is.null(object$X4)){      
  theta12.st <- object$VC$X4%*%t(bs[,(object$VC$X1.d2 + object$VC$X2.d2 + object$VC$X3.d2 + 1):(object$VC$X1.d2 + object$VC$X2.d2 + object$VC$X3.d2 + object$VC$X4.d2)])    
  theta13.st <- object$VC$X5%*%t(bs[,(object$VC$X1.d2 + object$VC$X2.d2 + object$VC$X3.d2 + object$VC$X4.d2 + 1):(object$VC$X1.d2 + object$VC$X2.d2 + object$VC$X3.d2 + object$VC$X4.d2 + object$VC$X5.d2)])    
  theta23.st <- object$VC$X6%*%t(bs[,(object$VC$X1.d2 + object$VC$X2.d2 + object$VC$X3.d2 + object$VC$X4.d2 + object$VC$X5.d2 + 1):(object$VC$X1.d2 + object$VC$X2.d2 + object$VC$X3.d2 + object$VC$X4.d2 + object$VC$X5.d2 + object$VC$X6.d2)])  
}



if(object$VC$Chol == FALSE){
  
 epd12s <- teta.tr(object$VC, bs[, lf-2])$teta  
 epd13s <- teta.tr(object$VC, bs[, lf-1])$teta  
 epd23s <- teta.tr(object$VC, bs[, lf]  )$teta  
  
}
  
  
if(object$VC$Chol == TRUE){ 
  
  if(is.null(object$X4)){  
 
     for(i in 1:n.sim){ 

      Sigma <- PosDefCor(Sigma = 1, Chol = TRUE, theta12.st = bs[i, lf-2], theta13.st = bs[i, lf-1], theta23.st = bs[i, lf])
      
      epd12s[i] <- Sigma[1,2]
      epd13s[i] <- Sigma[1,3]
      epd23s[i] <- Sigma[2,3]  
    }
}    
    
    
   if(!is.null(object$X4)){  
      
    for(i in 1:n){    

       for(j in 1:n.sim){
       
      Sigma <- PosDefCor(Sigma = 1, Chol = TRUE, theta12.st = theta12.st[i, j], theta13.st = theta13.st[i, j], theta23.st = theta23.st[i, j])
     
      epd12s[i,j] <- Sigma[1,2]
      epd13s[i,j] <- Sigma[1,3]
      epd23s[i,j] <- Sigma[2,3] 
      
      }}
      
    }    
    
}
               
  
if(is.null(object$X4)){   
  CI12s <- rowQuantiles(t(epd12s), probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE) 
  CI13s <- rowQuantiles(t(epd13s), probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE) 
  CI23s <- rowQuantiles(t(epd23s), probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE) 
}

if(!is.null(object$X4)){   
  CI12s <- rowQuantiles(epd12s, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE) 
  CI13s <- rowQuantiles(epd13s, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE) 
  CI23s <- rowQuantiles(epd23s, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE) 
}



  if(object$VC$gc.l == TRUE) gc()
  
  susuR <- susu(object, SE, Vb)
  
  tableN <- susuR$tableN
  table  <- susuR$table  
  
  
rm(bs, SE, Vb) 
 
  res <- list(formula = object$formula, tableP1=table[[1]], tableP2=table[[2]], tableP3=table[[3]], 
              tableP4=table[[4]], tableP5=table[[5]], tableP6=table[[6]],
              tableP7=table[[7]],tableP8=table[[8]],
              tableNP1=tableN[[1]], tableNP2=tableN[[2]], tableNP3=tableN[[3]], 
              tableNP4=tableN[[4]], tableNP5=tableN[[5]], tableNP6=tableN[[6]], 
              tableNP7=tableN[[7]],tableNP8=tableN[[8]],
              n=n, n.sel1 = object$n.sel1, n.sel2 = object$n.sel2,  
              theta12.a=object$theta12.a, theta12=object$theta12, Model = object$Model,
              theta13.a=object$theta13.a, theta13=object$theta13,
              theta23.a=object$theta23.a, theta23=object$theta23,
              formula1=object$gam1$formula, formula2=object$gam2$formula, formula3=object$gam3$formula,
              formula4=object$gam4$formula, formula5=object$gam5$formula, formula6=object$gam6$formula,
              formula7=object$gam7$formula,formula8=object$gam8$formula,
              t.edf=object$t.edf, CI12s=CI12s, CI13s=CI13s, CI23s=CI23s,
              l.sp1 = object$l.sp1, l.sp2 = object$l.sp2, l.sp3 = object$l.sp3, 
              l.sp4 = object$l.sp4, l.sp5 = object$l.sp5, l.sp6 = object$l.sp6,
              l.sp7 = object$l.sp7, l.sp8 = object$l.sp8,
              univar.gamlss = FALSE, margins = object$margins, K1 = NULL
              )
  class(res) <- "summary.SemiParTRIV"
      
                                        

res

}

