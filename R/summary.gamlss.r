summary.gamlss <- function(object, n.sim = 100, prob.lev = 0.05, ...){

  bs <- SE <- Vb <- epds <- sigma2.st <- sigma2 <- nu.st <- nu <- est.RHOb <- XX <- Xt <- V <- 1
  
cont1par  <- c(object$VC$m1d, object$VC$bl)   
cont2par  <- c(object$VC$m2,object$VC$m2d) 
cont3par  <- c(object$VC$m3,object$VC$m3d) 
  
  n <- object$n
  epsilon <- 0.0000001; max.p   <- 0.9999999

  lf <- length(object$coefficients)
  
  if(object$VC$surv.flex == FALSE) Vb <- object$Vb else Vb <- object$Vb.t  
  
  SE <- sqrt(diag(Vb)) 
  bs <- rMVN(n.sim, mean = object$coefficients, sigma=Vb) # this is not correct when surv.flex == TRUE but actually I need no sampling here for output so ok 

#######
# CIs
#######

   susutsnR <- susutsn(object, bs, lf, cont1par, cont2par, cont3par, prob.lev, type = "gamls") 
   CIsig2   <- susutsnR$CIsig21
   CInu     <- susutsnR$CInu1

#########################


  if(object$VC$gc.l == TRUE) gc()

  susuR <- susu(object, SE, Vb)
  
  tableN <- susuR$tableN
  table  <- susuR$table

#########################
 
rm(bs, SE, Vb, XX, Xt, V) 
 
  res <- list(tableP1=table[[1]], tableP2=table[[2]], tableP3=table[[3]], 
              tableP4=table[[4]], tableP5=table[[5]], tableP6=table[[6]], tableP7=table[[7]], tableP8=table[[8]],
              tableNP1=tableN[[1]], tableNP2=tableN[[2]], tableNP3=tableN[[3]], 
              tableNP4=tableN[[4]], tableNP5=tableN[[5]], tableNP6=tableN[[6]], tableNP7=tableN[[7]], tableNP8=tableN[[8]], 
              n=n, 
              sigma2=object$sigma2,  
              nu=object$nu,  
              sigma2.a=object$sigma2.a, 
              nu.a=object$nu.a, 
              formula = object$formula,
              formula1=object$gam1$formula, formula2=object$gam2$formula, formula3=object$gam3$formula,
              formula4=object$gam4$formula, formula5=object$gam5$formula, formula6=object$gam6$formula, 
              formula7=object$gam7$formula, formula8=object$gam8$formula,
              t.edf=object$t.edf, 
              CIsig2=CIsig2, CInu=CInu,
              margins = object$margins, 
              l.sp1 = object$l.sp1, l.sp2 = object$l.sp2, l.sp3 = object$l.sp3, 
              l.sp4 = object$l.sp4, l.sp5 = object$l.sp5, l.sp6 = object$l.sp6, 
              l.sp7 = object$l.sp7, l.sp8 = object$l.sp8,
              X2.null = is.null(object$X2), univar.gamlss = TRUE, surv.flex = object$surv.flex
              )
              
              
              
  class(res) <- "summary.gamlss"
      
                                        

res

}

