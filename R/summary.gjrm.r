summary.gjrm <- function(object, n.sim = 100, prob.lev = 0.05, ...){

  bs <- SE <- Vb <- epds <- sigma2.st <- sigma2 <- nu.st <- nu <- est.RHOb <- XX <- Xt <- V <- 1 
  
  cont1par  <- c(object$VC$m1d,object$VC$bl)   
  cont2par  <- c(object$VC$m2,object$VC$m2d) 
  cont3par  <- c(object$VC$m3,object$VC$m3d) 
  n <- object$n  

  lf <- length(object$coefficients)
  Vb <- object$Vb 
  SE <- sqrt(diag(Vb)) 
  bs <- rMVN(n.sim, mean = object$coefficients, sigma=Vb)  

#######
# CIs
#######

   susutsnR <- susutsn(object, bs, lf, cont1par, cont2par, cont3par, prob.lev)
   CIrs     <- susutsnR$CIrs
   CIkt     <- susutsnR$CIkt
   CIsig21  <- susutsnR$CIsig21
   CIsig22  <- susutsnR$CIsig22 
   CInu1    <- susutsnR$CInu1
   CInu2    <- susutsnR$CInu2
   CIdof    <- susutsnR$CIdof
   est.RHOb <- susutsnR$est.RHOb
   
   if(sum(CIsig21) == 0) CIsig21 <- NULL 
   if(sum(CIsig22) == 0) CIsig22 <- NULL
   
#########################
         

  if(object$VC$gc.l == TRUE) gc()

  susuR <- susu(object, SE, Vb)
  
  tableN <- susuR$tableN
  table  <- susuR$table
  
#########################  

if(!is.null(CIrs))    { if(dim(CIrs)[1]    > 1)  dimnames(CIrs)[[1]]    <- dimnames(object$theta)[[1]]  else dimnames(CIrs)[[1]]    <- ""  }
if(!is.null(CIsig21)) { if(dim(CIsig21)[1] > 1)  dimnames(CIsig21)[[1]] <- dimnames(object$sigma1)[[1]] else dimnames(CIsig21)[[1]] <- ""  }
if(!is.null(CIsig22)) { if(dim(CIsig22)[1] > 1)  dimnames(CIsig22)[[1]] <- dimnames(object$sigma2)[[1]] else dimnames(CIsig22)[[1]] <- ""  }
if(!is.null(CInu1))   { if(dim(CInu1)[1]   > 1)  dimnames(CInu1)[[1]]   <- dimnames(object$nu1)[[1]]    else dimnames(CInu1)[[1]]   <- ""  }
if(!is.null(CInu2))   { if(dim(CInu2)[1]   > 1)  dimnames(CInu2)[[1]]   <- dimnames(object$nu2)[[1]]    else dimnames(CInu2)[[1]]   <- ""  }
if(!is.null(CIdof))   { if(dim(CIdof)[1]   > 1)  dimnames(CIdof)[[1]]   <- dimnames(object$dof)[[1]]    else dimnames(CIdof)[[1]]   <- ""  }

#########################

rm(bs, SE, Vb, XX, Xt, V) 
 
  res <- list(tableP1=table[[1]], tableP2=table[[2]], tableP3=table[[3]], est.RHOb = est.RHOb,
              tableP4=table[[4]], tableP5=table[[5]], tableP6=table[[6]], 
              tableP7=table[[7]], tableP8=table[[8]],
              tableNP1=tableN[[1]], tableNP2=tableN[[2]], tableNP3=tableN[[3]], 
              tableNP4=tableN[[4]], tableNP5=tableN[[5]], tableNP6=tableN[[6]], tableNP7=tableN[[7]], 
              tableNP8=tableN[[8]], Model = object$Model,
              n=n, theta=object$theta, theta.a=object$theta.a,
              dof=object$dof, dof.a=object$dof.a,
              sigma21=object$sigma21, sigma22=object$sigma22, 
              sigma1=object$sigma1, sigma2=object$sigma2, 
              nu1=object$nu1, nu2=object$nu2, #tau=object$tau, 
              sigma21.a=object$sigma21.a, sigma22.a=object$sigma22.a, 
              sigma1.a=object$sigma1.a, sigma2.a=object$sigma2.a, 
              nu1.a=object$nu1.a, nu2.a=object$nu2.a, 
              #tau.a=object$tau.a, 
              formula = object$formula,
              formula1=object$gam1$formula, formula2=object$gam2$formula, formula3=object$gam3$formula,
              formula4=object$gam4$formula, formula5=object$gam5$formula, formula6=object$gam6$formula, 
              formula7=object$gam7$formula, formula8=object$gam8$formula,
              t.edf=object$t.edf, 
              CItheta=CIrs, CIsigma1=CIsig21, CIsigma2=CIsig22, CInu1=CInu1, CInu2=CInu2, #CItau = CIkt,
              CIdof = CIdof, BivD=object$BivD, margins = object$margins, 
              l.sp1 = object$l.sp1, l.sp2 = object$l.sp2, l.sp3 = object$l.sp3, 
              l.sp4 = object$l.sp4, l.sp5 = object$l.sp5, l.sp6 = object$l.sp6, 
              l.sp7 = object$l.sp7, l.sp8 = object$l.sp8,
              X3.null = is.null(object$X3), univar.gamlss = FALSE, m2 = object$VC$m2, m3 = object$VC$m3, surv = object$surv,
              surv.flex = object$surv.flex, K1 = NULL, end.surv = object$end.surv, n.sel = object$n.sel, model = object$model)
              
              
              
  class(res) <- "summary.gjrm"
      
                                        

res

}

