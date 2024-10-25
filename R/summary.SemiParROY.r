summary.SemiParROY <- function(object, n.sim = 100, prob.lev = 0.05, ...){

  bs <- SE <- Vb <- epds <- sigma2.st <- sigma2 <- nu.st <- nu <- et1s <- et2s <- p1s <- p2s <- p11s <- p10s <- p00s <- p01s <- ORs <- GMs <- 1
  CIor <- CIgm <- NULL 
  
  cont1par  <- c(object$VC$m1d)  
  cont2par  <- c(object$VC$m2,object$VC$m2d) 
  cont3par  <- object$VC$m3  
  bin.link  <- object$VC$bl  
  
  n      <- object$n
  n.sel0 <- object$n.se0
  n.sel1 <- object$n.se1

  
  lf <- length(object$coefficients)
  Vb <- object$Vb 
  SE <- sqrt(diag(Vb)) 
  
  
#######
# CIs
#######

susutsnR <- susutsn(object, bs = NULL, lf, cont1par, cont2par, cont3par, prob.lev, type = "ROY", bin.link = bin.link, n.sim = n.sim, K1 = object$VC$K1)

CIsig1 <- susutsnR$CIsig1
CInu1  <- susutsnR$CInu1
CIsig2 <- susutsnR$CIsig2
CInu2  <- susutsnR$CInu2

CIkt1  <- susutsnR$CIkt1
CIkt2  <- susutsnR$CIkt2
CIrs1  <- susutsnR$CIrs1
CIrs2  <- susutsnR$CIrs2

########################   
   
if(object$VC$gc.l == TRUE) gc()

susuR <- susu(object, SE, Vb, K1 = object$VC$K1)
  
tableN <- susuR$tableN
table  <- susuR$table
  
#########################

if(!is.null(CIrs1))  { if(dim(CIrs1)[1]  > 1)  dimnames(CIrs1)[[1]]  <- dimnames(object$theta12)[[1]]  else dimnames(CIrs1)[[1]]  <- ""  }
if(!is.null(CIrs2))  { if(dim(CIrs2)[1]  > 1)  dimnames(CIrs2)[[1]]  <- dimnames(object$theta13)[[1]]  else dimnames(CIrs2)[[1]]  <- ""  }
if(!is.null(CIsig1)) { if(dim(CIsig1)[1] > 1)  dimnames(CIsig1)[[1]] <- dimnames(object$sigma2)[[1]]   else dimnames(CIsig1)[[1]] <- ""  }
if(!is.null(CIsig2)) { if(dim(CIsig2)[1] > 1)  dimnames(CIsig2)[[1]] <- dimnames(object$sigma3)[[1]]   else dimnames(CIsig2)[[1]] <- ""  }
if(!is.null(CInu1))  { if(dim(CInu1)[1]  > 1)  dimnames(CInu1)[[1]]  <- dimnames(object$nu2)[[1]]      else dimnames(CInu1)[[1]]  <- ""  }
if(!is.null(CInu2))  { if(dim(CInu2)[1]  > 1)  dimnames(CInu2)[[1]]  <- dimnames(object$nu3)[[1]]      else dimnames(CInu2)[[1]]  <- ""  }

###########################


rm(bs, SE, Vb, et1s, et2s, p1s, p2s, p11s, p10s, p00s, p01s, ORs, GMs) 
 
  res <- list(formula = object$formula,              
              formula1 = object$gam1$formula, formula2 = object$gam2$formula, formula3 = object$gam3$formula,
              formula4 = object$gam4$formula, formula5 = object$gam5$formula, formula6 = object$gam6$formula,
              formula7 = object$gam7$formula, formula8 = object$gam8$formula, formula9 = object$gam9$formula,
              tableP1 = table[[1]], tableP2 = table[[2]], tableP3 = table[[3]], 
              tableP4 = table[[4]], tableP5 = table[[5]], tableP6 = table[[6]],
              tableP7 = table[[7]], tableP8 = table[[8]], tableP9 = table[[9]],
              tableNP1 = tableN[[1]], tableNP2 = tableN[[2]], tableNP3 = tableN[[3]], 
              tableNP4 = tableN[[4]], tableNP5 = tableN[[5]], tableNP6 = tableN[[6]], 
              tableNP7 = tableN[[7]], tableNP8 = tableN[[8]], tableNP9 = tableN[[9]],
              n = n, n.sel0 = n.sel0, n.sel1 = n.sel1, 
              l.sp1 = object$l.sp1, l.sp2 = object$l.sp2, l.sp3 = object$l.sp3, 
              l.sp4 = object$l.sp4, l.sp5 = object$l.sp5, l.sp6 = object$l.sp6, 
              l.sp7 = object$l.sp7, l.sp8 = object$l.sp8, l.sp9 = object$l.sp9,              
              BivD1 = object$BivD1, BivD2 = object$BivD2, margins = object$margins, bin.link = bin.link, bl = bin.link, t.edf=object$t.edf,
              Model = object$Model, univar.gamlss = FALSE,
              dof12 = object$dof12, dof13 = object$dof13, dof12.a = object$dof12.a, dof13.a = object$dof13.a,   
              m2 = object$VC$m2, m3 = object$VC$m3, m1 = object$VC$m1,
              K1 = object$VC$K1,
              theta12 = object$theta12, theta13 = object$theta13, theta12.a = object$theta12.a, theta13.a = object$theta13.a,
              #tau12 = object$tau12, tau12.a = object$tau12.a, tau13 = object$tau13, tau13.a = object$tau13.a,
              sigma2.a = object$sigma2.a, sigma3.a = object$sigma3.a, sigma2 = object$sigma2, sigma3 = object$sigma3, 
              nu2.a = object$nu2.a, nu3.a = object$nu3.a, nu2 = object$nu2, nu3 = object$nu3, 
              CIsigma2 = CIsig1, CInu2  = CInu1, CIsigma3 = CIsig2, CInu3  = CInu2,
              CItheta12 = CIrs1, CItheta13 = CIrs2)#, CItau12 = CIkt1, CItau13 = CIkt2)

class(res) <- "summary.SemiParROY"
                                             
res

}

