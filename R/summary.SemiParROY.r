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
CInu2  <- susutsnR$CInu1


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
              dof1 = object$dof1, dof2 = object$dof2, dof1.a = object$dof1.a, dof2.a = object$dof2.a,   
              m2 = object$VC$m2, m3 = object$VC$m3, m1 = object$VC$m1,
              K1 = object$VC$K1,
              theta1 = object$theta1, theta2 = object$theta2, theta1.a = object$theta1.a, theta2.a = object$theta2.a,
              tau1 = object$tau1, tau1.a = object$tau1.a, tau2 = object$tau2, tau2.a = object$tau2.a,
              sigma1.a = object$sigma1.a, sigma2.a = object$sigma2.a, sigma1 = object$sigma1, sigma2 = object$sigma2, 
              nu1.a = object$nu1.a, nu2.a = object$nu2.a, nu1 = object$nu1, nu2 = object$nu2, 
              CIsig1 = CIsig1, CInu1  = CInu1, CIsig2 = CIsig2, CInu2  = CInu2,
              CItheta1 = CIrs1, CItheta2 = CIrs2, CItau1 = CIkt1, CItau2 = CIkt2)

class(res) <- "summary.SemiParROY"
                                             
res

}

