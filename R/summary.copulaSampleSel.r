summary.copulaSampleSel <- function(object, n.sim = 100, prob.lev = 0.05, ...){

bs <- SE <- Vb <- epds <- sigma2.st <- sigma2 <- nu.st <- nu <- est.RHOb <- et1s <- et2s <- p1s <- p2s <- p11s <- p10s <- p00s <- p01s <- ORs <- GMs <- XX <- Xt <- V <- 1
  
cont1par  <- c(object$VC$m1d)  
cont2par  <- c(object$VC$m2,object$VC$m2d) 
cont3par  <- object$VC$m3  
bin.link  <- object$VC$bl  
n <- object$n; n.sel <- object$n.sel

lf <- length(object$coefficients)
Vb <- object$Vb 
SE <- sqrt(diag(Vb)) 
bs <- rMVN(n.sim, mean = object$coefficients, sigma=Vb)  
  
#######
# CIs
#######

   susutsnR <- susutsn(object, bs, lf, cont1par, cont2par, cont3par, prob.lev, type = "copSS")
   CIrs     <- susutsnR$CIrs
   CIkt     <- susutsnR$CIkt
   CIsig21  <- susutsnR$CIsig21
   CIsig22  <- susutsnR$CIsig22 
   CInu1    <- susutsnR$CInu1
   CInu2    <- susutsnR$CInu2
   CIsig2   <- susutsnR$CIsig2 
   CInu     <- susutsnR$CInu   
   
#######

  if(object$VC$gc.l == TRUE) gc()

  susuR <- susu(object, SE, Vb)
  
  tableN <- susuR$tableN
  table  <- susuR$table

#######  
 
 
 
rm(bs, SE, Vb, epds, sigma2.st, sigma2, est.RHOb, et1s, et2s, p1s, p2s, p11s, p10s, p00s, p01s, ORs, GMs, XX, Xt, V) 
 
  res <- list(tableP1=table[[1]], tableP2=table[[2]], tableP3=table[[3]], formula = object$formula,
              tableP4=table[[4]], tableP5=table[[5]], tableP6=table[[6]],
              tableP7=table[[7]], tableP8=table[[8]],
              tableNP1=tableN[[1]], tableNP2=tableN[[2]], tableNP3=tableN[[3]], 
              tableNP4=tableN[[4]], tableNP5=tableN[[5]], tableNP6=tableN[[6]], 
              tableNP7=tableN[[7]], tableNP8=tableN[[8]],
              n=n, theta.a=object$theta.a, sigma2.a=object$sigma2.a, sigma.a=object$sigma2.a, nu.a=object$nu.a, 
              theta=object$theta, sigma2=object$sigma2, sigma=object$sigma2, nu=object$nu, 
              formula1=object$gam1$formula, formula2=object$gam2$formula, formula3=object$gam3$formula,
              formula4=object$gam4$formula, formula5=object$gam5$formula, formula6=object$gam6$formula,
              formula7=object$gam7$formula, formula8=object$gam8$formula,
              t.edf=object$t.edf, CItheta=CIrs, CIsig=CIsig2, CInu=CInu,  
              n.sel=n.sel, tau=object$tau, tau.a=object$tau.a, CItau = CIkt, 
              BivD=object$BivD, margins = object$margins, bin.link = bin.link, 
              l.sp1 = object$l.sp1, l.sp2 = object$l.sp2, l.sp3 = object$l.sp3, 
              l.sp4 = object$l.sp4, l.sp5 = object$l.sp5, l.sp6 = object$l.sp6, 
              l.sp7 = object$l.sp7, l.sp8 = object$l.sp8, univar.gamlss = FALSE,
              bl = bin.link, dof=object$dof, dof.a=object$dof, K1 = NULL)
  class(res) <- "summary.copulaSampleSel"
      
                                        

res

}
