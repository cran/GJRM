summary.SemiParBIV <- function(object, n.sim = 100, prob.lev = 0.05, gm = FALSE, ...){

  bs <- SE <- Vb <- epds <- sigma2.st <- sigma2 <- nu.st <- nu <- et1s <- et2s <- p1s <- p2s <- p11s <- p10s <- p00s <- p01s <- ORs <- GMs <- 1
  CIor <- CIgm <- NULL 
  cont1par  <- c(object$VC$m1d)  
  cont2par  <- c(object$VC$m2,object$VC$m2d) 
  cont3par  <- object$VC$m3  
  bin.link  <- object$VC$bl  
  n <- object$n; n.sel <- object$n.sel

  
  lf <- length(object$coefficients)
  Vb <- object$Vb 
  SE <- sqrt(diag(Vb)) 
  
  
#######
# CIs
#######

   susutsnR <- susutsn(object, bs = NULL, lf, cont1par, cont2par, cont3par, prob.lev, type = "biv", bin.link = bin.link, n.sim = n.sim, K1 = object$VC$K1)
   CIrs     <- susutsnR$CIrs
   CIkt     <- susutsnR$CIkt
   CIsig2   <- susutsnR$CIsig2
   CInu     <- susutsnR$CInu
   est.RHOb <- susutsnR$est.RHOb
       
########################   
   
   
if (gm == TRUE && is.null(object$VC$theta.fx)){   

bs <- susutsnR$bs

   ####
   # for OR and GM
   ##   
   
if(object$VC$margins[2] %in% bin.link && object$VC$Model != "BPO0"){   
     
   et1s <- object$X1%*%t(bs[,1:object$X1.d2])    
   et2s <- object$X2%*%t(bs[,(object$X1.d2+1):(object$X1.d2+object$X2.d2)]) 
     
   p1s <- probm(et1s, object$margins[1], min.dn = object$VC$min.dn, min.pr = object$VC$min.pr, max.pr = object$VC$max.pr)$pr 
   p2s <- probm(et2s, object$margins[2], min.dn = object$VC$min.dn, min.pr = object$VC$min.pr, max.pr = object$VC$max.pr)$pr
   
   p11s <- matrix(NA,dim(p1s)[1],dim(p1s)[2])
     
   if( !is.null(object$X3) ) { for(i in 1:n.sim) p11s[,i] <- mm(BiCDF(p1s[,i], p2s[,i], object$nC, est.RHOb[,i], object$VC$dof), min.pr = object$VC$min.pr, max.pr = object$VC$max.pr  ) }
   if(  is.null(object$X3) ) { for(i in 1:n.sim) p11s[,i] <- mm(BiCDF(p1s[,i], p2s[,i], object$nC, est.RHOb[i], object$VC$dof), min.pr = object$VC$min.pr, max.pr = object$VC$max.pr  )  }
    
 p11s <- mm(p11s,                       min.pr = object$VC$min.pr, max.pr = object$VC$max.pr)              # does not check if sum = 1 but not longer used anyway
 p11s <- mm(p11s,                       min.pr = object$VC$min.pr, max.pr = object$VC$max.pr)  
 p10s <- mm(p1s - p11s,                 min.pr = object$VC$min.pr, max.pr = object$VC$max.pr) 
 p00s <- mm((1 - p2s) - ( p1s - p11s ), min.pr = object$VC$min.pr, max.pr = object$VC$max.pr)
 p01s <- mm(p2s - p11s,                 min.pr = object$VC$min.pr, max.pr = object$VC$max.pr)
 

 ORs <- (p00s*p11s)/(p01s*p10s)

 ORs  <- ifelse(ORs  ==  Inf,  8.218407e+307, ORs ) 
 ORs  <- ifelse(ORs  == -Inf, -8.218407e+307, ORs ) 


 GMs <- colMeans((ORs - 1)/(ORs + 1))
 ORs <- colMeans(ORs)


 CIor <- as.numeric(quantile(ORs,c(prob.lev/2,1-prob.lev/2),na.rm=TRUE))
 CIgm <- as.numeric(quantile(GMs,c(prob.lev/2,1-prob.lev/2),na.rm=TRUE))


}


if(object$VC$gc.l == TRUE) gc()

}
 
#########################
                 
  susuR <- susu(object, SE, Vb, K1 = object$VC$K1)
  
  tableN <- susuR$tableN
  table  <- susuR$table
  
#########################
 
rm(bs, SE, Vb, et1s, et2s, p1s, p2s, p11s, p10s, p00s, p01s, ORs, GMs) 
 
  res <- list(tableP1=table[[1]], tableP2=table[[2]], tableP3=table[[3]], formula = object$formula,
              tableP4=table[[4]], tableP5=table[[5]], tableP6=table[[6]],
              tableP7=table[[7]], tableP8=table[[8]],
              tableNP1=tableN[[1]], tableNP2=tableN[[2]], tableNP3=tableN[[3]], 
              tableNP4=tableN[[4]], tableNP5=tableN[[5]], tableNP6=tableN[[6]], 
              tableNP7=tableN[[7]], tableNP8=tableN[[8]],
              n=n, theta.a=object$theta.a, sigma2.a=object$sigma2.a, sigma.a=object$sigma2.a, nu.a=object$nu.a, 
              theta=object$theta, sigma2=object$sigma2, sigma=object$sigma2, nu=object$nu,
              OR = object$OR, GM = object$GM, 
              formula1=object$gam1$formula, formula2=object$gam2$formula, formula3=object$gam3$formula,
              formula4=object$gam4$formula, formula5=object$gam5$formula, formula6=object$gam6$formula,
              formula7=object$gam7$formula, formula8=object$gam8$formula,
              t.edf=object$t.edf, CItheta=CIrs, CIsig=CIsig2, CInu=CInu,  
              n.sel=n.sel, CIor = CIor, CIgm = CIgm, CItau = CIkt, tau=object$tau, tau.a=object$tau.a,  
              BivD=object$BivD, margins = object$margins, bin.link = bin.link, 
              Model=object$Model, m2 = object$VC$m2, m3 = object$VC$m3, 
              l.sp1 = object$l.sp1, l.sp2 = object$l.sp2, l.sp3 = object$l.sp3, 
              l.sp4 = object$l.sp4, l.sp5 = object$l.sp5, l.sp6 = object$l.sp6, 
              l.sp7 = object$l.sp7, l.sp8 = object$l.sp8, univar.gamlss = FALSE,
              bl = bin.link, K1 = object$VC$K1, 
              dof=object$dof, dof.a=object$dof)

class(res) <- "summary.SemiParBIV"
      
                                        

res

}

