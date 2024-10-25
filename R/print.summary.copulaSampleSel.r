print.summary.copulaSampleSel <- function(x, digits = max(3, getOption("digits") - 3),
                                             signif.stars = getOption("show.signif.stars"),...){

 ppR <- pp(x)  
 
 cont1par <- ppR$cont1par
 cont2par <- ppR$cont2par
 cont3par <- ppR$cont3par
 cop      <- ppR$cop
 lind     <- ppR$lind
 m1l      <- ppR$m1l
 m2l      <- ppR$m2l      
 bin.link <- x$bin.link  
    
 as.p <- x$theta.a  
 main.t <- "\nCOPULA:"; cp <- "  theta = "     

 cat(main.t,cop) 
 
 pscr0(x, type = "copSS")   


 pscr(x, lind, m1l, m2l, cont1par, cont2par, cont3par, type = "copSS", digits, signif.stars, ...)

  #kt.p <- x$tau.a 
  #CIkt <- colMeans(x$CItau, na.rm = TRUE)
  CIrs <- colMeans(x$CItheta, na.rm = TRUE)
  if(x$margins[2] %in% cont2par)  CIsig2 <- colMeans(x$CIsigma, na.rm = TRUE)
  if(x$margins[2] %in% cont3par){ CIsig2 <- colMeans(x$CIsigma, na.rm = TRUE); CInu <- colMeans(x$CInu, na.rm = TRUE)}

  nodi <- 3
  
  if(x$margins[2] %in% cont1par ) cat("\ntheta = ",format(as.p,digits=nodi),"(",format(CIrs[1],digits=nodi),",",format(CIrs[2],digits=nodi),")","\nn = ",x$n,"  n.sel = ", x$n.sel,"\ntotal edf = ",format(x$t.edf,digits=nodi),"\n\n", sep="")  
  
  if(x$margins[2] %in% cont2par ) cat("\nsigma2 = ",format(x$sigma2.a,digits=nodi),"(",format(CIsig2[1],digits=nodi),",",format(CIsig2[2],digits=nodi),")","\ntheta = ",format(as.p,digits=nodi),"(",format(CIrs[1],digits=nodi),",",format(CIrs[2],digits=nodi),")","\nn = ",x$n,"  n.sel = ", x$n.sel,"\ntotal edf = ",format(x$t.edf,digits=nodi),"\n\n", sep="")  
       
  if(x$margins[2] %in% cont3par ) cat("\nsigma2 = ",format(x$sigma2.a,digits=nodi),"(",format(CIsig2[1],digits=nodi),",",format(CIsig2[2],digits=nodi),")","  nu2 = ",format(x$nu.a,digits = nodi),"(",format(CInu[1],digits=nodi),",",format(CInu[2],digits=nodi),")","\ntheta = ",format(as.p,digits=nodi),"(",format(CIrs[1],digits=nodi),",",format(CIrs[2],digits=nodi),")","\nn = ",x$n,"  n.sel = ", x$n.sel,"\ntotal edf = ",format(x$t.edf,digits=nodi),"\n\n", sep="")  
       
       
invisible(x)
                
}


