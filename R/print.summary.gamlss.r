print.summary.gamlss <- function(x, digits = max(3, getOption("digits") - 3),
                                             signif.stars = getOption("show.signif.stars"),...){

   ppR <- pp(x)  
   
   cont1par <- ppR$cont1par
   cont2par <- ppR$cont2par
   cont3par <- ppR$cont3par
   m1l      <- ppR$m1l
 
   s1 <- "sigma2 = "; s1.p <- x$sigma2.a
   n1 <- "nu = "; n1.p <- x$nu.a
   
   pscr0(x, type = "gamlss")    
  
   pscr(x, lind = NULL, m1l, m2l = NULL, cont1par, cont2par, cont3par, type = "gamls", digits, signif.stars, ...)
  
  
  
  if(x$margins[1] %in% c(cont2par,cont3par) )  CIsig2 <- colMeans(x$CIsig2, na.rm = TRUE)
  if(x$margins[1] %in% cont3par)               CInu   <- colMeans(x$CInu, na.rm = TRUE)

  nodi <- 3
  
  
  if( x$margins[1] %in% cont2par) cat(s1,format(s1.p,digits=nodi),"(",format(CIsig2[1],digits=nodi),",",format(CIsig2[2],digits=nodi),")",                                                               
                                      "\nn = ",x$n, "  total edf = ",format(x$t.edf,digits=nodi),"\n\n", sep="")  

  if( x$margins[1] %in% cont1par) cat("n = ",x$n, "  total edf = ",format(x$t.edf,digits=nodi),"\n\n", sep="")  
  
  if( x$margins[1] %in% cont3par) cat(s1,format(s1.p,digits=nodi),"(",format(CIsig2[1],digits=nodi),",",format(CIsig2[2],digits=nodi),")",
                                     "  ",n1,format(n1.p,digits=nodi),"(",format(CInu[1],digits=nodi),",",format(CInu[2],digits=nodi),")",
                                     "\nn = ",x$n,"  total edf = ",format(x$t.edf,digits=nodi),"\n\n", sep="")  

     
invisible(x)
                
}

