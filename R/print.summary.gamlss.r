print.summary.gamlss <- function(x, digits = max(3, getOption("digits") - 3),
                                             signif.stars = getOption("show.signif.stars"),...){

   ppR <- pp(x)  
   
     
   paretos <- c("GP","GPII","GPo","DGP","DGPII")
   
   cont1par <- ppR$cont1par
   cont2par <- ppR$cont2par
   cont3par <- ppR$cont3par
   m1l      <- ppR$m1l
 
   s1 <- "sigma = "; s1.p <- x$sigma2.a
   n1 <- "nu = "; n1.p <- x$nu.a
   
   pscr0(x, type = "gamlss")    

  #   if(x$robust == TRUE) cat("\nROBUST Fit")


  
   pscr(x, lind = NULL, m1l, m2l = NULL, cont1par, cont2par, cont3par, type = "gamls", digits, signif.stars, ...)
  
  
  
  if(x$margins[1] %in% c(cont2par,cont3par) )  CIsig2 <- colMeans(x$CIsigma, na.rm = TRUE)
  if(x$margins[1] %in% cont3par)               CInu   <- colMeans(x$CInu, na.rm = TRUE)
  
  
  if(x$margins[1] %in% paretos){
  
    mu   <- mean(x$mu, na.rm = TRUE)
    CImu <- colMeans(x$CImu, na.rm = TRUE)
  
    } 

  
  
  
  

  nodi <- 3
  
  
  if( x$margins[1] %in% paretos) cat(  "mu = ",format(mu,digits=nodi),"(",format(CImu[1],digits=nodi),",",format(CImu[2],digits=nodi),")  ",
                                       s1,format(s1.p,digits=nodi),"(",format(CIsig2[1],digits=nodi),",",format(CIsig2[2],digits=nodi),")",                                                               
                                      "\nn = ",x$n, "  total edf = ",format(x$t.edf,digits=nodi),"\n\n", sep="")    
  
  
  
  if( x$margins[1] %in% cont2par && !(x$margins[1] %in% paretos)) cat(s1,format(s1.p,digits=nodi),"(",format(CIsig2[1],digits=nodi),",",format(CIsig2[2],digits=nodi),")",                                                               
                                      "\nn = ",x$n, "  total edf = ",format(x$t.edf,digits=nodi),"\n\n", sep="")  
                                      

  if( x$margins[1] %in% cont1par) cat("n = ",x$n, "  total edf = ",format(x$t.edf,digits=nodi),"\n\n", sep="")  
  
  if( x$margins[1] %in% cont3par) cat(s1,format(s1.p,digits=nodi),"(",format(CIsig2[1],digits=nodi),",",format(CIsig2[2],digits=nodi),")",
                                     "  ",n1,format(n1.p,digits=nodi),"(",format(CInu[1],digits=nodi),",",format(CInu[2],digits=nodi),")",
                                     "\nn = ",x$n,"  total edf = ",format(x$t.edf,digits=nodi),"\n\n", sep="")  

  #if(x$robust == TRUE) warning("Interval(s) reported at the bottom of the summary output and\np-values for smooth functions (if present) should be used\nwith caution. See ?summary.gamlss for more details.", call. = FALSE) 

  #if(x$robust == TRUE) warning("P-values for smooth functions (if present) should be used\nwith caution. See ?summary.gamlss for more details.", call. = FALSE) 



  if(x$margins[1] %in% paretos){
  
  
  if( sum(as.numeric(x$indx == FALSE)) != 0) warning("Condition violated. You need to check the model.", call. = FALSE) 


  } 
   
invisible(x)
                
}

