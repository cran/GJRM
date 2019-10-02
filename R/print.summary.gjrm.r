print.summary.gjrm <- function(x, digits = max(3, getOption("digits") - 3),
                                             signif.stars = getOption("show.signif.stars"),...){


 ppR <- pp(x)  
 
 cont1par <- ppR$cont1par
 cont2par <- ppR$cont2par
 cont3par <- ppR$cont3par
 cop      <- ppR$cop
 lind     <- ppR$lind
 m1l      <- ppR$m1l
 m2l      <- ppR$m2l      
     
     
 main.t <- "\nCOPULA:  "     
 cp <- "  theta = "; as.p <- x$theta.a; dof <- x$dof.a
 ct <- "  tau = "; kt.p <- x$tau.a
 s1 <- "sigma.1 = "; s1.p <- x$sigma21.a
 s2 <- "sigma.2 = "; s2.p <- x$sigma22.a 
 n1 <- "nu.1 = "; n1.p <- x$nu1.a
 n2 <- "nu.2 = "; n2.p <- x$nu2.a 
   
   
 cat(main.t,cop) 
  
 pscr0(x)  
 
 pscr(x, lind, m1l, m2l, cont1par, cont2par, cont3par, type = "copR", digits, signif.stars, ...)
  
 CIrs    <- colMeans(x$CItheta, na.rm = TRUE)
 CIkt    <- colMeans(x$CItau, na.rm = TRUE)
 if( x$margins[1] %in% c(cont2par,cont3par) ) CIsig21 <- colMeans(x$CIsig1, na.rm = TRUE)
 if( x$margins[2] %in% c(cont2par,cont3par) ) CIsig22 <- colMeans(x$CIsig2, na.rm = TRUE)
  
 if(x$margins[1] %in% cont3par)  CInu1 <- colMeans(x$CInu1, na.rm = TRUE)
 if(x$margins[2] %in% cont3par)  CInu2 <- colMeans(x$CInu2, na.rm = TRUE)  
  
  
 BivD <- x$BivD 
 if(x$surv == TRUE) BivD <- "N" 
  
 if(BivD == "T" && x$margins[1] %in% c(x$m2,x$m3) && x$margins[2] %in% c(x$m2,x$m3)) CIdof <- colMeans(x$CIdof, na.rm = TRUE)


  nodi <- 3
 
 
 
 if(BivD == "T" && x$margins[1] %in% c(x$m2,x$m3) && x$margins[2] %in% c(x$m2,x$m3)){
 
 
   if( x$margins[1] %in% cont2par && x$margins[2] %in% cont2par) cat(s1,format(s1.p,digits=nodi),"(",format(CIsig21[1],digits=nodi),",",format(CIsig21[2],digits=nodi),")",
                                                                      "  ",s2,format(s2.p,digits=nodi),"(",format(CIsig22[1],digits=nodi),",",format(CIsig22[2],digits=nodi),")",
                                                                      "\ndof = ",format(dof,digits=nodi),"(",format(CIdof[1],digits=nodi),",",format(CIdof[2],digits=nodi),")",                                                                     
                                                                      "\ntheta = ",format(as.p,digits=nodi),"(",format(CIrs[1],digits=nodi),",",format(CIrs[2],digits=nodi),")",
                                                                      ct,format(kt.p,digits=nodi),"(",format(CIkt[1],digits=nodi),",",format(CIkt[2],digits=nodi),")",
                                                                      "\nn = ",x$n, "  total edf = ",format(x$t.edf,digits=nodi),"\n\n", sep="")  
 

 
   if( x$margins[1] %in% cont3par && x$margins[2] %in% cont3par) cat(s1,format(s1.p,digits=nodi),"(",format(CIsig21[1],digits=nodi),",",format(CIsig21[2],digits=nodi),")",
                                                                      "  ",s2,format(s2.p,digits=nodi),"(",format(CIsig22[1],digits=nodi),",",format(CIsig22[2],digits=nodi),")",
                                                                      "\n",n1,format(n1.p,digits=nodi),"(",format(CInu1[1],digits=nodi),",",format(CInu1[2],digits=nodi),")",
                                                                      "  ",n2,format(n2.p,digits=nodi),"(",format(CInu2[1],digits=nodi),",",format(CInu2[2],digits=nodi),")",   
                                                                      "\ndof = ",format(dof,digits=nodi),"(",format(CIdof[1],digits=nodi),",",format(CIdof[2],digits=nodi),")",
                                                                      "\ntheta = ",format(as.p,digits=nodi),"(",format(CIrs[1],digits=nodi),",",format(CIrs[2],digits=nodi),")",
                                                                      ct,format(kt.p,digits=nodi),"(",format(CIkt[1],digits=nodi),",",format(CIkt[2],digits=nodi),")",
                                                                      "\nn = ",x$n,"  total edf = ",format(x$t.edf,digits=nodi),"\n\n", sep="")  
 
 
 
 if( x$margins[1] %in% cont2par && x$margins[2] %in% cont3par) cat(  s1,format(s1.p,digits=nodi),"(",format(CIsig21[1],digits=nodi),",",format(CIsig21[2],digits=nodi),")",
                                                                      "  ",s2,format(s2.p,digits=nodi),"(",format(CIsig22[1],digits=nodi),",",format(CIsig22[2],digits=nodi),")",
                                                                      "\n",n2,format(n2.p,digits=nodi),"(",format(CInu2[1],digits=nodi),",",format(CInu2[2],digits=nodi),")", 
                                                                      "\ndof = ",format(dof,digits=nodi),"(",format(CIdof[1],digits=nodi),",",format(CIdof[2],digits=nodi),")",
                                                                      "\ntheta = ",format(as.p,digits=nodi),"(",format(CIrs[1],digits=nodi),",",format(CIrs[2],digits=nodi),")",
                                                                      ct,format(kt.p,digits=nodi),"(",format(CIkt[1],digits=nodi),",",format(CIkt[2],digits=nodi),")",
                                                                      "\nn = ",x$n, "  total edf = ",format(x$t.edf,digits=nodi),"\n\n", sep="")  
                           
 
   if( x$margins[1] %in% cont3par && x$margins[2] %in% cont2par) cat(s1,format(s1.p,digits=nodi),"(",format(CIsig21[1],digits=nodi),",",format(CIsig21[2],digits=nodi),")",
                                                                      "  ",s2,format(s2.p,digits=nodi),"(",format(CIsig22[1],digits=nodi),",",format(CIsig22[2],digits=nodi),")",
                                                                      "\n",n1,format(n1.p,digits=nodi),"(",format(CInu1[1],digits=nodi),",",format(CInu1[2],digits=nodi),")",
                                                                      "\ndof = ",format(dof,digits=nodi),"(",format(CIdof[1],digits=nodi),",",format(CIdof[2],digits=nodi),")",
                                                                      "\ntheta = ",format(as.p,digits=nodi),"(",format(CIrs[1],digits=nodi),",",format(CIrs[2],digits=nodi),")",
                                                                      ct,format(kt.p,digits=nodi),"(",format(CIkt[1],digits=nodi),",",format(CIkt[2],digits=nodi),")",
                                                                      "\nn = ",x$n, "  total edf = ",format(x$t.edf,digits=nodi),"\n\n", sep="")  
         
 
 
 }else{
 
 
   if( x$margins[1] %in% cont2par && x$margins[2] %in% cont1par && x$surv == TRUE  ) cat(s1,format(s1.p,digits=nodi),"(",format(CIsig21[1],digits=nodi),",",format(CIsig21[2],digits=nodi),")",
                                                                      "\ntheta = ",format(as.p,digits=nodi),"(",format(CIrs[1],digits=nodi),",",format(CIrs[2],digits=nodi),")",
                                                                      ct,format(kt.p,digits=nodi),"(",format(CIkt[1],digits=nodi),",",format(CIkt[2],digits=nodi),")",
                                                                      "\nn = ",x$n, "  total edf = ",format(x$t.edf,digits=nodi),"\n\n", sep="")  

                        
  if( x$margins[1] %in% cont3par && x$margins[2] %in% cont1par && x$surv == TRUE) cat(s1,format(s1.p,digits=nodi),"(",format(CIsig21[1],digits=nodi),",",format(CIsig21[2],digits=nodi),")",
                                                                     "\n",n1,format(n1.p,digits=nodi),"(",format(CInu1[1],digits=nodi),",",format(CInu1[2],digits=nodi),")",
                                                                     "\ntheta = ",format(as.p,digits=nodi),"(",format(CIrs[1],digits=nodi),",",format(CIrs[2],digits=nodi),")",
                                                                     ct,format(kt.p,digits=nodi),"(",format(CIkt[1],digits=nodi),",",format(CIkt[2],digits=nodi),")",
                                                                     "\nn = ",x$n, "  total edf = ",format(x$t.edf,digits=nodi),"\n\n", sep="")  
       

  
  if( x$margins[1] %in% cont2par && x$margins[2] %in% cont2par) cat(s1,format(s1.p,digits=nodi),"(",format(CIsig21[1],digits=nodi),",",format(CIsig21[2],digits=nodi),")",
                                                                     "  ",s2,format(s2.p,digits=nodi),"(",format(CIsig22[1],digits=nodi),",",format(CIsig22[2],digits=nodi),")",
                                                                     "\ntheta = ",format(as.p,digits=nodi),"(",format(CIrs[1],digits=nodi),",",format(CIrs[2],digits=nodi),")",
                                                                     ct,format(kt.p,digits=nodi),"(",format(CIkt[1],digits=nodi),",",format(CIkt[2],digits=nodi),")",
                                                                     "\nn = ",x$n, "  total edf = ",format(x$t.edf,digits=nodi),"\n\n", sep="")  


  if( x$margins[1] %in% cont1par && x$margins[2] %in% cont2par) cat(s2,format(s2.p,digits=nodi),"(",format(CIsig22[1],digits=nodi),",",format(CIsig22[2],digits=nodi),")",
                                                                     "\ntheta = ",format(as.p,digits=nodi),"(",format(CIrs[1],digits=nodi),",",format(CIrs[2],digits=nodi),")",
                                                                     ct,format(kt.p,digits=nodi),"(",format(CIkt[1],digits=nodi),",",format(CIkt[2],digits=nodi),")",
                                                                     "\nn = ",x$n, "  total edf = ",format(x$t.edf,digits=nodi),"\n\n", sep="")  

  if( x$margins[1] %in% cont1par && x$margins[2] %in% cont1par) cat("theta = ",format(as.p,digits=nodi),"(",format(CIrs[1],digits=nodi),",",format(CIrs[2],digits=nodi),")",
                                                                     ct,format(kt.p,digits=nodi),"(",format(CIkt[1],digits=nodi),",",format(CIkt[2],digits=nodi),")",
                                                                     "\nn = ",x$n, "  total edf = ",format(x$t.edf,digits=nodi),"\n\n", sep="")  

if( x$margins[1] %in% cont1par && x$margins[2] %in% cont3par) cat(s2,format(s2.p,digits=nodi),"(",format(CIsig22[1],digits=nodi),",",format(CIsig22[2],digits=nodi),")",
                                                                     "\n",n2,format(n2.p,digits=nodi),"(",format(CInu2[1],digits=nodi),",",format(CInu2[2],digits=nodi),")", 
                                                                     "\ntheta = ",format(as.p,digits=nodi),"(",format(CIrs[1],digits=nodi),",",format(CIrs[2],digits=nodi),")",
                                                                     ct,format(kt.p,digits=nodi),"(",format(CIkt[1],digits=nodi),",",format(CIkt[2],digits=nodi),")",
                                                                     "\nn = ",x$n, "  total edf = ",format(x$t.edf,digits=nodi),"\n\n", sep="")                           

  if( x$margins[1] %in% cont3par && x$margins[2] %in% cont3par) cat(s1,format(s1.p,digits=nodi),"(",format(CIsig21[1],digits=nodi),",",format(CIsig21[2],digits=nodi),")",
                                                                     "  ",s2,format(s2.p,digits=nodi),"(",format(CIsig22[1],digits=nodi),",",format(CIsig22[2],digits=nodi),")",
                                                                     "\n",n1,format(n1.p,digits=nodi),"(",format(CInu1[1],digits=nodi),",",format(CInu1[2],digits=nodi),")",
                                                                     "  ",n2,format(n2.p,digits=nodi),"(",format(CInu2[1],digits=nodi),",",format(CInu2[2],digits=nodi),")",   
                                                                     "\ntheta = ",format(as.p,digits=nodi),"(",format(CIrs[1],digits=nodi),",",format(CIrs[2],digits=nodi),")",
                                                                     ct,format(kt.p,digits=nodi),"(",format(CIkt[1],digits=nodi),",",format(CIkt[2],digits=nodi),")",
                                                                     "\nn = ",x$n,"  total edf = ",format(x$t.edf,digits=nodi),"\n\n", sep="")  

if( x$margins[1] %in% cont2par && x$margins[2] %in% cont3par) cat(  s1,format(s1.p,digits=nodi),"(",format(CIsig21[1],digits=nodi),",",format(CIsig21[2],digits=nodi),")",
                                                                     "  ",s2,format(s2.p,digits=nodi),"(",format(CIsig22[1],digits=nodi),",",format(CIsig22[2],digits=nodi),")",
                                                                     "\n",n2,format(n2.p,digits=nodi),"(",format(CInu2[1],digits=nodi),",",format(CInu2[2],digits=nodi),")", 
                                                                     "\ntheta = ",format(as.p,digits=nodi),"(",format(CIrs[1],digits=nodi),",",format(CIrs[2],digits=nodi),")",
                                                                     ct,format(kt.p,digits=nodi),"(",format(CIkt[1],digits=nodi),",",format(CIkt[2],digits=nodi),")",
                                                                     "\nn = ",x$n, "  total edf = ",format(x$t.edf,digits=nodi),"\n\n", sep="")  

                        
  if( x$margins[1] %in% cont3par && x$margins[2] %in% cont2par) cat(s1,format(s1.p,digits=nodi),"(",format(CIsig21[1],digits=nodi),",",format(CIsig21[2],digits=nodi),")",
                                                                     "  ",s2,format(s2.p,digits=nodi),"(",format(CIsig22[1],digits=nodi),",",format(CIsig22[2],digits=nodi),")",
                                                                     "\n",n1,format(n1.p,digits=nodi),"(",format(CInu1[1],digits=nodi),",",format(CInu1[2],digits=nodi),")",
                                                                     "\ntheta = ",format(as.p,digits=nodi),"(",format(CIrs[1],digits=nodi),",",format(CIrs[2],digits=nodi),")",
                                                                     ct,format(kt.p,digits=nodi),"(",format(CIkt[1],digits=nodi),",",format(CIkt[2],digits=nodi),")",
                                                                     "\nn = ",x$n, "  total edf = ",format(x$t.edf,digits=nodi),"\n\n", sep="")  
       

         

}  
  
invisible(x)
                
}

