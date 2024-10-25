print.summary.SemiParTRIV <- function(x, digits = max(3, getOption("digits") - 3),
                                             signif.stars = getOption("show.signif.stars"),...){

nodi <- 3
cop  <- "Gaussian"
lind <- "atanh"
  
as.p12 <- x$theta12.a
as.p13 <- x$theta13.a
as.p23 <- x$theta23.a
  
main.t <- "\nCOPULA:"     

  if(x$margins[1]=="probit")  m1l <- "probit"
  if(x$margins[1]=="logit")   m1l <- "logit"
  if(x$margins[1]=="cloglog") m1l <- "cloglog"
  
  if(x$margins[2]=="probit")  m2l <- "probit"
  if(x$margins[2]=="logit")   m2l <- "logit"
  if(x$margins[2]=="cloglog") m2l <- "cloglog"
 
  if(x$margins[3]=="probit")  m3l <- "probit"
  if(x$margins[3]=="logit")   m3l <- "logit"
  if(x$margins[3]=="cloglog") m3l <- "cloglog" 
      
cat(main.t,cop) 
cat("\nMARGIN 1: Bernoulli")  
cat("\nMARGIN 2: Bernoulli")
cat("\nMARGIN 3: Bernoulli")   


 pscr(x, lind, m1l, m2l, cont1par = NULL, cont2par = NULL, cont3par = NULL, type = "triv", digits, signif.stars, m3l = m3l, ...)

    
    
if(length(x$CItheta12) == 2){    
  CI12 <- x$CItheta12
  CI13 <- x$CItheta13
  CI23 <- x$CItheta23
}else{
  CI12 <- colMeans(x$CItheta12)
  CI13 <- colMeans(x$CItheta13)
  CI23 <- colMeans(x$CItheta23)
}


# tess can be improved with n.sel
  
if(x$Model == "T")   cat(
                 	"\ntheta12 = ",format(as.p12,digits=nodi),"(",format(CI12[1],digits=nodi),",",format(CI12[2],digits=nodi),")",
                 	"\ntheta13 = ",format(as.p13,digits=nodi),"(",format(CI13[1],digits=nodi),",",format(CI13[2],digits=nodi),")",
                 	"\ntheta23 = ",format(as.p23,digits=nodi),"(",format(CI23[1],digits=nodi),",",format(CI23[2],digits=nodi),")",   
                 	"\nn = ",x$n, "  total edf = ",format(x$t.edf,digits=nodi),
                 	"\n\n", sep="")  


if(x$Model == "TSS") cat(
                 	"\ntheta12 = ",format(as.p12,digits=nodi),"(",format(CI12[1],digits=nodi),",",format(CI12[2],digits=nodi),")",
                 	"\ntheta13 = ",format(as.p13,digits=nodi),"(",format(CI13[1],digits=nodi),",",format(CI13[2],digits=nodi),")",
                 	"\ntheta23 = ",format(as.p23,digits=nodi),"(",format(CI23[1],digits=nodi),",",format(CI23[2],digits=nodi),")",    
                 	"\nn = ",x$n,"  n.sel1 = ",x$n.sel1,"  n.sel2 = ",x$n.sel2,"\ntotal edf = ",format(x$t.edf,digits=nodi),
                 	"\n\n", sep="")  

if(x$Model == "TESS") cat(
                 	"\ntheta12 = ",format(as.p12,digits=nodi),"(",format(CI12[1],digits=nodi),",",format(CI12[2],digits=nodi),")",
                 	"\ntheta13 = ",format(as.p13,digits=nodi),"(",format(CI13[1],digits=nodi),",",format(CI13[2],digits=nodi),")",
                 	"\ntheta23 = ",format(as.p23,digits=nodi),"(",format(CI23[1],digits=nodi),",",format(CI23[2],digits=nodi),")",  
                 	"\nn = ",x$n,"  n.sel1 = ",x$n.sel1,"  total edf = ",format(x$t.edf,digits=nodi),
                 	"\n\n", sep="")  


       
invisible(x)
                
}


