print.summary.SemiParBIV <- function(x, digits = max(3, getOption("digits") - 3),
                                             signif.stars = getOption("show.signif.stars"),...){

 ppR <- pp(x)  
 
 cont1par <- ppR$cont1par
 cont2par <- ppR$cont2par
 cont3par <- ppR$cont3par
 
 cont1par <- cont1par[1:2] # test this
 
 cop      <- ppR$cop
 lind     <- ppR$lind
 
 m1l      <- ppR$m1l
 m2l      <- ppR$m2l  
   
 bin.link <- x$bin.link  
 main.t <- "\nCOPULA:  "     
 cp <- "  theta = "; as.p <- x$theta.a
 ct <- "  tau = "; kt.p <- x$tau.a  

 cat(main.t,cop) 
 
 pscr0(x, type = "copSS")  
 
 
 pscr(x, lind, m1l, m2l, cont1par, cont2par, cont3par, type = "biv", digits, signif.stars, ...)


  CIrs <- colMeans(x$CItheta, na.rm = TRUE)
  if(x$margins[2] %in% cont2par)  CIsig2 <- colMeans(x$CIsig2, na.rm = TRUE)
  if(x$margins[2] %in% cont3par) {CIsig2 <- colMeans(x$CIsig2, na.rm = TRUE) ; CInu <- colMeans(x$CInu, na.rm = TRUE)}
  CIkt    <- colMeans(x$CItau, na.rm = TRUE)


  nodi <- 3

  if( (x$Model=="B" || x$Model=="BPO") && x$margins[2] %in% bin.link) cat("\nn = ",x$n,cp,format(as.p,digits=nodi),"(",format(CIrs[1],digits=nodi),",",format(CIrs[2],digits=nodi),")",ct,format(kt.p,digits=nodi),"(",format(CIkt[1],digits=nodi),",",format(CIkt[2],digits=nodi),")","\ntotal edf = ",format(x$t.edf,digits=nodi),"\n\n", sep="")  
  if( x$Model == "BPO0" ) cat("\nn = ",x$n,"  total edf = ",format(x$t.edf,digits=nodi),"\n\n", sep="")  
  if(x$Model=="BSS") cat("\nn = ",x$n,"  n.sel = ",x$n.sel,cp,format(as.p,digits=nodi),"(",format(CIrs[1],digits=nodi),",",format(CIrs[2],digits=nodi),")","\ntau = ",format(kt.p,digits=nodi),"(",format(CIkt[1],digits=nodi),",",format(CIkt[2],digits=nodi),")","  total edf = ",format(x$t.edf,digits=nodi),"\n\n", sep="") 
     

  if(x$Model=="B" && x$margins[2] %in% cont1par ) cat("\nn = ",x$n,"\ntheta = ",format(as.p,digits=nodi),"(",format(CIrs[1],digits=nodi),",",format(CIrs[2],digits=nodi),")",ct,format(kt.p,digits=nodi),"(",format(CIkt[1],digits=nodi),",",format(CIkt[2],digits=nodi),")","\ntotal edf = ",format(x$t.edf,digits=nodi),"\n\n", sep="")  
     
  if(x$Model=="B" && x$margins[2] %in% cont2par ) cat("\nn = ",x$n,"  sigma2 = ",format(x$sigma2.a,digits=nodi),"(",format(CIsig2[1],digits=nodi),",",format(CIsig2[2],digits=nodi),")","\ntheta = ",format(as.p,digits=nodi),"(",format(CIrs[1],digits=nodi),",",format(CIrs[2],digits=nodi),")",ct,format(kt.p,digits=nodi),"(",format(CIkt[1],digits=nodi),",",format(CIkt[2],digits=nodi),")","\ntotal edf = ",format(x$t.edf,digits=nodi),"\n\n", sep="")  
       
  if(x$Model=="B" && x$margins[2] %in% cont3par ) cat("\nn = ",x$n,"  sigma2 = ",format(x$sigma2.a,digits=nodi),"(",format(CIsig2[1],digits=nodi),",",format(CIsig2[2],digits=nodi),")","  nu = ",format(x$nu.a,digits = nodi),"(",format(CInu[1],digits=nodi),",",format(CInu[2],digits=nodi),")","\ntheta = ",format(as.p,digits=nodi),"(",format(CIrs[1],digits=nodi),",",format(CIrs[2],digits=nodi),")",ct,format(kt.p,digits=nodi),"(",format(CIkt[1],digits=nodi),",",format(CIkt[2],digits=nodi),")","\ntotal edf = ",format(x$t.edf,digits=nodi),"\n\n", sep="")  
      





invisible(x)
                
}


