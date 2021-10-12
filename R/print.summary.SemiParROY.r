print.summary.SemiParROY <- function(x, digits = max(3, getOption("digits") - 3),
                                             signif.stars = getOption("show.signif.stars"),...){
  
 nodi <- 3

 ppR <- ppROY(x)  
 
 cont1par <- ppR$cont1par
 cont2par <- ppR$cont2par
 cont3par <- ppR$cont3par
 
 cop1      <- ppR$cop1
 lind1     <- ppR$lind1

 cop2      <- ppR$cop2
 lind2     <- ppR$lind2 
 
 m1l      <- ppR$m1l
 m2l      <- ppR$m2l 
 m3l      <- ppR$m3l 
 
 

 doff1 <- doff2 <- "log(\u00B7 - 2)"




 bin.link <- x$bl
 
 cp1 <- "\ntheta.1 = "; as.p1 <- x$theta1.a
 cp2 <- "  theta.2 = "; as.p2 <- x$theta2.a
 
 ct1 <- "\ntau.1 = "; kt.p1 <- x$tau1.a
 ct2 <- "  tau.2 = "; kt.p2 <- x$tau2.a

 cs1 <- "\nsigma.1 = "; cs.p1 <- x$sigma1.a
 cs2 <- "  sigma.2 = "; cs.p2 <- x$sigma2.a
 
 cn1 <- "\nnu.1 = "; cn.p1 <- x$nu1.a
 cn2 <- "  nu.2 = "; cn.p2 <- x$nu2.a 

 
 
 main.t1 <- "\nCOPULA 1:"; main.t2 <- "\nCOPULA 2:"     
 cat(main.t1,cop1); cat(main.t2,cop2)  
  
 pscr0(x, type = "ROY")  
 
 
 pscr(x, lind1, m1l, m2l, cont1par, cont2par, cont3par, type = "ROY", digits, signif.stars, m3l, lind2, ...)
 


  CIrs1 <- colMeans(x$CItheta1, na.rm = TRUE)
  CIrs2 <- colMeans(x$CItheta2, na.rm = TRUE)
  CIkt1 <- colMeans(x$CItau1, na.rm = TRUE)
  CIkt2 <- colMeans(x$CItau2, na.rm = TRUE)

  if( x$margins[2] %in% c(cont2par,cont3par) && x$margins[3] %in% c(cont2par,cont3par) ){ CIsi1 <- colMeans(x$CIsig1, na.rm = TRUE); CIsi2 <- colMeans(x$CIsig2, na.rm = TRUE) } 
  if( x$margins[2] %in% cont3par && x$margins[3] %in% cont3par)                         { CIn1  <- colMeans(x$CInu1, na.rm = TRUE);  CIn2  <- colMeans(x$CInu2, na.rm = TRUE)  }


  cat("\n")


  if(x$margins[2] %in% cont1par && x$margins[3] %in% cont1par) cat("n = ",x$n,"  n.sel0 = ",x$n.sel0,"  n.sel1 = ",x$n.sel1, 
                                                                    cp1, format(as.p1, digits=3), "(", format(CIrs1[1],digits=nodi),",",format(CIrs1[2],digits=nodi),")",
                                                                    cp2, format(as.p2, digits=3), "(", format(CIrs2[1],digits=nodi),",",format(CIrs2[2],digits=nodi),")", 
                                                                    ct1, format(kt.p1, digits=3), "(", format(CIkt1[1],digits=nodi),",",format(CIkt1[2],digits=nodi),")",
                                                                    ct2, format(kt.p2, digits=3), "(", format(CIkt2[1],digits=nodi),",",format(CIkt2[2],digits=nodi),")",                                                                     
                                                                    "\ntotal edf = ",format(x$t.edf, digits=3),"\n\n", sep="")

  if(x$margins[2] %in% cont2par && x$margins[3] %in% cont2par) cat("n = ",x$n,"  n.sel0 = ",x$n.sel0,"  n.sel1 = ",x$n.sel1, 
                                                                    cs1, format(cs.p1, digits=3), "(", format(CIsi1[1],digits=nodi),",",format(CIsi1[2],digits=nodi),")",
                                                                    cs2, format(cs.p2, digits=3), "(", format(CIsi2[1],digits=nodi),",",format(CIsi2[2],digits=nodi),")",                                                                  
                                                                    cp1, format(as.p1, digits=3), "(", format(CIrs1[1],digits=nodi),",",format(CIrs1[2],digits=nodi),")",
                                                                    cp2, format(as.p2, digits=3), "(", format(CIrs2[1],digits=nodi),",",format(CIrs2[2],digits=nodi),")", 
                                                                    ct1, format(kt.p1, digits=3), "(", format(CIkt1[1],digits=nodi),",",format(CIkt1[2],digits=nodi),")",
                                                                    ct2, format(kt.p2, digits=3), "(", format(CIkt2[1],digits=nodi),",",format(CIkt2[2],digits=nodi),")",                                                                     
                                                                    "\ntotal edf = ",format(x$t.edf, digits=3),"\n\n", sep="")

  if(x$margins[2] %in% cont3par && x$margins[3] %in% cont3par) cat("n = ",x$n,"  n.sel0 = ",x$n.sel0,"  n.sel1 = ",x$n.sel1, 
                                                                    cs1, format(cs.p1, digits=3), "(", format(CIsi1[1],digits=nodi),",",format(CIsi1[2],digits=nodi),")",
                                                                    cs2, format(cs.p2, digits=3), "(", format(CIsi2[1],digits=nodi),",",format(CIsi2[2],digits=nodi),")",  
                                                                    cn1, format(cn.p1, digits=3), "(", format(CIn1[1],digits=nodi),",",format(CIn1[2],digits=nodi),")",
                                                                    cn2, format(cn.p2, digits=3), "(", format(CIn2[1],digits=nodi),",",format(CIn2[2],digits=nodi),")",                                                                     
                                                                    cp1, format(as.p1, digits=3), "(", format(CIrs1[1],digits=nodi),",",format(CIrs1[2],digits=nodi),")",
                                                                    cp2, format(as.p2, digits=3), "(", format(CIrs2[1],digits=nodi),",",format(CIrs2[2],digits=nodi),")", 
                                                                    ct1, format(kt.p1, digits=3), "(", format(CIkt1[1],digits=nodi),",",format(CIkt1[2],digits=nodi),")",
                                                                    ct2, format(kt.p2, digits=3), "(", format(CIkt2[1],digits=nodi),",",format(CIkt2[2],digits=nodi),")",                                                                     
                                                                    "\ntotal edf = ",format(x$t.edf, digits=3),"\n\n", sep="")




invisible(x)
                
}


