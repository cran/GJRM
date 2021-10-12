print.SemiParROY <- function(x, ...){

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

 
 
 main.t1 <- "\nCOPULA 1:"; main.t2 <- "\nCOPULA 2:"     
 cat(main.t1,cop1); cat(main.t2,cop2)  
  
 pscr0(x, type = "ROY")  
  

  cat("\n\nEQUATION 1")
  cat("\nLink function for mu.1:",m1l,"\n")
  cat("Formula: "); print(x$formula[[1]])

  cat("\nEQUATION 2")
  
  cat("\nLink function for mu.2:",m2l,"\n")
  cat("Formula: "); print(x$formula[[2]])
  

  cat("\nEQUATION 3")
  
  cat("\nLink function for mu.3:",m3l,"\n")
  cat("Formula: "); print(x$formula[[3]])



if(length(x$formula) > 3 && x$margins[2] %in% c(cont1par) && x$margins[3] %in% c(cont1par)){

  cat("\nEQUATION 4")
  cat("\nLink function for theta.1:",lind1,"\n") 
  cat("Formula: "); print(x$formula[[4]]) 



  cat("\nEQUATION 5")
  cat("\nLink function for theta.2:",lind2,"\n") 
  cat("Formula: "); print(x$formula[[5]]) 


}



if(length(x$formula) > 3 && x$margins[2] %in% c(cont2par) && x$margins[3] %in% c(cont2par)){

  cat("\nEQUATION 4")
  if(x$margins[2] != "BE") cat("\nLink function for sigma.1:","log","\n") else cat("\nLink function for sigma.1:","qlogis","\n") 
  cat("Formula: "); print(x$formula[[4]]) 



  cat("\nEQUATION 5")
  if(x$margins[3] != "BE") cat("\nLink function for sigma.2:","log","\n") else cat("\nLink function for sigma.2:","qlogis","\n")
  cat("Formula: "); print(x$formula[[5]]) 
  

  cat("\nEQUATION 6")
  cat("\nLink function for theta.1:",lind1,"\n") 
  cat("Formula: "); print(x$formula[[6]]) 



  cat("\nEQUATION 7")
  cat("\nLink function for theta.2:",lind2,"\n") 
  cat("Formula: "); print(x$formula[[7]])  


}



if(length(x$formula) > 3 && x$margins[2] %in% c(cont3par) && x$margins[3] %in% c(cont3par)){

  cat("\nEQUATION 4")
  cat("\nLink function for sigma.1:","log","\n") 
  cat("Formula: "); print(x$formula[[4]]) 



  cat("\nEQUATION 5")
  cat("\nLink function for sigma.2:","log","\n")
  cat("Formula: "); print(x$formula[[5]]) 
  

  cat("\nEQUATION 6")
  cat("\nLink function for nu.1:","log","\n")  
  cat("Formula: "); print(x$formula[[6]]) 



  cat("\nEQUATION 7")
  cat("\nLink function for nu.2:","log","\n") 
  cat("Formula: "); print(x$formula[[7]])  
  
  
  cat("\nEQUATION 8")
  cat("\nLink function for theta.1:",lind1,"\n") 
  cat("Formula: "); print(x$formula[[8]]) 



  cat("\nEQUATION 9")
  cat("\nLink function for theta.2:",lind2,"\n") 
  cat("Formula: "); print(x$formula[[9]])    
  
  


}





  
  


  cat("\n")
  
  
  



  if(x$margins[2] %in% cont1par && x$margins[3] %in% cont1par) cat("n = ",x$n,"  n.sel0 = ",x$n.se0,"  n.sel1 = ",x$n.se1, cp1, format(as.p1, digits=3),cp2,format(as.p2, digits=3), "  total edf = ",format(x$t.edf, digits=3),"\n\n", sep="")
  if(x$margins[2] %in% cont2par && x$margins[3] %in% cont2par) cat("n = ",x$n,"  n.sel0 = ",x$n.se0,"  n.sel1 = ",x$n.se1, "\nsigma.1 = ",x$sigma1.a, "  sigma.2 = ",x$sigma2.a, cp1, format(as.p1, digits=3),cp2,format(as.p2, digits=3), "  total edf = ",format(x$t.edf, digits=3),"\n\n", sep="")
  if(x$margins[2] %in% cont3par && x$margins[3] %in% cont3par) cat("n = ",x$n,"  n.sel0 = ",x$n.se0,"  n.sel1 = ",x$n.se1, "\nsigma.1 = ",x$sigma1.a, "  sigma.2 = ",x$sigma2.a,"\nnu.1 = ",x$nu1.a, "  nu.2 = ",x$nu2.a, cp1, format(as.p1, digits=3),cp2,format(as.p2, digits=3), "  total edf = ",format(x$t.edf, digits=3),"\n\n", sep="")
 
  
  
  




invisible(x)

}

