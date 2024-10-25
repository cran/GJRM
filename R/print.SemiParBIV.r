print.SemiParBIV <- function(x, ...){

 ppR <- pp(x)  
 
 cont1par <- ppR$cont1par
 cont2par <- ppR$cont2par
 cont3par <- ppR$cont3par
 
 cop      <- ppR$cop
 lind     <- ppR$lind
 
 m1l      <- ppR$m1l
 m2l      <- ppR$m2l 
 
 doff <- "log(\u00B7 - 2)"

 bin.link <- x$bl
 cp <- "theta = "; as.p <- x$theta.a
 main.t <- "\nCOPULA:"     
 cat(main.t,cop) 
  
 pscr0(x, type = "copSS")  
  

  cat("\n\nEQUATION 1")
  cat("\nLink function for mu1:",m1l,"\n")
  cat("Formula: "); print(x$formula[[1]])

  cat("\nEQUATION 2")
  
  cat("\nLink function for mu2:",m2l,"\n")
  cat("Formula: "); print(x$formula[[2]])
  
  
  
  if(!is.null(x$X3) && is.null(x$X4) ){
  
  cat("\nEQUATION 3")
  cat("\nLink function for theta:",lind,"\n") 
  cat("Formula: "); print(x$formula[[3]]) 
  
  }
  
  if(!is.null(x$X3) && !is.null(x$X4) &&  is.null(x$X5)){
  

  cat("\nEQUATION 3")
  if(x$margins[2] != "BE") cat("\nLink function for sigma2:","log","\n") else cat("\nLink function for sigma2:","qlogis","\n") 
  cat("Formula: ");  print(x$formula[[3]]) 
  

  cat("\nEQUATION 4")
  cat("\nLink function for theta:",lind,"\n") 
  cat("Formula: ");  print(x$formula[[4]]) 
  
  } 
  
  
  if(!is.null(x$X3) && !is.null(x$X4) && !is.null(x$X5)){
  
  cat("\nEQUATION 3")
  
  #if(!(x$margins[2] %in% c("TW"))) 
  cat("\nLink function for sigma2:","log","\n") #else cat("\nLink function for sigma2:","qlogis","\n")  

  
  #cat("\nLink function for sigma:","log","\n") 
  cat("Formula: "); print(x$formula[[3]])  
  
  
  cat("\nEQUATION 4")
  if(x$margins[2] %in% c("DAGUM","SM")) cat("\nLink function for nu2:","log","\n")  
  if(x$margins[2] %in% c("TW")) cat("\nLink function for nu2:","qlogis","\n")  
  
  cat("Formula: "); print(x$formula[[4]]) 
    

  cat("\nEQUATION 5")
  cat("\nLink function for theta:",lind,"\n") 
  cat("Formula: "); print(x$formula[[5]]) 
  
  }  

  
  


  cat("\n")
  
  
  
  if(x$Model %in% c("B","BPO") && x$margins[2] %in% cont1par) cat(cp,format(as.p, digits=3),"\nn = ",x$n,"  total edf = ",format(x$t.edf, digits=3),"\n\n", sep="")
  if(x$Model == "BPO0")                                       cat("n = ",x$n,"  total edf = ",format(x$t.edf, digits=3),"\n\n", sep="")
  if(x$Model == "BSS")                                        cat(cp,format(as.p, digits=3),"\nn = ",x$n,"  n.sel = ",x$n.sel,"\ntotal edf = ",format(x$t.edf, digits=3),"\n\n", sep="")
    
  if(x$Model=="B" && x$margins[2] %in% cont2par ) cat("sigma2 = ",x$sigma2.a, "\ntheta = ", format(as.p, digits=3),"\nn = ",x$n,"  total edf = ",format(x$t.edf, digits=3),"\n\n", sep="")
  if(x$Model=="B" && x$margins[2] %in% cont3par ) cat("sigma2 = ",x$sigma2.a, "  nu2 = ",x$nu.a, "\ntheta = ", format(as.p, digits=3),"\nn = ",x$n,"  total edf = ",format(x$t.edf, digits=3),"\n\n", sep="")





invisible(x)

}

