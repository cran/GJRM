print.copulaSampleSel <- function(x, ...){


 ppR <- pp(x)  
 
 cont1par <- ppR$cont1par
 cont2par <- ppR$cont2par
 cont3par <- ppR$cont3par
 cop      <- ppR$cop
 lind     <- ppR$lind
 m1l      <- ppR$m1l
 m2l      <- ppR$m2l 
 bin.link <- x$bl
 
 as.p <- x$theta.a
 main.t <- "\nCOPULA:  "   
 
 cat(main.t,cop) 
        
 pscr0(x, type = "copSS")    






  cat("\n\nEQUATION 1")
  cat("\nLink function for mu.1:",m1l,"\n")
  cat("Formula: "); print(x$formula[[1]])



  cat("\nEQUATION 2")
  cat("\nLink function for mu.2:",m2l,"\n")
  cat("Formula: "); print(x$formula[[2]])
  
  
  
  if(!is.null(x$X3) && is.null(x$X4) ){
  
  cat("\nEQUATION 3")
  cat("\nLink function for theta:",lind,"\n") 
  cat("Formula: "); print(x$formula[[3]]) 
  
  }
  
  
  
  
  
  
  if(!is.null(x$X3) && !is.null(x$X4) &&  is.null(x$X5)){
  

  cat("\nEQUATION 3")
  if(x$margins[2] != "BE") cat("\nLink function for sigma:","log","\n") else cat("\nLink function for sigma:","qlogis","\n") 
  cat("Formula: ");  print(x$formula[[3]]) 
  

  cat("\nEQUATION 4")
  cat("\nLink function for theta:",lind,"\n") 
  cat("Formula: ");  print(x$formula[[4]]) 
  
  } 
  
  
  
  
  
  
  if(!is.null(x$X3) && !is.null(x$X4) && !is.null(x$X5)){
  
  cat("\nEQUATION 3")
  # if(!(x$margins[2] %in% c("TW"))) 
  cat("\nLink function for sigma:","log","\n") #else cat("\nLink function for sigma:","qlogis","\n")  
  cat("Formula: "); print(x$formula[[3]]) 
  
  
  cat("\nEQUATION 4")
  if(x$margins[2] %in% c("DAGUM","SM")) cat("\nLink function for nu:","log","\n")  
  if(x$margins[2] %in% c("TW")) cat("\nLink function for nu:","qlogis","\n")  

  cat("Formula: "); print(x$formula[[4]]) 
    

  cat("\nEQUATION 5")
  cat("\nLink function for theta:",lind,"\n") 
  cat("Formula: "); print(x$formula[[5]]) 
  
  }  
  
  
 
  cat("\n")
  
  if(x$margins[2] %in% cont1par ) cat("n = ",x$n,"  n.sel = ", x$n.sel,"\ntheta = ", format(as.p, digits=3),"  total edf = ",format(x$t.edf, digits=3),"\n\n", sep="")
  if(x$margins[2] %in% cont2par ) cat("n = ",x$n,"  n.sel = ", x$n.sel,"\nsigma = ",x$sigma2.a, "\ntheta = ", format(as.p, digits=3),"  total edf = ",format(x$t.edf, digits=3),"\n\n", sep="")
  if(x$margins[2] %in% cont3par ) cat("n = ",x$n,"  n.sel = ", x$n.sel,"\nsigma = ",x$sigma2.a, "  nu = ",x$nu.a, "\ntheta = ", format(as.p, digits=3),"  total edf = ",format(x$t.edf, digits=3),"\n\n", sep="")

invisible(x)

}

