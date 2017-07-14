print.SemiParTRIV <- function(x, ...){

  cop <- "Gaussian"; lind <- "atanh"
  
  as.p12 <- x$theta12.a
  as.p13 <- x$theta13.a
  as.p23 <- x$theta23.a
  
  main.t <- "\nCOPULA:  "     
      
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
      
  cat("\n\nEQUATION 1")
  cat("\nLink function for mu.1:",m1l,"\n")
  cat("Formula: "); print(x$formula[[1]])

  cat("\nEQUATION 2")
  cat("\nLink function for mu.2:",m2l,"\n")
  cat("Formula: "); print(x$formula[[2]])
  
  cat("\nEQUATION 3")
  cat("\nLink function for mu.3:",m3l,"\n")
  cat("Formula: "); print(x$formula[[3]])  
 
if(!is.null(x$X4)){

  cat("\nEQUATION 4\n")
  cat("Formula: "); print(x$formula[[4]])  

  cat("\nEQUATION 5\n")
  cat("Formula: "); print(x$formula[[5]])  
  
  cat("\nEQUATION 6\n")
  cat("Formula: "); print(x$formula[[6]])    

}



 
 
  cat("\n")
            
if(x$Model == "T" || x$Model == "TESS")   cat("n = ",x$n,"\ntheta12 = ", format(as.p12, digits=3),"  theta13 = ", format(as.p13, digits=3),"  theta23 = ", format(as.p23, digits=3),"\ntotal edf = ",format(x$t.edf, digits=3),"\n\n", sep="")
if(x$Model == "TSS") cat("n = ",x$n,"  n.sel1 = ",x$n.sel1,"  n.sel2 = ",x$n.sel2,"\ntheta12 = ", format(as.p12, digits=3),"  theta13 = ", format(as.p13, digits=3),"  theta23 = ", format(as.p23, digits=3),"\ntotal edf = ",format(x$t.edf, digits=3),"\n\n", sep="")

# tess can be improved with n.sel

invisible(x)

}

