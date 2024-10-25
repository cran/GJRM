pscr <- function(x, lind, m1l, m2l, cont1par, cont2par, cont3par, type = "copR", digits, signif.stars, m3l = NULL, lind2 = NULL, ...){


if(type == "ROY"){


  cat("\n\nEQUATION 1 - Switching Mechanism")
  cat("\nLink function for mu1:",m1l,"\n")
  cat("Formula: "); print(x$formula[[1]])
  cat("\n") 
  cat("Parametric coefficients:\n")
  printCoefmat(x$tableP1,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp1!=0){ 
    cat("Approximate significance of smooth terms:\n")
    printCoefmat(x$tableNP1,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }



  cat("\nEQUATION 2 - Regime 0")
  cat("\nLink function for mu2:",m2l,"\n")
  cat("Formula: "); print(x$formula[[2]])
  cat("\n")
  cat("Parametric coefficients:\n")
  printCoefmat(x$tableP2,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp2!=0){
    cat("Approximate significance of smooth terms:\n")
    printCoefmat(x$tableNP2,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }  
    
    

  cat("\nEQUATION 3 - Regime 1")
  cat("\nLink function for mu3:",m3l,"\n")
  cat("Formula: "); print(x$formula[[3]])
  cat("\n")
  cat("Parametric coefficients:\n")
  printCoefmat(x$tableP3,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp3!=0){
    cat("Approximate significance of smooth terms:\n")
    printCoefmat(x$tableNP3,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }   



#***********
if(length(x$formula) > 3 && x$margins[2] %in% cont1par && x$margins[3] %in% cont1par){

  cat("\nEQUATION 4 - Regime 0")
  cat("\nLink function for theta12:",lind,"\n") 
  cat("Formula: "); print(x$formula[[4]]) 
  cat("\n")
  cat("Parametric coefficients:\n")
  printCoefmat(x$tableP4,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp4!=0){
    cat("Approximate significance of smooth terms:\n")
    printCoefmat(x$tableNP4,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }    


  cat("\nEQUATION 5 - Regime 1")
  cat("\nLink function for theta13:",lind2,"\n") 
  cat("Formula: "); print(x$formula[[5]]) 
  cat("\n")
  cat("Parametric coefficients:\n")
  printCoefmat(x$tableP5,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp5!=0){
    cat("Approximate significance of smooth terms:\n")
    printCoefmat(x$tableNP5,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }    

}


#***********
if(length(x$formula) > 3 && x$margins[2] %in% c(cont2par) && x$margins[3] %in% c(cont2par)){


  cat("\nEQUATION 4 - Regime 0")
  if(x$margins[2] != "BE") cat("\nLink function for sigma2:","log","\n") else cat("\nLink function for sigma2:","qlogis","\n") 
  cat("Formula: "); print(x$formula[[4]]) 
  cat("\n")
  cat("Parametric coefficients:\n")
  printCoefmat(x$tableP4,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp4!=0){
    cat("Approximate significance of smooth terms:\n")
    printCoefmat(x$tableNP4,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }    


  cat("\nEQUATION 5 - Regime 1")
  if(x$margins[3] != "BE") cat("\nLink function for sigma3:","log","\n") else cat("\nLink function for sigma3:","qlogis","\n")
  cat("Formula: "); print(x$formula[[5]]) 
  cat("\n")
  cat("Parametric coefficients:\n")
  printCoefmat(x$tableP5,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp5!=0){
    cat("Approximate significance of smooth terms:\n")
    printCoefmat(x$tableNP5,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    } 


  cat("\nEQUATION 6 - Regime 0")
  cat("\nLink function for theta12:",lind,"\n") 
  cat("Formula: "); print(x$formula[[6]]) 
  cat("\n")
  cat("Parametric coefficients:\n")
  printCoefmat(x$tableP6,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp6!=0){
    cat("Approximate significance of smooth terms:\n")
    printCoefmat(x$tableNP6,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }    


  cat("\nEQUATION 7 - Regime 1")
  cat("\nLink function for theta13:",lind2,"\n") 
  cat("Formula: "); print(x$formula[[7]]) 
  cat("\n")
  cat("Parametric coefficients:\n")
  printCoefmat(x$tableP7,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp7!=0){
    cat("Approximate significance of smooth terms:\n")
    printCoefmat(x$tableNP7,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }    

}




#***********
if(length(x$formula) > 3 && x$margins[2] %in% c(cont3par) && x$margins[3] %in% c(cont3par)){


  cat("\nEQUATION 4 - Regime 0")
  cat("\nLink function for sigma2:","log","\n") 
  cat("Formula: "); print(x$formula[[4]]) 
  cat("\n")
  cat("Parametric coefficients:\n")
  printCoefmat(x$tableP4,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp4!=0){
    cat("Approximate significance of smooth terms:\n")
    printCoefmat(x$tableNP4,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }    


  cat("\nEQUATION 5 - Regime 1")
  cat("\nLink function for sigma3:","log","\n")
  cat("Formula: "); print(x$formula[[5]]) 
  cat("\n")
  cat("Parametric coefficients:\n")
  printCoefmat(x$tableP5,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp5!=0){
    cat("Approximate significance of smooth terms:\n")
    printCoefmat(x$tableNP5,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    } 
    
    
    
    
    
  cat("\nEQUATION 6 - Regime 0")
  cat("\nLink function for nu2:","log","\n")  
  cat("Formula: "); print(x$formula[[6]]) 
  cat("\n")
  cat("Parametric coefficients:\n")
  printCoefmat(x$tableP6,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp6!=0){
    cat("Approximate significance of smooth terms:\n")
    printCoefmat(x$tableNP6,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }    


  cat("\nEQUATION 7 - Regime 1")
  cat("\nLink function for nu3:","log","\n") 
  cat("Formula: "); print(x$formula[[7]]) 
  cat("\n")
  cat("Parametric coefficients:\n")
  printCoefmat(x$tableP7,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp7!=0){
    cat("Approximate significance of smooth terms:\n")
    printCoefmat(x$tableNP7,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    } 
    
    


  cat("\nEQUATION 8 - Regime 0")
  cat("\nLink function for theta12:",lind,"\n") 
  cat("Formula: "); print(x$formula[[8]]) 
  cat("\n")
  cat("Parametric coefficients:\n")
  printCoefmat(x$tableP8,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp8!=0){
    cat("Approximate significance of smooth terms:\n")
    printCoefmat(x$tableNP8,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }    


  cat("\nEQUATION 9 - Regime 1")
  cat("\nLink function for theta13:",lind2,"\n") 
  cat("Formula: "); print(x$formula[[9]]) 
  cat("\n")
  cat("Parametric coefficients:\n")
  printCoefmat(x$tableP9,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp9!=0){
    cat("Approximate significance of smooth terms:\n")
    printCoefmat(x$tableNP9,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }    

}






}








if(type != "ROY"){








if(type == "biv"){


  cat("\n\nEQUATION 1")  
  cat("\nLink function for mu1:",m1l,"\n")
  cat("Formula: "); print(x$formula[[1]])  
  cat("\n") 
  
  if(is.null(x$K1) || (!is.null(x$K1) && x$K1ns > 1) ){ 
  
  cat("Parametric coefficients:\n")
  printCoefmat(x$tableP1,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")
                                          } 

    if(x$l.sp1!=0){ 
    cat("Approximate significance of smooth terms:\n")
    printCoefmat(x$tableNP1,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }
    
    
  cat("\nEQUATION 2")
  cat("\nLink function for mu2:",m2l,"\n")
  cat("Formula: "); print(x$formula[[2]])   
  cat("\n")
  
  
  if(is.null(x$K2) || (!is.null(x$K2) && x$K2ns > 1) ){
  
  cat("Parametric coefficients:\n")
  printCoefmat(x$tableP2,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

                                        }
                                        
    if(x$l.sp2!=0){
    cat("Approximate significance of smooth terms:\n")
    printCoefmat(x$tableNP2,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }
    


if(!is.null(x$tableP3) && is.null(x$tableP4)  ){

  cat("\nEQUATION 3")
  cat("\nLink function for theta:",lind,"\n") 
  cat("Formula: "); print(x$formula[[3]])
  cat("\n")
  cat("Parametric coefficients:\n")
  printCoefmat(x$tableP3,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp3!=0){
    cat("Approximate significance of smooth terms:\n")
    printCoefmat(x$tableNP3,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }    
    
} 


if(!is.null(x$tableP3) && !is.null(x$tableP4) && is.null(x$tableP5)  ){


  cat("\nEQUATION 3")
  if(x$margins[2]!="BE") cat("\nLink function for sigma2:","log","\n") else cat("\nLink function for sigma2:","qlogis","\n")
  cat("Formula: "); print(x$formula[[3]])
  cat("\n")
  cat("Parametric coefficients:\n")
  printCoefmat(x$tableP3,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp3!=0){
    cat("Approximate significance of smooth terms:\n")
    printCoefmat(x$tableNP3,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }  


  cat("\nEQUATION 4")
  cat("\nLink function for theta:",lind,"\n") 
  cat("Formula: "); print(x$formula[[4]])
  cat("\n")
  cat("Parametric coefficients:\n")
  printCoefmat(x$tableP4,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp4!=0){
    cat("Approximate significance of smooth terms:\n")
    printCoefmat(x$tableNP4,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }    
    
} 


if(!is.null(x$tableP3) && !is.null(x$tableP4) && !is.null(x$tableP5)  ){


  cat("\nEQUATION 3")
  #cat("\nLink function for sigma:","log","\n") 
  
  #if(!(x$margins[2] %in% c("TW"))) 
  cat("\nLink function for sigma2:","log","\n") #else cat("\nLink function for sigma2:","qlogis","\n")  

  
  cat("Formula: "); print(x$formula[[3]])
  cat("\n")
  cat("Parametric coefficients:\n")
  printCoefmat(x$tableP3,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp3!=0){
    cat("Approximate significance of smooth terms:\n")
    printCoefmat(x$tableNP3,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }  
    
  cat("\nEQUATION 4")
  if(x$margins[2] %in% c("DAGUM","SM")) cat("\nLink function for nu2:","log","\n") 
  if(x$margins[2] %in% c("TW")) cat("\nLink function for nu2:","qlogis","\n") 
  cat("Formula: "); print(x$formula[[4]])
  cat("\n")
  cat("Parametric coefficients:\n")
  printCoefmat(x$tableP4,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp4!=0){
    cat("Approximate significance of smooth terms:\n")
    printCoefmat(x$tableNP4,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }     


  cat("\nEQUATION 5")
  cat("\nLink function for theta:",lind,"\n") 
  cat("Formula: "); print(x$formula[[5]])
  cat("\n")
  cat("Parametric coefficients:\n")
  printCoefmat(x$tableP5,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp5!=0){
    cat("Approximate significance of smooth terms:\n")
    printCoefmat(x$tableNP5,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }    
    
} 







}














if(type == "triv"){

  cat("\n\nEQUATION 1")  
  cat("\nLink function for mu1:",m1l,"\n")
  cat("Formula: "); print(x$formula[[1]]) 
  cat("\n") 
  cat("Parametric coefficients:\n")
  printCoefmat(x$tableP1,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp1!=0){ 
    cat("Approximate significance of smooth terms:\n")
    printCoefmat(x$tableNP1,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }
    
    
  cat("\nEQUATION 2")
  cat("\nLink function for mu2:",m2l,"\n")
  cat("Formula: "); print(x$formula[[2]])   
  cat("\n")
  cat("Parametric coefficients:\n")
  printCoefmat(x$tableP2,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp2!=0){
    cat("Approximate significance of smooth terms:\n")
    printCoefmat(x$tableNP2,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }
    
  cat("\nEQUATION 3")
  cat("\nLink function for mu3:",m3l,"\n")
  cat("Formula: "); print(x$formula[[3]])   
  cat("\n")
  cat("Parametric coefficients:\n")
  printCoefmat(x$tableP3,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp3!=0){
    cat("Approximate significance of smooth terms:\n")
    printCoefmat(x$tableNP3,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }


if(!is.null(x$tableP4)){


  cat("\nEQUATION 4\n")
  cat("Formula: "); print(x$formula[[4]])   
  cat("\n")
  cat("Parametric coefficients:\n")
  printCoefmat(x$tableP4,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp4!=0){
    cat("Approximate significance of smooth terms:\n")
    printCoefmat(x$tableNP4,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }

  cat("\nEQUATION 5\n")
  cat("Formula: "); print(x$formula[[5]])   
  cat("\n")
  cat("Parametric coefficients:\n")
  printCoefmat(x$tableP5,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp5!=0){
    cat("Approximate significance of smooth terms:\n")
    printCoefmat(x$tableNP5,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }


  cat("\nEQUATION 6\n")
  cat("Formula: "); print(x$formula[[6]])   
  cat("\n")
  cat("Parametric coefficients:\n")
  printCoefmat(x$tableP6,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp6!=0){
    cat("Approximate significance of smooth terms:\n")
    printCoefmat(x$tableNP6,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }






}





}








if(type == "copR"){


  doff <- "log(\u00B7 - 2)"



  cat("\n\nEQUATION 1")    
  if(x$surv.flex == FALSE || (x$surv.flex == TRUE && x$margins[1] %in% c(x$m2,x$m3)) ) cat("\nLink function for mu1:",m1l,"\n")
  if(x$surv.flex == TRUE && !(x$surv.flex == TRUE && x$margins[1] %in% c(x$m2,x$m3))) cat("\n")  

  cat("Formula: "); print(x$formula[[1]]) 
  cat("\n") 
  cat("Parametric coefficients:\n")
  printCoefmat(x$tableP1,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp1!=0){ 
    cat("Approximate significance of smooth terms:\n")
    printCoefmat(x$tableNP1,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }
    
    
  cat("\nEQUATION 2") 
  if(x$surv.flex == FALSE) cat("\nLink function for mu2:",m2l,"\n")
  if(x$surv.flex == TRUE) cat("\n")  

  cat("Formula: "); print(x$formula[[2]])   
  cat("\n")
  cat("Parametric coefficients:\n")
  printCoefmat(x$tableP2,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp2!=0){
    cat("Approximate significance of smooth terms:\n")
    printCoefmat(x$tableNP2,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }







  if( x$X3.null == FALSE ){ ###



if( x$margins[1] %in% c(x$m2,x$m3) && x$margins[2] %in% c(x$m2,x$m3) && x$BivD == "T"){



if( x$margins[1] %in% cont2par && x$margins[2] %in% cont2par){


  cat("\nEQUATION 3")
  if(x$margins[1] !="BE") cat("\nLink function for sigma1:","log","\n") else cat("\nLink function for sigma1:","qlogis","\n") 
  cat("Formula: "); print(x$formula[[3]]) 
  cat("\n")
    cat("Parametric coefficients:\n")
    printCoefmat(x$tableP3,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
    cat("\n")
  
      if(x$l.sp3!=0){
      cat("Approximate significance of smooth terms:\n")
      printCoefmat(x$tableNP3,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
      cat("\n")
      }
    
     
  cat("\nEQUATION 4")
  if(x$margins[2] !="BE") cat("\nLink function for sigma2:","log","\n") else cat("\nLink function for sigma2:","qlogis","\n")  
  cat("Formula: "); print(x$formula[[4]]) 
    cat("\n")
      cat("Parametric coefficients:\n")
      printCoefmat(x$tableP4,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
      cat("\n")
    
        if(x$l.sp4!=0){
        cat("Approximate significance of smooth terms:\n")
        printCoefmat(x$tableNP4,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
        cat("\n")
      }
  
  cat("\nEQUATION 5")
  cat("\nLink function for dof:",doff,"\n") 
  cat("Formula: "); print(x$formula[[5]])
    cat("\n")
      cat("Parametric coefficients:\n")
      printCoefmat(x$tableP5,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
      cat("\n")
    
        if(x$l.sp5!=0){
        cat("Approximate significance of smooth terms:\n")
        printCoefmat(x$tableNP5,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
        cat("\n")
      }  
      
      
  cat("\nEQUATION 6")
  cat("\nLink function for theta:",lind,"\n") 
  cat("Formula: "); print(x$formula[[6]])
    cat("\n")
      cat("Parametric coefficients:\n")
      printCoefmat(x$tableP6,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
      cat("\n")
    
        if(x$l.sp6!=0){
        cat("Approximate significance of smooth terms:\n")
        printCoefmat(x$tableNP6,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
        cat("\n")
      }      
      


}






     if( x$margins[1] %in% cont3par && x$margins[2] %in% cont3par){


  cat("\nEQUATION 3")
  cat("\nLink function for sigma1:","log","\n") 
  cat("Formula: "); print(x$formula[[3]])  
    cat("\n")
      cat("Parametric coefficients:\n")
      printCoefmat(x$tableP3,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
      cat("\n")
    
        if(x$l.sp3!=0){
        cat("Approximate significance of smooth terms:\n")
        printCoefmat(x$tableNP3,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
        cat("\n")
      }
     
  cat("\nEQUATION 4")
  cat("\nLink function for sigma2:","log","\n") 
  cat("Formula: "); print(x$formula[[4]])  
      cat("\n")
        cat("Parametric coefficients:\n")
        printCoefmat(x$tableP4,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
        cat("\n")
      
          if(x$l.sp4!=0){
          cat("Approximate significance of smooth terms:\n")
          printCoefmat(x$tableNP4,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
          cat("\n")
      }
  
  cat("\nEQUATION 5")
  cat("\nLink function for nu1:","log","\n") 
  cat("Formula: "); print(x$formula[[5]])  
      cat("\n")
        cat("Parametric coefficients:\n")
        printCoefmat(x$tableP5,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
        cat("\n")
      
          if(x$l.sp5!=0){
          cat("Approximate significance of smooth terms:\n")
          printCoefmat(x$tableNP5,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
          cat("\n")
      }  
     
  cat("\nEQUATION 6")
  cat("\nLink function for nu2:","log","\n") 
  cat("Formula: "); print(x$formula[[6]])  
      cat("\n")
        cat("Parametric coefficients:\n")
        printCoefmat(x$tableP6,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
        cat("\n")
      
          if(x$l.sp6!=0){
          cat("Approximate significance of smooth terms:\n")
          printCoefmat(x$tableNP6,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
          cat("\n")
      }  
  
  cat("\nEQUATION 7")
  cat("\nLink function for dof:",doff,"\n") 
  cat("Formula: "); print(x$formula[[7]])
      cat("\n")
        cat("Parametric coefficients:\n")
        printCoefmat(x$tableP7,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
        cat("\n")
      
          if(x$l.sp7!=0){
          cat("Approximate significance of smooth terms:\n")
          printCoefmat(x$tableNP7,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
          cat("\n")
      }  
      
  cat("\nEQUATION 8")
  cat("\nLink function for theta:",lind,"\n") 
  cat("Formula: "); print(x$formula[[8]])
      cat("\n")
        cat("Parametric coefficients:\n")
        printCoefmat(x$tableP8,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
        cat("\n")
      
          if(x$l.sp8!=0){
          cat("Approximate significance of smooth terms:\n")
          printCoefmat(x$tableNP8,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
          cat("\n")
      }      



}






     if( x$margins[1] %in% cont2par && x$margins[2] %in% cont3par){


  cat("\nEQUATION 3")
  if(x$margins[1] !="BE") cat("\nLink function for sigma1:","log","\n") else cat("\nLink function for sigma1:","qlogis","\n")
  cat("Formula: "); print(x$formula[[3]])  
    cat("\n")
      cat("Parametric coefficients:\n")
      printCoefmat(x$tableP3,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
      cat("\n")
    
        if(x$l.sp3!=0){
        cat("Approximate significance of smooth terms:\n")
        printCoefmat(x$tableNP3,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
        cat("\n")
      }
     
  cat("\nEQUATION 4")
  cat("\nLink function for sigma2:","log","\n") 
  cat("Formula: "); print(x$formula[[4]])  
      cat("\n")
        cat("Parametric coefficients:\n")
        printCoefmat(x$tableP4,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
        cat("\n")
      
          if(x$l.sp4!=0){
          cat("Approximate significance of smooth terms:\n")
          printCoefmat(x$tableNP4,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
          cat("\n")
      }
  
 
     
  cat("\nEQUATION 5")
  cat("\nLink function for nu2:","log","\n") 
  cat("Formula: "); print(x$formula[[5]])  
      cat("\n")
        cat("Parametric coefficients:\n")
        printCoefmat(x$tableP5,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
        cat("\n")
      
          if(x$l.sp5!=0){
          cat("Approximate significance of smooth terms:\n")
          printCoefmat(x$tableNP5,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
          cat("\n")
      }  
  
  cat("\nEQUATION 6")
  cat("\nLink function for dof:",doff,"\n") 
  cat("Formula: "); print(x$formula[[6]])
      cat("\n")
        cat("Parametric coefficients:\n")
        printCoefmat(x$tableP6,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
        cat("\n")
      
          if(x$l.sp6!=0){
          cat("Approximate significance of smooth terms:\n")
          printCoefmat(x$tableNP6,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
          cat("\n")
      }  
      
      
  cat("\nEQUATION 7")
  cat("\nLink function for theta:",lind,"\n") 
  cat("Formula: "); print(x$formula[[7]])
      cat("\n")
        cat("Parametric coefficients:\n")
        printCoefmat(x$tableP7,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
        cat("\n")
      
          if(x$l.sp7!=0){
          cat("Approximate significance of smooth terms:\n")
          printCoefmat(x$tableNP7,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
          cat("\n")
      }       



}









     if( x$margins[1] %in% cont3par && x$margins[2] %in% cont2par){


  cat("\nEQUATION 3")
  cat("\nLink function for sigma1:","log","\n") 
  cat("Formula: "); print(x$formula[[3]])  
    cat("\n")
      cat("Parametric coefficients:\n")
      printCoefmat(x$tableP3,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
      cat("\n")
    
        if(x$l.sp3!=0){
        cat("Approximate significance of smooth terms:\n")
        printCoefmat(x$tableNP3,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
        cat("\n")
      }
     
  cat("\nEQUATION 4")
  if(x$margins[2] !="BE") cat("\nLink function for sigma2:","log","\n") else cat("\nLink function for sigma2:","qlogis","\n")
  cat("Formula: "); print(x$formula[[4]])  
      cat("\n")
        cat("Parametric coefficients:\n")
        printCoefmat(x$tableP4,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
        cat("\n")
      
          if(x$l.sp4!=0){
          cat("Approximate significance of smooth terms:\n")
          printCoefmat(x$tableNP4,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
          cat("\n")
      }
  
  cat("\nEQUATION 5")
  cat("\nLink function for nu1:","log","\n") 
  cat("Formula: "); print(x$formula[[5]])  
      cat("\n")
        cat("Parametric coefficients:\n")
        printCoefmat(x$tableP5,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
        cat("\n")
      
          if(x$l.sp5!=0){
          cat("Approximate significance of smooth terms:\n")
          printCoefmat(x$tableNP5,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
          cat("\n")
      }  
     
 
  
  cat("\nEQUATION 6")
  cat("\nLink function for dof:",doff,"\n") 
  cat("Formula: "); print(x$formula[[6]])
      cat("\n")
        cat("Parametric coefficients:\n")
        printCoefmat(x$tableP6,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
        cat("\n")
      
          if(x$l.sp6!=0){
          cat("Approximate significance of smooth terms:\n")
          printCoefmat(x$tableNP6,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
          cat("\n")
      }  



  cat("\nEQUATION 7")
  cat("\nLink function for theta:",lind,"\n") 
  cat("Formula: "); print(x$formula[[7]])
      cat("\n")
        cat("Parametric coefficients:\n")
        printCoefmat(x$tableP7,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
        cat("\n")
      
          if(x$l.sp7!=0){
          cat("Approximate significance of smooth terms:\n")
          printCoefmat(x$tableNP7,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
          cat("\n")
      }  



}###







}else{



if( x$margins[1] %in% cont2par && x$margins[2] %in% cont1par && x$surv.flex == TRUE){


  cat("\nEQUATION 3")
  if(!(x$margins[1] %in% c("BE"))) cat("\nLink function for sigma1:","log","\n") else cat("\nLink function for sigma1:","qlogis","\n") 
  cat("Formula: "); print(x$formula[[3]]) 
  cat("\n")
    cat("Parametric coefficients:\n")
    printCoefmat(x$tableP3,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
    cat("\n")
  
      if(x$l.sp3!=0){
      cat("Approximate significance of smooth terms:\n")
      printCoefmat(x$tableNP3,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
      cat("\n")
      }
    
  
  cat("\nEQUATION 4")
  cat("\nLink function for theta:",lind,"\n") 
  cat("Formula: "); print(x$formula[[4]])
    cat("\n")
      cat("Parametric coefficients:\n")
      printCoefmat(x$tableP4,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
      cat("\n")
    
        if(x$l.sp4!=0){
        cat("Approximate significance of smooth terms:\n")
        printCoefmat(x$tableNP4,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
        cat("\n")
      }  


}








if( x$margins[1] %in% cont3par && x$margins[2] %in% cont1par && x$surv.flex == TRUE){


  cat("\nEQUATION 3")
  cat("\nLink function for sigma1:","log","\n") 
  cat("Formula: "); print(x$formula[[3]])  
    cat("\n")
      cat("Parametric coefficients:\n")
      printCoefmat(x$tableP3,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
      cat("\n")
    
        if(x$l.sp3!=0){
        cat("Approximate significance of smooth terms:\n")
        printCoefmat(x$tableNP3,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
        cat("\n")
      }
     

  cat("\nEQUATION 4")
  cat("\nLink function for nu1:","log","\n") 
  cat("Formula: "); print(x$formula[[4]])  
      cat("\n")
        cat("Parametric coefficients:\n")
        printCoefmat(x$tableP4,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
        cat("\n")
      
          if(x$l.sp4!=0){
          cat("Approximate significance of smooth terms:\n")
          printCoefmat(x$tableNP4,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
          cat("\n")
      }  
     

  cat("\nEQUATION 5")
  cat("\nLink function for theta:",lind,"\n") 
  cat("Formula: "); print(x$formula[[5]])
      cat("\n")
        cat("Parametric coefficients:\n")
        printCoefmat(x$tableP5,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
        cat("\n")
      
          if(x$l.sp5!=0){
          cat("Approximate significance of smooth terms:\n")
          printCoefmat(x$tableNP5,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
          cat("\n")
      }  

}







if( x$margins[1] %in% cont2par && x$margins[2] %in% cont2par){


  cat("\nEQUATION 3")
  if(x$margins[1] !="BE") cat("\nLink function for sigma1:","log","\n") else cat("\nLink function for sigma1:","qlogis","\n") 
  cat("Formula: "); print(x$formula[[3]]) 
  cat("\n")
    cat("Parametric coefficients:\n")
    printCoefmat(x$tableP3,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
    cat("\n")
  
      if(x$l.sp3!=0){
      cat("Approximate significance of smooth terms:\n")
      printCoefmat(x$tableNP3,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
      cat("\n")
      }
    
     
  cat("\nEQUATION 4")
  if(x$margins[2] !="BE") cat("\nLink function for sigma2:","log","\n") else cat("\nLink function for sigma2:","qlogis","\n")  
  cat("Formula: "); print(x$formula[[4]]) 
    cat("\n")
      cat("Parametric coefficients:\n")
      printCoefmat(x$tableP4,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
      cat("\n")
    
        if(x$l.sp4!=0){
        cat("Approximate significance of smooth terms:\n")
        printCoefmat(x$tableNP4,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
        cat("\n")
      }
  
  cat("\nEQUATION 5")
  cat("\nLink function for theta:",lind,"\n") 
  cat("Formula: "); print(x$formula[[5]])
    cat("\n")
      cat("Parametric coefficients:\n")
      printCoefmat(x$tableP5,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
      cat("\n")
    
        if(x$l.sp5!=0){
        cat("Approximate significance of smooth terms:\n")
        printCoefmat(x$tableNP5,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
        cat("\n")
      }  


}





    if( x$margins[1] %in% cont1par && x$margins[2] %in% cont1par ){

 
  
  cat("\nEQUATION 3")
  cat("\nLink function for theta:",lind,"\n") 
  cat("Formula: "); print(x$formula[[3]])
    cat("\n")
      cat("Parametric coefficients:\n")
      printCoefmat(x$tableP3,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
      cat("\n")
    
        if(x$l.sp3!=0){
        cat("Approximate significance of smooth terms:\n")
        printCoefmat(x$tableNP3,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
        cat("\n")
      }  


}










     if( x$margins[1] %in% cont1par && x$margins[2] %in% cont2par ){

 
  cat("\nEQUATION 3")
  if(x$margins[2] !="BE") cat("\nLink function for sigma2:","log","\n") else cat("\nLink function for sigma2:","qlogis","\n")  
  cat("Formula: "); print(x$formula[[3]]) 
    cat("\n")
      cat("Parametric coefficients:\n")
      printCoefmat(x$tableP3,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
      cat("\n")
    
        if(x$l.sp3!=0){
        cat("Approximate significance of smooth terms:\n")
        printCoefmat(x$tableNP3,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
        cat("\n")
      }
  
  cat("\nEQUATION 4")
  cat("\nLink function for theta:",lind,"\n") 
  cat("Formula: "); print(x$formula[[4]])
    cat("\n")
      cat("Parametric coefficients:\n")
      printCoefmat(x$tableP4,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
      cat("\n")
    
        if(x$l.sp4!=0){
        cat("Approximate significance of smooth terms:\n")
        printCoefmat(x$tableNP4,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
        cat("\n")
      }  


}













     if( x$margins[1] %in% cont1par && x$margins[2] %in% cont3par ){


  cat("\nEQUATION 3")
  cat("\nLink function for sigma2:","log","\n") 
  cat("Formula: "); print(x$formula[[3]])  
      cat("\n")
        cat("Parametric coefficients:\n")
        printCoefmat(x$tableP3,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
        cat("\n")
      
          if(x$l.sp3!=0){
          cat("Approximate significance of smooth terms:\n")
          printCoefmat(x$tableNP3,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
          cat("\n")
      }
  
 
     
  cat("\nEQUATION 4")
  cat("\nLink function for nu2:","log","\n") 
  cat("Formula: "); print(x$formula[[4]])  
      cat("\n")
        cat("Parametric coefficients:\n")
        printCoefmat(x$tableP4,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
        cat("\n")
      
          if(x$l.sp4!=0){
          cat("Approximate significance of smooth terms:\n")
          printCoefmat(x$tableNP4,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
          cat("\n")
      }  
  
  cat("\nEQUATION 5")
  cat("\nLink function for theta:",lind,"\n") 
  cat("Formula: "); print(x$formula[[5]])
      cat("\n")
        cat("Parametric coefficients:\n")
        printCoefmat(x$tableP5,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
        cat("\n")
      
          if(x$l.sp5!=0){
          cat("Approximate significance of smooth terms:\n")
          printCoefmat(x$tableNP5,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
          cat("\n")
      }  



}







     if( x$margins[1] %in% cont3par && x$margins[2] %in% cont3par){


  cat("\nEQUATION 3")
  cat("\nLink function for sigma1:","log","\n") 
  cat("Formula: "); print(x$formula[[3]])  
    cat("\n")
      cat("Parametric coefficients:\n")
      printCoefmat(x$tableP3,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
      cat("\n")
    
        if(x$l.sp3!=0){
        cat("Approximate significance of smooth terms:\n")
        printCoefmat(x$tableNP3,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
        cat("\n")
      }
     
  cat("\nEQUATION 4")
  cat("\nLink function for sigma2:","log","\n") 
  cat("Formula: "); print(x$formula[[4]])  
      cat("\n")
        cat("Parametric coefficients:\n")
        printCoefmat(x$tableP4,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
        cat("\n")
      
          if(x$l.sp4!=0){
          cat("Approximate significance of smooth terms:\n")
          printCoefmat(x$tableNP4,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
          cat("\n")
      }
  
  cat("\nEQUATION 5")
  cat("\nLink function for nu1:","log","\n") 
  cat("Formula: "); print(x$formula[[5]])  
      cat("\n")
        cat("Parametric coefficients:\n")
        printCoefmat(x$tableP5,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
        cat("\n")
      
          if(x$l.sp5!=0){
          cat("Approximate significance of smooth terms:\n")
          printCoefmat(x$tableNP5,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
          cat("\n")
      }  
     
  cat("\nEQUATION 6")
  cat("\nLink function for nu2:","log","\n") 
  cat("Formula: "); print(x$formula[[6]])  
      cat("\n")
        cat("Parametric coefficients:\n")
        printCoefmat(x$tableP6,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
        cat("\n")
      
          if(x$l.sp6!=0){
          cat("Approximate significance of smooth terms:\n")
          printCoefmat(x$tableNP6,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
          cat("\n")
      }  
  
  cat("\nEQUATION 7")
  cat("\nLink function for theta:",lind,"\n") 
  cat("Formula: "); print(x$formula[[7]])
      cat("\n")
        cat("Parametric coefficients:\n")
        printCoefmat(x$tableP7,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
        cat("\n")
      
          if(x$l.sp7!=0){
          cat("Approximate significance of smooth terms:\n")
          printCoefmat(x$tableNP7,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
          cat("\n")
      }  



}









     if( x$margins[1] %in% cont2par && x$margins[2] %in% cont3par){


  cat("\nEQUATION 3")
  if(x$margins[1] !="BE") cat("\nLink function for sigma1:","log","\n") else cat("\nLink function for sigma1:","qlogis","\n")
  cat("Formula: "); print(x$formula[[3]])  
    cat("\n")
      cat("Parametric coefficients:\n")
      printCoefmat(x$tableP3,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
      cat("\n")
    
        if(x$l.sp3!=0){
        cat("Approximate significance of smooth terms:\n")
        printCoefmat(x$tableNP3,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
        cat("\n")
      }
     
  cat("\nEQUATION 4")
  cat("\nLink function for sigma2:","log","\n") 
  cat("Formula: "); print(x$formula[[4]])  
      cat("\n")
        cat("Parametric coefficients:\n")
        printCoefmat(x$tableP4,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
        cat("\n")
      
          if(x$l.sp4!=0){
          cat("Approximate significance of smooth terms:\n")
          printCoefmat(x$tableNP4,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
          cat("\n")
      }
  
 
     
  cat("\nEQUATION 5")
  cat("\nLink function for nu2:","log","\n") 
  cat("Formula: "); print(x$formula[[5]])  
      cat("\n")
        cat("Parametric coefficients:\n")
        printCoefmat(x$tableP5,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
        cat("\n")
      
          if(x$l.sp5!=0){
          cat("Approximate significance of smooth terms:\n")
          printCoefmat(x$tableNP5,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
          cat("\n")
      }  
  
  cat("\nEQUATION 6")
  cat("\nLink function for theta:",lind,"\n") 
  cat("Formula: "); print(x$formula[[6]])
      cat("\n")
        cat("Parametric coefficients:\n")
        printCoefmat(x$tableP6,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
        cat("\n")
      
          if(x$l.sp6!=0){
          cat("Approximate significance of smooth terms:\n")
          printCoefmat(x$tableNP6,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
          cat("\n")
      }  



}









     if( x$margins[1] %in% cont3par && x$margins[2] %in% cont2par){


  cat("\nEQUATION 3")
  cat("\nLink function for sigma1:","log","\n") 
  cat("Formula: "); print(x$formula[[3]])  
    cat("\n")
      cat("Parametric coefficients:\n")
      printCoefmat(x$tableP3,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
      cat("\n")
    
        if(x$l.sp3!=0){
        cat("Approximate significance of smooth terms:\n")
        printCoefmat(x$tableNP3,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
        cat("\n")
      }
     
  cat("\nEQUATION 4")
  if(x$margins[2] !="BE") cat("\nLink function for sigma2:","log","\n") else cat("\nLink function for sigma2:","qlogis","\n")
  cat("Formula: "); print(x$formula[[4]])  
      cat("\n")
        cat("Parametric coefficients:\n")
        printCoefmat(x$tableP4,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
        cat("\n")
      
          if(x$l.sp4!=0){
          cat("Approximate significance of smooth terms:\n")
          printCoefmat(x$tableNP4,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
          cat("\n")
      }
  
  cat("\nEQUATION 5")
  cat("\nLink function for nu1:","log","\n") 
  cat("Formula: "); print(x$formula[[5]])  
      cat("\n")
        cat("Parametric coefficients:\n")
        printCoefmat(x$tableP5,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
        cat("\n")
      
          if(x$l.sp5!=0){
          cat("Approximate significance of smooth terms:\n")
          printCoefmat(x$tableNP5,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
          cat("\n")
      }  
     
 
  
  cat("\nEQUATION 6")
  cat("\nLink function for theta:",lind,"\n") 
  cat("Formula: "); print(x$formula[[6]])
      cat("\n")
        cat("Parametric coefficients:\n")
        printCoefmat(x$tableP6,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
        cat("\n")
      
          if(x$l.sp6!=0){
          cat("Approximate significance of smooth terms:\n")
          printCoefmat(x$tableNP6,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
          cat("\n")
      }  



}






} # null



} # t





} # cReg













if(type == "gamls"){




  cat("\n\nEQUATION 1")    
  if(x$surv.flex == FALSE && x$margins[1] != "DGP0") cat("\nLink function for mu:",m1l,"\n")
  if(x$surv.flex == FALSE && x$margins[1] == "DGP0") cat("\nLink function for sigma:",m1l,"\n")
  
  if(x$surv.flex == TRUE) cat("\n")

  cat("Formula: "); print(x$formula[[1]])  
  cat("\n") 
  cat("Parametric coefficients:\n")
  printCoefmat(x$tableP1,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp1!=0){ 
    cat("Approximate significance of smooth terms:\n")
    printCoefmat(x$tableNP1,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }
    
    

  if( x$X2.null == FALSE ){
  

  if( x$margins[1] %in% cont2par){
  
  cat("\nEQUATION 2")
  
  if( !(x$margins[1] %in%c("BE")) ) cat("\nLink function for sigma:","log","\n") else cat("\nLink function for sigma:","qlogis","\n") 
  
  #if(x$margins[1] !="BE" && !(x$margins[1] %in% c("GP","DGP","DGPII"))) cat("\nLink function for sigma:","log","\n")
  
  #if(x$margins[1] =="BE") cat("\nLink function for sigma:","qlogis","\n") 
  
  #if(x$margins[1] %in% c("GP","DGP","DGPII")) cat("\nLink function for sigma:","identity","\n")   
  

  cat("Formula: "); print(x$formula[[2]]) 
  cat("\n")
    cat("Parametric coefficients:\n")
    printCoefmat(x$tableP2,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
    cat("\n")
  
      if(x$l.sp2!=0){
      cat("Approximate significance of smooth terms:\n")
      printCoefmat(x$tableNP2,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
      cat("\n")
                    }
                                  }
                                  
                                  
                                  

if( x$margins[1] %in% cont3par){

  cat("\nEQUATION 2")
  
  
  #cat("\nLink function for sigma:","log","\n") 
  
  if( !(x$margins[1] %in%c("BE")) ) cat("\nLink function for sigma:","log","\n") else cat("\nLink function for sigma:","qlogis","\n") 

  
  
  
  
  cat("Formula: "); print(x$formula[[2]])  
    cat("\n")
      cat("Parametric coefficients:\n")
      printCoefmat(x$tableP2,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
      cat("\n")
    
        if(x$l.sp2!=0){
        cat("Approximate significance of smooth terms:\n")
        printCoefmat(x$tableNP2,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
        cat("\n")
                      }
 
  cat("\nEQUATION 3")
  
  
  if( !(x$margins[1] %in% c("TW")) ) cat("\nLink function for nu:","log","\n") else cat("\nLink function for nu:","qlogis","\n")  
  
  
  cat("Formula: "); print(x$formula[[3]])  
      cat("\n")
        cat("Parametric coefficients:\n")
        printCoefmat(x$tableP3,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
        cat("\n")
      
          if(x$l.sp3!=0){
          cat("Approximate significance of smooth terms:\n")
          printCoefmat(x$tableNP3,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
          cat("\n")
                        }  
                        
                        
                                  }




}                                  








}







if(type == "copSS"){





  cat("\n\nEQUATION 1")  
  cat("\nLink function for mu1:",m1l,"\n")
  cat("Formula: "); print(x$formula[[1]]) 
  cat("\n") 
  cat("Parametric coefficients:\n")
  printCoefmat(x$tableP1,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp1!=0){ 
    cat("Approximate significance of smooth terms:\n")
    printCoefmat(x$tableNP1,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }
    
    
  cat("\nEQUATION 2")
  cat("\nLink function for mu2:",m2l,"\n")
  cat("Formula: "); print(x$formula[[2]])    
  cat("\n")
  cat("Parametric coefficients:\n")
  printCoefmat(x$tableP2,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp2!=0){
    cat("Approximate significance of smooth terms:\n")
    printCoefmat(x$tableNP2,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }
    


if(!is.null(x$tableP3) && is.null(x$tableP4)  ){

  cat("\nEQUATION 3")
  cat("\nLink function for theta:",lind,"\n") 
  cat("Formula: "); print(x$formula[[3]])
  cat("\n")
  cat("Parametric coefficients:\n")
  printCoefmat(x$tableP3,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp3!=0){
    cat("Approximate significance of smooth terms:\n")
    printCoefmat(x$tableNP3,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }    
    
} 


if(!is.null(x$tableP3) && !is.null(x$tableP4) && is.null(x$tableP5)  ){


  cat("\nEQUATION 3")
  if(!(x$margins[2]%in%c("BE"))) cat("\nLink function for sigma2:","log","\n") else cat("\nLink function for sigma2:","qlogis","\n")
  cat("Formula: "); print(x$formula[[3]])
  cat("\n")
  cat("Parametric coefficients:\n")
  printCoefmat(x$tableP3,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp3!=0){
    cat("Approximate significance of smooth terms:\n")
    printCoefmat(x$tableNP3,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }  


  cat("\nEQUATION 4")
  cat("\nLink function for theta:",lind,"\n") 
  cat("Formula: "); print(x$formula[[4]])
  cat("\n")
  cat("Parametric coefficients:\n")
  printCoefmat(x$tableP4,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp4!=0){
    cat("Approximate significance of smooth terms:\n")
    printCoefmat(x$tableNP4,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }    
    
} 


if(!is.null(x$tableP3) && !is.null(x$tableP4) && !is.null(x$tableP5)  ){


  cat("\nEQUATION 3")
 # cat("\nLink function for sigma2:","log","\n") 
  
    #if(!(x$margins[2] %in% c("TW"))) 
    cat("\nLink function for sigma2:","log","\n") #else cat("\nLink function for sigma2:","qlogis","\n")  

  
  cat("Formula: "); print(x$formula[[3]])
  cat("\n")
  cat("Parametric coefficients:\n")
  printCoefmat(x$tableP3,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp3!=0){
    cat("Approximate significance of smooth terms:\n")
    printCoefmat(x$tableNP3,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }  
    
  cat("\nEQUATION 4")
  if(x$margins[2] %in% c("DAGUM","SM")) cat("\nLink function for nu2:","log","\n")  
  if(x$margins[2] %in% c("TW"))         cat("\nLink function for nu2:","qlogis","\n")  
  
  cat("Formula: "); print(x$formula[[4]])
  cat("\n")
  cat("Parametric coefficients:\n")
  printCoefmat(x$tableP4,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp4!=0){
    cat("Approximate significance of smooth terms:\n")
    printCoefmat(x$tableNP4,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }     


  cat("\nEQUATION 5")
  cat("\nLink function for theta:",lind,"\n") 
  cat("Formula: "); print(x$formula[[5]])
  cat("\n")
  cat("Parametric coefficients:\n")
  printCoefmat(x$tableP5,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp5!=0){
    cat("Approximate significance of smooth terms:\n")
    printCoefmat(x$tableNP5,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }    
    
} 















}


} # !ROY

                
}

