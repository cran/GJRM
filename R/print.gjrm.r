print.gjrm <- function(x, ...){

   
 ppR <- pp(x)  
 
 cont1par <- ppR$cont1par
 cont2par <- ppR$cont2par
 cont3par <- ppR$cont3par
 cop      <- ppR$cop
 lind     <- ppR$lind
 m1l      <- ppR$m1l
 m2l      <- ppR$m2l 

  cp <- "theta = "; as.p <- x$theta.a
  s1 <- "sigma.1 = "; s1.p <- x$sigma21.a
  s2 <- "sigma.2 = "; s2.p <- x$sigma22.a  
  n1 <- "nu.1 = "; n1.p <- x$nu1.a
  n2 <- "nu.2 = "; n2.p <- x$nu2.a   
  
  doff <- "log(\u00B7 - 2)"
  

  cat("\nCOPULA:  ",cop)
  
  
  pscr0(x, type = "copR")    
  



  cat("\n\nEQUATION 1")
  if(x$surv.flex == FALSE )                                             cat("\nLink function for mu.1:",m1l,"\n")
  if(x$surv.flex == TRUE && x$VC$margins[1] %in% c(x$VC$m2,x$VC$m3) )   cat("\nLink function for mu.1:",m1l,"\n")
  if(x$surv.flex == TRUE && !(x$VC$margins[1] %in% c(x$VC$m2,x$VC$m3))) cat("\n")

  cat("Formula: "); print(x$gam1$formula) 
  
  cat("\nEQUATION 2") 
  if(x$surv.flex == FALSE) cat("\nLink function for mu.2:",m2l,"\n")
  if(x$surv.flex == TRUE)  cat("\n")  
  cat("Formula: "); print(x$gam2$formula) 
  
  
  
  if(x$margins[1] %in% c("TW")  ||  x$margins[2] %in% c("TW")) stop("Amend function to accommodate Tweedie.")  
  
  
  
  if( !is.null(x$X3) ){
  
 

 
 
if(x$BivD == "T" && x$VC$margins[1] %in% c(x$VC$m2,x$VC$m3) && x$VC$margins[2] %in% c(x$VC$m2,x$VC$m3)){


     
  if( x$margins[1] %in% cont2par && x$margins[2] %in% cont2par){
  

  cat("\nEQUATION 3")
  if(x$margins[1] !="BE") cat("\nLink function for sigma.1:","log","\n") else cat("\nLink function for sigma.1:","qlogis","\n") 
  cat("Formula: "); print(x$formula[[3]])   
  
  cat("\nEQUATION 4")
  if(x$margins[2] !="BE") cat("\nLink function for sigma.2:","log","\n") else cat("\nLink function for sigma.2:","qlogis","\n") 
  cat("Formula: "); print(x$formula[[4]])   
  
  cat("\nEQUATION 5")
  cat("\nLink function for dof:",doff,"\n") 
  cat("Formula: "); print(x$formula[[5]]) 
  
  cat("\nEQUATION 6")
  cat("\nLink function for theta:",lind,"\n") 
  cat("Formula: "); print(x$formula[[6]])   
     
     }
          
  

     
if( x$margins[1] %in% cont3par && x$margins[2] %in% cont3par){
     
  cat("\nEQUATION 3")
  cat("\nLink function for sigma.1:","log","\n") 
  cat("Formula: "); print(x$formula[[3]]) 
     
  cat("\nEQUATION 4")
  cat("\nLink function for sigma.2:","log","\n") 
  cat("Formula: "); print(x$formula[[4]])  
  
  cat("\nEQUATION 5")
  cat("\nLink function for nu.1:","log","\n") 
  cat("Formula: "); print(x$formula[[5]])  
     
  cat("\nEQUATION 6")
  cat("\nLink function for nu.2:","log","\n") 
  cat("Formula: "); print(x$formula[[6]])    
  
  cat("\nEQUATION 7")
  cat("\nLink function for dof:",doff,"\n") 
  cat("Formula: "); print(x$formula[[7]]) 
  
  cat("\nEQUATION 8")
  cat("\nLink function for theta:",lind,"\n") 
  cat("Formula: "); print(x$formula[[8]]) 
  
     
     }     
    
  

  
 
     
     
     if( x$margins[1] %in% cont3par && x$margins[2] %in% cont2par){
     
  cat("\nEQUATION 3")
  cat("\nLink function for sigma.1:","log","\n") 
  cat("Formula: "); print(x$formula[[3]])   
     
  cat("\nEQUATION 4")
  if(x$margins[2] !="BE") cat("\nLink function for sigma.2:","log","\n") else cat("\nLink function for sigma.2:","qlogis","\n")
  cat("Formula: "); print(x$formula[[4]])  
  
  cat("\nEQUATION 5")
  cat("\nLink function for nu.1:","log","\n") 
  cat("Formula: "); print(x$formula[[5]]) 
         
  cat("\nEQUATION 6")
  cat("\nLink function for dof:",doff,"\n") 
  cat("Formula: "); print(x$formula[[6]]) 
  
    cat("\nEQUATION 7")
    cat("\nLink function for theta:",lind,"\n") 
  cat("Formula: "); print(x$formula[[7]]) 
     
     }     
      
  
  
  
  
    

       
       
 if( x$margins[1] %in% cont2par && x$margins[2] %in% cont3par){
       
       
    cat("\nEQUATION 3")
    if(x$margins[1] !="BE") cat("\nLink function for sigma.1:","log","\n") else cat("\nLink function for sigma.1:","qlogis","\n")
    cat("Formula: "); print(x$formula[[3]])  
       
    cat("\nEQUATION 4")
    cat("\nLink function for sigma.2:","log","\n") 
    cat("Formula: "); print(x$formula[[4]]) 
       
    cat("\nEQUATION 5")
    cat("\nLink function for nu.2:","log","\n") 
    cat("Formula: "); print(x$formula[[5]])      
    
    cat("\nEQUATION 6")
    cat("\nLink function for dof:",doff,"\n") 
    cat("Formula: "); print(x$formula[[6]]) 
    
    cat("\nEQUATION 7")
    cat("\nLink function for theta:",lind,"\n") 
    cat("Formula: "); print(x$formula[[7]])     
       
       }       
    
  
##########################

}else{
 
 
  if( x$margins[1] %in% cont1par && x$margins[2] %in% cont1par ){
  

  cat("\nEQUATION 3")
  cat("\nLink function for theta:",lind,"\n") 
  cat("Formula: "); print(x$formula[[3]]) 
     
     } 
 
 

  if( x$margins[1] %in% cont1par && x$margins[2] %in% cont2par ){
  
  
  cat("\nEQUATION 3")
  if(x$margins[2] !="BE") cat("\nLink function for sigma.2:","log","\n") else cat("\nLink function for sigma.2:","qlogis","\n") 
  cat("Formula: "); print(x$formula[[3]])  
  
 
  cat("\nEQUATION 4")
  cat("\nLink function for theta:",lind,"\n") 
  cat("Formula: "); print(x$formula[[4]]) 
     
     } 
     
     
     
     

    
 if( x$margins[1] %in% cont1par && x$margins[2] %in% cont3par ){
       

    cat("\nEQUATION 3")
    cat("\nLink function for sigma.2:","log","\n") 
    cat("Formula: "); print(x$formula[[3]]) 
       
    cat("\nEQUATION 4")
    cat("\nLink function for nu.2:","log","\n") 
    cat("Formula: "); print(x$formula[[4]])     
    
    cat("\nEQUATION 5")
    cat("\nLink function for theta:",lind,"\n") 
    cat("Formula: "); print(x$formula[[5]]) 
       
       }
          
 
 
 
 

  if( x$margins[1] %in% cont2par && x$margins[2] %in% cont2par){
  

  cat("\nEQUATION 3")
  if(x$margins[1] !="BE") cat("\nLink function for sigma.1:","log","\n") else cat("\nLink function for sigma.1:","qlogis","\n") 
  cat("Formula: "); print(x$formula[[3]])   
  
  cat("\nEQUATION 4")
  if(x$margins[2] !="BE") cat("\nLink function for sigma.2:","log","\n") else cat("\nLink function for sigma.2:","qlogis","\n") 
  cat("Formula: "); print(x$formula[[4]])   
  
  cat("\nEQUATION 5")
  cat("\nLink function for theta:",lind,"\n") 
  cat("Formula: "); print(x$formula[[5]]) 
     
     }
     
     
  if( x$margins[1] %in% cont2par && x$margins[2] %in% cont1par && x$surv.flex == TRUE){
  

  cat("\nEQUATION 3")
  if(x$margins[1] !="BE") cat("\nLink function for sigma.1:","log","\n") else cat("\nLink function for sigma.1:","qlogis","\n") 
  cat("Formula: "); print(x$formula[[3]])   
    
  cat("\nEQUATION 4")
  cat("\nLink function for theta:",lind,"\n") 
  cat("Formula: "); print(x$formula[[4]]) 
     
     }     
     
  if( x$margins[1] %in% cont3par && x$margins[2] %in% cont1par && x$surv.flex == TRUE){
  

  cat("\nEQUATION 3")
  if(x$margins[1] !="BE") cat("\nLink function for sigma.1:","log","\n") else cat("\nLink function for sigma.1:","qlogis","\n") 
  cat("Formula: "); print(x$formula[[3]]) 
  
  cat("\nEQUATION 4")
  cat("\nLink function for nu.1:","log","\n") 
  cat("Formula: "); print(x$formula[[4]])  
    
  cat("\nEQUATION 5")
  cat("\nLink function for theta:",lind,"\n") 
  cat("Formula: "); print(x$formula[[5]]) 
     
     }      
     
     
     
          
  
if( x$margins[1] %in% cont3par && x$margins[2] %in% cont3par){
     
  cat("\nEQUATION 3")
  cat("\nLink function for sigma.1:","log","\n") 
  cat("Formula: "); print(x$formula[[3]])  
     
  cat("\nEQUATION 4")
  cat("\nLink function for sigma.2:","log","\n") 
  cat("Formula: "); print(x$formula[[4]]) 
  
  cat("\nEQUATION 5")
  cat("\nLink function for nu.1:","log","\n") 
  cat("Formula: "); print(x$formula[[5]])  
     
  cat("\nEQUATION 6")
  cat("\nLink function for nu.2:","log","\n") 
  cat("Formula: "); print(x$formula[[6]])    
  
  cat("\nEQUATION 7")
  cat("\nLink function for theta:",lind,"\n") 
  cat("Formula: "); print(x$formula[[7]]) 
     
     }
     
    
  

  
     if( x$margins[1] %in% cont3par && x$margins[2] %in% cont2par){
     
  cat("\nEQUATION 3")
  cat("\nLink function for sigma.1:","log","\n") 
  cat("Formula: "); print(x$formula[[3]])   
     
  cat("\nEQUATION 4")
  if(x$margins[2] !="BE") cat("\nLink function for sigma.2:","log","\n") else cat("\nLink function for sigma.2:","qlogis","\n")
  cat("Formula: "); print(x$formula[[4]])  
  
  cat("\nEQUATION 5")
  cat("\nLink function for nu.1:","log","\n") 
  cat("Formula: "); print(x$formula[[5]]) 
         
  cat("\nEQUATION 6")
  cat("\nLink function for theta:",lind,"\n") 
  cat("Formula: "); print(x$formula[[6]]) 
     
     }
     
     
      
  
  
  
  
    
 if( x$margins[1] %in% cont2par && x$margins[2] %in% cont3par){
       
       
    cat("\nEQUATION 3")
    if(x$margins[1] !="BE") cat("\nLink function for sigma.1:","log","\n") else cat("\nLink function for sigma.1:","qlogis","\n")
    cat("Formula: "); print(x$formula[[3]])  
       
    cat("\nEQUATION 4")
    cat("\nLink function for sigma.2:","log","\n") 
    cat("Formula: "); print(x$formula[[4]]) 
       
    cat("\nEQUATION 5")
    cat("\nLink function for nu.2:","log","\n") 
    cat("Formula: "); print(x$formula[[5]])      
    
    cat("\nEQUATION 6")
    cat("\nLink function for theta:",lind,"\n") 
    cat("Formula: "); print(x$formula[[6]]) 
       
       }
       
        
    
  
  
 }#else 
  
  
  
}# is.null
  
  
  
  
  
  cat("\n")
  
  
  
if(x$BivD == "T" && x$VC$margins[1] %in% c(x$VC$m2,x$VC$m3) && x$VC$margins[2] %in% c(x$VC$m2,x$VC$m3)){

                                                                
  if(x$margins[1] %in% cont2par && x$margins[2] %in% cont2par) cat( s1, format(s1.p, digits=3), "  ",s2, format(s2.p, digits=3),
                                                                     "\n",cp, format(as.p, digits=3),"  dof = ", format(x$dof.a, digits=3),
                                                                     "\nn = ",x$n,"  total edf = ",format(x$t.edf, digits=3),"\n\n", sep="")                                                                     
                                                                     

                                                                     
                                                                     
  if(x$margins[1] %in% cont3par && x$margins[2] %in% cont3par) cat( s1, format(s1.p, digits=3), "  ",s2, format(s2.p, digits=3),
                                                                     "\n",n1, format(n1.p, digits=3), "  ",n2, format(n2.p, digits=3),
                                                                     "\n",cp, format(as.p, digits=3),"  dof = ", format(x$dof.a, digits=3),
                                                                     "\nn = ",x$n,"  total edf = ",format(x$t.edf, digits=3),"\n\n", sep="")                                                                      
                                                                  

                                                                     
  if(x$margins[1] %in% cont2par && x$margins[2] %in% cont3par) cat( s1, format(s1.p, digits=3), "  ",s2, format(s2.p, digits=3),
                                                                     "\n",n2, format(n2.p, digits=3),
                                                                     "\n",cp, format(as.p, digits=3),"  dof = ", format(x$dof.a, digits=3),
                                                                     "\nn = ",x$n,"  total edf = ",format(x$t.edf, digits=3),"\n\n", sep="")                                                                      

 

  if(x$margins[1] %in% cont3par && x$margins[2] %in% cont2par) cat( s1, format(s1.p, digits=3), "  ",s2, format(s2.p, digits=3),
                                                                     "\n",n1, format(n1.p, digits=3), 
                                                                     "\n",cp, format(as.p, digits=3),"  dof = ", format(x$dof.a, digits=3),
                                                                     "\nn = ",x$n,"  total edf = ",format(x$t.edf, digits=3),"\n\n", sep="") 

}else{




  
            
  if(x$margins[1] %in% cont2par && x$margins[2] %in% cont2par) cat( s1, format(s1.p, digits=3), "  ",s2, format(s2.p, digits=3),
                                                                     "\n",cp, format(as.p, digits=3),
                                                                     "\nn = ",x$n,"  total edf = ",format(x$t.edf, digits=3),"\n\n", sep="")
                                                                     
                                                                 
  if(x$margins[1] %in% cont3par && x$margins[2] %in% cont3par) cat( s1, format(s1.p, digits=3), "  ",s2, format(s2.p, digits=3),
                                                                     "\n",n1, format(n1.p, digits=3), "  ",n2, format(n2.p, digits=3),
                                                                     "\n",cp, format(as.p, digits=3),
                                                                     "\nn = ",x$n,"  total edf = ",format(x$t.edf, digits=3),"\n\n", sep="")   
                                                                     
                                                                                                                              
  if(x$margins[1] %in% cont1par && x$margins[2] %in% cont2par ) cat( s2, format(s2.p, digits=3),
                                                                     "\n",cp, format(as.p, digits=3),
                                                                     "\nn = ",x$n,"  total edf = ",format(x$t.edf, digits=3),"\n\n", sep="")
                                                                     
  if(x$margins[1] %in% cont1par && x$margins[2] %in% cont1par ) cat( cp, format(as.p, digits=3),
                                                                     "\nn = ",x$n,"  total edf = ",format(x$t.edf, digits=3),"\n\n", sep="")                                                                     
                                                                     
  if(x$margins[1] %in% cont1par && x$margins[2] %in% cont3par ) cat( s2, format(s2.p, digits=3),
                                                                     "\n",n2, format(n2.p, digits=3),
                                                                     "\n",cp, format(as.p, digits=3),
                                                                     "\nn = ",x$n,"  total edf = ",format(x$t.edf, digits=3),"\n\n", sep="")  
                                                                     

  if(x$margins[1] %in% cont2par && x$margins[2] %in% cont3par) cat( s1, format(s1.p, digits=3), "  ",s2, format(s2.p, digits=3),
                                                                     "\n",n2, format(n2.p, digits=3),
                                                                     "\n",cp, format(as.p, digits=3),
                                                                     "\nn = ",x$n,"  total edf = ",format(x$t.edf, digits=3),"\n\n", sep="")   
                                                                     
                                                                     

  if(x$margins[1] %in% cont3par && x$margins[2] %in% cont2par) cat( s1, format(s1.p, digits=3), "  ",s2, format(s2.p, digits=3),
                                                                     "\n",n1, format(n1.p, digits=3), 
                                                                     "\n",cp, format(as.p, digits=3),
                                                                     "\nn = ",x$n,"  total edf = ",format(x$t.edf, digits=3),"\n\n", sep="")   

 

}


invisible(x)

}

