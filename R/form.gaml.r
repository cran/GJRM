form.gaml <- function(formula, l.flist, M, type = "copR"){
  
formula.gamlss1 <- formula.gamlss2 <- formula.gamlss3 <- NULL



if(type == "biv"){


if(l.flist == 2) formula.gamlss2 <- list(formula[[2]])   # good for all cases with no distrib stuff
if(l.flist == 3) formula.gamlss2 <- list(formula[[2]])   # good for distrib stuff, PO
if(l.flist == 4) formula.gamlss2 <- list(formula[[2]],formula[[3]]) # distr stuff, N etc
if(l.flist == 5) formula.gamlss2 <- list(formula[[2]],formula[[3]],formula[[4]]) # distr stuff, DAGUM etc


}






if(type == "copR"){


if(l.flist == 2){

formula.gamlss1 <- list(formula[[1]])
formula.gamlss2 <- list(formula[[2]])

}



if(l.flist > 2){##


if(M$margins[1] %in% c(M$m2) && M$margins[2] %in% c(M$bl)){
    
formula.gamlss1 <- list(formula[[1]], formula[[3]])
formula.gamlss2 <- list(formula[[2]])
    
}

if(M$margins[1] %in% c(M$m3) && M$margins[2] %in% c(M$bl)){
    
formula.gamlss1 <- list(formula[[1]], formula[[3]], formula[[4]])
formula.gamlss2 <- list(formula[[2]])
    
}


if(M$margins[1] %in% c(M$m1d,M$bl) && M$margins[2] %in% c(M$m1d,M$bl)){
    
formula.gamlss1 <- list(formula[[1]])
formula.gamlss2 <- list(formula[[2]])
    
}   
  
if(M$margins[1] %in% c(M$m1d) && M$margins[2] %in% c(M$m2,M$m2d)){
    
formula.gamlss1 <- list(formula[[1]])
formula.gamlss2 <- list(formula[[2]],formula[[3]]) 
    
}  
  
  
if(M$margins[1] %in% c(M$m1d) && M$margins[2] %in% c(M$m3,M$m3d)){
      
formula.gamlss1 <- list(formula[[1]])
formula.gamlss2 <- list(formula[[2]],formula[[3]],formula[[4]])   
   
}    
    
    
if(M$margins[1] %in% c(M$m2,M$m2d) && M$margins[2] %in% c(M$m2,M$m2d)){
    
formula.gamlss1 <- list(formula[[1]],formula[[3]])
formula.gamlss2 <- list(formula[[2]],formula[[4]])      
   
}   
  
  
if(M$margins[1] %in% c(M$m3,M$m3d) && M$margins[2] %in% c(M$m3,M$m3d)){
    
formula.gamlss1 <- list(formula[[1]],formula[[3]],formula[[5]])
formula.gamlss2 <- list(formula[[2]],formula[[4]],formula[[6]])     

}    
  
      
  
if(M$margins[1] %in% c(M$m2,M$m2d) && M$margins[2] %in% c(M$m3,M$m3d)){
    
formula.gamlss1 <- list(formula[[1]],formula[[3]])
formula.gamlss2 <- list(formula[[2]],formula[[4]],formula[[5]])      
}    
    
if(M$margins[1] %in% c(M$m3,M$m3d) && M$margins[2] %in% c(M$m2,M$m2d)){
    
formula.gamlss1 <- list(formula[[1]],formula[[3]],formula[[5]])
formula.gamlss2 <- list(formula[[2]],formula[[4]])   
}    
      
  
} ##



}









if(type == "ROY"){


if(l.flist == 3){

     formula.gamlss2 <- list(formula[[2]])
     formula.gamlss3 <- list(formula[[3]])

}



if(l.flist > 3){##


      if(M$margins[2] %in% c(M$bl, M$m1d) && M$margins[2] %in% c(M$m1d, M$bl)){
    
        formula.gamlss2 <- list(formula[[2]])
        formula.gamlss3 <- list(formula[[3]])
    
                                                                  }
                                                                  
      if(M$margins[2] %in% c(M$m2d,M$m2) && M$margins[2] %in% c(M$m2d,M$m2)){
    
        formula.gamlss2 <- list(formula[[2]],formula[[4]])
        formula.gamlss3 <- list(formula[[3]],formula[[5]])
    
                                                                            }       
                                                                            
      if(M$margins[2] %in% c(M$bl, M$m1d) && M$margins[2] %in% c(M$m2d)){
    
        formula.gamlss2 <- list(formula[[2]])
        formula.gamlss3 <- list(formula[[3]],formula[[4]])
    
                                                                   }    
                                                                   
                                                                            
      if(M$margins[2] %in% c(M$m2d) && M$margins[2] %in% c(M$bl, M$m1d)){
    
        formula.gamlss2 <- list(formula[[2]],formula[[4]])
        formula.gamlss3 <- list(formula[[3]])
    
                                                                   }                                                                     
 
      if(M$margins[2] %in% c(M$m3) && M$margins[2] %in% c(M$m3)){
    
        formula.gamlss2 <- list(formula[[2]],formula[[4]],formula[[6]])
        formula.gamlss3 <- list(formula[[3]],formula[[5]],formula[[7]])
    
                                                                }    
 
      if(M$margins[2] %in% c(M$m2) && M$margins[2] %in% c(M$m3)){
    
        formula.gamlss2 <- list(formula[[2]],formula[[4]])
        formula.gamlss3 <- list(formula[[3]],formula[[5]],formula[[6]])
    
                                                                }   
  
      if(M$margins[2] %in% c(M$m3) && M$margins[2] %in% c(M$m2)){
    
        formula.gamlss2 <- list(formula[[2]],formula[[4]],formula[[6]])
        formula.gamlss3 <- list(formula[[3]],formula[[5]])
    
                                                                }   
  
  
} ##



}



list(formula.gamlss1 = formula.gamlss1, formula.gamlss2 = formula.gamlss2, formula.gamlss3 = formula.gamlss3)



}

