pred.var <- function(formula, l.flist, gaml = FALSE, triv = FALSE, informative = "no", ROY = FALSE){

ig <- interpret.gam(formula)
v3 <- v2 <- NULL
  
  
if(ROY == TRUE){  

    or1 <- as.character(formula[[1]][2])
    or2 <- as.character(formula[[2]][2])
    or3 <- as.character(formula[[3]][2])

    v1 <- all.vars(as.formula(formula[[1]]))[1]
    v1 <- c(v1, ig[[1]]$pred.names)
    v2 <- all.vars(as.formula(formula[[2]]))[1]
    v2 <- c(v2, ig[[2]]$pred.names)
    v3 <- all.vars(as.formula(formula[[3]]))[1]
    v3 <- c(v3, ig[[3]]$pred.names)    
    
    pred.n <- union(v1, c(v2, v3, or1, or2, or3))
                       
                    
  if( l.flist == 5 ){  
      
    v4 <- ig[[4]]$pred.names
    v5 <- ig[[5]]$pred.names 
    
    pred.n <- union(pred.n, c(v4,v5))
    
                    } 
                    
  if( l.flist == 6 ){  
      
    v4 <- ig[[4]]$pred.names
    v5 <- ig[[5]]$pred.names 
    v6 <- ig[[6]]$pred.names 
    
    
    pred.n <- union(pred.n, c(v4,v5,v6))
    
                    }
                    
  if( l.flist == 7 ){  
      
    v4 <- ig[[4]]$pred.names
    v5 <- ig[[5]]$pred.names 
    v6 <- ig[[6]]$pred.names 
    v7 <- ig[[7]]$pred.names 
    
    pred.n <- union(pred.n, c(v4,v5,v6,v7))
    
                    }
                    
  if( l.flist == 8 ){  
      
    v4 <- ig[[4]]$pred.names
    v5 <- ig[[5]]$pred.names 
    v6 <- ig[[6]]$pred.names 
    v7 <- ig[[7]]$pred.names 
    v8 <- ig[[8]]$pred.names 
    
    pred.n <- union(pred.n, c(v4,v5,v6,v7,v8))
    
                    }
                    
  if( l.flist == 9 ){  
      
    v4 <- ig[[4]]$pred.names
    v5 <- ig[[5]]$pred.names 
    v6 <- ig[[6]]$pred.names 
    v7 <- ig[[7]]$pred.names 
    v8 <- ig[[8]]$pred.names  
    v9 <- ig[[9]]$pred.names     
    
    
    pred.n <- union(pred.n, c(v4,v5,v6,v7,v8,v9))
    
                    }                    
                    
                      
} 
 
 
 
 
 
 
 
 if(ROY == FALSE){  # ROY

 
 
 
 
if(triv == TRUE){    

    v1 <- all.vars(as.formula(formula[[1]]))[1]
    v1 <- c(v1, ig[[1]]$pred.names)
    v2 <- all.vars(as.formula(formula[[2]]))[1]
    v2 <- c(v2, ig[[2]]$pred.names)
    v3 <- all.vars(as.formula(formula[[3]]))[1]
    v3 <- c(v3, ig[[3]]$pred.names)    
    
    pred.n <- union(v1, c(v2, v3))
                       
                    
  if( l.flist == 6 ){  

    v4 <- ig[[4]]$pred.names  
    v5 <- ig[[5]]$pred.names  
    v6 <- ig[[6]]$pred.names  
   
    pred.n <- union(v1, c(v2, v3, v4, v5, v6))
    
                    }                     
                      
}
    
    
if(gaml == FALSE && triv == FALSE){   # not very efficient but ok for now 

    or1 <- as.character(formula[[1]][2])
    or2 <- as.character(formula[[2]][2])

  if( l.flist == 2 ){  
    v1 <- all.vars(as.formula(formula[[1]]))[1]
    v1 <- c(v1, ig[[1]]$pred.names)
    v2 <- all.vars(as.formula(formula[[2]]))[1]
    v2 <- c(v2, ig[[2]]$pred.names)
    pred.n <- union(v1,c(v2,or1,or2))
                    }
  
  if( l.flist == 3 ){  
    v1 <- all.vars(as.formula(formula[[1]]))[1]
    v1 <- c(v1, ig[[1]]$pred.names)
    v2 <- all.vars(as.formula(formula[[2]]))[1]
    v2 <- c(v2, ig[[2]]$pred.names)
    v3 <- ig[[3]]$pred.names 
    
    pred.n <- union(v1,c(v2,v3,or1,or2))
                    } 
  
  if( l.flist == 4 ){  
    v1 <- all.vars(as.formula(formula[[1]]))[1]
    v1 <- c(v1, ig[[1]]$pred.names)
    v2 <- all.vars(as.formula(formula[[2]]))[1]
    v2 <- c(v2, ig[[2]]$pred.names)
    v3 <- ig[[3]]$pred.names  
    v4 <- ig[[4]]$pred.names  
    pred.n <- union(v1,c(v2,v3,v4,or1,or2))
                    } 
  
  if( l.flist == 5 ){  
    v1 <- all.vars(as.formula(formula[[1]]))[1]
    v1 <- c(v1, ig[[1]]$pred.names)
    v2 <- all.vars(as.formula(formula[[2]]))[1]
    v2 <- c(v2, ig[[2]]$pred.names)
    v3 <- ig[[3]]$pred.names 
    v4 <- ig[[4]]$pred.names
    v5 <- ig[[5]]$pred.names 
    pred.n <- union(v1,c(v2,v3,v4,v5,or1,or2))
                    }   
  
  if( l.flist == 6 ){  
    v1 <- all.vars(as.formula(formula[[1]]))[1]
    v1 <- c(v1, ig[[1]]$pred.names)
    v2 <- all.vars(as.formula(formula[[2]]))[1]
    v2 <- c(v2, ig[[2]]$pred.names)
    v3 <- ig[[3]]$pred.names
    v4 <- ig[[4]]$pred.names
    v5 <- ig[[5]]$pred.names  
    v6 <- ig[[6]]$pred.names
    pred.n <- union(v1,c(v2,v3,v4,v5,v6,or1,or2))
                    }  

  if( l.flist == 7 ){  
    v1 <- all.vars(as.formula(formula[[1]]))[1]
    v1 <- c(v1, ig[[1]]$pred.names)
    v2 <- all.vars(as.formula(formula[[2]]))[1]
    v2 <- c(v2, ig[[2]]$pred.names)
    v3 <- ig[[3]]$pred.names
    v4 <- ig[[4]]$pred.names
    v5 <- ig[[5]]$pred.names  
    v6 <- ig[[6]]$pred.names
    v7 <- ig[[7]]$pred.names
    pred.n <- union(v1,c(v2,v3,v4,v5,v6,v7,or1,or2))
                    }   
                    
  if( l.flist == 8 ){  
    v1 <- all.vars(as.formula(formula[[1]]))[1]
    v1 <- c(v1, ig[[1]]$pred.names)
    v2 <- all.vars(as.formula(formula[[2]]))[1]
    v2 <- c(v2, ig[[2]]$pred.names)
    v3 <- ig[[3]]$pred.names
    v4 <- ig[[4]]$pred.names
    v5 <- ig[[5]]$pred.names  
    v6 <- ig[[6]]$pred.names
    v7 <- ig[[7]]$pred.names
    v8 <- ig[[8]]$pred.names
    pred.n <- union(v1,c(v2,v3,v4,v5,v6,v7,v8,or1,or2))
                    }                      
                

}



if(gaml == TRUE && triv == FALSE){  # not very efficient but ok for now  


 or1 <- as.character(formula[[1]][2])


  if( l.flist == 1 ){  
    v1 <- all.vars(as.formula(formula[[1]]))[1]
    v1 <- c(v1, ig[[1]]$pred.names)
    pred.n <- union(v1, or1)
                    }

  if( l.flist == 2 ){  
    v1 <- all.vars(as.formula(formula[[1]]))[1]
    v1 <- c(v1, ig[[1]]$pred.names)
    v2 <- ig[[2]]$pred.names
    pred.n <- union(v1,c(v2, or1))
                    }
  
  if( l.flist == 3 ){  
    v1 <- all.vars(as.formula(formula[[1]]))[1]
    v1 <- c(v1, ig[[1]]$pred.names)
    v2 <- ig[[2]]$pred.names
    v3 <- ig[[3]]$pred.names 
    
    pred.n <- union(v1,c(v2, v3, or1))
                    } 

    if(informative == "yes") v2 <- all.vars(as.formula(formula[[2]]))[1]


}


} # ROY FALSE
                
  list(v1 = v1, v2 = v2, v3 = v3, pred.n = pred.n)      
      
}