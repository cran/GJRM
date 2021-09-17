form.check <- function(formula, l.flist, gamlss = FALSE, ROY = FALSE){
 
    if(l.flist > 1 && gamlss == TRUE){
    
    f2t <- try(formula[[2]][[3]], silent = TRUE)  
    if(class(f2t)!="try-error") stop("The second equation does not require a response.") 
 
    } 
 
 
 
    if(l.flist > 2){
    
    
    if(ROY == FALSE){
    
    f3t <- try(formula[[3]][[3]], silent = TRUE)  
    if(class(f3t)!="try-error") stop("The third equation does not require a response.")
    
    }
    
    	if(l.flist > 3){
    	
    f4t <- try(formula[[4]][[3]], silent = TRUE)  
    if(class(f4t)!="try-error") stop("The fourth equation does not require a response.")  
    
    		if(l.flist > 4){
    		
    f5t <- try(formula[[5]][[3]], silent = TRUE)  
    if(class(f5t)!="try-error") stop("The fifth equation does not require a response.")    
    
    			if(l.flist > 5){
    			
    f6t <- try(formula[[6]][[3]], silent = TRUE)  
    if(class(f6t)!="try-error") stop("The sixth equation does not require a response.")      			
    			
    
        		      if(l.flist > 6){
        			
    f7t <- try(formula[[7]][[3]], silent = TRUE)  
    if(class(f7t)!="try-error") stop("The seventh equation does not require a response.")   
    			
    			
          		      if(l.flist > 7){
        			
    f8t <- try(formula[[8]][[3]], silent = TRUE)  
    if(class(f8t)!="try-error") stop("The eighth equation does not require a response.")   
    			
    		  
    		
              		      if(l.flist > 8){
        			
    f9t <- try(formula[[9]][[3]], silent = TRUE)  
    if(class(f9t)!="try-error") stop("The ninth equation does not require a response.")   
    			
    		}  		
    			
    			
    			
    			                               }
    			
    			
    			                        }    
    					}
    				}    				
                            }

}


}

