pscr0 <- function(x, type = "copR"){

if(type == "copR"){

  if(x$margins[1]%in%c("N"))              cat("\nMARGIN 1: Gaussian")  
  if(x$margins[1]%in%c("tN"))             cat("\nMARGIN 1: Truncated Gaussian")  
  
  if(x$margins[1]=="GU")                  cat("\nMARGIN 1: Gumbel")    
  if(x$margins[1]=="GP")                  cat("\nMARGIN 1: Generalised Pareto")   
  if(x$margins[1]=="GPII")                cat("\nMARGIN 1: Generalised Pareto II")    
  if(x$margins[1]=="GPo")                 cat("\nMARGIN 1: Generalised Pareto (Orthogonal Parameterisation)")    
  if(x$margins[1]=="DGP")                 cat("\nMARGIN 1: Discrete Generalised Pareto")    
  if(x$margins[1]=="DGPII")               cat("\nMARGIN 1: Discrete Generalised Pareto II")  
  if(x$margins[1]=="DGP0")                cat("\nMARGIN 1: Discrete Generalised Pareto (shape = 0)")    

  
  
  if(x$margins[1]=="rGU")                 cat("\nMARGIN 1: Reverse Gumbel")  
  if(x$margins[1]=="LO")                  cat("\nMARGIN 1: Logistic")   
  if(x$margins[1]=="LN")                  cat("\nMARGIN 1: Log-normal") 
  if(x$margins[1]=="WEI")                 cat("\nMARGIN 1: Weibull") 
  #if(x$margins[1]=="GO")                  cat("\nMARGIN 1: Gompertz")   
  if(x$margins[1]=="IG")                  cat("\nMARGIN 1: Inverse Gaussian") 
  if(x$margins[1]%in%c("GA","GAi"))       cat("\nMARGIN 1: Gamma")  
  #if(x$margins[1]%in%c("GGA"))            cat("\nMARGIN 1: Generalised Gamma")   
  #if(x$margins[1]%in%c("GA2"))            cat("\nMARGIN 1: Gamma (parametrisation 2)")   
  if(x$margins[1]=="BE")                  cat("\nMARGIN 1: Beta")    
  if(x$margins[1]=="DAGUM")               cat("\nMARGIN 1: Dagum")  
  if(x$margins[1]=="TW")                  cat("\nMARGIN 1: Tweedie")  
  
  if(x$margins[1]=="SM")                  cat("\nMARGIN 1: Singh-Maddala") 
  if(x$margins[1]=="FISK")                cat("\nMARGIN 1: Fisk") 
  if(x$margins[1]%in%c("NBI"))            cat("\nMARGIN 1: Negative Binomial - Type I") 
  if(x$margins[1]%in%c("NBII"))           cat("\nMARGIN 1: Negative Binomial - Type II")  
  if(x$margins[1]=="PIG")                 cat("\nMARGIN 1: Poisson Inverse Gaussian")
  if(x$margins[1]%in%c("tNBI"))          cat("\nMARGIN 1: Truncated Negative Binomial - Type I") 
  if(x$margins[1]%in%c("tNBII"))         cat("\nMARGIN 1: Truncated Negative Binomial - Type II")  
  if(x$margins[1]=="tPIG")               cat("\nMARGIN 1: Truncated Poisson Inverse Gaussian")  
  if(x$margins[1]=="P")                   cat("\nMARGIN 1: Poisson")   
  if(x$margins[1]=="tP")                 cat("\nMARGIN 1: Truncated Poisson")    

if(x$surv.flex == FALSE){ 

  if(x$margins[1]=="probit")              cat("\nMARGIN 1: probit")    
  if(x$margins[1]=="logit")               cat("\nMARGIN 1: logit")    
  if(x$margins[1]=="cloglog")             cat("\nMARGIN 1: cloglog")    

}
  
  
  if(x$margins[2]%in%c("N"))              cat("\nMARGIN 2: Gaussian")  
  if(x$margins[2]%in%c("tN"))             cat("\nMARGIN 2: Truncated Gaussian")  
  
  if(x$margins[2]=="GP")                  cat("\nMARGIN 2: Generalised Pareto")  
  if(x$margins[2]=="GPII")                cat("\nMARGIN 2: Generalised Pareto II")  
  if(x$margins[2]=="GPo")                 cat("\nMARGIN 2: Generalised Pareto (Orthogonal Parameterisation)")    
  

  if(x$margins[2]=="DGP")                 cat("\nMARGIN 2: Discrete Generalised Pareto")    
  if(x$margins[2]=="DGP0")                cat("\nMARGIN 2: Discrete Generalised Pareto (shape = 0)")    
  
  if(x$margins[2]=="DGPII")               cat("\nMARGIN 2: Discrete Generalised Pareto II")    
  if(x$margins[2]=="GU")     		  cat("\nMARGIN 2: Gumbel")    
  if(x$margins[2]=="rGU")    		  cat("\nMARGIN 2: Reverse Gumbel")  
  if(x$margins[2]=="LO")     		  cat("\nMARGIN 2: Logistic")   
  if(x$margins[2]=="LN")     		  cat("\nMARGIN 2: Log-normal") 
  if(x$margins[2]=="WEI")    		  cat("\nMARGIN 2: Weibull") 
  if(x$margins[2]=="GO")                  cat("\nMARGIN 2: Gompertz")   
  if(x$margins[2]=="IG")     		  cat("\nMARGIN 2: Inverse Gaussian") 
  if(x$margins[2]%in%c("GA","GAi")) 	  cat("\nMARGIN 2: Gamma")   
  #if(x$margins[2]%in%c("GGA"))            cat("\nMARGIN 2: Generalised Gamma")   
  #if(x$margins[2]%in%c("GA2")) 	          cat("\nMARGIN 2: Gamma (parametrisation 2)")   
  if(x$margins[2]=="BE")     		  cat("\nMARGIN 2: Beta")    
  if(x$margins[2]=="DAGUM")  		  cat("\nMARGIN 2: Dagum")
  if(x$margins[2]=="TW")  		  cat("\nMARGIN 2: Tweedie")
  
  
  if(x$margins[2]=="SM")     		  cat("\nMARGIN 2: Singh-Maddala") 
  if(x$margins[2]=="FISK")   		  cat("\nMARGIN 2: Fisk") 
  if(x$margins[2]%in%c("NBI"))            cat("\nMARGIN 2: Negative Binomial - Type I") 
  if(x$margins[2]%in%c("NBII"))           cat("\nMARGIN 2: Negative Binomial - Type II")  
  if(x$margins[2]=="PIG")                 cat("\nMARGIN 2: Poisson Inverse Gaussian") 
  if(x$margins[2]%in%c("tNBI"))          cat("\nMARGIN 2: Truncated Negative Binomial - Type I") 
  if(x$margins[2]%in%c("tNBII"))         cat("\nMARGIN 2: Truncated Negative Binomial - Type II")  
  if(x$margins[2]=="tPIG")               cat("\nMARGIN 2: Truncated Poisson Inverse Gaussian")    
  if(x$margins[2]=="P")                   cat("\nMARGIN 2: Poisson")   
  if(x$margins[2]=="tP")                 cat("\nMARGIN 2: Truncated Poisson")    

if(x$surv.flex == FALSE){ 
 
  if(x$margins[2]=="probit")              cat("\nMARGIN 2: probit")    
  if(x$margins[2]=="logit")               cat("\nMARGIN 2: logit")    
  if(x$margins[2]=="cloglog")             cat("\nMARGIN 2: cloglog")  
  
}  
  
if(x$surv.flex == TRUE){


 if(x$end.surv == FALSE){  
  if(x$margins[1]=="probit")              cat("\nMARGIN 1: Survival with -probit link")    
  if(x$margins[1]=="logit")               cat("\nMARGIN 1: Survival with -logit link")    
  if(x$margins[1]=="cloglog")             cat("\nMARGIN 1: Survival with -cloglog link")  
                        }
 
 if(x$end.surv == TRUE){
  if(x$margins[1]=="probit")              cat("\nMARGIN 1: Bernoulli with probit link")    
  if(x$margins[1]=="logit")               cat("\nMARGIN 1: Bernoulli with logit link")    
  if(x$margins[1]=="cloglog")             cat("\nMARGIN 1: Bernoulli with cloglog link")  
                       } 

  if(x$margins[2]=="probit")              cat("\nMARGIN 2: Survival with -probit link") 
  if(x$margins[2]=="logit")               cat("\nMARGIN 2: Survival with -logit link")
  if(x$margins[2]=="cloglog")             cat("\nMARGIN 2: Survival with -cloglog link")

}
  
    
}


if(type == "gamlss"){

  if(x$margins[1]=="GP")               cat("\nDistribution: Generalised Pareto")   
  if(x$margins[1]=="GPII")             cat("\nDistribution: Generalised Pareto II")  
  if(x$margins[1]=="GPo")              cat("\nDistribution: Generalised Pareto (Orthogonal Parameterisation)")    

  if(x$margins[1]=="DGP")              cat("\nDistribution: Discrete Generalised Pareto")  
  if(x$margins[1]=="DGPII")            cat("\nDistribution: Discrete Generalised Pareto II")  
  if(x$margins[1]=="DGP0")             cat("\nDistribution: Discrete Generalised Pareto (shape = 0)")  
  
  
  if(x$margins[1]%in%c("N","N2"))      cat("\nDistribution: Gaussian") 
  if(x$margins[1]%in%c("tN"))          cat("\nDistribution: Truncated Gaussian")  
  if(x$margins[1]=="GU")               cat("\nDistribution: Gumbel")    
  if(x$margins[1]=="rGU")              cat("\nDistribution: Reverse Gumbel")  
  if(x$margins[1]=="LO")               cat("\nDistribution: Logistic")   
  if(x$margins[1]=="LN")               cat("\nDistribution: Log-normal") 
  if(x$margins[1]=="WEI")              cat("\nDistribution: Weibull") 
  if(x$margins[1]=="GO")               cat("\nDistribution: Gompertz")   
  if(x$margins[1]=="IG")               cat("\nDistribution: Inverse Gaussian") 
  if(x$margins[1]%in%c("GA","GAi"))    cat("\nDistribution: Gamma")  
  if(x$margins[1]%in%c("GGA"))         cat("\nDistribution: Generalised Gamma")  
  if(x$margins[1]%in%c("GA2"))         cat("\nDistribution: Gamma (parametrisation 2)")    
  if(x$margins[1]=="FISK")             cat("\nDistribution: Fisk") 
  if(x$margins[1]=="BE")               cat("\nDistribution: Beta")    
  if(x$margins[1]=="DAGUM")            cat("\nDistribution: Dagum")  
  if(x$margins[1]=="TW")               cat("\nDistribution: Tweedie")  
  
  
  
  if(x$margins[1]=="SM")               cat("\nDistribution: Singh-Maddala") 
  if(x$margins[1]%in%c("NBI"))         cat("\nDistribution: Negative Binomial - Type I") 
  if(x$margins[1]%in%c("NBII"))        cat("\nDistribution: Negative Binomial - Type II")
  if(x$margins[1]=="PIG")              cat("\nDistribution: Poisson inverse Gaussian")
  if(x$margins[1]%in%c("tNBI"))       cat("\nDistribution: Truncated Negative Binomial - Type I") 
  if(x$margins[1]%in%c("tNBII"))      cat("\nDistribution: Truncated Negative Binomial - Type II")
  if(x$margins[1]=="tPIG")            cat("\nDistribution: Truncated Poisson Inverse Gaussian")   
  if(x$margins[1]=="P")                cat("\nDistribution: Poisson")   
  if(x$margins[1]=="tP")              cat("\nDistribution: Truncated Poisson")    
  if(x$margins[1]=="GEVlink")          cat("\nDistribution: Bernoulli") 
  if(x$margins[1]=="probit")           cat("\nFlexible survival model with -probit link")    
  if(x$margins[1]=="logit")            cat("\nFlexible survival model with -logit link")    
  if(x$margins[1]=="cloglog")          cat("\nFlexible survival model with loglog link")    

  

}



if(type == "copSS"){

if(x$margins[1] %in% x$bl && is.null(x$K1)) cat("\nMARGIN 1: Bernoulli") 
if(x$margins[1] %in% x$bl && !is.null(x$K1)) cat("\nMARGIN 1: Ordinal") 


if(x$margins[2] %in% x$bl && is.null(x$K2)) cat("\nMARGIN 2: Bernoulli") 
if(x$margins[2] %in% x$bl && !is.null(x$K2)) cat("\nMARGIN 2: Ordinal") 


if(is.null(x$K2)){

if(x$margins[2]%in%c("N"))           cat("\nMARGIN 2: Gaussian")  
if(x$margins[2]%in%c("tN"))          cat("\nMARGIN 2: Truncated Gaussian")  

if(x$margins[2]=="GP")               cat("\nMARGIN 2: Generalised Pareto")  
if(x$margins[2]=="GPII")             cat("\nMARGIN 2: Generalised Pareto II")  
if(x$margins[2]=="GPo")              cat("\nMARGIN 2: Generalised Pareto (Orthogonal Parameterisation)")    

if(x$margins[2]=="DGP")              cat("\nMARGIN 2: Discrete Generalised Pareto")  
if(x$margins[2]=="DGPII")            cat("\nMARGIN 2: Discrete Generalised Pareto II")  
if(x$margins[2]=="DGP0")            cat("\nMARGIN 2: Discrete Generalised Pareto (shape = 0)")  



if(x$margins[2]=="GU")                 cat("\nMARGIN 2: Gumbel")    
if(x$margins[2]=="rGU")                cat("\nMARGIN 2: Reverse Gumbel")  
if(x$margins[2]=="LO")                 cat("\nMARGIN 2: Logistic")   
if(x$margins[2]=="LN")                 cat("\nMARGIN 2: Log-normal") 
if(x$margins[2]=="WEI")                cat("\nMARGIN 2: Weibull") 
if(x$margins[2]=="GO")                 cat("\nMARGIN 2: Gompertz") 
if(x$margins[2]=="IG")                 cat("\nMARGIN 2: Inverse Gaussian") 
if(x$margins[2]%in%c("GA","GAi"))      cat("\nMARGIN 2: Gamma")
#if(x$margins[2]%in%c("GGA"))          cat("\nMARGIN 2: Generalised Gamma")  
#if(x$margins[2]%in%c("GA2"))          cat("\nMARGIN 2: Gamma (parametrisation 2)")   
if(x$margins[2]=="BE")                 cat("\nMARGIN 2: Beta")    
if(x$margins[2]=="DAGUM")              cat("\nMARGIN 2: Dagum")
if(x$margins[2]=="TW")                 cat("\nMARGIN 2: Tweedie")


if(x$margins[2]=="SM")                 cat("\nMARGIN 2: Singh-Maddala") 
if(x$margins[2]=="FISK")               cat("\nMARGIN 2: Fisk") 
if(x$margins[2]=="FISK2")              cat("\nMARGIN 2: Fisk2") 
if(x$margins[2] %in% c("NBI"))         cat("\nMARGIN 2: Negative Binomial - Type I") 
if(x$margins[2]%in% c("NBII"))         cat("\nMARGIN 2: Negative Binomial - Type II")
if(x$margins[2]=="PIG")    	       cat("\nMARGIN 2: Poisson inverse Gaussian") 
if(x$margins[2] %in% c("tNBI"))        cat("\nMARGIN 2: Truncated Negative Binomial - Type I") 
if(x$margins[2]%in% c("tNBII"))        cat("\nMARGIN 2: Truncated Negative Binomial - Type II")
if(x$margins[2]=="tPIG")    	       cat("\nMARGIN 2: Truncated Poisson Inverse Gaussian") 
if(x$margins[2]=="P")     	       cat("\nMARGIN 2: Poisson")   
if(x$margins[2]=="tP")    	       cat("\nMARGIN 2: Truncated Poisson")  

}


}




if(type == "ROY"){

if(x$margins[1] %in% x$bl &&  is.null(x$K1)) cat("\nMARGIN 1: Switching Mechanism - Bernoulli") 
if(x$margins[1] %in% x$bl && !is.null(x$K1)) cat("\nMARGIN 1: Switching Mechanism - Ordinal") 
if(x$margins[2] %in% x$bl && x$surv == FALSE) cat("\nMARGIN 2: Regime 0 - Bernoulli") 
if(x$margins[3] %in% x$bl && x$surv == FALSE) cat("\nMARGIN 3: Regime 1 - Bernoulli") 


if(x$surv == TRUE){
 
 if(x$margins[2] == "probit")  cat("\nMARGIN 2: Regime 0 - Survival with -probit link") 
 if(x$margins[2] == "logit")   cat("\nMARGIN 2: Regime 0 - Survival with -logit link") 
 if(x$margins[2] == "cloglog") cat("\nMARGIN 2: Regime 0 - Survival with -cloglog link") 
 
 if(x$margins[3] == "probit")  cat("\nMARGIN 3: Regime 1 - Survival with -probit link") 
 if(x$margins[3] == "logit")   cat("\nMARGIN 3: Regime 1 - Survival with -logit link") 
 if(x$margins[3] == "cloglog") cat("\nMARGIN 3: Regime 1 - Survival with -cloglog link") 
 
}


if(x$margins[2] %in% c("N"))            cat("\nMARGIN 2: Regime 0 - Gaussian")  
if(x$margins[2] %in% c("tN"))           cat("\nMARGIN 2: Regime 0 - Truncated Gaussian")  
if(x$margins[2] == "GP")                cat("\nMARGIN 2: Regime 0 - Generalised Pareto")  
if(x$margins[2] == "GPII")              cat("\nMARGIN 2: Regime 0 - Generalised Pareto II")  
if(x$margins[2] == "GPo")               cat("\nMARGIN 2: Regime 0 - Generalised Pareto (Orthogonal Parameterisation)")    
if(x$margins[2] == "DGP")               cat("\nMARGIN 2: Regime 0 - Discrete Generalised Pareto")  
if(x$margins[2] == "DGPII")             cat("\nMARGIN 2: Regime 0 - Discrete Generalised Pareto II")  
if(x$margins[2] == "DGP0")              cat("\nMARGIN 2: Regime 0 - Discrete Generalised Pareto (shape = 0)")  
if(x$margins[2] == "GU")                cat("\nMARGIN 2: Regime 0 - Gumbel")    
if(x$margins[2] == "rGU")               cat("\nMARGIN 2: Regime 0 - Reverse Gumbel")  
if(x$margins[2] == "LO")                cat("\nMARGIN 2: Regime 0 - Logistic")   
if(x$margins[2] == "LN")                cat("\nMARGIN 2: Regime 0 - Log-normal") 
if(x$margins[2] == "WEI")               cat("\nMARGIN 2: Regime 0 - Weibull") 
if(x$margins[2] == "GO")                cat("\nMARGIN 2: Regime 0 - Gompertz") 
if(x$margins[2] == "IG")                cat("\nMARGIN 2: Regime 0 - Inverse Gaussian") 
if(x$margins[2] %in% c("GA","GAi"))     cat("\nMARGIN 2: Regime 0 - Gamma")  
if(x$margins[2] == "BE")                cat("\nMARGIN 2: Regime 0 - Beta")    
if(x$margins[2] == "DAGUM")             cat("\nMARGIN 2: Regime 0 - Dagum")
if(x$margins[2] == "TW")                cat("\nMARGIN 2: Regime 0 - Tweedie")
if(x$margins[2] == "SM")                cat("\nMARGIN 2: Regime 0 - Singh-Maddala") 
if(x$margins[2] == "FISK")              cat("\nMARGIN 2: Regime 0 - Fisk") 
if(x$margins[2] == "FISK2")             cat("\nMARGIN 2: Regime 0 - Fisk2") 
if(x$margins[2] %in% c("NBI"))          cat("\nMARGIN 2: Regime 0 - Negative Binomial - Type I") 
if(x$margins[2] %in% c("NBII"))         cat("\nMARGIN 2: Regime 0 - Negative Binomial - Type II")
if(x$margins[2] == "PIG")    	        cat("\nMARGIN 2: Regime 0 - Poisson Inverse Gaussian") 
if(x$margins[2] %in% c("tNBI"))        cat("\nMARGIN 2: Regime 0 - Truncated Negative Binomial - Type I") 
if(x$margins[2] %in% c("tNBII"))       cat("\nMARGIN 2: Regime 0 - Truncated Negative Binomial - Type II")
if(x$margins[2] == "tPIG")    	        cat("\nMARGIN 2: Regime 0 - Truncated Poisson Inverse Gaussian") 
if(x$margins[2] == "P")     	        cat("\nMARGIN 2: Regime 0 - Poisson")   
if(x$margins[2] == "tP")    	        cat("\nMARGIN 2: Regime 0 - Truncated Poisson")  


if(x$margins[3] %in% c("N"))            cat("\nMARGIN 3: Regime 1 - Gaussian")  
if(x$margins[3] %in% c("tN"))           cat("\nMARGIN 3: Regime 1 - Truncated Gaussian")  
if(x$margins[3] == "GP")                cat("\nMARGIN 3: Regime 1 - Generalised Pareto")  
if(x$margins[3] == "GPII")              cat("\nMARGIN 3: Regime 1 - Generalised Pareto II")  
if(x$margins[3] == "GPo")               cat("\nMARGIN 3: Regime 1 - Generalised Pareto (Orthogonal Parameterisation)")    
if(x$margins[3] == "DGP")               cat("\nMARGIN 3: Regime 1 - Discrete Generalised Pareto")  
if(x$margins[3] == "DGPII")             cat("\nMARGIN 3: Regime 1 - Discrete Generalised Pareto II")  
if(x$margins[3] == "DGP0")              cat("\nMARGIN 3: Regime 1 - Discrete Generalised Pareto (shape = 0)")  
if(x$margins[3] == "GU")                cat("\nMARGIN 3: Regime 1 - Gumbel")    
if(x$margins[3] == "rGU")               cat("\nMARGIN 3: Regime 1 - Reverse Gumbel")  
if(x$margins[3] == "LO")                cat("\nMARGIN 3: Regime 1 - Logistic")   
if(x$margins[3] == "LN")                cat("\nMARGIN 3: Regime 1 - Log-normal") 
if(x$margins[3] == "WEI")               cat("\nMARGIN 3: Regime 1 - Weibull") 
if(x$margins[3] == "GO")                cat("\nMARGIN 3: Regime 1 - Gompertz") 
if(x$margins[3] == "IG")                cat("\nMARGIN 3: Regime 1 - Inverse Gaussian") 
if(x$margins[3] %in% c("GA","GAi"))     cat("\nMARGIN 3: Regime 1 - Gamma")  
if(x$margins[3] == "BE")                cat("\nMARGIN 3: Regime 1 - Beta")    
if(x$margins[3] == "DAGUM")             cat("\nMARGIN 3: Regime 1 - Dagum")
if(x$margins[3] == "TW")                cat("\nMARGIN 3: Regime 1 - Tweedie")
if(x$margins[3] == "SM")                cat("\nMARGIN 3: Regime 1 - Singh-Maddala") 
if(x$margins[3] == "FISK")              cat("\nMARGIN 3: Regime 1 - Fisk") 
if(x$margins[3] == "FISK2")             cat("\nMARGIN 3: Regime 1 - Fisk2") 
if(x$margins[3] %in% c("NBI"))          cat("\nMARGIN 3: Regime 1 - Negative Binomial - Type I") 
if(x$margins[3] %in% c("NBII"))         cat("\nMARGIN 3: Regime 1 - Negative Binomial - Type II")
if(x$margins[3] == "PIG")    	        cat("\nMARGIN 3: Regime 1 - Poisson Inverse Gaussian") 
if(x$margins[3] %in% c("tNBI"))        cat("\nMARGIN 3: Regime 1 - Truncated Negative Binomial - Type I") 
if(x$margins[3] %in% c("tNBII"))       cat("\nMARGIN 3: Regime 1 - Truncated Negative Binomial - Type II")
if(x$margins[3] == "tPIG")    	        cat("\nMARGIN 3: Regime 1 - Truncated Poisson Inverse Gaussian") 
if(x$margins[3] == "P")     	        cat("\nMARGIN 3: Regime 1 - Poisson")   
if(x$margins[3] == "tP")    	        cat("\nMARGIN 3: Regime 1 - Truncated Poisson")  




}











         
}

