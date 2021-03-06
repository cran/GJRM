pscr0 <- function(x, type = "copR"){

if(type == "copR"){

  if(x$margins[1]%in%c("N"))              cat("\nMARGIN 1: Gaussian")  
  if(x$margins[1]=="GU")                  cat("\nMARGIN 1: Gumbel")    
  if(x$margins[1]=="GP")                  cat("\nMARGIN 1: generalised Pareto")   
  if(x$margins[1]=="GPII")                cat("\nMARGIN 1: generalised Pareto II")    
  if(x$margins[1]=="GPo")                 cat("\nMARGIN 1: generalised Pareto (Orthogonal Parameterisation)")    
  if(x$margins[1]=="DGP")                 cat("\nMARGIN 1: discrete generalised Pareto")    
  if(x$margins[1]=="DGPII")               cat("\nMARGIN 1: discrete generalised Pareto II")    
  
  
  if(x$margins[1]=="rGU")                 cat("\nMARGIN 1: reverse Gumbel")  
  if(x$margins[1]=="LO")                  cat("\nMARGIN 1: logistic")   
  if(x$margins[1]=="LN")                  cat("\nMARGIN 1: log-normal") 
  if(x$margins[1]=="WEI")                 cat("\nMARGIN 1: Weibull") 
  #if(x$margins[1]=="GO")                  cat("\nMARGIN 1: Gompertz")   
  if(x$margins[1]=="iG")                  cat("\nMARGIN 1: inverse Gaussian") 
  if(x$margins[1]%in%c("GA","GAi"))       cat("\nMARGIN 1: gamma")  
  #if(x$margins[1]%in%c("GGA"))            cat("\nMARGIN 1: generalised gamma")   
  #if(x$margins[1]%in%c("GA2"))            cat("\nMARGIN 1: gamma (parametrisation 2)")   
  if(x$margins[1]=="BE")                  cat("\nMARGIN 1: beta")    
  if(x$margins[1]=="DAGUM")               cat("\nMARGIN 1: Dagum")  
  if(x$margins[1]=="TW")                  cat("\nMARGIN 1: Tweedie")  
  
  if(x$margins[1]=="SM")                  cat("\nMARGIN 1: Singh-Maddala") 
  if(x$margins[1]=="FISK")                cat("\nMARGIN 1: Fisk") 
  if(x$margins[1]%in%c("NBI","NBIa"))     cat("\nMARGIN 1: Negative Binomial - Type I") 
  if(x$margins[1]%in%c("NBII","NBIIa"))   cat("\nMARGIN 1: Negative Binomial - Type II")  
  if(x$margins[1]=="PIG")                 cat("\nMARGIN 1: Poisson inverse Gaussian")
  if(x$margins[1]=="PO")                  cat("\nMARGIN 1: Poisson")   
  if(x$margins[1]=="ZTP")                 cat("\nMARGIN 1: Zero Truncated Poisson")    

if(x$surv.flex == FALSE){ 

  if(x$margins[1]=="probit")              cat("\nMARGIN 1: probit")    
  if(x$margins[1]=="logit")               cat("\nMARGIN 1: logit")    
  if(x$margins[1]=="cloglog")             cat("\nMARGIN 1: cloglog")    

}
  
  
  if(x$margins[2]%in%c("N"))              cat("\nMARGIN 2: Gaussian")  
  if(x$margins[2]=="GP")                  cat("\nMARGIN 2: generalised Pareto")  
  if(x$margins[2]=="GPII")                cat("\nMARGIN 2: generalised Pareto II")  
  if(x$margins[2]=="GPo")                 cat("\nMARGIN 2: generalised Pareto (Orthogonal Parameterisation)")    
  

  if(x$margins[2]=="DGP")                 cat("\nMARGIN 2: discrete generalised Pareto")    
  if(x$margins[2]=="DGPII")               cat("\nMARGIN 2: discrete generalised Pareto II")    
  if(x$margins[2]=="GU")     		  cat("\nMARGIN 2: Gumbel")    
  if(x$margins[2]=="rGU")    		  cat("\nMARGIN 2: reverse Gumbel")  
  if(x$margins[2]=="LO")     		  cat("\nMARGIN 2: logistic")   
  if(x$margins[2]=="LN")     		  cat("\nMARGIN 2: log-normal") 
  if(x$margins[2]=="WEI")    		  cat("\nMARGIN 2: Weibull") 
  if(x$margins[2]=="GO")                  cat("\nMARGIN 2: Gompertz")   
  if(x$margins[2]=="iG")     		  cat("\nMARGIN 2: inverse Gaussian") 
  if(x$margins[2]%in%c("GA","GAi")) 	  cat("\nMARGIN 2: gamma")   
  #if(x$margins[2]%in%c("GGA"))            cat("\nMARGIN 2: generalised gamma")   
  #if(x$margins[2]%in%c("GA2")) 	          cat("\nMARGIN 2: gamma (parametrisation 2)")   
  if(x$margins[2]=="BE")     		  cat("\nMARGIN 2: beta")    
  if(x$margins[2]=="DAGUM")  		  cat("\nMARGIN 2: Dagum")
  if(x$margins[2]=="TW")  		  cat("\nMARGIN 2: Tweedie")
  
  
  if(x$margins[2]=="SM")     		  cat("\nMARGIN 2: Singh-Maddala") 
  if(x$margins[2]=="FISK")   		  cat("\nMARGIN 2: Fisk") 
  if(x$margins[2]%in%c("NBI","NBIa"))     cat("\nMARGIN 2: Negative Binomial - Type I") 
  if(x$margins[2]%in%c("NBII","NBIIa"))   cat("\nMARGIN 2: Negative Binomial - Type II")  
  if(x$margins[2]=="PIG")                 cat("\nMARGIN 2: Poisson inverse Gaussian")  
  if(x$margins[2]=="PO")                  cat("\nMARGIN 2: Poisson")   
  if(x$margins[2]=="ZTP")                 cat("\nMARGIN 2: Zero Truncated Poisson")    

if(x$surv.flex == FALSE){ 
 
  if(x$margins[2]=="probit")              cat("\nMARGIN 2: probit")    
  if(x$margins[2]=="logit")               cat("\nMARGIN 2: logit")    
  if(x$margins[2]=="cloglog")             cat("\nMARGIN 2: cloglog")  
  
}  
  
if(x$surv.flex == TRUE){

  if(x$margins[1]=="probit")              cat("\nMARGIN 1: survival with -probit link")    
  if(x$margins[1]=="logit")               cat("\nMARGIN 1: survival with -logit link")    
  if(x$margins[1]=="cloglog")             cat("\nMARGIN 1: survival with loglog link")  

  if(x$margins[2]=="probit")              cat("\nMARGIN 2: survival with -probit link") 
  if(x$margins[2]=="logit")               cat("\nMARGIN 2: survival with -logit link")
  if(x$margins[2]=="cloglog")             cat("\nMARGIN 2: survival with loglog link")

}
  
    
}


if(type == "gamlss"){

  if(x$margins[1]=="GP")               cat("\nDistribution: generalised Pareto")   
  if(x$margins[1]=="GPII")             cat("\nDistribution: generalised Pareto II")  
  if(x$margins[1]=="GPo")              cat("\nDistribution: generalised Pareto (Orthogonal Parameterisation)")    

  if(x$margins[1]=="DGP")              cat("\nDistribution: discrete generalised Pareto")  
  if(x$margins[1]=="DGPII")            cat("\nDistribution: discrete generalised Pareto II")  
  
  if(x$margins[1]%in%c("N","N2"))      cat("\nDistribution: Gaussian")  
  if(x$margins[1]=="GU")               cat("\nDistribution: Gumbel")    
  if(x$margins[1]=="rGU")              cat("\nDistribution: reverse Gumbel")  
  if(x$margins[1]=="LO")               cat("\nDistribution: logistic")   
  if(x$margins[1]=="LN")               cat("\nDistribution: log-normal") 
  if(x$margins[1]=="WEI")              cat("\nDistribution: Weibull") 
  if(x$margins[1]=="GO")               cat("\nDistribution: Gompertz")   
  if(x$margins[1]=="iG")               cat("\nDistribution: inverse Gaussian") 
  if(x$margins[1]%in%c("GA","GAi"))    cat("\nDistribution: gamma")  
  if(x$margins[1]%in%c("GGA"))         cat("\nDistribution: generalised gamma")  
  if(x$margins[1]%in%c("GA2"))         cat("\nDistribution: gamma (parametrisation 2)")    
  if(x$margins[1]=="FISK")             cat("\nDistribution: Fisk") 
  if(x$margins[1]=="BE")               cat("\nDistribution: beta")    
  if(x$margins[1]=="DAGUM")            cat("\nDistribution: Dagum")  
  if(x$margins[1]=="TW")               cat("\nDistribution: Tweedie")  
  
  
  
  if(x$margins[1]=="SM")               cat("\nDistribution: Singh-Maddala") 
  if(x$margins[1]%in%c("NBI","NBIa"))  cat("\nDistribution: Negative Binomial - Type I") 
  if(x$margins[1]%in%c("NBII","NBIIa"))cat("\nDistribution: Negative Binomial - Type II")
  if(x$margins[1]=="PIG")              cat("\nDistribution: Poisson inverse Gaussian") 
  if(x$margins[1]=="PO")               cat("\nDistribution: Poisson")   
  if(x$margins[1]=="ZTP")              cat("\nDistribution: Zero Truncated Poisson")    
  if(x$margins[1]=="GEVlink")          cat("\nDistribution: Bernoulli") 
  if(x$margins[1]=="probit")           cat("\nFlexible survival model with -probit link")    
  if(x$margins[1]=="logit")            cat("\nFlexible survival model with -logit link")    
  if(x$margins[1]=="cloglog")          cat("\nFlexible survival model with loglog link")    

  

}



if(type == "copSS"){

if(x$margins[1] %in% x$bl && is.null(x$K1)) cat("\nMARGIN 1: Bernoulli") 
if(x$margins[1] %in% x$bl && !is.null(x$K1)) cat("\nMARGIN 1: categorical") 


if(x$margins[2] %in% x$bl) cat("\nMARGIN 2: Bernoulli") 


if(x$margins[2]%in%c("N"))           cat("\nMARGIN 2: Gaussian")  
if(x$margins[2]=="GP")               cat("\nMARGIN 2: generalised Pareto")  
if(x$margins[2]=="GPII")             cat("\nMARGIN 2: generalised Pareto II")  
if(x$margins[2]=="GPo")              cat("\nMARGIN 2: generalised Pareto (Orthogonal Parameterisation)")    

if(x$margins[2]=="DGP")              cat("\nMARGIN 2: discrete generalised Pareto")  
if(x$margins[2]=="DGPII")            cat("\nMARGIN 2: discrete generalised Pareto II")  


if(x$margins[2]=="GU")                 cat("\nMARGIN 2: Gumbel")    
if(x$margins[2]=="rGU")                cat("\nMARGIN 2: reverse Gumbel")  
if(x$margins[2]=="LO")                 cat("\nMARGIN 2: logistic")   
if(x$margins[2]=="LN")                 cat("\nMARGIN 2: log-normal") 
if(x$margins[2]=="WEI")                cat("\nMARGIN 2: Weibull") 
if(x$margins[2]=="GO")                 cat("\nMARGIN 2: Gompertz") 
if(x$margins[2]=="iG")                 cat("\nMARGIN 2: inverse Gaussian") 
if(x$margins[2]%in%c("GA","GAi"))      cat("\nMARGIN 2: gamma")
#if(x$margins[2]%in%c("GGA"))          cat("\nMARGIN 2: generalised gamma")  
#if(x$margins[2]%in%c("GA2"))          cat("\nMARGIN 2: gamma (parametrisation 2)")   
if(x$margins[2]=="BE")                 cat("\nMARGIN 2: beta")    
if(x$margins[2]=="DAGUM")              cat("\nMARGIN 2: Dagum")
if(x$margins[2]=="TW")                 cat("\nMARGIN 2: Tweedie")


if(x$margins[2]=="SM")                 cat("\nMARGIN 2: Singh-Maddala") 
if(x$margins[2]=="FISK")               cat("\nMARGIN 2: Fisk") 
if(x$margins[2]=="FISK2")              cat("\nMARGIN 2: Fisk2") 
if(x$margins[2] %in% c("NBI","NBIa"))  cat("\nMARGIN 2: Negative Binomial - Type I") 
if(x$margins[2]%in% c("NBII","NBIIa")) cat("\nMARGIN 2: Negative Binomial - Type II")
if(x$margins[2]=="PIG")    	       cat("\nMARGIN 2: Poisson inverse Gaussian") 
if(x$margins[2]=="PO")     	       cat("\nMARGIN 2: Poisson")   
if(x$margins[2]=="ZTP")    	       cat("\nMARGIN 2: Zero Truncated Poisson")  

}






         
}

