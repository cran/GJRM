lmc <- function(y, X){

################################################################
# Transformation functions
################################################################

unconstr <- function(x){log(x[-1]/(1-sum(x[-1])))}
constr   <- function(x){delta <- exp(c(0,x)); delta/sum(delta) }

################################################################
# starting values
################################################################

s.v <- lm(y ~ -1 + X) 
s.v <- c( abs(coef(s.v))/sum(abs(coef(s.v))) )
beta.st <- unconstr(s.v) 


###############################################################
# Likelihood
###############################################################

func.l <- function(beta.st, X, y){

  x <- beta.st
  M <- X
  xpar <- constr(x)
  Mh   <- t(M)%*%M
  
  t(y)%*%y + t(xpar)%*%Mh%*%xpar - 2*t(y)%*%M%*%xpar

}

###############################################################
# Gradient
###############################################################

func.g <- function(beta.st, X, y){

  x <- beta.st
  M <- X

  xpar <- constr(x)
  
  der <- rbind(-exp(x)/(1+sum(exp(x)))^2, 
              -exp(x)%*%t(exp(x))/(1+sum(exp(x)))^2 +diag(exp(x)/(1+sum(exp(x)))))
                
  t(der)%*%(2*t(M)%*%M%*%xpar - 2*t(M)%*%y)
   
}

###############################################################
# Hessian
###############################################################

func.h <- function(beta.st, X, y){

  x <- beta.st
  M <- X


  xpar = constr(x)
  m    = length(x)
  MM   = 2*t(M)%*%M

  h    = matrix(0,m,m)
  der  = rbind(-exp(x)/(1+sum(exp(x)))^2, 
               -exp(x)%*%t(exp(x))/(1+sum(exp(x)))^2 +diag(exp(x)/(1+sum(exp(x))))  )
               
  for (s in 1:m){
    for (i in 1:m){
      if (s==i){
        #d/dt_i^2
        der2  = rbind(2*exp(2*x[i])/(1+sum(exp(x)))^3-exp(x[i])/(1+sum(exp(x)))^2,
                      cbind(2*exp(2*x[i]+x)/(1+sum(exp(x)))^3-exp(x[i]+x)/(1+sum(exp(x)))^2))
        der2[i+1]=  der2[i+1]+ exp(x[i])/(1+sum(exp(x)))-2*exp(2*x[i])/(1+sum(exp(x)))^2
        h[i,i] =t(der[,i])%*%(MM)%*%der[,i] + 
          + t(der2)%*%(MM%*%xpar - 2*t(M)%*%y)
      } else{
        #d/dt_st_i
        der2  = rbind(2*exp(x[s])*t(exp(x[i]))/(1+sum(exp(x)))^3,
                      cbind( 2*exp(x)*exp(x[s])*exp(x[i])/(1+sum(exp(x)))^3))
        der2[s+1] = der2[s+1]-exp(x[s])*exp(x[i])/(1+sum(exp(x)))^2
        der2[i+1] = der2[i+1]-exp(x[s])*exp(x[i])/(1+sum(exp(x)))^2
        h[s,i] =t(der[,s])%*%(MM)%*%der[,i]+ 
          + t(der2)%*%(MM%*%xpar - 2*t(M)%*%y)
      }  }  }
      
  h
  
}


###############################################################
# Optimisation 
###############################################################

main.f <- function(beta.st, X, y){

    f <- func.l(beta.st, X, y)
    g <- func.g(beta.st, X, y)
    h <- func.h(beta.st, X, y)
      
    list(value = f, gradient = g, hessian = h)
}

fit <- trust(main.f, beta.st, X = X, y = y, rinit = 1, rmax = 100, blather = TRUE)

L <- list(coefficients = fit$argument, c.coefficients = constr(fit$argument), fit = fit, hess = TRUE,
          fp = FALSE, l.sp1 = 0, l.sp2 = 0, l.sp3 = 0, l.sp4 = 0, l.sp5 = 0, l.sp6 = 0, l.sp7 = 0, l.sp8 = 0,
          iter.if = fit$iterations)

class(L) <- c("lmc")

L

}










