lmc <- function(y, X, start.v = NULL, lambda = 1, pen = "none", gamma = 1, a = 3.7){

#stop("Function not ready yet.")

if( !(pen %in% c("ridge", "lasso", "alasso", "scad", "none")) ) stop("pen should be one of: none, ridge, lasso, alasso, scad")

n <- length(y)


################################################################
# Transformation functions
################################################################

unconstr <- function(x, eps = .Machine$double.eps){
                              
               ums <- ( 1 - sum(x[-1]) )
               ums <- ifelse(ums < eps, eps, ums)
               
               log( x[-1]/ums )
               
               }




constr   <- function(x){delta <- exp(c(0,x)); delta/sum(delta) }

check.be <- function(beta, eps = .Machine$double.eps){

       eps <- .Machine$double.eps

    if( any(beta < eps) == TRUE ){		
        beta <- ifelse(beta < eps, eps, beta)
        beta <- beta/sum(beta) # safety measure
		}
beta

}


################################################################
# starting values
################################################################

if(is.null(start.v)){
   s.v      <- lm(y ~ -1 + X)
   nm.c     <- names(coef(s.v))
   w.alasso <- check.be( c( abs(coef(s.v))/sum(abs(coef(s.v))) ) )

}else{

  w.alasso <- check.be( c( abs(start.v)/sum(abs(start.v)) ) )   
  nm.c     <- names(w.alasso)


}



beta.st  <- unconstr(w.alasso) 
D        <- diag(rep(0, length(w.alasso)))
   
################################################################
################################################################
   
func.l <- function(beta.st, X, y, D){

  x <- beta.st
  M <- X
  xpar <- constr(x)
  Mh   <- t(M)%*%M
  
  t(y)%*%y + t(xpar)%*%Mh%*%xpar - 2*t(y)%*%M%*%xpar + t(xpar)%*%(D)%*%xpar
   
}

####

func.g <- function(beta.st, X, y, D){

  x <- beta.st
  M <- X

  xpar <- constr(x)
  
  der <- rbind(-exp(x)/(1+sum(exp(x)))^2, 
              -exp(x)%*%t(exp(x))/(1+sum(exp(x)))^2 +diag(exp(x)/(1+sum(exp(x)))))
                
  t(der)%*%(2*t(M)%*%M%*%xpar - 2*t(M)%*%y +2*D%*%xpar)
 
}

####

func.h <- function(beta.st, X, y, D){

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

        der2  = rbind(2*exp(2*x[i])/(1+sum(exp(x)))^3-exp(x[i])/(1+sum(exp(x)))^2,
                      cbind(2*exp(2*x[i]+x)/(1+sum(exp(x)))^3-exp(x[i]+x)/(1+sum(exp(x)))^2))
        der2[i+1]=  der2[i+1]+ exp(x[i])/(1+sum(exp(x)))-2*exp(2*x[i])/(1+sum(exp(x)))^2
        h[i,i] =t(der[,i])%*%(MM + 2*D)%*%der[,i] + 
          + t(der2)%*%(MM%*%xpar - 2*t(M)%*%y + 2*D%*%xpar)
      } else{

        der2  = rbind(2*exp(x[s])*t(exp(x[i]))/(1+sum(exp(x)))^3,
                      cbind( 2*exp(x)*exp(x[s])*exp(x[i])/(1+sum(exp(x)))^3))
        der2[s+1] = der2[s+1]-exp(x[s])*exp(x[i])/(1+sum(exp(x)))^2
        der2[i+1] = der2[i+1]-exp(x[s])*exp(x[i])/(1+sum(exp(x)))^2
        h[s,i] =t(der[,s])%*%(MM + 2*D)%*%der[,i]+ 
          + t(der2)%*%(MM%*%xpar - 2*t(M)%*%y + 2*D%*%xpar)
      }  }  }
  h 
  
 }




###############################################################
# Optimisation 
###############################################################
# note: lambda inside Dpens and final penalty multiplied by n/2


main.f <- function(beta.st, X, y, lambda, w.alasso, gamma, a){

 beta    <- check.be( constr(beta.st) )
 beta.st <- unconstr(beta) 


 if(pen != "none"){
 
   params <- constr(beta.st) 
   D      <- n/2*Dpens(params, type = pen, lambda, w.alasso, gamma, a)
   
 }

    f <- func.l(beta.st, X, y, D)
    g <- func.g(beta.st, X, y, D)
    h <- func.h(beta.st, X, y, D)
      
    list(value = f, gradient = g, hessian = h)

}


fit <- trust(main.f, beta.st, X = X, y = y, rinit = 1, rmax = 100, blather = TRUE, lambda = lambda, 
             w.alasso = w.alasso, gamma = gamma, a = a)
             
c.coefficients <- constr(fit$argument)
names(c.coefficients) <- nm.c 

###############################################################
###############################################################

if(pen != "none"){ # re-check this theoretically but should be ok
                   # because of simplification 

D  <- n/2*Dpens(constr(fit$argument), type = pen, lambda, w.alasso, gamma, a)
XX <- t(X)%*%X 
Vb <- solve(XX + 2*D)

edf   <- diag(Vb%*%XX) 
t.edf <- sum(edf)

}else{edf <- rep(1, length(fit$argument) + 1 ); t.edf <- sum(edf)}


###############################################################
###############################################################


L <- list(coefficients = fit$argument, c.coefficients = c.coefficients, fit = fit, hess = TRUE,
          fp = FALSE, l.sp1 = 0, l.sp2 = 0, l.sp3 = 0, l.sp4 = 0, l.sp5 = 0, l.sp6 = 0, l.sp7 = 0, l.sp8 = 0,
          iter.if = fit$iterations, logLik = fit$value, n = n, t.edf = t.edf, edf = edf)

class(L) <- c("lmc")

L

}










