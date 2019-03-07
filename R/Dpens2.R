Dpens2 <- function(params, type = "lasso", lambda = 1, w.alasso = NULL, gamma = 1, a = 3.7, eps = 1e-10){ 

# "lasso", "alasso", "scad" 

if(type == "lasso"){


#S <- lambda * diag( 1/sqrt(params^2 + eps) )*as.integer(params != 0)


f1 <- function(x) sqrt(x^2 + eps)
f2 <- function(x) x/sqrt(x^2 + eps) 
f3 <- function(x) 0.5 * (2 * (x^2 + eps)^-0.5 - 2 * x * ((x^2 + eps)^-(0.5 + 1) * (0.5 * (2 * x))))


S2 <- lambda*f1(params)
S1 <- lambda*f2(params)
S  <- lambda*f3(params)


}


#if(type == "alasso"){
#
#    if( is.null(w.alasso) ) w.alasso <- 1
#      
#  w.al <- 1/abs(w.alasso)^gamma
#  S    <- lambda * diag( w.al/sqrt(params^2 + eps) )*as.integer(params != 0)
#
#}
#
#
#
#if(type == "scad"){
#
#
#    theta <- abs(params) 
#    p     <- length(params)
#    f1    <- sapply(theta, function(theta) {max (a * lambda - theta, 0) / ((a - 1) * lambda)})
#    f.d   <- lambda * ((theta <= lambda) + f1 * (theta > lambda))
#
#    S <- diag( f.d / ( sqrt(params^2 + eps) + 1e-06 ) )
#
#}


L <- list(S = S, S1 = S1, S2 = S2)
L

}


