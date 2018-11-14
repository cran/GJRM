penAprL <- function(params, epsilon = 1e-08){ # , type = "lasso1", alpha = 10000){




#if(type = "lasso1"){

pen0 <- params^2/sqrt(params^2 + epsilon)
pen1 <- (params^3 + 2*params*epsilon)/(params^2 + epsilon)^1.5  
pen2 <- epsilon*(2*epsilon - params^2)/(params^2 + epsilon)^2.5 

#}



#if(type = "lasso2"){
#1/alpha*( log( 1 + exp(-alpha*params) ) + log( 1 + exp(alpha*params) ) )
#pen0 <- 1/alpha*( log( 1 + exp(-alpha*params) ) + log( 1 + exp(alpha*params) ) )
#pen1 <- 1/alpha * (exp(alpha * params) * alpha/(1 + exp(alpha * params)) - exp(-alpha * params) * alpha/(1 + exp(-alpha * params)))
#pen2 <- 2*alpha*exp(alpha*params) / ( exp(alpha*params) + 1 )^2
#}



L <- list(pen0 = pen0, pen1 = pen1, pen2 = pen2)

L

}


