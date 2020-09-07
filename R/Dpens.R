Dpens <- function(params, type = "lasso", lambda = 1, w.alasso = NULL, gamma = 1, a = 3.7, eps = 1e-08){ 



# "lasso", "alasso", "scad" 



if(type == "ridge"){

 S <- lambda * diag( rep(1, length(params) ) )
 
 
}


if(type == "lasso"){

 S <- lambda * diag( 1/sqrt(params^2 + eps) )#*as.integer(params != 0)

}


if(type == "alasso"){

    if( is.null(w.alasso) ) w.alasso <- 1
      
  w.al <- 1/abs(w.alasso)^gamma
  S    <- lambda * diag( w.al/sqrt(params^2 + eps) )#*as.integer(params != 0)

}



if(type == "scad"){


    theta <- abs(params) 
    p     <- length(params)
    f1    <- sapply(theta, function(theta) {max (a * lambda - theta, 0) / ((a - 1) * lambda)})
    f.d   <- lambda * ((theta <= lambda) + f1 * (theta > lambda))

    S <- diag( f.d / ( sqrt(params^2 + eps) + 1e-06 ) )

}


S

}


