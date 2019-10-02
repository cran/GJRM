ggmtrust <- function(s, n, data = NULL, lambda = 1, pen = "lasso", params = NULL, method = "BHHH", w.alasso = NULL, gamma = 1, a = 3.7){


my.env <- new.env(); my.env$countPD <- 0 


if( !(pen %in% c("ridge", "lasso", "alasso", "scad")) ) stop("pen should be one of: ridge, lasso, alasso, scad")


if( !(method %in% c("H", "BHHH", "TB")) ) stop("method should be one of: H, BHHH, TB")


fun.sk <- function(x, sigma, params, idx){
  a <- ((sigma-x)+t(sigma-x)-diag(diag(sigma-x)))/2
  a <- a[lower.tri(a, diag = TRUE)]
  a[idx] <- a[idx]*params[idx]
  return(a)
}



VC <- list(my.env = my.env, method = method, dat = data, fun.sk = fun.sk, sk = vector(mode = "list"), NLM = FALSE )

 
p   <- ncol(s)
k   <- 0:(p-1) 
idx <- ( k*p - (k-1)*k/2 ) + 1  


if(is.null(params)){

 omega.ini <- diag(1/diag(s)) 
 params    <- omega.ini[lower.tri(omega.ini, diag = TRUE)]
 
}


# make a check on the length

 params[idx] <- log(params[idx])  



parsc <- rep(1, length(params))

miter <- 100000

if( method %in% c("H", "BHHH") ) fit <- trust(ggm.Deriv, params, rinit = 1, rmax = 100, parscale = parsc, s = s, n = n, idx = idx, lambda = lambda, pen = pen, 
                                              VC = VC, blather = TRUE, iterlim = miter, w.alasso = w.alasso, gamma = gamma, a = a)  

if( method %in% c("TB") ) fit <- trust.optim(params, fn = ggm.DerivOPT1, gr = ggm.DerivOPT2, method = "BFGS", control = list(maxit = miter, function.scale.factor = 1, report.level = -1, report.freq = -1), s = s, n = n, idx = idx, 
                                             lambda = lambda, pen = pen, VC = VC, w.alasso = w.alasso, gamma = gamma, a = a)


gc()  


  ## post-estimation ##

if( method %in% c("H") ){ 
   
    He      <- fit$hessian
    S       <- fit$S
    logLik  <- -fit$l
    est.par <- fit$argument
    fCount  <- fit$countPD
    iters   <- fit$iterations

}



if( method %in% c("BHHH") ){ 
   
    est.par <- fit$argument

    VC$method <- "H"
    VC$dat    <- data
    ggmobj <- ggm.Deriv(est.par, s, n, idx, lambda = lambda, pen = pen, VC, w.alasso, gamma, a)

    He      <- ggmobj$hessian
    S       <- ggmobj$S
    logLik  <- -fit$l
    fCount  <- fit$countPD
    iters   <- fit$iterations
    VC$method <- method

}





if( method %in% c("TB") ){ 
   
    est.par <- fit$solution

    VC$method <- "H"
    ggmobj <- ggm.Deriv(est.par, s, n, idx, lambda = lambda, pen = pen, VC, w.alasso, gamma, a)

    He      <- ggmobj$hessian
    S       <- ggmobj$S
    logLik  <- -fit$fval
    fCount  <- NA
    iters   <- fit$iterations
    VC$method <- method

}

gc()  

    
   
    #Vb <- PDef(He)$res.inv  # this option or the alternative below
    
    
    Vb <- chol2inv(chol(He))
    
    
    #He.eig <- eigen(He, symmetric = TRUE)
    #if(min(He.eig$values) < sqrt(.Machine$double.eps) && sign( min( sign(He.eig$values) ) ) == -1) He.eig$values <- abs(He.eig$values)  
    #if(min(He.eig$values) < sqrt(.Machine$double.eps) ) { pep <- which(He.eig$values < sqrt(.Machine$double.eps)); He.eig$values[pep] <- epsilon }
    #Vb <- He.eig$vectors%*%tcrossprod(diag(1/He.eig$values, nrow = length(He.eig$values), ncol =  length(He.eig$values)),He.eig$vectors)     
        
    
    
    
    Vb <- (Vb + t(Vb) )/2 
    HeSh <- He - S 
    
    edf  <- diag(Vb%*%HeSh)
    Ve <- Vb%*%HeSh%*%Vb   # frequentist covariance matrix
    Ve <- (Ve + t(Ve) )/2
    
    t.edf <- sum(edf)

    omega <- matrix(0, p, p)
    
    est.par[idx] <- esp.tr(est.par[idx], "N")$vrb # exp taken here
    
    omega[lower.tri(omega, diag = TRUE)] <- est.par
    omega <- t(omega) + omega - diag(diag(omega))
    


# main outputs here are omega and t.edf, so all good for now



L <- list(coefficients = est.par, fit = fit, omega = omega, edf = edf, t.edf = t.edf, Vb = Vb, Ve = Ve, logLik = logLik, n = n, p = p,
          hess = TRUE, countPD = fCount, idx = idx,
          l.sp1 = 0, l.sp2 = 0, l.sp3 = 0, l.sp4 = 0, l.sp5 = 0, l.sp6 = 0, l.sp7 = 0, fp = FALSE,
          iter.if = iters)

class(L) <- c("ggmtrust")

L

}










