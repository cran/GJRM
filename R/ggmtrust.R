ggmtrust <- function(s, n, lambda = 1, pen = "lasso", params = NULL){

my.env <- new.env(); my.env$countPD <- 0 

VC <- list(my.env = my.env)

 
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


fit <- trust(ggm.Deriv, params, rinit = 1, rmax = 100, parscale = parsc, s = s, n = n, idx = idx, 
             lambda = lambda, pen = pen, VC = VC, blather = TRUE, iterlim = 100)  

  ## post-estimation ##

    
    He <- fit$hessian
    Vb <- PDef(He)$res.inv    
    Vb <- (Vb + t(Vb) )/2 
    HeSh <- He - fit$S # S is always here assuming that a penalty is always employed
    
    edf  <- diag(Vb%*%HeSh)
    Ve <- Vb%*%HeSh%*%Vb   # frequentist covariance matrix
    Ve <- (Ve + t(Ve) )/2
    
    t.edf <- sum(edf)
    logLik <- -fit$l
    p <- fit$p
    
    omega <- matrix(0, p, p)
    
    est.par <- fit$argument
    est.par[idx] <- esp.tr(est.par[idx], "N")$vrb 
    
    omega[lower.tri(omega, diag = TRUE)] <- est.par
    omega <- t(omega) + omega - diag(diag(omega))
    

# main outputs here are omega and t.edf, so al god for now

L <- list(coefficients = fit$argument, fit = fit, omega = omega, edf = edf, t.edf = t.edf, Vb = Vb, Ve = Ve, logLik = logLik, n = n, p = p,
          hess = TRUE, countPD = fit$countPD, idx = idx,
          l.sp1 = 0, l.sp2 = 0, l.sp3 = 0, l.sp4 = 0, l.sp5 = 0, l.sp6 = 0, l.sp7 = 0, fp = FALSE,
          iter.if = fit$iterations)

class(L) <- c("ggmtrust")

L

}










