vuong.test <- function(obj1, obj2, sig.lev = 0.05){

    l1 <- obj1$logLik
    l2 <- obj2$logLik
    n  <- obj1$n
    n1 <- obj2$n
    if(n != n1)  stop("The two competing models have different sample sizes.")   
    if(l1 == l2) stop("The two competing models have identical log-likelihoods!")
    
    p1 <- obj1$t.edf
    p2 <- obj2$t.edf
    
    li12 <- obj1$fit$l.par - obj2$fit$l.par
    w <- sqrt(var(li12) * (n - 1))
    
    vt <- (l1 - l2 - (p1 - p2)/2 * log(n))/w
    
    if(abs(vt) <= qnorm(1 - sig.lev/2)) dvt <- "Neither model is significantly preferred"
    if(    vt   > qnorm(1 - sig.lev/2)) dvt <- "Model 1 is preferred"
    if(    vt  < -qnorm(1 - sig.lev/2)) dvt <- "Model 2 is preferred"

    cat("\n", dvt, "\n\n", sep = "")
    
    return(invisible(as.numeric(vt)))
    
}

