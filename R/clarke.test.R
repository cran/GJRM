clarke.test <- function(obj1, obj2, sig.lev = 0.05){

    l1 <- obj1$logLik
    l2 <- obj2$logLik
    n  <- obj1$n
    n1 <- obj2$n
    p1 <- obj1$t.edf
    p2 <- obj2$t.edf    
    
    if(n != n1)  stop("The two competing models have different sample sizes.")   
    if(l1 == l2) stop("The two competing models have identical log-likelihoods!")
         
    li12  <- obj1$fit$l.par - obj2$fit$l.par
    
    li12b <- li12 - (p1 - p2)/(2 * n) * log(n)
    
    b <- sum(li12b > 0)
    
    db <- "Neither model is significantly preferred"
    
    if(b >= n/2){
        pvalue <- 2 * (1 - pbinom(b - 1, n, 0.5))
        if(pvalue <= sig.lev) db <- "Model 1 is preferred"
    }
    if(b < n/2){
        pvalue <- 2 * (pbinom(b, n, 0.5))
        if(pvalue <= sig.lev) db <- "Model 2 is preferred"
               }
                              
    cat("\n", db, "\n\n", sep = "")
    
    return(invisible(c(b, pvalue)))

}
