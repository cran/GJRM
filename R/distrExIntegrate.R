distrExIntegrate <- function (f, lower, upper, subdivisions = 100, rel.tol = .Machine$double.eps^0.25, 
    abs.tol = rel.tol, stop.on.error = TRUE, distr, order = 5000, ..., diagnostic = FALSE){


# this function is from the distrEx R package, Ruckdeschel P, Kohl M, Stabla T, Camphausen F (2006). S4 Classes for Distributions. R News, 6(2), 2-6.
# some amendments have been carried out as not all functionalities are required here

distrExOpt <- list(GLIntegrateTruncQuantile = 2.220446e-15,GLIntegrateOrder = 5000, MCIterations = 1e+05, ElowerTruncQuantile = 1e-07, EupperTruncQuantile = 1e-07, ErelativeTolerance = 0.0001220703, m1dfLowerTruncQuantile = 0, m1dfRelativeTolerance = 0.0001220703, m2dfLowerTruncQuantile = 0, m2dfRelativeTolerance = 0.0001220703, nDiscretize = 100, hSmooth = 0.05, IQR.fac = 15, propagate.names.functionals = TRUE)


GLIntegrate <- function(f, lower, upper, order = 500, ...){
    #if (order %in% c(50, 100, 400, 500, 800, 1000, 4000, 5000, 
    #    8000, 10000, 40000, 50000, 80000, 1e+05)) 
    
    AW <- getFromNamespace(paste(".AW", as.character(order), sep = "."), ns = "distrEx")
    
    #else AW <- .GLaw(order)
    
    xl <- (upper - lower)/2
    W <- xl * AW[, 2]
    A <- xl * AW[, 1] + (lower + upper)/2
    res <- W * c(f(A, ...))
    sum(res)
}




.filterFunargs <- function (dots, fun, neg = FALSE){
    if (length(dots) == 0) 
        return(NULL)
    formFunNames <- names(formals(fun))
    if (neg) 
        return(dots[!names(dots) %in% formFunNames])
    else return(dots[names(dots) %in% formFunNames])
}

    mc <- match.call()
    time <- proc.time()
    ppt <- function(y) {
        if (!is.na(y[4L])) 
            y[1L] <- y[1L] + y[4L]
        if (!is.na(y[5L])) 
            y[2L] <- y[2L] + y[5L]
        paste(formatC(y[1L:3L]), collapse = " ")
    }
    dots <- list(...)
    dotsFun <- .filterFunargs(dots, f)
    funwD <- function(x) do.call(f, c(list(x), dotsFun))
    on.exit(message("Timing stopped at: ", ppt(proc.time() - 
        time)))
    dotsInt <- if (length(names(dots))) 
        dots[names(dots) %in% names(formals(integrate))]
    else NULL
    res <- try(do.call(integrate, c(list(funwD, lower = lower, 
        upper = upper, rel.tol = rel.tol, abs.tol = abs.tol, 
        stop.on.error = stop.on.error), dotsInt)), silent = TRUE)
    if (!is(res, "try-error")) {
        val <- res$value
        if (diagnostic) {
            diagn <- list(call = mc, method = "integrate", args = c(list(lower = lower, 
                upper = upper, rel.tol = rel.tol, abs.tol = abs.tol, 
                stop.on.error = stop.on.error), list(...)), result = res)
            res <- val
        }
        else res <- val
    }
    else {
        errmess <- res
        Zi <- 1
        if (lower >= upper) {
            lo <- lower
            lower <- upper
            upper <- lo
            Zi <- -1
        }
        
        
    #    if (!missing(distr)) {
    #        q.lDots <- NULL
    #        if (length(names(dots))) {
    #            q.lDots <- dots[names(dots) %in% names(formals(q.l(distr)))]
    #            q.lDots[["p"]] <- q.lDots[["lower.tail"]] <- NULL
    #        }
    #    }
        
        
        if(!is.finite(lower)) 
            
        if(missing(distr)) stop(res)
            
            #else {
            #    lower <- do.call(q.l(distr), c(list(distrExOpt$GLIntegrateTruncQuantile), 
            #      q.lDots))
            #}
            
        if(!is.finite(upper)) 
        
        if(missing(distr)) stop(res)
          
          #  else {
          #      q.lArgs <- if ("lower.tail" %in% names(formals(distr@q))) 
          #        list(p = distrExOpt$GLIntegrateTruncQuantile, 
          #          lower.tail = FALSE)
          #      else list(p = 1 - distrExOpt$GLIntegrateTruncQuantile)
          #      q.lArgs <- c(q.lArgs, q.lDots)
          #      upper <- do.call(q.l(distr), q.lArgs)
          #  }
          
        dotsGLInt <- NULL
        if (length(names(dots))) 
            dotsGLInt <- dots[names(dots) %in% names(formals(GLIntegrate))]
        res <- Zi * do.call(GLIntegrate, c(list(f = funwD, lower = lower, 
            upper = upper, order = order), dotsGLInt))
        if (diagnostic) {
            diagn <- list(call = mc, method = "GLIntegrate", 
                args = c(list(lower = lower, upper = upper, order = order), 
                  list(...)), result = list(GLIresult = res, 
                  errorMessage = errmess), distrExOptions = distrExOpt)
        }
    }
    new.time <- proc.time()
    on.exit()
    if (diagnostic) {
        diagn$time <- structure(new.time - time, class = "proc_time")
        attr(res, "diagnostic") <- diagn
        class(attr(res, "diagnostic")) <- "DiagnosticClass"
    }
    return(res)
}









































