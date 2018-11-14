postVb <- function(SemiParFit, VC){

epsilon <- 0.0000001 
Vb.t    <- coef.t <- NULL

He <- SemiParFit$fit$hessian
                                   
                                   
    # replace with new function PDef?                                  
                                   
    He.eig <- eigen(He, symmetric = TRUE)
    if(min(He.eig$values) < sqrt(.Machine$double.eps) && sign( min( sign(He.eig$values) ) ) == -1) He.eig$values <- abs(He.eig$values)  
    if(min(He.eig$values) < sqrt(.Machine$double.eps) ) { pep <- which(He.eig$values < sqrt(.Machine$double.eps)); He.eig$values[pep] <- epsilon }
    Vb <- He.eig$vectors%*%tcrossprod(diag(1/He.eig$values, nrow = length(He.eig$values), ncol =  length(He.eig$values)),He.eig$vectors)  # this could be taken from magic as well but we don't 
    Vb <- (Vb + t(Vb) ) / 2 
    
                                     
if( (VC$l.sp1!=0 || VC$l.sp2!=0 || VC$l.sp3!=0 || VC$l.sp4!=0 || VC$l.sp5!=0 || VC$l.sp6!=0 || VC$l.sp7!=0 || VC$l.sp8!=0) && VC$fp==FALSE){

    HeSh <- He - SemiParFit$fit$S.h
    
    
    if(VC$surv.flex == TRUE){

    est.c <- SemiParFit$fit$argument
    cm <- rep(1, length(est.c))
    
    
    
    
    if(VC$informative == "no"){
    
        if(!is.null(VC$mono.sm.pos)) mono.sm.pos <- VC$mono.sm.pos else mono.sm.pos <- c(VC$mono.sm.pos1, VC$mono.sm.pos2 + VC$X1.d2)  
    
    }
    
    if(VC$informative == "yes") mono.sm.pos <- which(names(est.c) %in% names(est.c[VC$mono.sm.pos])) 
    
    
    cm[mono.sm.pos] <- exp( est.c[mono.sm.pos] )    
    
    Vb.t <- t(cm*t(cm*Vb))  
    coef.t <- SemiParFit$fit$argument
    coef.t[mono.sm.pos] <- exp( coef.t[mono.sm.pos] )      
    
    #invHeSh <- solve(HeSh)
    #
    #invHeSh.t <- t(cm*t(cm*invHeSh))
    #
    #HeSh.t <- solve(invHeSh.t)
    #
    #F.t <- Vb.t%*%HeSh.t 
    #
    #Ve.t <- F.t%*%Vb.t 
    #
    # based on test above, edfs should be the same and 
    # those from magic should be just a bit better maybe
    # 2*diag(F.t) - rowSums(t(F.t)*F.t) gives same res using F and magic
    # should do Ve as well however we only use it
    # for random effects, which are never affected by monotonic smoothing  
    # Ve.t <- t(cm*t(cm*Ve))  
    
    }
    
    
    
    F  <- diag(SemiParFit$magpp$edf )             # Vb%*%HeSh #  this is not correct when fixing sp  
    F1 <- diag(SemiParFit$magpp$edf1)             # needed for testing
    R  <- SemiParFit$bs.mgfit$R                   # needed for testing especially for RE smoothers, it should not be affected by mono smooth
    Ve <- Vb%*%HeSh%*%Vb                          # diag(SemiParFit$magpp$Ve) and diag(SemiParFit$magpp$Vb) but need to be careful with dispersion parameters
    Ve <- (Ve + t(Ve) ) / 2
    
    
    if(VC$robust == TRUE){
    
    SemiParFitT <- SemiParFit 
    
    SemiParFitT$VC <- VC
    SemiParFitT$coefficients <- SemiParFitT$fit$argument 
    SemiParFitT$Vb <- diag(1, length(SemiParFitT$coefficients)) 
    
    Q <- rIC(SemiParFitT)$hbs
    
    Vb <- Vb%*%Q%*%Vb         # this is just freq and not Bayesian but it is ok for now. So Vb and Ve are the same 
    Vb <- (Vb + t(Vb) ) / 2
    
    Ve <- Vb
    
    }
    
    
    
    
    
    
}else{ 

HeSh <- He
Ve <- Vb

    if(VC$robust == TRUE){
    
    SemiParFitT <- SemiParFit 
    
    SemiParFitT$VC <- VC
    SemiParFitT$coefficients <- SemiParFitT$fit$argument 
    SemiParFitT$Vb <- diag(1, length(SemiParFitT$coefficients)) 
    
    Q <- rIC(SemiParFitT)$hbs
    
    Vb <- Ve%*%Q%*%Ve       
    Vb <- (Vb + t(Vb) ) / 2
    
    Ve <- Vb
    
    }

F <- F1 <- diag(rep(1,dim(Vb)[1]))
R <- SemiParFit$bs.mgfit$R

} 

t.edf <- sum(diag(F))

##################################################
# need this for poisson intercept only gamlss case
Vb11 <- F11 <- SemiParFit$fit$hessian
d1 <- dim(Vb11)[1]
d2 <- dim(Vb11)[2]
Vb11[1:d1, 1:d2] <- Vb[1:d1, 1:d2] 
 F11[1:d1, 1:d2] <-  F[1:d1, 1:d2] 
Vb <- Vb11
F <- F11
############

dimnames(SemiParFit$fit$hessian)[[1]] <- dimnames(SemiParFit$fit$hessian)[[2]] <- dimnames(Vb)[[1]] <- dimnames(Vb)[[2]] <- dimnames(HeSh)[[1]] <- dimnames(HeSh)[[2]] <- dimnames(F)[[1]] <- dimnames(F)[[2]] <- dimnames(He)[[1]] <- dimnames(He)[[2]] <- names(SemiParFit$fit$argument)   




list(He = He, Vb = Vb, Vb.t = Vb.t, HeSh = HeSh, F = F, F1 = F1, R = R, Ve = Ve, t.edf = t.edf, SemiParFit = SemiParFit, coef.t = coef.t)

}

