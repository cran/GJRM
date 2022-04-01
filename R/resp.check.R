resp.check <- function(y, margin = "N", 
                           main = "Histogram and Density of Response",
                           xlab = "Response", print.par = FALSE, plots = TRUE, 
                           loglik = FALSE, os = FALSE, 
                           intervals = FALSE, n.sim = 100, prob.lev = 0.05, 
                           i.f = FALSE, 
                           min.dn = 1e-40, min.pr = 1e-16, max.pr = 0.999999, ...){

m2 <- c("N","GU","rGU","LO","LN","WEI","iG","GA","GAi","BE","FISK","GP","GPII","GPo")
m3 <- c("DAGUM","SM","TW")
nu <- NULL
m1d  <- c("PO","ZTP","DGP0")
m2d  <- c("NBI", "NBII","NBIa", "NBIIa","PIG","DGP","DGPII")
m3d  <- c("DEL","SICHEL")

y1m <- NA
y <- as.numeric( na.omit(y) )

if(!(margin %in% c(m2,m3,m1d,m2d)) ) stop("Error in margin value. It should be one of:\nN, GU, rGU, LO, LN, WEI, iG, GA, DAGUM, TW, SM, BE, FISK, NBI, NBII, PIG, PO, ZTP, GP, GPII, GPo, DGP, DGPII, DGP0.") 

if(margin %in% c("LN","WEI","iG","GA","GAi","DAGUM","SM","FISK") && min(y, na.rm = TRUE) <= 0) stop("The response must be positive.")
if(margin %in% c("TW") && min(y, na.rm = TRUE) < 0) stop("The response must be > = 0.")

if(margin %in% c("BE") && (min(y, na.rm = TRUE) <= 0 || max(y, na.rm = TRUE) >= 1) ) stop("The response must be in the interval (0,1).")
   
    if(margin %in% c(m1d,m2d,m3d) && min(y, na.rm = TRUE) < 0) stop("The response must be positive.")
    if(margin %in% c(m1d,m2d,m3d)){
    
    is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
    if(sum(as.numeric(is.wholenumber(y))) != length(y)) stop("The response must be discrete.")     
    }
    if(margin %in% c("ZTP") && min(y, na.rm = TRUE) < 1) stop("The response must be greater than 0.") 
    





if(margin == "LN") y <- log(y)

margins <- c(margin, margin) # not important to chance probit here

VC <- list(X1 = matrix(1, nrow = length(y), ncol = 1), X1.d2 = 1,
           X2 = NULL, X2.d2 = 1,
           X3 = NULL, X3.d2 = 1, robust = FALSE,
           l.sp1 = 0, l.sp2 = 0, l.sp3 = 0, l.sp4 = 0, l.sp5 = 0, l.sp6 = 0, l.sp7 = 0, l.sp8 = 0, l.sp9 = 0, 
           weights = 1, m2 = m2, m3 = m3, m1d = m1d, m2d = m2d, m3d = m3d, 
           margins = margins, fp = TRUE,
           extra.regI = "t", Cont = "NO", ccss = "no", triv = FALSE, surv = FALSE, zero.tol = 1e-02, 
           min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)

ps <- list(S.h = 0, S.h1 = 0, S.h2 = 0)

respvec <- list(y1 = y)
           
if( margin %in% c("PO","ZTP") )         st.v <- c( log( mean((y + mean(y))/2) ) )           
if( margin %in% c("NBI","NBIa","PIG") ) st.v <- c( log(mean((y + mean(y))/2)), log( max( (var(y) - mean(y))/mean(y)^2, 0.1) ) )
if( margin %in% c("NBII","NBIIa") )     st.v <- c( log(mean((y + mean(y))/2)), log( max( (var(y)/mean(y)) - 1, 0.1) ) )    
if( margin %in% c("DEL") )              st.v <- c( log(mean((y + mean(y))/2)), log( max( (var(y) - mean(y))/mean(y)^2, 0.1) ), qlogis(0.5) )  
if( margin %in% c("SICHEL") )           st.v <- c( log(mean((y + mean(y))/2)), log( max( (var(y) - mean(y))/mean(y)^2, 0.1) ), -0.5 )    
if( margin %in% c("N","LN") )           st.v <- c( mean((y + mean(y))/2) ,           log( var(y) ) )  
if( margin %in% c("N2") )               st.v <- c( mean((y + mean(y))/2) ,           log( sqrt(var(y)) ) )  
if( margin %in% c("LO") )               st.v <- c( mean((y + mean(y))/2) ,           log(  3*var(y)/pi^2 ) )  
if( margin %in% c("iG") )               st.v <- c( log( mean((y + mean(y))/2) ) , log( var(y)/mean(y)^3)  )    
if( margin %in% c("GU") )               st.v <- c( mean(y) + 0.57722*sqrt(var(y)/1.64493) ,  log(6*var(y)/pi^2) )    
if( margin %in% c("rGU") )              st.v <- c( mean(y) - 0.57722*sqrt(var(y)/1.64493) ,  log(6*var(y)/pi^2) )   
if( margin %in% c("WEI") )              st.v <- c( log( mean( exp(log(y) + 0.5772/(1.283/sqrt(var(log(y))))) )  ) , log( ( 1.283/sqrt(var(log(y))) )^2 ) ) 
if( margin %in% c("GA") )               st.v <- c( log(mean((y + mean(y))/2)), log(var(y)/mean(y)^2)  ) # log( 1^2 )             
if( margin %in% c("GAi") )              st.v <- c( mean((y + mean(y))/2), log(var(y)/mean(y)^2)  ) # log( 1^2 )
if( margin %in% c("TW") )              {st.v <- st.vTW <- coef(gam(list(y ~ 1, ~ 1, ~ 1), family = twlss()))
                                        st.vTW[2] <- st.v[3]; st.vTW[3] <- st.v[2]
                                        st.v <- st.vTW 
                                        }




#if( margin %in% c("GP","DGP") )         st.v <- c( log(10) ,   0.0025 )  
#if( margin %in% c("GP","DGP") )         st.v <- c( 0.05, log(100)  )  
#if( margin %in% c("GA2") )              st.v <- c( mean(y)^2/var(y), log(mean(y)/var(y))  ) # log( 1^2 )
#if( margin %in% c("GGA") )              st.v <- c( log(mean(y)), log(1), log(1) )



if( margin %in% c("GP","GPII","GPo","DGP","DGPII","DGP0") ){

 est.ob <- try(  gpd.fit(y, threshold = 0, siglink = exp, show = FALSE), silent = TRUE) 
 
 if( inherits(est.ob, "try-error")   ) est.ob <- try( gpd.fit(y, threshold = 0, siglink = exp, show = FALSE, siginit = mean(y + mean(y))/2, shinit = 0.0025), silent = TRUE)  
 if( inherits(est.ob, "try-error")   ) st.v <- c( 0.0025, log( mean((y + mean(y))/2)  )  ) 
 if( !inherits(est.ob, "try-error")  ) st.v <- c( est.ob$mle[2], est.ob$mle[1]  ) 
  
 if( margin %in% c("GPo","GPII","DGPII") ) { if(st.v[1] < 0) st.v[1] <- 0.001  }

 if( margin %in% c("GPII","GPo") ) st.v[1] <- log(st.v[1] + 0.5) 
 if( margin %in% c("DGPII") )      st.v[1] <- log(st.v[1]) # sqrt(st.v[1])
 if( margin %in% c("DGP0") )       st.v    <- st.v[2]
 
 

}                                                              

if( margin %in% c("DAGUM","SM") )       st.v <- c( log(mean((y + mean(y))/2)), log(sqrt(2)), log(1) )   # log(0.01) #  log(sqrt(2))       # 0.1    
if( margin %in% c("FISK") )             st.v <- c( log(mean((y + mean(y))/2)), log(sqrt(2)))    
if( margin %in% c("BE") )               st.v <- c( qlogis(mean((y + mean(y))/2)), qlogis( var(y)/( mean(y)*(1-mean(y)) )  )  )              

if( margin %in% c(m1d) )    names(st.v) <- c("mu.star")
if( margin %in% c(m2,m2d) ) names(st.v) <- c("mu.star", "sigma.star")
if( margin %in% c(m3,m3d) ) names(st.v) <- c("mu.star", "sigma.star", "nu.star")


if(print.par == TRUE && os == TRUE){ 

if( margin %in% c(m1d) )    names(st.v) <- c("mu")
if( margin %in% c(m2,m2d) ) names(st.v) <- c("mu", "sigma")
if( margin %in% c(m3,m3d) ) names(st.v) <- c("mu", "sigma", "nu")

}



if(margin %in% c(m1d,m2,m2d)) univfit <-  try(trust(bprobgHsContUniv, st.v, rinit = 1, rmax = 100, respvec = respvec, 
                                            VC = VC, ps = ps, blather = TRUE), silent = TRUE)
                                            
if(margin %in% c(m3,m3d)) univfit <-  try(trust(bprobgHsContUniv3, st.v, rinit = 1, rmax = 100, respvec = respvec, 
                                            VC = VC, ps = ps, blather = TRUE), silent = TRUE)                                            
                   
if(  inherits(univfit, "try-error")   ) stop("The parameters of the chosen distribution could not be estimated. Try a different distribution.")                  
 
 
    
if(plots == TRUE){ ##       
             

if(margin == "LN") y <- exp(y)


    if(margin %in% c("ZTP","DGP","DGPII","DGP0")){
     
    ly1 <- length(y)
    y1m <- list()
    my1 <- max(y)
    for(i in 1:ly1){ y1m[[i]] <- seq(1, y[i]); length(y1m[[i]]) <- my1} 
    y1m <- do.call(rbind, y1m)  
    
    if(max(y) > 170) y1m <- mpfr( y1m, pmax(53, getPrec(y))) 

     
    }

    


if(margin %in% m2)   pp <-      distrHsAT(y, univfit$argument[1], esp.tr(univfit$argument[2], margin)$vrb, 1, margin2 = margin, min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)
if(margin %in% m3)   pp <-      distrHsAT(y, univfit$argument[1], esp.tr(univfit$argument[2], margin)$vrb, enu.tr(univfit$argument[3], margin)$vrb, margin2 = margin, min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)
if(margin %in% m1d)  pp <- distrHsATDiscr(y, univfit$argument[1], 1, 1, margin2 = margin, y2m = y1m, robust = FALSE, min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)
if(margin %in% m2d)  pp <- distrHsATDiscr(y, univfit$argument[1], esp.tr(univfit$argument[2], margin)$vrb, 1, margin2 = margin, y2m = y1m, robust = FALSE, min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)


p <- pp$p2
d <- pp$pdf2

if(margin %in% c(m1d,m2d)) p <- runif(y, p - d, p)


par(mfrow = c(1, 2))
hist(y, freq = FALSE, ylim=c(0, max(d, hist(y, plot = FALSE)$density) ),
     main=main,
     xlab=xlab, ...)

lines(sort(y),d[order(y)],lwd=2)


if(any(is.na(p)) == TRUE) stop("It is not possible to produce a QQ-plot.\nThe chosen distribution (unconditional on covariates)\nis not probably a good fit.")


if(intervals == FALSE){qqnorm(qnorm(p)); abline(0, 1, col = "red")}
if(intervals == TRUE) {univfit$y2 <- y
                       univfit$Cont <- "NO" 
                       univfit$VC$ccss <- "no"
                       univfit$univar.gamlss <- FALSE
                       univfit$n <- length(y) 
                       int.postcheck(univfit, margin, n.rep = n.sim, prob.lev = prob.lev, y2m = y1m)}




#if(print.par == TRUE) print( univfit$argument )  


                 } ##
                 
                 
                 
                 
if(print.par == TRUE && os == FALSE && i.f == FALSE) print( univfit$argument )                  
         
         
if(print.par == TRUE && os == TRUE && i.f == FALSE){         
         
mu    <- eta.tr(univfit$argument[1], margin)
mupos <- c("LN","WEI","iG","GA","DAGUM","SM","TW","FISK",m1d,m2d,m3d)
mub   <- c("BE")
if(margin %in% mupos) mu <- exp(mu)
if(margin %in% mub)   mu <- plogis(mu)

if(!(margin %in% m1d)) sigma <- NULL else sigma <- esp.tr(univfit$argument[2], margin)$vrb

if(margin %in% m3)  nu <- enu.tr(univfit$argument[3], margin)$vrb   
if(margin %in% m3d) nu <- enu.tr(univfit$argument[3], margin)$vrb   


if(!(margin %in% m1d)) print( c(mu) ) else print( c(mu, sigma, nu) )
                 
}



#####################################################################################
# this is need for estimation in internal routines 
##

if(plots == FALSE && print.par == TRUE && loglik == FALSE && os == FALSE && i.f == TRUE) return( univfit$argument )  

##
#####################################################################################




if(loglik == TRUE){ ##

if(margin == "LN"){ 

if(plots == FALSE) d <- distrHsAT(exp(y), univfit$argument[1], esp.tr(univfit$argument[2], margin)$vrb, 1, margin, min.dn = min.dn, min.pr = min.pr, max.pr = max.pr)$pdf2
lk <- sum(log(d))

                  }



if(margin != "LN") lk <- -univfit$l

attr(lk, "nobs")     <- length(y)
attr(lk, "df")       <- length(st.v)
#attr(lk, "est.pars") <- univfit$argument

class(lk) <- "logLik"
lk

                  } ##
                  
                  



}

