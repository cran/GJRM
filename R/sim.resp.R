sim.resp <- function(margin, rsim, eta, sigma2, nu, setseed = TRUE){



#if(margin == "GO")    y <- rgompertz(rsim, shape = exp(eta), rate = sqrt(sigma2))
#if(margin == "GA2")   y <- rgamma(rsim, shape = sqrt(sigma2), rate = exp(eta) ) 
#if(margin == "GGA")   y <- exp( log(rgamma(rsim, shape = nu))/sqrt(sigma2) + log(exp(eta)) )


if(margin == "ZTP") rZTP <- function(n, mu) qpois(runif(n, dpois(0, mu), 1), mu)

if(margin == "DGP"){

       rDGP <- function(n, mu, sigma){

                minp <- distrHsATDiscr(0, mu, sigma^2, NULL, "DGP", NULL, robust = TRUE)$pdf2
                
                p    <- runif(rsim)
                
                indv <- p <= minp
                
                q <- ifelse(indv == TRUE, 0, ceiling( sigma/mu*( (1 - p)^(-mu) - 1   ) - 1 ))
                
                q
                
                                      }
}


if(setseed == TRUE) set.seed(1)

if(margin == "N")     y <- rNO(   rsim,    mu =     eta,     sigma = sqrt(sigma2)) 
if(margin == "N2")    y <- rNO(   rsim,    mu =     eta,     sigma =      sigma2) 
if(margin == "GU")    y <- rGU(   rsim,    mu =     eta,     sigma = sqrt(sigma2)) 
if(margin == "rGU")   y <- rRG(   rsim,    mu =     eta,     sigma = sqrt(sigma2)) 
if(margin == "LO")    y <- rLO(   rsim,    mu =     eta,     sigma = sqrt(sigma2)) 
if(margin == "LN")    y <- rLOGNO(rsim,    mu =     eta,     sigma = sqrt(sigma2)) 
if(margin == "WEI")   y <- rWEI(  rsim,    mu = exp(eta),    sigma = sqrt(sigma2)) 
if(margin == "iG")    y <- rIG(   rsim,    mu = exp(eta),    sigma = sqrt(sigma2)) 
if(margin == "GA")    y <- rGA(   rsim,    mu = exp(eta),    sigma = sqrt(sigma2)) 
if(margin == "GAi")   y <- rGA(   rsim,    mu =     eta,     sigma = sqrt(sigma2))


if(margin == "GP"){

             if(rsim == 1) y <- rgpd(  rsim, loc = 0, shape = eta, scale = sqrt(sigma2))
             
             if(rsim > 1){
                y <- NA
                if(length(sigma2) == 1) sigma2 <- rep(sigma2, rsim) 
                for(i in 1:rsim) y[i] <- rgpd(  1, loc = 0, shape = eta[i], scale = sqrt(sigma2[i]))
             }
}


if(margin == "DAGUM") y <- rGB2(  rsim,    mu = exp(eta),    sigma = sqrt(sigma2), nu = nu, tau = 1) 
if(margin == "SM")    y <- rGB2(  rsim,    mu = exp(eta),    sigma = sqrt(sigma2), nu = 1 , tau = nu) 
if(margin == "BE")    y <- rBE(   rsim,    mu = plogis(eta), sigma = sqrt(sigma2)) 
if(margin == "FISK")  y <- rGB2(  rsim,    mu = exp(eta),    sigma = sqrt(sigma2), nu = 1 , tau = 1)
if(margin == "NBI")   y <- rNBI(  rsim,    mu = exp(eta),    sigma = sqrt(sigma2)) 
if(margin == "NBII")  y <- rNBII( rsim,    mu = exp(eta),    sigma = sqrt(sigma2)) 
if(margin == "PIG")   y <- rPIG(  rsim,    mu = exp(eta),    sigma = sqrt(sigma2)) 
if(margin == "PO")    y <- rPO(   rsim,    mu = exp(eta)) 
if(margin == "ZTP")   y <- rZTP(  rsim,    mu = exp(eta)) 
if(margin == "DGP")   y <- rDGP( rsim,    mu = eta,         sigma = sqrt(sigma2)) 






if(margin %in% c("probit", "logit", "cloglog"))  y <- rbinom(rsim, 1, prob = probm(eta, margin)$pr )
                                                

if(setseed == TRUE) rm(list = ".Random.seed", envir = globalenv()) 

y

}    